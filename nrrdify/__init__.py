#!/usr/bin/env python

# ========================================================================
#  Copyright Het Nederlands Kanker Instituut - Antoni van Leeuwenhoek
#
#  Licensed under the 3-clause BSD License
# ========================================================================

import logging
import os

import pydicom
import SimpleITK as sitk
import tqdm

from . import dicomvolume

logger = logging.getLogger('nrrdify')
counter = 0
post_processing = None

if len(logger.handlers) == 0:
  print('Adding handler for logger')
  handler = logging.StreamHandler()
  formatter = logging.Formatter('[%(asctime)-.19s] %(levelname)-.1s: %(message)s')
  handler.setFormatter(formatter)
  handler.setLevel(logging.INFO)

  logger.addHandler(handler)
  logger.setLevel(logging.INFO)


def walk_folder(source,
                destination,
                filename=None,
                fileformat='nrrd',
                overwrite=False,
                just_check=False,
                process_per_folder=False,
                structure=None,
                output_writer=None,
                dump_protocol=False):
  global counter, logger, post_processing
  if not os.path.isdir(source):
    logger.error('Source directory (%s) does not exist! Exiting...', source)
    return
  if not os.path.isdir(destination):
    logger.error('Destination directory (%s) does not exist! Exiting...', destination)
    return
  counter = 0
  logger.info('Input (%s) and output (%s) valid, scanning input folder for DICOM files', source, destination)
  datasets = {}  # Holds the dicom files, sorted by series UID ({seriesUID: [files]})
  for curdir, dirnames, fnames in os.walk(source):
    if len(fnames) > 0:  # Only process folder if it contains files
      logger.info('Processing folder %s', curdir)

      with tqdm.tqdm(fnames, desc='Processing files') as bar:  # Progress reporting
        for fname in bar:  # for each file in current folder
          try:
            # Check if it contains a valid DICOM header (first 4 bytes = DICM)
            with open(os.path.join(curdir, fname), mode='rb') as openFile:
              openFile.seek(128)
              header = openFile.read(4)
              if header.decode() != 'DICM':
                # Not a valid DICOM file, skip to next
                continue  # Go to next file

            # Load dicom file using PyDicom (needed for name extraction, sorting of series and slices)
            dicfile = pydicom.read_file(os.path.join(curdir, fname), stop_before_pixels=True)

            imagetype = getattr(dicfile, 'ImageType', None)
            sop_class = getattr(dicfile, 'SOPClassUID', None)  # Check if it is a dicomfile containing an image
            series_uid = getattr(dicfile, 'SeriesInstanceUID', None)  # Series UID
            modality = getattr(dicfile, 'Modality', '')

            if series_uid is None:
              continue  # Error cannot sort, so skip and go To next file
            if sop_class is None:
              continue  # not image dicom file, so skip and go to next file

            if modality == 'RTDOSE':
              logger.debug('Dicom file %s is RT dose storage', os.path.join(curdir, fname))
            elif 'Image Storage' not in sop_class.name:  # Not RT Dose, is the sopclass valid?
                continue
            elif imagetype is None:  # Sop class is valid, does it contain image type?
                logger.debug("missing Image Type tag in dicom file %s", os.path.join(curdir, fname))
                continue  # Error cannot sort, so skip and go To next file

            if imagetype is None:
              imagetype = 'RTDOSE'
            else:
              imagetype = tuple(imagetype)

            if series_uid not in datasets:
              datasets[series_uid] = {}

            if imagetype not in datasets[series_uid]:
              datasets[series_uid][imagetype] = dicomvolume.DicomVolume(post_processing)

            datasets[series_uid][imagetype].addSlice(dicfile)
          except KeyboardInterrupt:
            return
          except:
            logger.error('DOH!! Something went wrong!', exc_info=True)
      if process_per_folder:
        if structure == 'source':
          dest = os.path.join(destination, os.path.relpath(curdir, source))
          if not os.path.isdir(dest):
            logger.debug('Creating output directory "%s"', dest)
            os.makedirs(dest)
        else:
          dest = destination

        _processResults(datasets, dest, filename, fileformat, overwrite, just_check, structure == 'dicom', output_writer, dump_protocol)
        datasets = {}

  if not process_per_folder:
    _processResults(datasets, destination, filename, fileformat, overwrite, just_check, structure == 'dicom', output_writer, dump_protocol)


def _processResults(datasets, destination, filename, fileformat, overwrite, just_check, mkdirs, output_writer, dump_protocol):
    global logger
    if just_check:
      for ds in datasets:
        for volume_idx, volume in enumerate(datasets[ds].values()):
          checkVolume(volume, ds, volume_idx)
    else:
      # Done scanning files, now make some NRRDs out of them!
      logger.info('Input folder scanned, found %d unique DICOM series', len(datasets))
      if len(datasets) > 1:  # If more than 1 series is found, a custom filename is not possible
        filename = None
      for ds in datasets:  # Multiple datasets, so generate name from DICOM
        for volume_idx, volume in enumerate(datasets[ds].values()):
          processVolume(volume, destination, filename, fileformat, overwrite, volume_idx, mkdirs, output_writer, dump_protocol)


def processVolume(volume,
                  destination,
                  filename=None,
                  fileformat='nrrd',
                  overwrite=False,
                  file_idx=None,
                  mkdirs=False,
                  output_writer=None,
                  dump_protocol=False):
  global counter, logger
  try:
    if len(volume.slices) == 0:  # No files for this series UID (maybe not image storage?)
      logger.debug('No files for this series...')
      return

    patient_name = str(getattr(volume[0], 'PatientName', '')).split('^')[0]
    study_date = getattr(volume[0], 'StudyDate', '19000101')
    series_description = getattr(volume[0], 'SeriesDescription', 'Unkn')
    series_number = getattr(volume[0], 'SeriesNumber', -1)

    if filename is None:  # Generate a filename from DICOM metadata
      filename = volume.build_filename()

    if file_idx is not None and file_idx > 0:
      filename = '%s (%d)' % (filename, file_idx)

    if mkdirs:
      # Remove invalid characters from dir names
      patient_dir = patient_name
      study_dir = study_date
      for c in r'[]/\;,><&*:%=+@!#^()|?^':
        patient_dir = patient_dir.replace(c, '')
        study_dir = study_dir.replace(c, '')

      destination = os.path.join(destination, patient_dir, study_dir)

    if volume.check_4D():
      if not process_4D(volume,
                        patient_name,
                        study_date,
                        series_number,
                        series_description,
                        destination,
                        filename,
                        fileformat,
                        overwrite,
                        output_writer):
        logger.warning("Volume is 4D, but was not able to split (patient %s, studydate %s series %d. %s), skipping...",
                       patient_name, study_date, series_number, series_description)
        return
    else:
      im = volume.getSimpleITKImage()
      if im is None:
        return

      logger.info('Generating NRRD for pt %s, studydate %s, series %s:%s (%i slices)',
                  patient_name, study_date, series_number, series_description, len(volume.slices))

      _store_image(im, destination, filename, fileformat, patient_name, study_date, len(volume.slices), overwrite, output_writer)

    if dump_protocol:
      protocol_fname = os.path.join(destination, filename) + '_protocol.txt'
      volume.writeProtocol(protocol_fname)
    counter += 1

  except KeyboardInterrupt:
    raise
  except:
    logger.error('Oh Oh... something went wrong...', exc_info=True)


def process_4D(volume,
               patient_name,
               study_date,
               series_number,
               series_description,
               destination,
               filename=None,
               fileformat='nrrd',
               overwrite=False,
               output_writer=None):
  if getattr(volume.slices[0], 'DiffusionBValue', None) is not None:
    # Volume is DWI
    logger.info('Volume is 4D, attempting DWI splitting on standard bvalue tag')
    b_count = 0
    for b_val, b_im, sliceCount in volume.getSimpleITK4DImage(max_value=2000):
      b_count += 1
      if b_im is None:
        continue

      logger.info('Generating NRRD for pt %s, studydate %s, series %s:%s, bvalue %s (%i slices)',
                  patient_name, study_date, series_number, series_description, b_val, sliceCount)

      b_fname = filename + '_b' + str(b_val)
      _store_image(b_im, destination, b_fname, fileformat, patient_name, study_date, sliceCount, overwrite,
                   output_writer)
    if b_count == 0:
      logger.warning('4D volume contains standard bvalue tag, but unable to perform correct split.')
    else:
      logger.debug('DWI volume processing complete, found %i b values', b_count)
      return True

  if 0x0019100c in volume.slices[0]:
    # Volume is DWI
    logger.info('Volume is 4D, attempting DWI splitting on private bvalue tag')
    b_count = 0
    for b_val, b_im, sliceCount in volume.getSimpleITK4DImage(max_value=2000, splitTag=0x0019100c):
      b_count += 1
      if b_im is None:
        continue

      logger.info('Generating NRRD for pt %s, studydate %s, series %s:%s, bvalue %s (%i slices)',
                  patient_name, study_date, series_number, series_description, b_val, sliceCount)

      b_fname = filename + '_b' + str(b_val)
      _store_image(b_im, destination, b_fname, fileformat, patient_name, study_date, sliceCount, overwrite,
                   output_writer)
    if b_count == 0:
      logger.warning('4D volume contains private bvalue tag, but unable to perform correct split.')
    else:
      logger.debug('DWI volume processing complete, found %i b values', b_count)
      return True

  if getattr(volume.slices[0], 'TemporalPositionIdentifier', None) is not None:
    # Volume is 4D on Temporal Position
    logger.info('Volume is 4D, attempting splitting on Temporal Position tag')
    pos_count = 0
    for pos, temp_im, sliceCount in volume.getSimpleITK4DImage(max_value=2000, splitTag=0x00200100):
      pos_count += 1
      if temp_im is None:
        continue

      logger.info('Generating NRRD for pt %s, studydate %s, series %s:%s, temporal position %s (%i slices)',
                  patient_name, study_date, series_number, series_description, pos, sliceCount)

      pos_fname = filename + '_' + str(pos)
      _store_image(temp_im, destination, pos_fname, fileformat, patient_name, study_date, sliceCount, overwrite,
                   output_writer)
    if pos_count == 0:
      logger.warning('4D volume contains Temporal Position tag, but unable to perform correct split.')
    else:
      logger.debug('4D volume processing complete, found %i temporal positions', pos_count)
      return True

  return False


def _store_image(im, destination, fname, fileformat, patient, studydate, slicecount, overwrite=False, output_writer=None):
  global logger
  target = os.path.join(destination, fname)

  nrrd_fname = target + '.' + fileformat

  if destination != '' and not os.path.isdir(destination):
    logger.debug('Creating study directory "%s"', destination)
    os.makedirs(destination)
  elif os.path.isfile(fname):
    if overwrite:
      logger.warning('file "%s" already exists, overwriting...', fname)
    else:
      logger.warning('file "%s" already exists, skipping...', fname)
      return

  logger.info('Storing in %s', nrrd_fname)
  sitk.WriteImage(im, str(nrrd_fname))

  if output_writer is not None:
    logger.debug('Storing location in CSV output')
    output_writer.writerow([counter, patient, studydate, nrrd_fname.replace(os.path.sep, '/'), slicecount])


def checkVolume(dicomVolume, uid, volume_idx=0):
  global logger
  try:
    if len(dicomVolume.slices) == 0:  # No files for this series UID (maybe not image storage?)
      logger.debug('No files for this series...')
      return

    dicomVolume.sortSlices()
    if dicomVolume.is_equidistant and dicomVolume.is_valid:
      if volume_idx > 0:
        logger.info('DicomVolume %s, (volume %d) is valid...', uid, volume_idx + 1)
      else:
        logger.info('DicomVolume %s is valid...', uid)
  except KeyboardInterrupt:
    raise
  except:
    logger.error('Oh Oh... something went wrong...', exc_info=True)
