#!/usr/bin/env python

# ========================================================================
#  Copyright Het Nederlands Kanker Instituut - Antoni van Leeuwenhoek
#
#  Licensed under the 3-clause BSD License
# ========================================================================

import logging
import os
import re
import sys

import pydicom
import SimpleITK as sitk
import tqdm

from . import dicomvolume

logger = logging.getLogger('nrrdify')
counter = 0
post_processing = None
splitter_4D = None

if len(logger.handlers) == 0:
  print('Adding handler for logger')
  handler = logging.StreamHandler()
  formatter = logging.Formatter('[%(asctime)-.19s] %(levelname)-.1s: %(message)s')
  handler.setFormatter(formatter)
  handler.setLevel(logging.INFO)

  logger.addHandler(handler)
  logger.setLevel(logging.INFO)


def walk_folder(source, destination, structure=None, process_per_folder=False, **kwargs):
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

      with tqdm.tqdm(fnames, desc='Processing files', file=sys.stdout) as bar:  # Progress reporting
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
                imagetype = "UNK"

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

        _processResults(datasets, dest, mkdirs=structure == 'dicom', **kwargs)
        datasets = {}

  if not process_per_folder:
    _processResults(datasets, destination, mkdirs=structure == 'dicom', **kwargs)


def _processResults(datasets, destination, mkdirs, **kwargs):
    global logger
    # Done scanning files, now make some NRRDs out of them!
    logger.info('Input folder scanned, found %d unique DICOM series', len(datasets))
    if len(datasets) > 1:  # If more than 1 series is found, a custom filename is not possible
      kwargs['filename'] = None

    for ds in datasets:  # Multiple datasets, so generate name from DICOM
      volume_idx = 0
      for volume in datasets[ds].values():
        sub_volumes = {}
        split_3D_keys = kwargs.get('split_3D', [])

        volume.sortSlices()
        if volume.is_valid and not volume.is_equidistant and len(split_3D_keys) > 0:
          logger.info('Splitting volume by keys: %s', split_3D_keys)
          for slice in volume.slices:
            set_id = []
            for k in split_3D_keys:
              set_id.append(str(getattr(slice, k, None)))
            if len(set_id) == 1:
              set_id = set_id[0]
            else:
              set_id = tuple(set_id)

            if set_id not in sub_volumes:
              sub_volumes[set_id] = dicomvolume.DicomVolume(post_processing)

            sub_volumes[set_id].addSlice(slice)
        else:
          sub_volumes[0] = volume

        for v in sub_volumes:
          if kwargs.get('just_check', False):
            checkVolume(sub_volumes[v], ds, volume_idx)
          else:
            processVolume(sub_volumes[v], destination, volume_idx, mkdirs, **kwargs)

          volume_idx += 1


def processVolume(volume, destination, file_idx=None, mkdirs=False, **kwargs):
  global counter, logger
  try:
    filename = kwargs.get('filename', None)
    fileformat = kwargs.get('fileformat', '.nrrd')
    overwrite = kwargs.get('overwrite', False)
    output_writer = kwargs.get('output_writer', None)
    dump_protocol = kwargs.get('dump_protocol', False)

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
      t_count = 0
      if splitter_4D:
        for ext_4D, im_4D, sliceCount in splitter_4D(volume):
          t_count += 1
          logger.info('Generating NRRD for pt %s, studydate %s, series %s:%s, temporal value %s (%i slices)',
                      patient_name, study_date, series_number, series_description, ext_4D, sliceCount)
          fname_4D = '%s_%s' % (filename, ext_4D)
          store_name = _store_image(im_4D, destination, fname_4D, fileformat, overwrite)
          if output_writer is not None and store_name is not None:
            logger.debug('Storing location in CSV output')
            output_writer.writerow([
              counter,
              patient_name,
              study_date,
              store_name.replace(os.path.sep, '/'),
              sliceCount,
              volume.getOrientation()
            ])
      if t_count > 0:
        logger.debug('4D volume processing complete, found %i temporal positions', t_count)
      elif hasattr(volume.slices[0], 'DiffusionBValue'):
        # Volume is DWI
        logger.debug('Volume is 4D, attempting DWI splitting on standard bvalue tag')
        b_count = 0
        for b_val, b_im, sliceCount in volume.getSimpleITK4DImage(max_value=2000):
          b_count += 1
          if b_im is None:
            continue

          logger.info('Generating NRRD for pt %s, studydate %s, series %s:%s, bvalue %s (%i slices)',
                      patient_name, study_date, series_number, series_description, b_val, sliceCount)

          b_fname = filename + '_b' + str(b_val)
          store_name = _store_image(b_im, destination, b_fname, fileformat, overwrite)
          if output_writer is not None and store_name is not None:
            logger.debug('Storing location in CSV output')
            output_writer.writerow([
              counter,
              patient_name,
              study_date,
              store_name.replace(os.path.sep, '/'),
              sliceCount,
              volume.getOrientation()
            ])
        if b_count == 0:
          logger.warning('4D volume contains bvalue tag, but unable to perform correct split.')
        else:
          logger.debug('DWI volume processing complete, found %i b values', b_count)
      elif re.match('\*ep_b\d{1,4}', getattr(volume.slices[0], 'SequenceName', '')):
        # Volume is DWI
        logger.debug('Volume is 4D, attempting splitting on sequence name tag')
        t_count = 0
        for t_val, t_im, sliceCount in volume.getSimpleITK4DImage(splitTag='SequenceName'):
          t_count += 1
          if t_im is None:
            continue
          logger.info('Generating NRRD for pt %s, studydate %s, series %s:%s, temporal value %s (%i slices)',
                      patient_name, study_date, series_number, series_description, t_val, sliceCount)
          t_val = re.sub(r'\*ep_b(\d{1,4})t?', '\1', t_val)

          t_fname = filename + '_' + dicomvolume.DicomVolume.get_safe_filename(t_val)
          store_name = _store_image(t_im, destination, t_fname, fileformat, overwrite)
          if output_writer is not None and store_name is not None:
            logger.debug('Storing location in CSV output')
            output_writer.writerow([
              counter,
              patient_name,
              study_date,
              store_name.replace(os.path.sep, '/'),
              sliceCount,
              volume.getOrientation()
            ])
        if t_count == 0:
          logger.warning('4D volume contains Sequence Name tag, but unable to perform correct split.')
        else:
          logger.debug('4D volume processing complete, found %i temporal positions', t_count)
      else:
        logger.warning("Volume is 4D, but could not split (patient %s, studydate %s series %d. %s), skipping...",
                       patient_name, study_date, series_number, series_description)
        return
    else:
      im = volume.getSimpleITKImage()
      if im is None:
        return

      logger.info('Generating NRRD for pt %s, studydate %s, series %s:%s (%i slices)',
                  patient_name, study_date, series_number, series_description, len(volume.slices))

      store_name = _store_image(im, destination, filename, fileformat, overwrite)
      if output_writer is not None and store_name is not None:
        logger.debug('Storing location in CSV output')
        output_writer.writerow([
          counter,
          patient_name,
          study_date,
          store_name.replace(os.path.sep, '/'),
          len(volume.slices),
          volume.getOrientation()
        ])

    if dump_protocol:
      protocol_fname = os.path.join(destination, filename) + '_protocol.txt'
      volume.writeProtocol(protocol_fname)
    counter += 1

  except KeyboardInterrupt:
    raise
  except:
    logger.error('Oh Oh... something went wrong...', exc_info=True)


def _store_image(im, destination, fname, fileformat, overwrite=False):
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

  return nrrd_fname


def checkVolume(dicomVolume, uid, volume_idx=0):
  global logger
  try:
    if len(dicomVolume.slices) == 0:  # No files for this series UID (maybe not image storage?)
      logger.debug('No files for this series...')
      return

    dicomVolume.sortSlices()
    if dicomVolume.is_equidistant and dicomVolume.is_valid:
      if volume_idx > 0:
        logger.info('DicomVolume %s %s, (volume %d) is valid...',
                    uid, getattr(dicomVolume.slices[0], 'SeriesDescription', 'N/A'), volume_idx + 1)
      else:
        logger.info('DicomVolume %s %s is valid...',
                    uid, getattr(dicomVolume.slices[0], 'SeriesDescription', 'N/A'))
  except KeyboardInterrupt:
    raise
  except:
    logger.error('Oh Oh... something went wrong...', exc_info=True)
