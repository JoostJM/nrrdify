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
                mkdirs=True,
                output_writer=None):
  global counter, logger, post_processing
  if not os.path.isdir(source):
    logger.error('Source directory (%s) does not exist! Exiting...', source)
    return
  if not os.path.isdir(destination):
    logger.error('Destination directory (%s) does not exist! Exiting...', destination)
    return
  counter = 0
  logger.info('Input and output valid, scanning input folder for DICOM files')
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
            if imagetype is None:
              logger.debug("missing Image Type tag in dicom file %s", os.path.join(curdir, fname))
              continue  # Error cannot sort, so skip and go To next file
            if series_uid is None:
              continue  # Error cannot sort, so skip and go To next file
            if sop_class is None or 'Image Storage' not in str(sop_class):
              continue  # not image dicom file, so skip and go to next file

            imagetype = tuple(imagetype)

            if series_uid not in datasets:
              datasets[series_uid] = {}

            if imagetype not in datasets[series_uid]:
              datasets[series_uid][imagetype] = dicomvolume.DicomVolume(post_processing)

            datasets[series_uid][imagetype].addSlice(dicfile)
          except:
            logger.error('DOH!! Something went wrong!', exc_info=True)
      if process_per_folder:
        dest = os.path.join(destination, curdir)
        if not os.path.isdir(dest):
          logger.debug('Creating output directory "%s"', dest)
          os.makedirs(dest)

        _processResults(datasets, dest, filename, fileformat, overwrite, just_check, mkdirs, output_writer)
        datasets = {}

  if not process_per_folder:
    _processResults(datasets, destination, filename, fileformat, overwrite, just_check, mkdirs, output_writer)


def _processResults(datasets, destination, filename, fileformat, overwrite, just_check, mkdirs, output_writer):
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
          processVolume(volume, destination, filename, fileformat, overwrite, volume_idx, mkdirs, output_writer)


def processVolume(volume,
                  destination,
                  filename=None,
                  fileformat='nrrd',
                  overwrite=False,
                  file_idx=None,
                  mkdirs=False,
                  output_writer=None):
  global counter, logger
  try:
    if len(volume.slices) == 0:  # No files for this series UID (maybe not image storage?)
      logger.debug('No files for this series...')
      return

    patient_name = str(getattr(volume[0], 'PatientName', '')).split('^')[0]
    study_date = getattr(volume[0], 'StudyDate', '19000101')
    series_description = getattr(volume[0], 'SeriesDescription', 'Unkn')
    series_number = getattr(volume[0], 'SeriesNumber', -1)

    if volume.check_4D():
      logger.warning("Volume is 4D (patient %s, studydate %s series %d. %s), skipping...",
                     patient_name, study_date, series_number, series_description)
      return

    im = volume.getSimpleITKImage()
    if im is None:
      return

    logger.info('Generating NRRD for pt %s, studydate %s, series %s:%s' %
                (patient_name, study_date, series_number, series_description))

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

    filename = os.path.join(destination, filename)
    filename += '.' + fileformat

    if destination != '' and not os.path.isdir(destination):
      logger.debug('Creating study directory "%s"', destination)
      os.makedirs(destination)
    elif os.path.isfile(filename):
      if overwrite:
        logger.warning('file "%s" already exists, overwriting...', filename)
      else:
        logger.info('file "%s" already exists, skipping...', filename)
        return

    logger.info('Image file series read (%d files), storing in %s', len(volume.slices), filename)

    sitk.WriteImage(im, filename)

    counter += 1

    if output_writer is not None:
      logger.debug('Storing location in CSV output')
      output_writer.writerow([counter, patient_name, study_date, filename, len(volume.slices)])
  except:
    logger.error('Oh Oh... something went wrong...', exc_info=True)


def checkVolume(dicomVolume, uid, volume_idx=0):
  global logger
  try:
    if len(dicomVolume.dicFiles()) == 0:  # No files for this series UID (maybe not image storage?)
      logger.debug('No files for this series...')
      return

    dicomVolume.sortSlices()
    if dicomVolume.is_equidistant and dicomVolume.is_valid:
      if volume_idx > 0:
        logger.info('DicomVolume %s, (volume %d) is valid...', uid, volume_idx + 1)
      else:
        logger.info('DicomVolume %s is valid...', uid)
  except:
    logger.error('Oh Oh... something went wrong...', exc_info=True)
