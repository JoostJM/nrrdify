#!/usr/bin/env python

# =========================================================================
#  Copyright The Netherlands Cancer Institute - Antoni van Leeuwenhoek
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0.txt
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# ========================================================================

# Authors: Stefano Trebeschi, Joost van Griethuysen

# Python Requirements:
# SimpleITK
# PyDicom
# tqdm
# Numpy

import logging
import os
import argparse
import struct
import traceback

import dicom
import numpy as np
import SimpleITK as sitk
import tqdm

logger = logging.getLogger('nrrdify')

if len(logger.handlers) == 0:
  print('Adding handler for logger')
  handler = logging.StreamHandler()
  formatter = logging.Formatter('[%(asctime)-.19s] %(levelname)-.1s: %(message)s')
  handler.setFormatter(formatter)
  handler.setLevel(logging.INFO)

  logger.addHandler(handler)
  logger.setLevel(logging.INFO)


class DicomVolume:

  def __init__(self):
    self.slices = []
    self.slices4D = None
    self.logger = logging.getLogger('nrrdify.DicomVolume')

    self.is_valid = True
    self.is_equidistant = True
    self.is_sorted = False

    self.is_4D = False
    self.is_sorted4D = False

  def __getitem__(self, item):
    return self.slices[item]

  def addSlice(self, dicFile):
    self.slices += [dicFile]
    self.is_sorted = False
    self.is_sorted4D = False

  def _check_valid(self):
    required_tags = ['ImagePositionPatient']
    identical_tags = ['ImageOrientationPatient']

    for tag in required_tags:
      for dfile in self.slices:
        if getattr(dfile, tag, None) is None:
          self.logger.error('No value found for tag %s in file %s, invalid series!', tag, dfile.filename)
          return False
    for tag in identical_tags:
      val = getattr(self.slices[0], tag, None)
      if val is None:
        self.logger.error('No value found for tag %s in file %s, invalid series!', tag, self.slices[0].filename)
        return False
      for dfile in self.slices[1:]:
        val2 = getattr(dfile, tag, None)
        if val2 is None:
          self.logger.error('No value found for tag %s in file %s, invalid series!', tag, dfile.filename)
          return False
        if not np.allclose(val, val2, rtol=1e-2):
          self.logger.error('Non-matching values found for tag %s between files %s and %s',
                            tag, self.slices[0].filename, dfile.filename)
          return False
    return True

  def _get_slice_locations(self):
    self.logger.debug('Calculation slice positions...')
    image_orientation = self.slices[0].ImageOrientationPatient
    xvector = image_orientation[:3]
    yvector = image_orientation[3:]
    zvector = np.cross(xvector, yvector)  # This function assumes that the Z axis of the image stack is orthogonal to the Y and X axis

    locations = [np.dot(dfile.ImagePositionPatient, zvector) for dfile in self.slices]  # Z locations in mm of each slice
    return locations

  def _getImage(self, slices):
    if len(slices) == 1:  # e.g. enhanced image format file, or 2D image
      self.logger.debug('Single File, attempting SimpleITK.ReadImage')
      im = sitk.ReadImage(slices[0].filename)
    else:  # 'classic' DICOM file (1 file / slice)
      reader = sitk.ImageSeriesReader()
      self.logger.debug('Setting filenames for image reader...')
      reader.SetFileNames([f.filename for f in slices])
      self.logger.debug('Getting the image (%d files)...', len(slices))
      im = reader.Execute()
    return im

  def build_filename(self):
    patient_name = str(getattr(self.slices[0], 'PatientName', '')).split('^')[0]
    study_date = getattr(self.slices[0], 'StudyDate', '19000101')
    series_description = getattr(self.slices[0], 'SeriesDescription', 'Unkn')
    series_number = getattr(self.slices[0], 'SeriesNumber', -1)

    filename = '%s-%s-%s. %s' % (patient_name, study_date, series_number, series_description)
    # Remove invalid characters from filename
    for c in r'[]/\;,><&*:%=+@!#^()|?^':
      filename = filename.replace(c, '')

    return filename

  def sortSlices(self):
    if self.is_sorted:
      return

    if len(self.slices) < 2:
      self.is_sorted = True
      self.is_valid = True
      self.is_equidistant = True
      return

    if not self._check_valid():
      self.is_valid = False
      return

    locations = self._get_slice_locations()

    # Check if all slices are equidistant
    delta_slices = np.diff(np.array(sorted(locations)))

    boundaries = delta_slices[delta_slices > 0.03]  # Exclude slice-distances that are 0 (indicating 4D volume)
    if not np.allclose(boundaries, boundaries[0], rtol=1e-2):
      self.logger.warning('Slices are not equidistant!')
      self.logger.debug('Slice distances:\n%s', boundaries)
      self.is_equidistant = False
      return

    if len(delta_slices) > len(boundaries):
      self.logger.debug('DicomVolume is 4D')
      self.is_4D = True

    self.slices = [f for (d, f) in sorted(zip(locations, self.slices), key=lambda s: s[0])]
    self.is_sorted = True

  def getSimpleITKImage(self):
    if not self.is_sorted:
      self.sortSlices()

    if not (self.is_valid and self.is_equidistant):
      return

    if self.is_4D:
      self.logger.error('Cannot generate 3D image, DicomVolume is 4D!')
      return

    im = self._getImage(self.slices)
    return im

  def check_4D(self):
    if not self.is_sorted:
      self.sortSlices()

    return self.is_4D

  def split4D(self, splitTag='DiffusionBValue', max_value=None):
    self.slices4D = {}
    for s in self.slices:
      temporal_position = getattr(s, splitTag, None)
      if temporal_position is None:
        self.logger.error('Split Tag %s not valid, missing value in file %s', splitTag, s.filename)
        return False

      if not str(temporal_position).isdigit():  # not a valid B value
        try:
          temporal_position = struct.unpack('d', temporal_position)[0]
        except Exception:
          self.logger.error('Error unpacking value for tag %s in file %s', splitTag, s.filename)
          return False

      if max_value is not None and temporal_position > max_value:
        self.logger.warning('File %s excluded (temporal position (%d) exceeded max value %d)', s.filename, temporal_position, max_value)

      if temporal_position not in self.slices4D:
        self.slices4D[temporal_position] = []

      self.slices4D[temporal_position].append(s)
    return True

  def sortSlices4D(self):
    if self.is_sorted4D:
      return True  # Slices already sorted, no re-sort needed

    if self.slices4D is None:
      logger.warning("DicomVolume needs to be split by temporal position before 4D sorting can occur")
      return False

    slice_count = None
    for t in self.slices4D:
      if slice_count is None:
        slice_count = len(self.slices4D[t])
      elif len(self.slices4D[t]) != slice_count:
        self.logger.error('Different number of slices between temporal positions!')
        return False

      locations = self._get_slice_locations()
      self.slices4D[t] = [f for (d, f) in sorted(zip(locations, self.slices4D[t]), key=lambda s: s[0])]

    return True

  def getSimpleITK4DImage(self, splitTag='DiffusionBValue', max_value=None):
    if not self.is_sorted:
      self.sortSlices()

    if not (self.is_valid and self.is_equidistant):
      return

    if not self.is_4D:
      self.logger.error('Cannot generate 4D image, DicomVolume is 3D!')
      return

    if self.slices4D is None:
      if not self.split4D(splitTag, max_value):
        return

    if not self.is_sorted4D:
      if not self.sortSlices4D():
        return

    ims = {}
    for t in self.slices4D:
      self.logger.debug('Processing temporal position %d', t)
      ims[t] = self._getImage(self.slices4D[t])

    return ims


def main(source, destination, filename=None, fileformat='nrrd', overwrite=False, just_check=False):
  global logger
  if os.path.isdir(source) and os.path.isdir(destination):
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
              dicfile = dicom.read_file(os.path.join(curdir, fname), stop_before_pixels=True)

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
                datasets[series_uid][imagetype] = DicomVolume()

              datasets[series_uid][imagetype].addSlice(dicfile)
            except:
              logger.error('DOH!! Something went wrong \n\n%s' % traceback.format_exc())
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
          processVolume(volume, destination, filename, fileformat, overwrite, volume_idx)


def processVolume(dicomVolume, destination, filename=None, fileformat='nrrd', overwrite=False, file_idx=None):
  global logger
  try:
    if len(dicomVolume.slices) == 0:  # No files for this series UID (maybe not image storage?)
      logger.debug('No files for this series...')
      return

    if dicomVolume.check_4D():
      logger.warning("Volume is 4D, skipping...")
      return

    patient_name = str(getattr(dicomVolume[0], 'PatientName', '')).split('^')[0]
    study_date = getattr(dicomVolume[0], 'StudyDate', '19000101')
    series_description = getattr(dicomVolume[0], 'SeriesDescription', 'Unkn')
    series_number = getattr(dicomVolume[0], 'SeriesNumber', -1)

    logger.info('Generating NRRD for pt %s, studydate %s, series %s:%s' %
                (patient_name, study_date, series_number, series_description))

    if filename is None:  # Generate a filename from DICOM metadata
      filename = dicomVolume.build_filename()

    if file_idx is not None and file_idx > 0:
      filename = '%s (%d)' % (filename, file_idx)

    filename = os.path.join(destination, filename)
    filename += '.' + fileformat

    if os.path.isfile(filename):
      if overwrite:
        logger.warning('file "%s" already exists, overwriting...', filename)
      else:
        logger.info('file "%s" already exists, skipping...', filename)
        return

    im = dicomVolume.getSimpleITKImage()
    logger.info('Image file series read (%d files), storing in %s', len(dicomVolume.dicFiles()), filename)

    sitk.WriteImage(im, filename)
  except:
    logger.error('Oh Oh... something went wrong...', exc_info=True)


def checkVolume(dicomVolume, uid, volume_idx=None):
  try:
    if len(dicomVolume.dicFiles()) == 0:  # No files for this series UID (maybe not image storage?)
      logger.debug('No files for this series...')
      return

    dicomVolume.sortSlices()
    if dicomVolume.is_equidistant and dicomVolume.is_valid:
      if volume_idx is not None:
        logger.info('DicomVolume %s, (volume %d) is valid...', uid, volume_idx + 1)
      else:
        logger.info('DicomVolume %s is valid...', uid)
  except:
    logger.error('Oh Oh... something went wrong...', exc_info=True)


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("inputFolder", metavar="In", type=str, help="Folder containing the DICOM file(s) to convert.")
  parser.add_argument("--out", "-o", help="Folder to store converted files in. If omitted, stores "
                                          "files in parent directory of In folder.")
  parser.add_argument("--name", "-n", help="Filename for the new file, without extension. If omitted, or more series "
                                           "are found, Filename is generated from DICOM tags: "
                                           "<PatientName>-<StudyDate>-<SeriesNumber>. <SeriesDescription>")
  parser.add_argument("--format", "-f", nargs="?", default="nrrd", choices=["nrrd", "nii", "nii.gz"],
                      help="Image format to convert to. Default is the 'nrrd' format")
  parser.add_argument('--logging-level', metavar='LEVEL',
                      choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                      default='WARNING', help='Set capture level for logging')
  parser.add_argument('--log-file', metavar='FILE', type=argparse.FileType('w'), default=None,
                      help='File to append logger output to')
  parser.add_argument('--overwrite', action='store_true', help='if this argument is specified, script will overwrite existing files, '
                                                               'otherwise, file write for already existing files is skipped.')
  parser.add_argument('--check', action='store_true', help='if this argument is specified, DICOMS are checked but not converted.')

  args = parser.parse_args()

  # if specified, set up logging to a file
  if args.log_file is not None:
    log_handler = logging.StreamHandler(args.log_file)
    log_formatter = logging.Formatter('[%(asctime)-.19s] %(levelname)-.1s: %(message)s')
    log_handler.setFormatter(log_formatter)
    logLevel = getattr(logging, args.logging_level)
    log_handler.setLevel(logLevel)

    logger.addHandler(log_handler)
    if logLevel < logging.INFO:  # Lower logger level if necessary
      logger.setLevel(logLevel)

  source_folder = args.inputFolder
  destination_folder = args.out
  if destination_folder is None:
    destination_folder = os.path.dirname(source_folder)

  main(source_folder, destination_folder, args.name, args.format, args.overwrite, just_check=args.check)
