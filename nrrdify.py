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


class DicomVolume:

  def __init__(self):
    self._slices = []
    self.logger = logging.getLogger('nrrdify.DicomVolume')

    self.is_valid = True
    self.is_equidistant = True
    self.is_sorted = False

  def __getitem__(self, item):
    return self._slices[item]

  def addSlice(self, dicFile):
    self._slices += [dicFile]
    self.is_sorted = False

  def _check_valid(self):
    required_tags = ['ImagePositionPatient']
    identical_tags = ['ImageOrientationPatient']

    for tag in required_tags:
      for dfile in self._slices:
        if getattr(dfile, tag, None) is None:
          self.logger.error('No value found for tag %s in file %s, invalid series!', tag, dfile.filename)
          return False
    for tag in identical_tags:
      val = getattr(self._slices[0], tag, None)
      if val is None:
        self.logger.error('No value found for tag %s in file %s, invalid series!', tag, dicom_files[0].filename)
        return False
      for dfile in self._slices[1:]:
        val2 = getattr(dfile, tag, None)
        if val2 is None:
          self.logger.error('No value found for tag %s in file %s, invalid series!', tag, dfile.filename)
          return False
        if not np.allclose(val, val2, rtol=1e-2):
          self.logger.error('Non-matching values found for tag %s between files %s and %s',
                       tag, self._slices[0].filename, dfile.filename)
          return False
    return True

  def _get_slice_locations(self):
    self.logger.debug('Calculation slice positions...')
    image_orientation = self._slices[0].ImageOrientationPatient
    xvector = image_orientation[:3]
    yvector = image_orientation[3:]
    zvector = np.cross(xvector, yvector)  # This function assumes that the Z axis of the image stack is orthogonal to the Y and X axis

    locations = [np.dot(dfile.ImagePositionPatient, zvector) for dfile in self._slices]  # Z locations in mm of each slice
    return locations

  def build_filename(self):
    patient_name = getattr(self._slices[0], 'PatientName', '').split('^')[0]
    study_date = getattr(self._slices[0], 'StudyDate', '19000101')
    series_description = getattr(self._slices[0], 'SeriesDescription', 'Unkn')
    series_number = getattr(self._slices[0], 'SeriesNumber', -1)

    filename = '%s-%s-%s. %s' % (patient_name, study_date, series_number, series_description)
    # Remove invalid characters from filename
    for c in r'[]/\;,><&*:%=+@!#^()|?^':
      filename = filename.replace(c, '')

    return filename

  def dicFiles(self):
    return self._slices

  def sortSlices(self):
    if len(self._slices) < 2:
      return

    if not self._check_valid():
      self.is_valid = False
      return

    locations = self._get_slice_locations()

    # Check if all slices are equidistant
    delta_slices = np.diff(np.array(sorted(locations)))
    if not np.allclose(delta_slices, delta_slices[0], rtol=1e-2):
      logger.warning('Slices are not equidistant!')
      logger.debug('Slice distances:\n%s', delta_slices)
      self.is_equidistant = False
      return

    self._slices = [f for (d, f) in sorted(zip(locations, self._slices), key=lambda s: s[0])]
    self.is_sorted = True

  def getSimpleITKImage(self):
    if len(self._slices) == 1:  # e.g. enhanced image format file, or 2D image
      self.logger.debug('Single File, attempting SimpleITK.ReadImage')
      im = sitk.ReadImage(self._slices[0].filename)
      self.logger.info('Single File read, storing in %s', self._slices[0].filename)
    else:  # 'classic' DICOM file (1 file / slice)
      if not self.is_sorted:
          self.sortSlices()

      if not (self.is_valid and self.is_equidistant):
        return

      reader = sitk.ImageSeriesReader()
      self.logger.debug('Setting filenames for image reader...')
      reader.SetFileNames([f.filename for f in self._slices])
      self.logger.debug('Getting the image (%d files)...', len(self._slices))
      im = reader.Execute()
    return im


def main(source_folder, destination_folder, filename=None, fileformat='nrrd', overwrite=False, max_bval=1500, just_check=False):
  global logger
  if os.path.isdir(source_folder) and os.path.isdir(destination_folder):
    logger.info('Input and output valid, scanning input folder for DICOM files')
    datasets = {}  # Holds the dicom files, sorted by series UID ({seriesUID: [files]})
    for curdir, dirnames, fnames in os.walk(source_folder):
      if len(fnames) > 0:  # Only process folder if it contains files
        logger.info('Processing folder %s', curdir)

        with tqdm.tqdm(fnames, desc='Processing files') as bar:  # Progress reporting
          for fname in bar:  # for each file in current folder
            try:
              # Check if it contains a valid DICOM header (first 4 bytes = DICM)
              with open(os.path.join(curdir, fname), mode='rb') as openFile:
                openFile.seek(128)
                header = openFile.read(4)
                if header != 'DICM':
                  # Not a valid DICOM file, skip to next
                  continue  # Go to next file

              # Load dicom file using PyDicom (needed for name extraction, sorting of series and slices)
              dicfile = dicom.read_file(os.path.join(curdir, fname), stop_before_pixels=True)
              sop_class = getattr(dicfile, 'SOPClassUID', None)  # Check if it is a dicomfile containing an image
              series_uid = getattr(dicfile, 'SeriesInstanceUID', None)  # Series UID
              if series_uid is None:
                continue  # Error cannot sort, so skip and go To next file
              if sop_class is None or 'Image Storage' not in str(sop_class):
                continue  # not image dicom file, so skip and go to next file

              #b_value = getattr(dicfile, 'DiffusionBValue', None)
              #if b_value is None:
              #  continue
              #if not str(b_value).isdigit():  # not a valid B value
              #  try:
              #    b_value = struct.unpack('d', b_value)
              #  except Exception:
              #    continue

              #if b_value > max_bval:  # Do not use b values larger than max_bval (exclude e.g. b2000)
              #  continue

              if series_uid not in datasets:
                datasets[series_uid] = DicomVolume()

              datasets[series_uid].addSlice(dicfile)
            except Exception as e:
              logger.error('DOH!! Something went wrong \n\n%s' % traceback.format_exc())
    if just_check:
      for ds in datasets:
        checkVolume(datasets[ds], ds)
    else:
      # Done scanning files, now make some NRRDs out of them!
      logger.info('Input folder scanned, found %d unique DICOM series', len(datasets))
      if len(datasets) == 1:  # If only 1 series is found, a custom filename is possible
        processVolume(datasets[datasets.keys()[0]], destination_folder, filename, fileformat, overwrite)
      else:
        for ds in datasets:  # Multiple datasets, so generate name from DICOM
          processVolume(datasets[ds], destination_folder, fileformat=fileformat, overwrite=overwrite)


def processVolume(dicomVolume, destination_folder, filename=None, fileformat='nrrd', overwrite=False):
  global logger
  try:
    if len(dicomVolume.dicFiles()) == 0:  # No files for this series UID (maybe not image storage?)
      logger.debug('No files for this series...')
      return

    patient_name = getattr(dicomVolume[0], 'PatientName', '').split('^')[0]
    study_date = getattr(dicomVolume[0], 'StudyDate', '19000101')
    series_description = getattr(dicomVolume[0], 'SeriesDescription', 'Unkn')
    series_number = getattr(dicomVolume[0], 'SeriesNumber', -1)

    logger.info('Generating NRRD for pt %s, studydate %s, series %s:%s' %
                (patient_name, study_date, series_number, series_description))

    if filename is None:  # Generate a filename from DICOM metadata
      filename = dicomVolume.build_filename()
    filename = os.path.join(destination_folder, filename)
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

def checkVolume(dicomVolume, uid):
  try:
    if len(dicomVolume.dicFiles()) == 0:  # No files for this series UID (maybe not image storage?)
      logger.debug('No files for this series...')
      return

    dicomVolume.sortSlices()
    if dicomVolume.is_equidistant and dicomVolume.is_valid:
    	logger.info('DicomVolume %s is valid...', uid)
  except:
    logger.error('Oh Oh... something went wrong...', exc_info=True)


if __name__ == '__main__':
  if len(logger.handlers) == 0:
    print('Adding handler for logger')
    handler = logging.StreamHandler()
    formatter = logging.Formatter('[%(asctime)-.19s] %(levelname)-.1s: %(message)s')
    handler.setFormatter(formatter)
    handler.setLevel(logging.INFO)

    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

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
