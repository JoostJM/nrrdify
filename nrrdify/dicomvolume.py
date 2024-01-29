#!/usr/bin/env python

# ========================================================================
#  Copyright Het Nederlands Kanker Instituut - Antoni van Leeuwenhoek
#
#  Licensed under the 3-clause BSD License
# ========================================================================

import logging
import struct

import numpy as np
import SimpleITK as sitk
import six

from . import getPostProcessing


class DicomVolume:

  def __init__(self):
    self.slices = []
    self.slices4D = None
    self.logger = logging.getLogger('nrrdify.DicomVolume')

    self._descriptor = None

    self.is_valid = True
    self.is_equidistant = True
    self.is_sorted = False

    self.z_indices = None
    self.n_slices = 0

    self.is_4D = False
    self.is_sorted4D = False
    self.split_tag = None
    self.n_pos = 1

    self.post_processing = getPostProcessing()

  def __getitem__(self, item):
    return self.slices[item]

  def addSlice(self, dicFile):
    self.slices += [dicFile]
    self.is_sorted = False
    self.is_sorted4D = False

  def _check_valid(self):
    patientName = getattr(self.slices[0], 'PatientName', None)
    studyDate = getattr(self.slices[0], 'StudyDate', None)

    series_description = getattr(self.slices[0], 'SeriesDescription', 'Unkn')
    series_number = getattr(self.slices[0], 'SeriesNumber', -1)

    if patientName is None or studyDate is None:
      self.logger.error('Missing patient name in file %s', self.slices[0].filename)
      return False
    if studyDate is None:
      self.logger.error('Missing study date in file %s', self.slices[0].filename)
      return False

    required_tags = ['ImagePositionPatient', 'PatientName', 'StudyDate']  # Check pt name and study date in other slices
    identical_tags = ['ImageOrientationPatient']

    for tag in required_tags:
      for dfile in self.slices:
        if getattr(dfile, tag, None) is None:
          self.logger.error('No value found for tag %s in file %s (patient %s, studydate %s, series %d. %s), invalid series!',
                            tag, dfile.filename, patientName, studyDate, series_number, series_description)
          return False
    for tag in identical_tags:
      val = getattr(self.slices[0], tag, None)
      if val is None:
        self.logger.error('No value found for tag %s in file %s (patient %s, studydate %s, series %d. %s), invalid series!',
                          tag, self.slices[0].filename, patientName, studyDate, series_number, series_description)
        return False

      for dfile in self.slices[1:]:
        val2 = getattr(dfile, tag, None)
        if val2 is None:
          self.logger.error('No value found for tag %s in file %s (patient %s, studydate %s, series %d. %s), invalid series!',
                            tag, dfile.filename, patientName, studyDate, series_number, series_description)
          return False
        if not np.allclose(val, val2, rtol=1e-2, atol=1e-6):
          self.logger.error('Non-matching values found for tag %s between files %s and %s (patient %s, studydate %s, series %d. %s)',
                            tag, self.slices[0].filename, dfile.filename, patientName, studyDate, series_number, series_description)
          return False
    return True

  def _get_slice_locations(self, slices):
    self.logger.debug('Calculation slice positions...')
    image_orientation = slices[0].ImageOrientationPatient
    xvector = image_orientation[:3]
    yvector = image_orientation[3:]
    # This function assumes that the Z axis of the image stack is orthogonal to the Y and X axis
    zvector = np.cross(xvector, yvector)
    locations = [np.dot(dfile.ImagePositionPatient, zvector) for dfile in slices]  # Z locations in mm of each slice
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

    if self.post_processing is not None:
      im = self.post_processing(im, slices)

    return im

  @property
  def descriptor(self):
    if self._descriptor is None:
      self._descriptor = Descriptor(self.slices[0])
    return self._descriptor

  def build_filename(self):
    filename = '%s-%s-%s. %s' % (
      self.descriptor.patient_name,
      self.descriptor.study_date,
      self.descriptor.series_number,
      self.descriptor.series_description
    )
    return self.normalize_filename(filename)

  @staticmethod
  def normalize_filename(filename):
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

    locations = self._get_slice_locations(self.slices)

    # Check if all slices are equidistant
    delta_slices = np.diff(np.array(sorted(locations)))

    boundaries = delta_slices[delta_slices > 0.03]  # Exclude slice-distances that are 0 (indicating 4D volume)
    if len(boundaries) > 0 and not np.allclose(boundaries, boundaries[0], rtol=5e-2):
      self.logger.warning('Slices are not equidistant!')
      self.logger.debug('Slice distances:\n%s', boundaries)
      self.is_equidistant = False
      return

    self.n_slices = len(boundaries) + 1

    if len(delta_slices) > len(boundaries):
      self.logger.debug('DicomVolume is 4D')
      self.is_4D = True
      # Check if the total slice count is equally divisible by n_slices. If not, some volumes in this 4D series
      # may be missing some slices...
      if len(self.slices) % self.n_slices > 0:
        self.logger.warning('Total slice count not equally divisible by n_slices, missing 4D slices?')
        self.is_equidistant = False
        return

    else:
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
    self.split_tag = splitTag
    for s in self.slices:
      if isinstance(splitTag, int):
        temporal_position = s.get(splitTag, None)
        if temporal_position is not None:
          temporal_position = temporal_position.value
      else:
        temporal_position = getattr(s, splitTag, None)
      if temporal_position is None:
        self.logger.error('Split Tag %s not valid, missing value in file %s', splitTag, s.filename)
        return False

      is_string_temporal = False
      if isinstance(temporal_position, float):
        if temporal_position.is_integer():
          temporal_position = int(temporal_position)
      elif isinstance(temporal_position, six.string_types):
        if temporal_position.strip().isdigit():
          temporal_position = int(temporal_position)
        else:
          is_string_temporal = True
      elif not isinstance(temporal_position, int):
        try:
          temporal_position = float(struct.unpack('d', temporal_position)[0])
          if temporal_position.is_integer():
            temporal_position = int(temporal_position)
        except KeyboardInterrupt:
          raise
        except Exception:
          self.logger.error('Error unpacking value for tag %s in file %s', splitTag, s.filename)
          return False

      if not is_string_temporal and max_value is not None and temporal_position > max_value:
        self.logger.warning('File %s excluded (temporal position (%d) exceeded max value %d)', s.filename, temporal_position, max_value)

      if temporal_position not in self.slices4D:
        self.slices4D[temporal_position] = []

      self.slices4D[temporal_position].append(s)
    return True

  def sortSlices4D(self):
    if self.is_sorted4D:
      return True  # Slices already sorted, no re-sort needed

    if self.slices4D is None:
      self.logger.warning("DicomVolume needs to be split by temporal position before 4D sorting can occur")
      return False

    slice_count = None
    for t in self.slices4D:
      if slice_count is None:
        slice_count = len(self.slices4D[t])
      elif len(self.slices4D[t]) != slice_count:
        self.logger.warning('Different number of slices between temporal positions! Split tag "%s" not valid')
        return False

      locations = self._get_slice_locations(self.slices4D[t])

      delta_slices = np.diff(np.array(sorted(locations)))
      if np.any(delta_slices < 0.03):
        self.logger.warning('Found overlapping slicer in timepoint %s! Split tag "%s" not valid.', t, self.split_tag)
        return False

      if not np.allclose(delta_slices, delta_slices[0], rtol=1e-2):
        self.logger.warning('Slices are not equidistant!')
        self.logger.debug('Slice distances:\n%s', delta_slices)
        self.is_equidistant = False
        return
      self.slices4D[t] = [f for (d, f) in sorted(zip(locations, self.slices4D[t]), key=lambda s: s[0])]

    self.n_pos = len(self.slices4D.keys())

    self.is_sorted4D = True
    return True

  def getSimpleITK4DImage(self, splitTag=None, max_value=None):
    if not self.is_sorted:
      self.sortSlices()

    if not (self.is_valid and self.is_equidistant):
      return

    if not self.is_4D:
      self.logger.error('Cannot generate 4D image, DicomVolume is 3D!')
      return

    if self.slices4D is None or (splitTag is not None and self.split_tag != splitTag):
      if splitTag is None:
        splitTag = 'DiffusionBValue'
      if not self.split4D(splitTag, max_value):
        return

    if not self.is_sorted4D:
      if not self.sortSlices4D():
        return

    for t in self.slices4D:
      self.logger.debug('Processing temporal position %d', t)
      yield t, self._getImage(self.slices4D[t]), len(self.slices4D[t])

  def _compute_z_indices(self, ):
    if self.z_indices is not None:
      return

    if not self.is_sorted:
      self.sortSlices()

    if self.is_4D:
      if not self.is_sorted4D:
        raise ValueError("Volume is 4D, first ensure volume is 4D sorted before continuing!")
      self.z_indices = {}
      for t in self.slices4D:
        for idx, s in enumerate(self.slices4D[t]):
          self.z_indices[s.get('SOPInstanceUID')] = (idx, t)
    else:
      self.z_indices = {s.get('SOPInstanceUID'): (idx,) for idx, s in enumerate(self.slices)}

  def get_z_index(self, sopinstanceuid):
    if self.z_indices is None:
        self._compute_z_indices()

    return self.z_indices[sopinstanceuid][0]

  def get_temporal_position(self, sopinstanceuid):
    if not self.is_4D:
      return 0

    if self.z_indices is None:
      self._compute_z_indices()

    return self.z_indices[sopinstanceuid][1]

  def getOrientation(self):
    assert len(self.slices) > 0, "No slices in this volume!"
    assert hasattr(self.slices[0], 'ImageOrientationPatient'),\
        'Invalid volume! does not contain ImageOrientationPatient tag'

    # First, get cosines and make them absolute
    image_orientation = np.abs(self.slices[0].ImageOrientationPatient)
    # Xx = 0
    # Yx = 1
    # Zx = 2
    # Xy = 3
    # Yy = 4
    # Zy = 5
    # Xz = 6
    # Yz = 7
    # Zz = 8

    if image_orientation[3] >= image_orientation[0]:  # |Xy| >= |Xx| (saggital)
      return 'sag'
    elif image_orientation[4] >= image_orientation[7]:  # |Yy| >= |Yz| (Transversal)
      return 'tra'
    else:  # Coronal
      return 'cor'

  def writeProtocol(self, proctocol_fname):
    if not self.is_sorted:
      self.sortSlices()

    divider = """
    \n
    \n
    -----------------------------------------------------------------------------
    \n
    \n
    """

    with open(proctocol_fname, mode='w') as protocol_fs:

      for dicfile in self.slices:
        if 0x00280013 in dicfile:
          del dicfile[0x00280013]
        protocol_fs.write(str(dicfile) + divider)


class Descriptor:
  def __init__(self, dicfile):
    self.patient_name = str(getattr(dicfile, 'PatientName', '')).split('^')[0]
    self.study_date = getattr(dicfile, 'StudyDate', '19000101')
    self.series_description = getattr(dicfile, 'SeriesDescription', 'Unkn')
    self.series_number = getattr(dicfile, 'SeriesNumber', -1)
