#!/usr/bin/env python

# ========================================================================
#  Copyright Het Nederlands Kanker Instituut - Antoni van Leeuwenhoek
#
#  Licensed under the 3-clause BSD License
# ========================================================================

import logging
import os
import sys

import pydicom
import SimpleITK as sitk
import tqdm

from . import dicomvolume


class Walker:
  def __init__(self, **kwargs):
    self.logger = logging.getLogger('nrrdify.walker')
    self.config = kwargs
    self.counter = 0

  def run(self, source, destination, **kwargs):
    process_per_folder = self.config.get('process_per_folder', False)
    structure = self.config.get('structure', None)

    if not os.path.isdir(source):
      self.logger.error('Source directory (%s) does not exist! Exiting...', source)
      return
    if not os.path.isdir(destination):
      self.logger.error('Destination directory (%s) does not exist! Exiting...', destination)
      return
    self.logger.info('Input (%s) and output (%s) valid, scanning input folder for DICOM files', source, destination)
    datasets = {}  # Holds the dicom files, sorted by series UID ({seriesUID: [files]})
    for curdir, dirnames, fnames in os.walk(source):
      if len(fnames) > 0:  # Only process folder if it contains files
        self.logger.info('Processing folder %s', curdir)

        it = fnames.__iter__()
        if self.config.get('use_pbar', True):  # Progress reporting
          it = tqdm.tqdm(fnames, desc='Processing files', file=self.config.get('pbar_out', sys.stdout))
        try:
          for fname in it:  # for each file in current folder
            try:
              # Check if it contains a valid DICOM header (first 4 bytes = DICM)
              with open(os.path.join(curdir, fname), mode='rb') as openFile:
                openFile.seek(128)
                header = openFile.read(4)
                if header != b'DICM':
                  # Not a valid DICOM file, skip to next
                  self.logger.debug('File %s has no valid DICOM header', os.path.join(curdir, fname))
                  continue  # Go to next file

              # Load dicom file using PyDicom (needed for name extraction, sorting of series and slices)
              dicfile = pydicom.read_file(os.path.join(curdir, fname), stop_before_pixels=True)

              imagetype = getattr(dicfile, 'ImageType', None)
              sop_class = getattr(dicfile, 'SOPClassUID', None)  # Check if it is a dicomfile containing an image
              series_uid = getattr(dicfile, 'SeriesInstanceUID', None)  # Series UID
              modality = getattr(dicfile, 'Modality', '')

              if series_uid is None:
                self.logger.debug('Dicom file %s has no series UID', os.path.join(curdir, fname))
                continue  # Error cannot sort, so skip and go To next file
              if sop_class is None:
                self.logger.debug('Dicom file %s has no SOP Class', os.path.join(curdir, fname))
                continue  # not image dicom file, so skip and go to next file

              if modality == 'RTDOSE':
                self.logger.debug('Dicom file %s is RT dose storage', os.path.join(curdir, fname))
              elif 'Image Storage' not in sop_class.name:  # Not RT Dose, is the sopclass valid?
                self.logger.debug('Dicom file %s has no valid SOP Class', os.path.join(curdir, fname))
                continue
              elif imagetype is None:  # Sop class is valid, does it contain image type?
                self.logger.debug("missing Image Type tag in dicom file %s", os.path.join(curdir, fname))
                imagetype = 'UNK'

              if imagetype is None:
                imagetype = 'RTDOSE'
              else:
                imagetype = tuple(imagetype)

              if series_uid not in datasets:
                datasets[series_uid] = {}

              if imagetype not in datasets[series_uid]:
                datasets[series_uid][imagetype] = dicomvolume.DicomVolume()
              else:
                datasets[series_uid][imagetype].addSlice(dicfile)
            except KeyboardInterrupt:
              return
            except:
              self.logger.error('DOH!! Something went wrong!', exc_info=True)
        finally:
          if self.config.get('use_pbar', True):
            it.close()
        if process_per_folder:
          if structure == 'source':
            dest = os.path.join(destination, os.path.relpath(curdir, source))
          else:
            dest = destination

          self._process_results(datasets, dest, **kwargs)
          datasets = {}

    if not process_per_folder:
      self._process_results(datasets, destination, **kwargs)

  def _process_results(self, datasets, destination, **kwargs):
    # Done scanning files, now make some NRRDs out of them!
    self.logger.info('Input folder scanned, found %d unique DICOM series', len(datasets))
    if len(datasets) > 1:  # If more than 1 series is found, a custom filename is not possible
      kwargs['filename'] = None

    for ds in datasets:  # Multiple datasets, so generate name from DICOM
      volume_idx = 0
      for volume in datasets[ds].values():
        sub_volumes = {}
        split_3d_keys = self.config.get('split_3D', [])

        volume.sortSlices()
        if volume.is_valid and not volume.is_equidistant and len(split_3d_keys) > 0:
          self.logger.info('Splitting volume by keys: %s', split_3d_keys)
          for dic_slice in volume.slices:
            set_id = []
            for k in split_3d_keys:
              set_id.append(str(getattr(dic_slice, k, None)))
            if len(set_id) == 1:
              set_id = set_id[0]
            else:
              set_id = tuple(set_id)

            if set_id not in sub_volumes:
              sub_volumes[set_id] = dicomvolume.DicomVolume()

            sub_volumes[set_id].addSlice(dic_slice)
        else:
          sub_volumes[0] = volume

        for v in sub_volumes:
          if kwargs.get('just_check', False):
            self.check_volume(sub_volumes[v], ds, volume_idx)
          else:
            self.process_volume(sub_volumes[v], volume_idx, destination, **kwargs)

          volume_idx += 1

  def process_volume(self, volume, volume_idx, destination, **kwargs):
    filename = kwargs.get('filename', None)
    structure = self.config.get('structure', None)
    dump_protocol = self.config.get('dump_protocol', False)

    try:
      if len(volume.slices) == 0:  # No files for this series UID (maybe not image storage?)
        self.logger.debug('No files for this series...')
        return

      if filename is None:  # Generate a filename from DICOM metadata
        filename = volume.build_filename()

      if volume_idx is not None and volume_idx > 0:
        filename = '%s (%d)' % (filename, volume_idx)

      if structure == 'dicom':
        # Remove invalid characters from dir names
        patient_dir = volume.descriptor.patient_name
        study_dir = volume.descriptor.study_date
        for c in r'[]/\;,><&*:%=+@!#^()|?^':
          patient_dir = patient_dir.replace(c, '')
          study_dir = study_dir.replace(c, '')

        destination = os.path.join(destination, patient_dir, study_dir)

      if volume.check_4D():
        if not self._process4d(volume, filename, destination):
          self.logger.warning(
            "Volume is 4D, but was not able to split (patient {patient_name}, studydate {study_date}, "
            "series {series_number:d}. {series_description}), skipping...".format(**volume.descriptor.__dict__))
          return
      else:
        im = volume.getSimpleITKImage()
        if im is None:
          return

        self.logger.info('Generating NRRD for patient {patient_name}, studydate {study_date}, series {series_number:d}.'
                         ' {series_description}, (%i slices)'.format(**volume.descriptor.__dict__), len(volume.slices))

        self._store_image(im, destination, filename, volume.descriptor, len(volume.slices), volume.getOrientation())

      if dump_protocol:
        protocol_fname = os.path.join(destination, filename) + '_protocol.txt'
        volume.writeProtocol(protocol_fname)
      self.counter += 1

    except KeyboardInterrupt:
      raise
    except:
      self.logger.error('Oh Oh... something went wrong...', exc_info=True)

  def _process4d(self, volume: dicomvolume.DicomVolume, filename, destination):
    if volume.split4D('DiffusionBValue', 2000) and volume.sortSlices4D():
      # Volume is DWI
      self.logger.info('Volume is 4D, splitting DWI on standard bvalue tag')
      prefix = '_b'
      split_tag = 'DiffusionBValue'
      split_unit = 'sec/mm2'
    elif volume.split4D(0x0019100c, 2000) and volume.sortSlices4D():
      # Volume is DWI
      self.logger.info('Volume is 4D, splitting DWI on private bvalue tag')
      prefix = '_b'
      split_tag = '(0019, 100c)'
      split_unit = 'sec/mm2'
    elif volume.split4D('TemporalPositionIdentifier') and volume.sortSlices4D():
      # Volume is 4D by timepoint
      self.logger.info('Volume is 4D, splitting on Temporal Position tag')
      prefix = '_'
      split_tag = 'TemporalPosition'
      split_unit = 'sec'
    else:
      return False

    orientation = volume.getOrientation()
    if self.config.get('combine4d', False):
      slice_count = 0
      ims = {pos: im for pos, im, sliceCount in volume.getSimpleITK4DImage()}
      positions = list(sorted(ims.keys()))
      im4d: sitk.Image = sitk.Compose([ims[t] for t in positions])
      im4d.SetMetaData('MultiVolume.FrameIdentifyingDICOMTagName', split_tag)
      im4d.SetMetaData('MultiVolume.FrameIdentifyingDICOMTagUnits', split_unit)
      im4d.SetMetaData('MultiVolume.FrameLabels', ','.join(str(p) for p in positions))
      im4d.SetMetaData('MultiVolume.NumberOfFrames', str(len(positions)))

      self.logger.info('Generating NRRD for patient {patient_name}, studydate {study_date}, series {series_number:d}. '
                       '{series_description}, (%i slices)'.format(**volume.descriptor.__dict__), slice_count)
      self._store_image(im4d, destination, filename, volume.descriptor, slice_count, orientation)
      return True
    else:
      positions = 0
      for position, im, slice_count in volume.getSimpleITK4DImage():
        positions += 1
        if im is None:
          continue

        self.logger.info('Generating NRRD for patient {patient_name}, studydate {study_date}, '
                         'series {series_number:d}. {series_description}, '
                         'temporal position %s (%i slices)'.format(**volume.descriptor.__dict__),
                         position, slice_count)

        pos_fname = filename + prefix + str(position)
        self._store_image(im, destination, pos_fname, volume.descriptor, slice_count, orientation)

      if positions == 0:
        self.logger.warning('4D volume contains standard bvalue tag, but unable to perform correct split.')
      else:
        self.logger.debug('4D volume processing complete, found %i positions', positions)
        return True

    return False

  def _store_image(self, im, destination, fname, descriptor, slicecount, orientation):
    fileformat = self.config.get('format', 'nrrd')
    overwrite = self.config.get('overwrite', False)
    output_writer = self.config.get('output_writer', None)

    target = os.path.join(destination, fname)

    nrrd_fname = target + '.' + fileformat

    if destination != '' and not os.path.isdir(destination):
      self.logger.debug('Creating study directory "%s"', destination)
      os.makedirs(destination)
    elif os.path.isfile(fname):
      if overwrite:
        self.logger.warning('file "%s" already exists, overwriting...', fname)
      else:
        self.logger.warning('file "%s" already exists, skipping...', fname)
        return

    self.logger.info('Storing in %s', nrrd_fname)
    sitk.WriteImage(im, str(nrrd_fname), useCompression=self.config.get('compress', False))

    if output_writer is not None:
      self.logger.debug('Storing location in CSV output')
      output_writer.writerow([
        self.counter,
        descriptor.patient_name,
        descriptor.study_date,
        nrrd_fname.replace(os.path.sep, '/'),
        slicecount,
        orientation
      ])

  def check_volume(self, volume, uid, volume_idx=0):
    try:
      if len(volume.slices) == 0:  # No files for this series UID (maybe not image storage?)
        self.logger.debug('No files for this series...')
        return

      volume.sortSlices()
      if volume.is_equidistant and volume.is_valid:
        if volume_idx > 0:
          self.logger.info('DicomVolume %s %s, (volume %d) is valid...',
                           uid, getattr(volume.slices[0], 'SeriesDescription', 'N/A'), volume_idx + 1)
        else:
          self.logger.info('DicomVolume %s %s is valid...',
                           uid, getattr(volume.slices[0], 'SeriesDescription', 'N/A'))
    except KeyboardInterrupt:
      raise
    except:
      self.logger.error('Oh Oh... something went wrong...', exc_info=True)
