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

from . import dicomvolume, post_processing


class Walker:
  def __init__(self, **kwargs):
    self.logger = logging.getLogger('nrrdify.walker')
    self.config = kwargs
    self.counter = 0
    self.post_processing = post_processing

  def run(self, source, destination, filename=None, just_check=False):
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
                self.logger.debug('Dicom file %s is RT dose storage', os.path.join(curdir, fname))
              elif 'Image Storage' not in sop_class.name:  # Not RT Dose, is the sopclass valid?
                continue
              elif imagetype is None:  # Sop class is valid, does it contain image type?
                self.logger.debug("missing Image Type tag in dicom file %s", os.path.join(curdir, fname))
                continue  # Error cannot sort, so skip and go To next file

              if imagetype is None:
                imagetype = 'RTDOSE'
              else:
                imagetype = tuple(imagetype)

              if series_uid not in datasets:
                datasets[series_uid] = {}

              if imagetype not in datasets[series_uid]:
                datasets[series_uid][imagetype] = dicomvolume.DicomVolume(self.post_processing)
              else:
                datasets[series_uid][imagetype].addSlice(dicfile)
            except KeyboardInterrupt:
              return
            except:
              self.logger.error('DOH!! Something went wrong!', exc_info=True)
        if process_per_folder:
          if structure == 'source':
            dest = os.path.join(destination, os.path.relpath(curdir, source))
            if not os.path.isdir(dest):
              self.logger.debug('Creating output directory "%s"', dest)
              os.makedirs(dest)
          else:
            dest = destination

          self._processResults(datasets, dest, filename, just_check)
          datasets = {}

    if not process_per_folder:
      self._processResults(datasets, destination, filename, just_check)

  def _processResults(self, datasets, destination, filename=None, just_check=False):
    if just_check:
      for ds in datasets:
        for volume_idx, volume in enumerate(datasets[ds].values()):
          self.checkVolume(volume, ds, volume_idx)
    else:
      # Done scanning files, now make some NRRDs out of them!
      self.logger.info('Input folder scanned, found %d unique DICOM series', len(datasets))
      if len(datasets) > 1:  # If more than 1 series is found, a custom filename is not possible
        filename = None
      for ds in datasets:  # Multiple datasets, so generate name from DICOM
        for volume_idx, volume in enumerate(datasets[ds].values()):
          self.processVolume(volume, volume_idx, destination, filename)

  def processVolume(self, volume, volume_idx, destination, filename):
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
        if not self._process4D(volume, filename, destination):
          self.logger.warning(
            "Volume is 4D, but was not able to split (patient {patient_name}, studydate {study_date), "
            "series {series_number:d}. {series_description}), skipping...".format(**volume.descriptor.__dict__))
          return
      else:
        im = volume.getSimpleITKImage()
        if im is None:
          return

        self.logger.info('Generating NRRD for patient {patient_name}, studydate {study_date), '
                         'series {series_number:d}. {series_description}, (%i slices)'.format(**volume.descriptor.__dict__),
                    len(volume.slices))

        self._store_image(im, destination, filename, volume.descriptor, len(volume.slices))

      if dump_protocol:
        protocol_fname = os.path.join(destination, filename) + '_protocol.txt'
        volume.writeProtocol(protocol_fname)
      self.counter += 1

    except KeyboardInterrupt:
      raise
    except:
      self.logger.error('Oh Oh... something went wrong...', exc_info=True)

  def _process4D(self, volume: dicomvolume.DicomVolume, filename, destination):
    if volume.split4D('DiffusionBValue', 2000) and volume.sortSlices4D():
      # Volume is DWI
      self.logger.info('Volume is 4D, splitting DWI on standard bvalue tag')
      prefix = '_b'
      splitTag = 'DiffusionBValue'
      splitUnit = 'sec/mm2'
    elif volume.split4D(0x0019100c, 2000) and volume.sortSlices4D():
      # Volume is DWI
      self.logger.info('Volume is 4D, splitting DWI on private bvalue tag')
      prefix = '_b'
      splitTag = '(0019, 100c)'
      splitUnit = 'sec/mm2'
    elif volume.split4D('TemporalPositionIdentifier') and volume.sortSlices4D():
      # Volume is 4D by timepoint
      self.logger.info('Volume is 4D, splitting on Temporal Position tag')
      prefix = '_'
      splitTag = 'TemporalPosition'
      splitUnit = 'sec'
    else:
      return False

    if self.config.get('combine4d', False):
      sliceCount = 0
      ims = {pos: im for pos, im, sliceCount in volume.getSimpleITK4DImage()}
      positions = list(sorted(ims.keys()))
      im4d: sitk.Image = sitk.Compose([ims[t] for t in positions])
      im4d.SetMetaData('MultiVolume.FrameIdentifyingDICOMTagName', splitTag)
      im4d.SetMetaData('MultiVolume.FrameIdentifyingDICOMTagUnits', splitUnit)
      im4d.SetMetaData('MultiVolume.FrameLabels', ','.join(str(p) for p in positions))
      im4d.SetMetaData('MultiVolume.NumberOfFrames', str(len(positions)))

      self.logger.info('Generating NRRD for patient {patient_name}, studydate {study_date}, series {series_number:d}. '
                       '{series_description}, (%i slices)'.format(**volume.descriptor.__dict__), sliceCount)
      self._store_image(im4d, destination, filename, volume.descriptor, sliceCount)
      return True
    else:
      positions = 0
      for position, im, sliceCount in volume.getSimpleITK4DImage():
        positions += 1
        if im is None:
          continue

        self.logger.info('Generating NRRD for patient {patient_name}, studydate {study_date}, '
                         'series {series_number:d}. {series_description}, '
                         'temporal position %s (%i slices)'.format(**volume.descriptor.__dict__),
                          position, sliceCount)

        pos_fname = filename + prefix + str(position)
        self._store_image(im, destination, pos_fname, volume.descriptor, sliceCount)

      if positions == 0:
        self.logger.warning('4D volume contains standard bvalue tag, but unable to perform correct split.')
      else:
        self.logger.debug('4D volume processing complete, found %i positions', positions)
        return True

    return False

  def _store_image(self, im, destination, fname, descriptor, slicecount, header_updates=None):
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

    if header_updates is not None and fileformat == 'nrrd':
      self.logger.info('Updating NRRD headers of %s', nrrd_fname)
      with open(str(nrrd_fname), mode='rb+') as fs:
        hdr_lines = []
        new_line = fs.readline()
        while new_line != b'\n':
          hdr_lines.append(new_line)
          new_line = fs.readline()
        data = fs.read()
        fs.seek(0)
        fs.writelines(hdr_lines)
        fs.writelines([('%s:=%s\n' % kv).encode('ascii') for kv in header_updates.items()])
        fs.write(b'\n')
        fs.write(data)

    if output_writer is not None:
      self.logger.debug('Storing location in CSV output')
      output_writer.writerow([
        self.counter,
        descriptor.patient_name,
        descriptor.study_date,
        nrrd_fname.replace(os.path.sep, '/'),
        slicecount
      ])

  def checkVolume(self, dicomVolume, uid, volume_idx=0):
    try:
      if len(dicomVolume.slices) == 0:  # No files for this series UID (maybe not image storage?)
        self.logger.debug('No files for this series...')
        return

      dicomVolume.sortSlices()
      if dicomVolume.is_equidistant and dicomVolume.is_valid:
        if volume_idx > 0:
          self.logger.info('DicomVolume %s, (volume %d) is valid...', uid, volume_idx + 1)
        else:
          self.logger.info('DicomVolume %s is valid...', uid)
    except KeyboardInterrupt:
      raise
    except:
      self.logger.error('Oh Oh... something went wrong...', exc_info=True)
