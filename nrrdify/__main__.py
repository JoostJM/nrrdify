#!/usr/bin/env python

# ========================================================================
#  Copyright Het Nederlands Kanker Instituut - Antoni van Leeuwenhoek
#
#  Licensed under the 3-clause BSD License
# ========================================================================

import argparse
import csv
import logging

import nrrdify.walker


def main(args=None):
  parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
  parser.add_argument("inputFolder", metavar="In", type=str, help="Folder containing the DICOM file(s) to convert.")
  parser.add_argument("--out", "-o", help="Folder to store converted files in. If omitted, stores "
                                          "files in parent directory of In folder.")
  parser.add_argument("--name", "-n", help="Filename for the new file, without extension. If omitted, or more series "
                                           "are found, Filename is generated from DICOM tags: "
                                           "<PatientName>-<StudyDate>-<SeriesNumber>. <SeriesDescription>")
  parser.add_argument('--combine4d', action='store_true',
                      help='If specified, 4D volumes are combined into a single output file')
  parser.add_argument("--format", "-f", nargs="?", default="nrrd", choices=["nrrd", "nii", "nii.gz"],
                      help="Image format to convert to. Default is the 'nrrd' format")
  parser.add_argument('--structure', '-s', choices=['none', 'source', 'dicom'], default='none',
                      help='directory structure to use in the output. "none" (Default): no structure, just store the '
                           'files, "source": copy the structure in the input. N.B. processes each folder independently.'
                           '"dicom": create tree based on dicomtags: "<patient>/<studydate>/<volume>"')
  parser.add_argument('--process-set', '-ps', action='store_false',
                      help='If specified, process the input as one set, otherwise, process set per folder. If structure'
                           'is "source", this setting has no effect (always processed per folder). '
                           'N.B. processing per set can be very memory intensive!')
  parser.add_argument('--compress', action='store_true', help='If specified, applies compression on the output')
  parser.add_argument('--csv-output', '-co', type=argparse.FileType('w'), default=None,
                      help='Specifies a new CSV-file to store the locations of the generated files, overwrites'
                           'existing files. If omitted, no CSV output is generated.')
  parser.add_argument('--dump-protocol', '-dp', action='store_true',
                      help='If specified, nrrdify will dump the DICOM metadata as a text file for each successfully '
                           'converted volume')
  parser.add_argument('--logging-level', metavar='LEVEL',
                      choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                      default='WARNING', help='Set capture level for logging')
  parser.add_argument('--log-file', metavar='FILE', type=argparse.FileType('a'), default=None,
                      help='File to append logger output to')
  parser.add_argument('--overwrite', action='store_true',
                      help='if this argument is specified, script will overwrite existing files, otherwise, file write '
                           'for already existing files is skipped.')
  parser.add_argument('--check', action='store_true',
                      help='if this argument is specified, DICOMS are checked but not converted.')
  parser.add_argument('--split-tag', default=None,
                      help='Dicom tag to split a 4D volume. If invalid, also tries DWI b-value and TemporalPositionIdentifier')
  parser.add_argument('--split3d', action='append', default=None)

  args = parser.parse_args(args)

  # if specified, set up logging to a file
  if args.log_file is not None:
    log_handler = logging.StreamHandler(args.log_file)
    log_formatter = logging.Formatter('[%(asctime)-.19s] %(levelname)-.1s: %(message)s')
    log_handler.setFormatter(log_formatter)
    log_level = getattr(logging, args.logging_level)
    log_handler.setLevel(log_level)

    nrrdify.logger.addHandler(log_handler)
    if log_level < logging.INFO:  # Lower logger level if necessary
      nrrdify.logger.setLevel(log_level)

  kwargs = {
    'format': args.format,
    'overwrite': args.overwrite,
    'structure': args.structure,
    'dump_protocol': args.dump_protocol,
    'combine4d': args.combine4d,
    'compress': args.compress,
    'split_tag': args.split_tag,
    'process_per_folder': args.process_set
  }

  source_folder = args.inputFolder
  destination_folder = args.out
  if destination_folder is None:
    nrrdify.logger.info('No destination specified, using source folder %s', source_folder)
    destination_folder = source_folder

  if args.structure == 'source':
    kwargs['process_per_folder'] = True
  else:  # none or dicom:
    kwargs['process_per_folder'] = args.process_set

  if args.csv_output is not None:
    kwargs['output_writer'] = csv.writer(args.csv_output, lineterminator='\n')
    kwargs['output_writer'].writerow(['idx', 'patient', 'studydate', 'image', 'numSlices', 'orientation'])

  if args.split3d is not None:
    kwargs['split_3D'] = args.split3d

  walker = nrrdify.walker.Walker(**kwargs)
  walker.run(source_folder, destination_folder, filename=args.name, just_check=args.check)


if __name__ == '__main__':
  main()
