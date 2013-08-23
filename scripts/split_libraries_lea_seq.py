#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from itertools import izip
from qcli import parse_command_line_parameters, make_option
from qiime.golay import decode_golay_12, get_invalid_golay_barcodes
from qiime.parse import MinimalFastqParser
from qiime.split_libraries import check_map
from qiime.split_libraries_fastq import correct_barcode, FastqParseError
from qiime.util import create_dir

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = []
script_info['script_usage'].append(("","",""))
script_info['output_description']= ""
script_info['required_options'] = [
    make_option('-i', '--sequence_read_fps', type='existing_filepaths',
                help='the forward and reverse sequence read fastq files '
                     '(comma-separated)'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='directory to store output files'),
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='metadata mapping file')
]
script_info['optional_options'] = [
    make_option('--barcode_type', type='string',
                help='the type of barcode used. This can be an integer, e.g. '
                     '6 for length 6 barcodes, or golay_12 for golay error-'
                     'correcting barcodes. Error correction will only be '
                     'applied for golay_12 barcodes [default: %default]',
                default='golay_12'),
    make_option('--max_barcode_errors', type='float',
                help='the maximum allowable number of errors in the barcode '
                     'if passing --barcode_type golay_12 [default: %default]',
                default=1.5),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    barcode_type = opts.barcode_type
    max_barcode_errors = opts.max_barcode_errors

    if max_barcode_errors < 0:
        option_parser.error("--max_barcode_errors must be greater than or "
                            "equal to zero. You provided %.4f." %
                            max_barcode_errors)

    if barcode_type == 'golay_12':
        barcode_correction_fn = decode_golay_12
        barcode_len = 12
    else:
        barcode_correction_fn = None

        try:
            barcode_len = int(barcode_type)
        except ValueError:
            option_parser.error("Invalid barcode type '%s'. The barcode type "
                                "must be either golay_12 or a positive "
                                "integer indicating the barcode length." %
                                barcode_type)

    if barcode_len < 1:
        option_parser.error("Invalid barcode length: %d. Must be greater "
                            "than zero." % barcode_len)

    seq_fps = opts.sequence_read_fps

    if len(seq_fps) != 2:
        option_parser.error("You must provide exactly two sequence read "
                            "filepaths, the first for forward reads and "
                            "second for reverse reads. You specified %d "
                            "filepaths." % len(seq_fps))

    create_dir(opts.output_dir)

    with open(opts.mapping_fp, 'U') as map_f:
        # Ensures that sample IDs and barcodes are unique, that barcodes are
        # all the same length, and that primers are present. Ensures barcodes
        # and primers only contain valid characters.
        _, _, bc_to_sid, _, _, bc_to_primers, _ = check_map(map_f, False)

    # Make sure our barcodes (which are guaranteed to be the same length at
    # this point) are the correct length that the user specified.
    barcode_len_in_map = len(bc_to_sid.keys()[0])
    if barcode_len_in_map != barcode_len:
        option_parser.error("Barcodes in mapping file are of length %d, but "
                            "expected barcodes of length %d." %
                            (barcode_len_in_map, barcode_len))

    if barcode_type == 'golay_12':
        invalid_golay_barcodes = get_invalid_golay_barcodes(bc_to_sid.keys())

        if invalid_golay_barcodes:
            option_parser.error("Some or all barcodes in the mapping file are "
                                "not valid golay codes. Do they need to be "
                                "reverse complemented? If these are not golay "
                                "barcodes pass --barcode_type 12 to disable "
                                "barcode error correction, or pass "
                                "--barcode_type # if the barcodes are not 12 "
                                "base pairs, where # is the size of the "
                                "barcodes.\n\nInvalid barcodes: %s" %
                                ' '.join(invalid_golay_barcodes))

    header_idx = 0
    seq_idx = 1
    qual_idx = 2
    fwd_read_f = open(seq_fps[0], 'U')
    rev_read_f = open(seq_fps[1], 'U')

    count_barcode_errors_exceed_max = 0

    for fwd_read, rev_read in izip(
            MinimalFastqParser(fwd_read_f, strict=False),
            MinimalFastqParser(rev_read_f, strict=False)):
        # Confirm match between read headers.
        if fwd_read[header_idx] != rev_read[header_idx]:
            raise PairedEndParseError("Headers of forward and reverse reads "
                                      "do not match. Confirm that the forward "
                                      "and reverse read fastq files that you "
                                      "provided have headers that match one "
                                      "another.")
        else:
            header = fwd_read[header_idx]

        # Grab the barcode sequence. It is always at the very end of the
        # forward read.
        barcode = fwd_read[seq_idx][-barcode_len:]

        # Correct the barcode (if applicable) and map to sample ID.
        num_barcode_errors, corrected_barcode, _, sample_id = correct_barcode(
                barcode, bc_to_sid, barcode_correction_fn)

        # Skip barcodes with too many errors.
        if num_barcode_errors > max_barcode_errors:
          count_barcode_errors_exceed_max += 1
          continue

    fwd_read_f.close()
    rev_read_f.close()

    print count_barcode_errors_exceed_max


class PairedEndParseError(FastqParseError):
    pass


if __name__ == "__main__":
    main()
