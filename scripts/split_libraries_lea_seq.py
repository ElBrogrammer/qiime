#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from qcli import parse_command_line_parameters, make_option
from qiime.golay import decode_golay_12
from qiime.split_libraries import check_map
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
                default='golay_12')
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    barcode_type = opts.barcode_type

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


if __name__ == "__main__":
    main()
