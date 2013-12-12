#!/usr/bin/env python
# File created on 21 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Justin Kuczynski",
               "Jose Carlos Clemente Litran", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from qiime.util import make_option
from qiime.util import parse_command_line_parameters
from qiime.transform_coordinate_matrices import transform_coordinate_matrices

script_info={}
script_info['brief_description']="""Transform two or more coordinate matrices"""
script_info['script_description']="""This script transforms two or more coordinate matrices (e.g., the output of principal_coordinates.py) using procrustes analysis to minimize the distances between corresponding points. The first coordinate matrix provided is treated as the reference, and all other coordinate matrices are transformed to minimize distances to the reference points. Monte Carlo simulations can additionally be performed (-r random trials are run) to estimate the probability of seeing an M^2 value as extreme as the actual M^2."""
script_info['script_usage']=[]
script_info['script_usage'].append(("Write the transformed procrustes matrices to file","","""%prog -i unweighted_unifrac_pc.txt,weighted_unifrac_pc.txt -o procrustes_output"""))

script_info['script_usage'].append(("Generate transformed procrustes matrices and monte carlo p-values for two principal coordinate matrices","","""%prog -i unweighted_unifrac_pc.txt,weighted_unifrac_pc.txt -o mc_procrustes_output_2 -r 1000""",))
script_info['script_usage'].append(("Generate transformed procrustes matrices and monte carlo p-values for four principal coordinate matrices","","""%prog -i unweighted_unifrac_pc.txt,weighted_unifrac_pc.txt,euclidean_pc.txt,bray_curtis_pc.txt -o mc_procrustes_output_4 -r 1000""",))
script_info['script_usage'].append(("Generate transformed procrustes matrices and monte carlo p-values for three principal coordinate matrices where the sample ids must be mapped between matrices","","""%prog -i s1_pc.txt,s2_pc.txt,s3_pc.txt -s s1_s2_map.txt,s1_s3_map.txt -o mc_procrustes_output_3 -r 1000""",))

script_info['output_description']="""Two transformed coordinate matrices corresponding to the two input coordinate matrices, and (if -r was specified) a text file summarizing the results of the Monte Carlo simulations."""
script_info['required_options']=[
 make_option('-i','--input_fps',type='existing_filepaths',
             help='comma-separated list of input coordinate matrices'),
 make_option('-o','--output_dir',type='new_dirpath',
             help='the output directory'),
]
script_info['optional_options']=[
 make_option('-r','--random_trials',type='int',
    help='Number of random permutations of matrix2 to perform. '+
    ' [default: (no Monte Carlo analysis performed)]',default=None),
 make_option('-d','--num_dimensions',type='int',default=3,
    help='Number of dimensions to include in output matrices'+
    ' [default: %default]'),
 make_option('-s','--sample_id_map_fps',
    type='existing_filepaths',
    help='If sample id maps are provided, there must be exactly one fewer files here than there are coordinate matrices (as each nth sample id map will provide the mapping from the first input coordinate matrix to the n+1th coordinate matrix) [default: %default]',
    default=None),
 make_option('--store_trial_details',
    help='Store PC matrices for individual trials [default: %default]',
    default=False,action='store_true'),
]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    transform_coordinate_matrices(opts.output_dir, opts.input_fps,
                                  opts.sample_id_map_fps, opts.num_dimensions,
                                  opts.random_trials, opts.store_trial_details)

if __name__ == "__main__":
    main()
