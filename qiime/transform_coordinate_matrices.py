#!/usr/bin/env python
# File created on 21 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Justin Kuczynski", "Kyle Patnode",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

import csv
from os import makedirs
from os.path import split, splitext, exists, join
from random import shuffle
from numpy import array, mean, append, zeros
from cogent.cluster.procrustes import procrustes
from cogent.util.dict2d import Dict2D
from qiime.util import create_dir
from qiime.parse import fields_to_dict, parse_coords
from qiime.format import format_coords, format_p_value_for_num_iters

SUMMARY_HEADER = (
        ('# FP1', 'FP2', 'Num included dimensions', 'Monte Carlo p-value',
         'Count better', 'M^2'),
        ('# Warning: p-values in this file are NOT currently adjusted for '
         'multiple comparisons.',)
)

def shuffle_full_matrix(m):
    """  """
    shape = m.shape
    data = m.flatten()
    shuffle(data)
    return array(data).reshape(shape)
    
def shuffle_within_rows(m):
    result = []
    for k in m:
        l = list(k)
        shuffle(l)
        result.append(l)
    return array(result)
    
def shuffle_col_order(m):
    order = range(m.shape[1])
    shuffle(order)
    result = []
    current_row = []
    for k in m:
        for o in order:
            current_row.append(k[o])
        result.append(current_row)
        current_row = []
    return array(result)

def shuffle_row_order(m):
    """ jk thinks this should be the one used for monte carlo analyses

    procrustes takes samples (rows) by dimensions(cols)"""
    order = range(m.shape[0])
    shuffle(order)
    return array(m[order]) # new object


def map_sample_ids(sample_ids,sample_id_map):
    """ Map sample ids to new values in sample_id_map """
    try:
        result = [sample_id_map[sample_id] for sample_id in sample_ids]
    except KeyError:
        raise KeyError("Unknown sample ID: %s" % sample_id)
        
    return result

def reorder_coords(coords,sample_ids,order):
    """ Arrange the rows in coords to correspond to order 
    
        Note: order is the master list here -- if a sample id is
        not included in order, that coord will not be included in
        the results. All sample ids in order however must be in
        sample_ids
    
    """
    try:
        result =  array([coords[sample_ids.index(sample_id)]
                         for sample_id in order])
    except ValueError:
        raise ValueError('Unknown sample ID: %s' % sample_id)
    return result

def filter_coords_matrix(coords,dimensions_to_keep):
    return coords[:,:dimensions_to_keep]

def pad_coords_matrix(coords,dimensions_to_add):
    if dimensions_to_add < 0:
        raise ValueError('Dimensions to add must be >= 0.')
    elif dimensions_to_add == 0:
        return coords
    else:
        result = []
        # find a more efficent way to do this (don't have
        # internet right now)
        return array([list(c) + [0.] * dimensions_to_add for c in coords])

def pad_coords_matrices(coords1,coords2):
    # Determine how many dimensions are in each vector, and the difference
    coords1_len = coords1.shape[1]
    coords2_len = coords2.shape[1]
    coords_len_diff = coords1_len - coords2_len
    # If the first vector is shorter, pad it with zeros
    if coords_len_diff < 0:
        coords1 = pad_coords_matrix(coords1,-1 * coords_len_diff)
    # If the second vector is shorter, pad it with zeros        
    elif coords_len_diff > 0:
        coords2 = pad_coords_matrix(coords2,coords_len_diff)
    else:
        pass
    return coords1, coords2
    
def get_mean_percent_variation(v1,v2):
    return [mean([p1,p2]) for p1,p2 in zip(v1,v2)]

def get_mean_eigenvalues(v1,v2):
    return [mean([p1,p2]) for p1,p2 in zip(v1,v2)]

def get_procrustes_results(coords_f1,coords_f2,sample_id_map=None,
    randomize=None,max_dimensions=None,
    get_eigenvalues=get_mean_eigenvalues,
    get_percent_variation_explained=get_mean_percent_variation):
    """ """
    # Parse the PCoA files
    sample_ids1, coords1, eigvals1, pct_var1 = parse_coords(coords_f1)
    sample_ids2, coords2, eigvals2, pct_var2 = parse_coords(coords_f2)
    if sample_id_map:
        sample_ids1 = map_sample_ids(sample_ids1,sample_id_map)
        sample_ids2 = map_sample_ids(sample_ids2,sample_id_map)
    # rearrange the order of coords in coords2 to correspond to 
    # the order of coords in coords1 
    order = list(set(sample_ids1) & set(sample_ids2)) 
    coords1 = reorder_coords(coords1,sample_ids1,order)
    coords2 = reorder_coords(coords2,sample_ids2,order)
    if len(order) == 0:
        raise ValueError('No overlapping samples in the two files')
        
    # If this is a random trial, apply the shuffling function passed as 
    # randomize()
    if randomize:
        coords2 = randomize(coords2)
        randomized_coords2 = format_coords(coord_header=order, coords=coords2,
                                           eigvals=eigvals2, pct_var=pct_var2)
    else:
        randomized_coords2 = None
        
        
    coords1, coords2 = pad_coords_matrices(coords1,coords2)
    if max_dimensions is None:
        if len(pct_var1)>len(pct_var2):
            pct_var2 = append(pct_var2,zeros(len(pct_var1)-len(pct_var2)))
            eigvals2 = append(eigvals2,zeros(len(eigvals1)-len(eigvals2)))
        elif len(pct_var1)<len(pct_var2):
            pct_var1 = append(pct_var1,zeros(len(pct_var2)-len(pct_var1)))
            eigvals1 = append(eigvals1,zeros(len(eigvals2)-len(eigvals1)))
    else:
        coords1 = filter_coords_matrix(coords1,max_dimensions)
        coords2 = filter_coords_matrix(coords2,max_dimensions)
        pct_var1 = pct_var1[:max_dimensions]
        pct_var2 = pct_var2[:max_dimensions]
        eigvals1 = eigvals1[:max_dimensions]
        eigvals2 = eigvals2[:max_dimensions]

    # Run the Procrustes analysis
    transformed_coords_m1, transformed_coords_m2, m_squared =\
     procrustes(coords1,coords2)
    
    eigvals = get_eigenvalues(eigvals1, eigvals2)
    pct_var = get_percent_variation_explained(pct_var1,pct_var2)
    
    transformed_coords1 = format_coords(coord_header=order,
                                        coords=transformed_coords_m1,
                                        eigvals=eigvals,
                                        pct_var=pct_var)
    transformed_coords2 = format_coords(coord_header=order,
                                        coords=transformed_coords_m2,
                                        eigvals=eigvals,
                                        pct_var=pct_var)

    # Return the results
    return (transformed_coords1, transformed_coords2, m_squared,
            randomized_coords2)

def procrustes_monte_carlo(coords_f1,
                           coords_f2,
                           trials=1000,
                           max_dimensions=None,
                           shuffle_f=shuffle_row_order,
                           sample_id_map=None,
                           trial_output_dir=None):
    """ Run procrustes analysis with random trials
    
        This analysis could be made more efficient, as the current version 
        just calls get_procrustes_results() random_trials times, which involves
        re-parsing each time, etc. It is very fast though, so that's low 
        priority.
    
    """
    # Get the M^2 for the actual data
    actual_m_squared = get_procrustes_results(
     coords_f1,
     coords_f2,
     sample_id_map=sample_id_map,
     randomize=None,
     max_dimensions=max_dimensions)[2]

    if max_dimensions is None:
        max_dimensions_str = 'alldim'
    else:
        max_dimensions_str = 'maxdim_%d' % max_dimensions 

    # prep for storing trial details if a trial output dir 
    # was provided
    if trial_output_dir:
        create_dir(trial_output_dir)
        trail_summary_fp = '%s/trial_summary_%s.txt' % (trial_output_dir,
                                                        max_dimensions_str)
        trial_summary_f = open(trail_summary_fp,'w')
        trial_summary_f.write('trial id\ttrial M^2\n')
        
    # Get the M^2 for the random trials, and count how many
    # are lower than or equal to the actual M^2
    trial_m_squareds = []
    count_better = 0
    for i in range(trials):
        # perform the procrustes analysis
        transformed_coords1, transformed_coords2, trial_m_squared, randomized_coords2 =\
          get_procrustes_results(
             coords_f1,
             coords_f2,
             sample_id_map=sample_id_map,
             randomize=shuffle_f,
             max_dimensions=max_dimensions)
        trial_m_squareds.append(trial_m_squared)
        
        if trial_m_squared <= actual_m_squared:
            count_better += 1
            
        # write the transformed coordinate matrices, if they're being
        # stored
        if trial_output_dir:
            trial_id = '%s_trial%d' % (max_dimensions_str,i)
            output_matrix1_fp = '%s/pc1_transformed_%s.txt' \
             % (trial_output_dir,trial_id)
            output_matrix2_fp = '%s/pc2_transformed_%s.txt' \
             % (trial_output_dir,trial_id)
            output_matrix3_fp = '%s/pc2_randomized_trial%s.txt' \
             % (trial_output_dir,i)
            output_matrix1_f = open(output_matrix1_fp,'w')
            output_matrix1_f.write(transformed_coords1)
            output_matrix1_f.close()
            output_matrix2_f = open(output_matrix2_fp,'w')
            output_matrix2_f.write(transformed_coords2)
            output_matrix2_f.close()
            output_matrix3_f = open(output_matrix3_fp,'w')
            output_matrix3_f.write(randomized_coords2)
            output_matrix3_f.close()
            trial_summary_f.write('%s\t%2.2f\n' % (trial_id,trial_m_squared))
    
    # Close the trial summary file, if one is being created
    if trial_output_dir:
        trial_summary_f.close()
    
    return actual_m_squared, trial_m_squareds, count_better, count_better/trials

def transform_coordinate_matrices(output_dir, input_fps, sid_map_fps=None,
                                  num_dims=None, random_trials=None,
                                  store_trial_details=False):
    if num_dims is not None and num_dims < 1:
        raise ValueError("Invalid number of dimensions: %d. Must be greater "
                         "than zero." % num_dims)

    if random_trials is not None and random_trials < 10:
        raise ValueError("Invalid number of trials: %d. Must perform 10 or "
                         "more trials for Monte Carlo analysis." %
                         random_trials)

    if len(input_fps) < 2:
        raise ValueError("Must provide at least two coordinate matrices as "
                         "input.")

    if sid_map_fps and (len(sid_map_fps) + 1) != len(input_fps):
        raise ValueError("If providing sample id maps, there must be exactly "
                         "one fewer sample id maps than input coordinate "
                         "matrices.")

    if num_dims is None:
        max_dims_str = 'alldim'
    else:
        max_dims_str = str(num_dims)

    if not exists(output_dir):
        makedirs(output_dir)

    reference_input_fp = input_fps[0]
    reference_input_fp_dir, input_fn1 = split(reference_input_fp)
    reference_input_fp_basename, reference_input_fp_ext = splitext(input_fn1)

    summary_lines = []
    for i, query_input_fp in enumerate(input_fps[1:]):
        query_input_fp_dir, query_input_fn = split(query_input_fp)
        query_input_fp_basename, query_input_fp_ext = splitext(query_input_fn)
        output_matrix1_fp = join(output_dir,
                '%s_transformed_reference.txt' % reference_input_fp_basename)
        output_matrix2_fp = join(output_dir,
                '%s_transformed_q%d.txt' % (query_input_fp_basename, i + 1))

        if sid_map_fps:
            sample_id_map = _build_sample_id_map(open(sid_map_fps[i], 'U'))
        else:
            sample_id_map = None

        trans_coords1, trans_coords2, m2, _ = get_procrustes_results(
                open(reference_input_fp, 'U'), open(query_input_fp, 'U'),
                sample_id_map=sample_id_map, randomize=False,
                max_dimensions=num_dims)

        with open(output_matrix1_fp, 'w') as output_matrix1_f:
            output_matrix1_f.write(trans_coords1)
        with open(output_matrix2_fp, 'w') as output_matrix2_f:
            output_matrix2_f.write(trans_coords2)

        if random_trials:
            if store_trial_details:
                trial_output_dir = join(output_dir, 'trial_details_%d' % i + 2)
            else:
                trial_output_dir = None

            coords_f1 = list(open(reference_input_fp, 'U'))
            coords_f2 = list(open(query_input_fp, 'U'))
            _, _, count_better, p_val = procrustes_monte_carlo(coords_f1,
                    coords_f2, trials=random_trials, max_dimensions=num_dims,
                    sample_id_map=sample_id_map,
                    trial_output_dir=trial_output_dir)

            # truncate the p-value to the correct number of significant digits
            p_val_str = format_p_value_for_num_iters(p_val, random_trials)
        else:
            p_val_str = 'NA'
            count_better = 'NA'

        summary_lines.append((reference_input_fp, query_input_fp,
                              max_dims_str, p_val_str, count_better,
                              '%1.3f' % m2))

    # Write a TSV summary to disk.
    output_summary_fp = join(output_dir, 'procrustes_results.txt')
    _write_summary_lines(output_summary_fp, summary_lines)

def _build_sample_id_map(sample_id_map_f):
    sample_id_map = {}

    for k, v in fields_to_dict(sample_id_map_f).items():
        sample_id_map[k] = v[0]

    return sample_id_map

def _write_summary_lines(output_fp, lines, delimiter='\t'):
    with open(output_fp, 'wb') as output_f:
        summary_writer = csv.writer(output_f, delimiter=delimiter)
        summary_writer.writerows(SUMMARY_HEADER)
        summary_writer.writerows(lines)
