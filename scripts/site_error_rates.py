'''
Copyright (c) 2016 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import os
import sys
import argparse
import random

import pandas

from denovoFilter.load_candidates import load_candidates
from denovoFilter.preliminary_filtering import preliminary_filtering
from readDepth.get_depth import get_read_depths
from readDepth.save_depths import save_depths

BAMS_PATH = '/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/bam_paths.txt'
DE_NOVOS_PATH = '/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/denovo_gear_trios_extracted_passed_variants_11.05.15.tsv'
FAMILIES_PATH = '/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt'

def get_options():
    """ get the command line options
    """
    
    parser = argparse.ArgumentParser(description='Identifies error rates at'
        'sites within a subset of unaffected parents.')
    parser.add_argument('--bams', default=BAMS_PATH,
        help='Path to file listing BAM locations.')
    parser.add_argument('--families', default=FAMILIES_PATH,
        help='Path to file listing BAM locations.')
    parser.add_argument('--de-novos', default=DE_NOVOS_PATH,
        help='Path to file listing BAM locations.')
    parser.add_argument('--count', default=500, type=int,
        help='number of unaffected parents to sample from.')
    parser.add_argument('--compress', action='store_true', default=False,
        help='whether to gzip compress the output. If the output writes to '
            'standard out, then this option is unused and the output is not '
            'compressed.')
    
    parser.add_argument('--output', default=sys.stdout,
        help='Path to file for output. Defaults to standard out.')
    
    args = parser.parse_args()
    
    return args

def get_parental_subset(path, count=500):
    ''' select a random subset of unaffected parents
    
    Args:
        path: path to family relationships file, containing sample IDs, family
            information, and whether the individuals are affected.
        counts: number of unaffected parents to select.
    
    Returns:
        pandas Dataframe, one row per individual.
    '''
    
    families = pandas.read_table(path)
    parents = families[(families['dad_id'] == '0') & (families['mum_id'] == '0') & \
        (families['affected'] == 1)]
    
    subset = parents.loc[random.sample(list(parents.index), count), ]
    
    return subset

def open_sites(path):
    ''' load a candidate sites for checking reads depths.
    
    Args:
        path: path to file containing sites to check. Requires 'chrom', 'pos',
            'ref' and 'alt' columns.
    
    Returns:
        sorted list of variant tuples (chrom, pos, ref and alt alleles.)
    '''
    
    de_novos = load_candidates(path)
    snvs = de_novos[(de_novos.ref.str.len() == 1) & (de_novos.alt.str.len() == 1)]
    snvs['end_pos'] = snvs['pos']
    if 'maf' in de_novos.columns:
        snvs = snvs[preliminary_filtering(snvs)]
    
    # get the coordinates of all the SNVs, so we can check them all at once in a BAM
    coords = snvs[['chrom', 'pos', 'end_pos', 'ref', 'alt']]
    variants = sorted(set([ tuple(x) for x in coords.values ]))
    
    return variants

def get_bam_path(bams, person_id):
    ''' identify a bam path for a person ID. Where there are multiple BVAMs,
    pick the most likely best sample.
    
    Args:
        bams: pandas DataFrame of bam paths, listed by person_stable_id
        person_id: string ID for an individual (e.g. 'DDDP1XXXXX')
    
    Returns:
        path to BAM file. Can be an IRODS path.
    '''
    
    person_bams = bams[bams['person_stable_id'] == person_id]
    
    if len(person_bams) > 1:
        person_bams = person_bams[person_bams['sample_type'] == 'SALIVA']
    
    bam_path = person_bams['archive_uri']
    bam_path = sorted(bam_path)[-1]
    
    return bam_path

def main():
    
    args = get_options()
    
    bams = pandas.read_table(args.bams)
    variants = open_sites(args.de_novos)
    subset = get_parental_subset(args.families, count=args.count)
    
    # clear the output file if it already exists
    try:
        if os.path.exists(args.output):
            os.remove(args.output)
    except TypeError:
        pass
    
    for i, row in subset.iterrows():
        person_id = row['individual_id']
        bam_path = get_bam_path(bams, person_id)
        
        try:
            depths = get_read_depths(bam_path, variants, min_quality=0, map_quality=0, by_strand=True)
        except ValueError:
            continue
        
        include_header = i == 0
        save_depths(depths, person_id, args.output, args.compress, include_header)

if __name__ == '__main__':
    main()
