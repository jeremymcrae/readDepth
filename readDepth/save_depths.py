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

import gzip

import pandas

def save_depths(depths, person_id, path, compress=False, include_header=True):
    ''' write the depths as a table to a file.
    
    Args:
        depths: list of dictionaries of allele depths by base.
        person_id: sample ID for the individual.
        path: path to output file.
        compress: whether to gzip compress the data in the file.
        include_header: whether to include the table column names in the output.
    '''
    
    d = []
    for key in sorted(depths.keys()):
        x = depths[key]
        tup = (x['ref']['forward'], x['ref']['reverse'], x['alt']['forward'], x['alt']['reverse'])
        d.append(tup)
    
    depths = pandas.concat([pandas.DataFrame({'person_id': [person_id] * len(d)}),
        pandas.DataFrame(sorted(depths.keys()),
            columns=['chrom', 'start_pos', 'end_pos', 'ref', 'alt']),
        pandas.DataFrame(d, columns=['ref_F', 'ref_R', 'alt_F', 'alt_R'])], axis=1)
    
    try:
        # if we want to compress the data, then open the file with gzip.open.
        # NOTE: this will compress each individuals data seperately, which is
        # less efficient (space-wise) than if we wrote all individuals together,
        # but then we'd have to keep a >10 million line dataset in memory.
        handle = open(path, 'a')
        if compress:
            handle = gzip.open(path, 'at')
    except TypeError:
        # catch if the path is a buffer, such as sys.stdout.
        # TODO: figure out how to compress the data to the buffer
        handle = path
        
    depths.to_csv(handle, sep='\t', index=False, na_rep='NA', header=include_header)
