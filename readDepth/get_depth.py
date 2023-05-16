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

from readDepth.allele_counter import AlleleCounter

def allele_depths(bam, chrom, pos, ref, alt, min_quality=0, map_quality=30,
        by_strand=True, stepper=None):
    counter = AlleleCounter(bam, min_qual=min_quality, map_qual=map_quality,
        by_strand=by_strand)
    try:
        return counter(chrom, pos, ref, alt, stepper=stepper)
    except AssertionError:
        if by_strand:
            return {'ref': {'forward': None, 'reverse': None},
                    'alt': {'forward': None, 'reverse': None}}
        else:
            return {'ref': None, 'alt': None}

def get_multiple_depths(bam, variants, min_quality=0, map_quality=30,
        by_strand=True, stepper=None):
    ''' get the read depths of the alleles at variable sites
    
    Args:
        bam: pysam.AlignmentFile
        variants: list of (chrom, start, end, ref, alt) tuples.
        min_quality: threshold for minimum base quality to include bases in the
            counts. Defaults to zero, ie includes all reads.
        map_quality: threshold for mapping quality. Defaults to 30.
        by_strand: whether to return read depths by strand.
    
    Returns:
        list of dictionaries of allele depths by base.
    '''
    counts = {}
    for chrom, pos, ref, alt in variants:
        key = (chrom, pos, ref, alt)
        counts[key] = allele_depths(bam, chrom, pos, ref, alt, min_quality, map_quality,
                                    by_strand, stepper)
    return counts
