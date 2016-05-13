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

from __future__ import print_function

import pysam

class AlleleCounter(object):
    
    def __init__(self, bam_path, max_coverage=1e10, min_qual=0, map_qual=0, by_strand=False):
        ''' counts reads with different base calls at a chrom position
        
        Args:
            bam_path: path to indexed bam file
            max_coverage: maximum coverage at which we stop tallying the bases
            min_qual: minimum quality score required to include the base call
            by_strand: whether to separate the ref and alt counts by forward and
                reverse reads.
        '''
        
        self.bam = pysam.AlignmentFile(bam_path)
        self.max_coverage = max_coverage
        self.min_qual = min_qual
        self.min_map_qual = map_qual
        self.by_strand = by_strand
    
    def check_variant(self, chrom, start_pos, end_pos=None, ref=None, alt=None):
        ''' counts reads with different base calls at a chrom position
        
        Args:
            chrom: chromosome to use eg 'chr1' or '1' depending on how the BAM is
                set up (specifically, an ID found in the BAMs sequence dictionary).
            start_pos: start_base position to count bases at.
            end_pos: end base position to count bases at, or None if single base SNV.
            ref: reference allele (e.g. 'T').
            alt: alternate allele (e.g. 'G').
        
        Returns:
            dictionary of read counts indexed by base calls
        '''
        
        # Make sure that if we have specified a ref allele, then we also have an alt
        # allele (and vice versa if we don't have the ref allele).
        assert type(ref) == type(alt)
        
        start_pos = int(start_pos)
        
        if end_pos is None:
            end_pos = start_pos + len(ref) - 1
        
        end_pos = int(end_pos)
        length = abs(end_pos - start_pos) + 1
        
        alleles = {'ref': 0, 'alt': 0}
        if self.by_strand:
            alleles = {'ref': {'forward': 0, 'reverse': 0},
                'alt': {'forward': 0, 'reverse': 0}}
        
        # count each base at the required site
        # for pileupcolumn in self.bam.pileup(chrom, start_pos - 1, end_pos + 1, truncate=True, stepper='nofilter'):
        for pileupcolumn in self.bam.pileup(chrom, start_pos - 1, end_pos + 1, truncate=True):
            if pileupcolumn.pos != start_pos - 1:
                continue
            
            for read in pileupcolumn.pileups:
                allele = self.check_read(read, length, ref, alt)
                
                # only include tallies for the ref and alt alleles, don't count
                # reads with other alleles
                if allele is not None:
                    if self.by_strand:
                        strands = {False: 'forward', True: 'reverse'}
                        strand = strands[read.alignment.is_reverse]
                        alleles[allele][strand] += 1
                    else:
                        alleles[allele] += 1
        
        return alleles
    
    def check_read(self, read, length, ref, alt):
        ''' gets the allele code for a variant site in a read
        
        Args:
            read: AlignmentRead from pileup
            length: length of variant in bp, from start pos to end_pos
            ref: sequence for reference allele
            alt: sequence for alternate alllele
        '''
        
        pos = read.query_position
        qual = self.min_qual + 1
        if pos is not None:
            qual = read.alignment.query_qualities[pos:pos + length][0]
        
        # ignore low qual reads
        if qual < self.min_qual:
            return None
        
        if read.alignment.mapping_quality < self.min_map_qual:
            return None
        
        return self.check_allele(read, ref, alt, length)
    
    def check_allele(self, read, ref, alt, length):
        ''' checks whether the sequence at a site is for a 'ref' or 'alt' allele
        
        Args:
            read: AlignmentRead from pileup
            length: length of variant in bp, from start pos to end_pos
            ref: sequence for reference allele
            alt: sequence for alternate alllele
        
        Returns:
            'ref' or 'alt' depending on whether the allele is for the reference
            sequence, or alternate to the reference sequence.
        '''
        
        indel_length = (len(alt) - len(ref))
        
        # figure out whether the read supports the reference of alternate
        # allele (ignoring reads that don't support either). We need to be
        # careful of indels.
        # TODO: I'm not sure if indel depths are accurate in every situation
        allele = None
        if len(ref) > 1 or len(alt) > 1:
            if read.indel == 0:
                allele = 'ref'
            elif read.indel == indel_length:
                allele = 'alt'
        elif read.query_position is not None:
            # if the variant is for a SNV, get the base call as a string
            seq = read.alignment.query_sequence[read.query_position:read.query_position + length]
            if seq == ref:
                allele = 'ref'
            elif seq == alt:
                allele = 'alt'
        
        return allele
