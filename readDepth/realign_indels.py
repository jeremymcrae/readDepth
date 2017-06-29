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
import subprocess
import tempfile
import shutil
import argparse

import pysam

JAVA_BIN = '/software/jdk1.7.0_25/bin/java'
GATK = '/software/hgi/pkglocal/gatk-protected-3.5/GenomeAnalysisTK.jar'
GENOME = '/lustre/scratch115/projects/ddd/vrpipe-ihtp_fy4/dataroot/hs37d5.fa'
OPTS = ['-Xmx300m']

def get_options():
    ''' parse the command line arguments
    '''
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--inbam', help='path to source BAM')
    parser.add_argument('--invcf', help='path to source VCF')
    parser.add_argument('--regions',
        help='comma-separated list of chrom:pos sites e.g. 1:1000,2:4000')
    parser.add_argument('--outbam', help='path to write realigned BAM to')
    
    parser.add_argument('--genome', default=GENOME,
        help='path to indexed genome fasta file')
    
    args = parser.parse_args()
    
    # parse and tidy up the defined regions
    args.regions = args.regions.split(',')
    args.regions = [ x.split(':') for x in args.regions ]
    args.regions = [ (x[0], int(x[1])) for x in args.regions ]
    
    return args

def sort_regions(regions):
    ''' sort the list of chromosome regions
    
    Args:
        regions: list of (chrom, position) tuples
    
    Returns:
        list of tuples, now sorted by successive genome positions.
    '''
    
    # define the chromosome order, so we can sort correctly
    strings = ( str(x) for x in list(range(1, 23)) + ['X', 'Y'] )
    nums = range(24)
    chroms = dict(zip(strings, nums))
    
    return sorted(regions, key=lambda x: (chroms[x[0]], x[1]))

def make_indel_vcf(regions, outfile):
    ''' prepare a vcf containing the indels at selected sites
    
    Args:
        vcf_path: path to source VCF
        regions: list of (chrom, position, end, ref, alt) tuples for required indels
        outfile: handle for output VCF
    '''
    
    # write minimal VCF header
    outfile.write('##fileformat=VCFv4.2\n'.encode('utf8'))
    outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tperson_id\n'.encode('utf8'))
    
    for chrom, pos, end, ref, alt in regions:
        line = '{}\t{}\t.\t{}\t{}\t.\t.\t.\t.\t.\n'.format(chrom, pos, ref, alt)
        outfile.write(line.encode('utf8'))
    
    outfile.flush()

def make_indel_bam(bam_path, regions, outfile, window=100):
    ''' create a bam for the reads surrounding the required indel regions.
    
    If we make a stripped down bam for only the regions we are interested in,
    this makes the indel realignment much faster and requires less memory.
    Creating the new bam should take less than 0.1 seconds.
    
    Args:
        bam_path: path to source bam
        regions: list of (chrom, position) tuples for required indels
        outfile: file handle for output bam
        window: window (in bp) to get reads surrounding the indel site
    '''
    
    unsorted = tempfile.NamedTemporaryFile(suffix='unsorted.bam')
    bam = pysam.AlignmentFile(bam_path)
    new_bam = pysam.AlignmentFile(unsorted.name, mode='wb', template=bam)
    
    # copy each read in the required regions to the new bam file
    for region in sort_regions(regions):
        for read in bam.fetch(region[0], region[1] - window, region[1] + window):
            new_bam.write(read)
    
    new_bam.close()
    
    # some bams might not be sorted, so sort before indexing
    pysam.sort(unsorted.name, "-o", outfile.name)
    pysam.index(outfile.name)
    
    unsorted.close()

def make_targets(regions, handle, window=50):
    ''' make a file of genome regions
    
    Args:
        regions: list of (chrom, position) tuples for required indels
        handle: file handle for output
        window: window (in bp) to define target regions surrounding the indel
    '''
    
    for region in sort_regions(regions):
        chrom, pos = region[0], region[1]
        handle.write('{}:{}-{}\n'.format(chrom, pos - window, pos + window).encode('utf8'))
    
    handle.flush()

def realign_indels(regions, inbam, outbam, java=JAVA_BIN, java_opts=OPTS,
        gatk=GATK, ref_genome=GENOME):
    ''' create bam realigned around specific indels
    
    This subsets the BAM to reads near the indels sites, for computational and
    storage efficiency.
    
    Args:
        vcf: path to source VCf
        regions: list of (chrom, position, end, ref, alt) tuples for required indels
        inbam: path to source bam
        outbam: path to write realigned bam to
        java: path to java binary
        java_opts: list of optins to pass to the java binary (mainly memory
            requirements)
        gatk: path to JATK jar file
        ref_genome: path to indexed genome fasta file
    '''
    
    # set up the temporary files
    tempdir = tempfile.mkdtemp()
    region_bam = tempfile.NamedTemporaryFile(suffix='.bam', dir=tempdir, delete=False)
    region_vcf = tempfile.NamedTemporaryFile(suffix='.vcf', dir=tempdir)
    targets = tempfile.NamedTemporaryFile(suffix='.intervals', dir=tempdir)
    
    # prepare the input files for indel realignment
    make_indel_bam(inbam, regions, region_bam)
    make_indel_vcf(regions, region_vcf)
    make_targets(regions, targets)
    
    realigner = [java] + java_opts + [
        '-jar', gatk,
        '-T', 'IndelRealigner',
        '-R', ref_genome,
        '-known', region_vcf.name,
        '--targetIntervals', targets.name,
        '-I', region_bam.name,
        '-o',  outbam ]
    
    subprocess.check_call(realigner)
    
    # clean up the temporary files
    for handle in [region_bam, region_vcf, targets]:
        handle.close()
    shutil.rmtree(tempdir)

def main():
    args = get_options()
    realign_indels(args.invcf, args.regions, args.inbam, args.outbam,
        ref_genome=args.genome)

if __name__ == '__main__':
    main()
