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

import subprocess

def get_full_bam(bam, outpath):
    ''' extracts BAM and BAI files for a participant from IRODS
    
    Args:
        bam: location of BAM file on IRODS
        outpath: path to extract the BAM file to (lustre recommended).
        attempts: keep track of how many times this function is recursively
            called, so we can quit after too many attempts
    '''
    
    # make sure we strip the irods:// prefix from the bam path. This should keep
    # a single '/' at the start of the path
    start = 'irods://'
    if bam.startswith(start):
        bam = bam[len(start):]
    
    bai = bam + '.bai'
    
    # if we have passed in a file handle for the destination, get the
    # corresponding path.
    try:
        outpath = outpath.name
    except AttributeError:
        pass
    
    outbai = outpath + '.bai'
    
    # and pull the BAM and BAI files out of IRODs, so that we can work on them
    subprocess.check_call(['iget', '-f', '-v', '-X', outpath, '-K', '--retries', '3', bam, outpath],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.check_call(['iget', '-f', '-v', '-X', outbai, '-K', '--retries', '3', bai, outbai],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
