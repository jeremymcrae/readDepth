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
from getpass import getpass

import pysam

from readDepth.constants import SAMTOOLS

def check_access(pwd_path=None):
    ''' get permission to access IRODS via kerebos authentication
    '''
    
    # don't re-authenticate if we already have access
    klist = subprocess.Popen(['klist' ,'-s'], stderr=subprocess.PIPE,
        stdout=subprocess.PIPE)
    if not klist.stderr.read().startswith('klist: No credentials'):
        return
    
    # you can either store the password in file (but make sure the file isn't
    # readable by any other user), or enter the password at the prompt.
    if pwd_path is not None and os.path.exists(pwd_path):
        pwd = open(pwd_path).read().strip()
    else:
        pwd = getpass('kerebos password required: ')
    
    kinit = subprocess.Popen(['kinit'], stdin=subprocess.PIPE,
        stdout=open(os.devnull, 'w'))
    stdout = kinit.communicate(input=pwd)

def get_bam_slice(bam_path, slice_bam, regions):
    ''' extract a section of a BAM file from IRODS, or lustre.
    
    It can be much quicker extract a small slice of a BAM from IRODS, or to
    iterate through a small slice of a BAM.
    
    Args:
        bam_path: path to BAM on IRODS, or lustre
        temp_bam: file handle for storing the BAM slice, or path to store the
            BAM slice.
        regions: list of (chrom, start, end) tuples for the regions to be
            extracted.
    '''
    
    # format the list of variant regions for samtools
    regions = ['{0}:{1}-{2}'.format(x[0], x[1], x[2]) for x in regions]
    
    # get the slice of the BAM
    temp = tempfile.NamedTemporaryFile(mode='w')
    command = [SAMTOOLS, 'view', '-b', bam_path] + regions
    p = subprocess.Popen(command, stdout=temp, stderr=subprocess.PIPE)
    (outs, errs) = p.communicate()
    print(p.returncode, outs, errs)
    
    if p.returncode == 1:
        if errs != '[main_samview] random alignment retrieval only works for indexed BAM or CRAM files.':
            raise subprocess.CalledProcessError
    
    if type(slice_bam) == str:
        slice_bam = open(slice_bam, 'wb')
    
    slice_bam.write(pysam.sort(temp.name))
    slice_bam.flush()
    pysam.index(slice_bam.name)
    
    # extracting BAMs slices from IRODS with samtools also extracts the
    # corresponding bam index to the working directory. Since we have generated
    # an index file with pysam (in the matching directory to the temporary bam
    # slice), we can delete the extracted BAM index file.
    index = os.path.basename(bam_path.replace('irods://', '')) + '.bai'
    if os.path.exists(index):
        os.remove(index)
