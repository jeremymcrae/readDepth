
import os
import subprocess
import shutil
from getpass import getpass

import pysam

from readDepth.constants import SAMTOOLS

def check_access(pwd_path=None):
    ''' get permission to access IRODS via kerebos authentication
    '''
    
    # don't re-authenticate if we already have access
    klist = subprocess.Popen(["klist"], stderr=subprocess.PIPE,
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
    """ extract a section of a BAM file from IRODS, or lustre
    
    Args:
        bam_path: path to BAM on IRODS, or lustre
        temp_bam: file handle for storing the BAM slice, or path to store the
            BAM slice.
        regions: list of (chrom, start, end) tuples for the regions to be
            extracted.
    """
    
    # format the list of variant regions for samtools
    regions = ["{0}:{1}-{2}".format(x[0], x[1], x[2]) for x in regions]
    
    # get the slice of the BAM
    temp = tempfile.NamedTemporaryFile(mode="w")
    command = [SAMTOOLS, "view", "-b", bam_path] + regions
    code = subprocess.check_call(command, stdout=temp, stderr=open(os.devnull, "w"))
    
    if type(slice_bam) == str:
        slice_bam = open(slice_bam, 'wb')
    
    slice_bam.write(pysam.sort(temp.name))
    slice_bam.flush()
    pysam.index(slice_bam.name)
    
    # extracting BAMs slices from IRODS with samtools also extracts the
    # corresponding bam index to the working directory. Since we have generated
    # an index file with pysam (in the matching directory to the temporary bam
    # slice), we can delete the extracted BAM index file.
    index = os.path.basename(bam_path.replace("irods://", "")) + ".bai"
    if os.path.exists(index):
        os.remove(index)
