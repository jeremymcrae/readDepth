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

import os
import tempfile

from readDepth.get_slice import get_bam_slice, check_access

def get_read_depths(bam_path, variants, min_quality=0, by_strand=True, keep_bams=False):
    """ get the read depths of the alleles at variable sites
    
    Args:
        bam_path: path to BAM file (can be on a standard filesystem (e.g.
        /lustre/scratch114/projects/ddd/release-main_fy3/20140916/sample_improved_bams_hgi_2/DDD_MAIN5605962.bam,
        or on an IRODS filesystem e.g. irods:///humgen/projects/ddd/20150316/1866STDY5139786.bam)
        variants: list of (chromosome, start position, end position) tuples.
    
    Returns:
        list of dictionaries of allele depths by base.
    """
    
    pwd_path = os.path.expanduser('~/.kinit')
    check_access(pwd_path)
    
    temp_bam = tempfile.NamedTemporaryFile(mode="wb", delete=not(keep_bams))
    get_bam_slice(bam_path, temp_bam, variants)
    
    counts = {}
    counter = AlleleCounter(temp_bam.name, min_qual=min_quality, by_strand=by_strand)
    for variant in variants:
        try:
            counts[variant] = counter.check_variant(*variant)
        except AssertionError:
            counts[variant] = {"ref": None, "alt": None}
            if by_strand:
                counts[variant] = {"ref": {"forward": None, "reverse": None},
                    "alt": {"forward": None, "reverse": None}}
    
    if not keep_bams:
        # clean up the bam index (the BAM file will be automatically removed as it
        # is a tmpfile).
        if os.path.exists("{}.bai".format(temp_bam.name)):
            os.remove("{}.bai".format(temp_bam.name))
    
    return counts, temp_bam.name
