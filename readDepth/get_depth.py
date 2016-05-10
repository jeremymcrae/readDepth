
from __future__ import print_function

import os
import tempfile

from readDepth.get_slice import get_bam_slice, check_access

def get_read_depths(bam_path, variants, min_quality, by_strand=True, keep_bams=False):
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
