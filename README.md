### readDepth

#### Install
```sh
pip install git+https://github.com/jeremymcrae/readDepth.git

# Alternatively:
git clone https://github.com/jeremymcrae/readDepth.git
cd readDepth
python setup.py install --user
```

#### Usage (in python):
```python
from pysam import AlignmentFile
from readDepth.get_depth import allele_depths, get_multiple_depths

bam_path = '/PATH/TO/BAM'
variants = [('1', 1000, 'A', 'G'), ('2', 2000, 'C', 'T')]

bam = AlignmentFile(bam_path)
for chrom, start, end, ref, alt in variants:
    depths = allele_depths(bam, chrom, pos, ref, alt)

# by default this will return depths by strand, but you can turn this off with
depths = allele_depths(bam, chrom, pos, ref, alt, by_strand=False)

# aqnd you can get multiple variants at once with
depths = get_multiple_depths(bam, variants)
```
