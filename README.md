### readDepth: de novo mutation recurrence significance testing
program to calculate the significance of seeing N DNMs of a specific
combination of functional types in a particular gene in M trios

#### Install
```sh
pip install git+git://github.com/jeremymcrae/readDepth.git --user

# Alternatively:
git clone https://github.com/jeremymcrae/readDepth.git
cd readDepth
python setup.py install --user
```

#### Usage (in python):
```python
from readDepth.get_depth import get_read_depths

bam_path = '/PATH/TO/BAM'
variants = [('1', 1000, 1000, 'A', 'G'), ('2', 2000, 2000, 'C', 'T')]

depths = get_read_depths(bam_path, variants)

# by default this will return depths by strand, but you can turn this off with
depths = get_read_depths(bam_path, variants, by_strand=False)
```
