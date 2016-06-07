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


Script to estimate site-specific error rates from counts of reference and
alternate alleles aggregated across individuals.
'''

import argparse
import sys

import numpy
import pandas

import matplotlib
matplotlib.use('Agg')

from scipy.stats import gaussian_kde
from matplotlib import pyplot
import seaborn

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def get_options():
    """ get the command line options
    """
    
    parser = argparse.ArgumentParser(description='Calculate site-specific error '
        'rates.')
    parser.add_argument('--input', default=sys.stdin,
        help='Path to table of read counts per individual.')
    parser.add_argument('--output', default=sys.stdout,
        help='Path to file for output. Defaults to standard out.')
    
    args = parser.parse_args()
    
    return args

def get_error_rates(depths):
    ''' calculate site-specific error rates from aggregated read depths
    '''
    
    columns = ['chrom', 'start_pos', 'ref', 'alt', 'rate']
    error_rates = pandas.DataFrame(columns=columns)
    for key, group in depths.groupby(['chrom', 'start_pos']):
        ref = group[['ref_F', 'ref_R']].sum().sum()
        alt = group[['alt_F', 'alt_R']].sum().sum()
        
        # if the rate is zero, in order to get a good upper estimate of the error
        # rate, assume the next read would be an alt
        if alt == 0:
            alt += 1
        
        temp = pandas.DataFrame({'chrom': [key[0]], 'start_pos': [key[1]],
            'ref': [ref], 'alt': [alt], 'rate': [float(alt)/(alt + ref)]})
        
        error_rates = error_rates.append(temp, ignore_index=True)
    
    return error_rates[columns]

def plot_density(values, xlabel='values', path='plot.pdf'):
    ''' plot the density of variants at different error rates.
    '''
    
    density = gaussian_kde(values)
    
    x_values = numpy.arange(min(values), max(values), 0.01)
    y_values = density(x_values)
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    e = ax.plot(x_values, y_values)
    
    e = ax.set_xlabel(xlabel)
    e = ax.set_ylabel('Density')
    
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['right'].set_visible(False)
    
    e = ax.spines['bottom'].set_position(('outward', 10))
    e = ax.spines['left'].set_position(('outward', 10))
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    fig.savefig(path, format='pdf', transparent=True, bbox_inches='tight',
        pad_inches=0.0)

def main():
    args = get_options()
    depths = pandas.read_table(args.input, compression='gzip')
    
    error_rates = get_error_rates(depths)
    error_rates.to_csv(args.output, sep='\t', index=False)
    
    rates = numpy.log10(list(error_rates['rate']))
    plot_density(rates, xlabel='log10(error rate)', path='error_density.pdf')
    
    depth = error_rates[['alt', 'ref']].sum(axis=1)
    depth = numpy.log10(depth)
    plot_density(depth, xlabel='log10(aggregate depth)', path='depth_density.pdf')

if __name__ == '__main__':
    main()
