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

Quick script to check what sort of depth we would need to distingush between
error rates of 0.001 and 0.01. Assumes we can aggregate depths across unaffected
parents.
'''

from __future__ import division, print_function

import math

from scipy.stats import fisher_exact, chisquare, norm
import numpy

import scipy
if scipy.__version__ < '0.12.0':
    raise ValueError('needs scipy > 0.9.0, due to a bug in fisher_exact')

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot
import seaborn

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def check_difference(low, high, x, z, p_values):
    ''' calculate p-values for observing a difference in error rates.
    
    Args:
        low: low sequencing error rate (e.g. 0.001)
        high: high sequencing error rate (e.g. 0.01)
        x: number of aggregate reads checked
        z: quantile of a standard normal distribution for the confidence interval
        p_values: dictionary of p-values, with lists for 'lower', 'mid', and 'upper'
    '''
    
    low_delta = z * math.sqrt((low * (1 - low))/x)
    high_delta = z * math.sqrt((high * (1 - high))/x)
    
    # calculate the number of alts at the low and high error rates
    low_alts = low * x
    high_alts = high * x
    
    # calculate the alt counts at the lower confidence interval
    lo_low_alts = (low + low_delta) * x
    lo_high_alts = (high - high_delta) * x
    
    # calculate the alt counts at the upper confidence interval
    hi_low_alts = (low - low_delta) * x
    hi_high_alts = (high + high_delta) * x
    
    # create lists of alt and ref counts
    mid = [[low_alts, (x - low_alts)], [high_alts, (x - high_alts)]]
    lo = [[lo_low_alts, (x - lo_low_alts)], [lo_high_alts, (x - lo_high_alts)]]
    hi = [[hi_low_alts, (x - hi_low_alts)], [hi_high_alts, (x - hi_high_alts)]]
    
    values = {'lower': lo, 'mid': mid, 'upper': hi}
    
    for key in values:
        data = values[key]
        try:
            oddsratio, fisher = fisher_exact(data)
        except ValueError:
            fisher = 1
        p_values[key].append(abs(fisher))

def plot_depth_vs_p(sizes, mid, lower, upper, path='aggregate_depth_vs_p_value.pdf'):
    '''plot p-value by aggregate depth
    '''
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    e = ax.plot(sizes, -numpy.log10(mid), color='gray')
    e = ax.fill_between(sizes, -numpy.log10(lower), -numpy.log10(upper), alpha=0.5)
    
    e = ax.set_xlabel('aggregate depth (n)')
    e = ax.set_ylabel('-log10(P)')
    
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['right'].set_visible(False)
    
    e = ax.spines['bottom'].set_position(('outward', 10))
    e = ax.spines['left'].set_position(('outward', 10))
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    fig.savefig(path, format='pdf', transparent=True, bbox_inches='tight',
        pad_inches=0.0)

def main():
    low = 0.001
    high = 0.01
    ci = 0.95
    
    increment = 20
    z = norm.ppf(1 - ((1 - ci)/2))
    
    sizes = range(increment, 5000, increment)
    p_values = {'lower': [], 'mid': [], 'upper': []}
    for x in sizes:
        print(x)
        check_difference(low, high, x, z, p_values)
    
    plot_depth_vs_p(sizes, p_values['mid'], p_values['lower'], p_values['upper'])

if __name__ == '__main__':
    main()
