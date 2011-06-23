#!/usr/bin/python
# -*- coding: utf-8 -*-

#    Copyright (C) 2011  Mehmet Atakan Gürkan
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3 as
#    published by the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program (probably in a file named COPYING).
#    If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt

import sys
import argparse

from trail_length_calc import trail_length_1d, trail_length_3d

parser = argparse.ArgumentParser(description='Analyzes a given trail')

parser.add_argument('-i', '--input-file',
                    metavar='<input file>',
                    type=argparse.FileType('r'), dest='infile',
                    default=None,
                    help='name(s) of the input file(s) (use \'-\' for stdin)')
# XXX accepting multiple file names is not implemented yet.
# (will use nargs?)
parser.add_argument('-t',
                    type=float, metavar='<float>', default=1.0,
                    help='duration of motion (default: 1.0)')
parser.add_argument('--first-column-time', dest='firstcol',
                    action='store_true', 
                    help='first column in data file is time (overrides \'-t\'; default: time interval is determined by subdividing duration uniformly)')
parser.add_argument('--numpy', dest='inputformat', action='store_const',
                    const='numpy', default='undecided',
                    help='input in NumPy format (default: NumPy)')
parser.add_argument('--ascii', dest='inputformat', action='store_const',
                    const='numpy', default='undecided',
                    help='input in ASCII format (default: NumPy)')


args = parser.parse_args()

def br_pl(a1, a2, m1, m2, m3, b):
    '''A function that returns a function that makes a broken powerlaw.

a1, a2 : x coordinates of the break points.
b : y intersect of the first power law (x<=a1)
m1, m2, m3: slopes for x<a1, a1<x<a2 and a2<x .'''   
    def ifelse(x, y, z) :
        if x: return y
        else: return z
    k1 = a1*(m1-m2) + b
    k2 = a2*(m2-m3) + k1
    return lambda x: ifelse(x<a1, m1*x +b,
                                  ifelse(x<a2, m2*x+k1,
                                               m3*x+k2))

def set_ruler_lengths(rl_min, rl_max, tend) :
    '''A function that creates an array of ruler lengths/sampling intervals.

All values returned are in the closed interval of [rl_min, rl_max]. They
are exact divisors of tend, which is the absolute maximum.'''
    dummy_rl = [tend/1.0,  tend/2.0,  tend/3.0,  tend/4.0,
                tend/5.0,  tend/6.0,  tend/7.0,  tend/8.0,
                tend/10.0, tend/12.0, tend/14.0, tend/17.0,
                tend/20.0, tend/24.0, tend/28.0, tend/33.0,
                tend/40.0, tend/48.0, tend/56.0, tend/67.0]
    for i in range(100) :
        dummy_rl.append(dummy_rl[-4]/2.0)
    rl = []
    for drl in dummy_rl :
        if drl <= rl_max and drl >= rl_min :
            rl.append(drl)
    rl.reverse()
    return np.array(rl)

dn = np.loadtxt(args.infile)
if args.firstcol==True :
    if np.size(np.shape(dn))==2 and np.shape(dn)[1]==4 :
        tt = dn[:,0]
        dd = dn[:,1:]
    elif np.size(np.shape(dn))==2 and np.shape(dn)[1]==2 :
        tt = dn[:,0]
        dd = dn[:,1]
    else :
        print('input file is not 1D or 3D')
        print(np.shape(dn))
        sys.exit()
else :
    tt = np.linspace(0, args.t, np.shape(dn)[0])
    dd = dn


tend = tt[-1]
period = 2.3e-4
rl_min = period/5.0
rl_max = tend/2.0
ruler_lengths = set_ruler_lengths(rl_min, rl_max, tend)

if np.size(np.shape(dd))==2 and np.shape(dd)[1]==3 :
#    print('3d')
    trail_lengths = trail_length_3d(ruler_lengths, tt, dd)
elif np.size(np.shape(dd))==1 :
#    print('1d')
    trail_lengths = trail_length_1d(ruler_lengths, tt, dd)
else :
    print('input file is not 1D or 3D')
    print(np.shape(dd))
    sys.exit()

for i, rl in enumerate(ruler_lengths) :
    print(rl, trail_lengths[i])
