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
from numpy.random import random, seed
import argparse, sys

parser = argparse.ArgumentParser(description='Creates a random walk trail')

parser.add_argument('-d',
                    type=int, default=3,
                    help='topological dimension of the trail (default: 3)')
parser.add_argument('-N',
                    type=int, default=1000,
                    help='number of steps in the trail (default: 1000)')
parser.add_argument('-i',
                    type=float, default=1.0,
                    help='size of increments (will get normalized by 1/sqrt(N))')
parser.add_argument('-r',
                    type=int, default=0,
                    help='repository size (default: 0)')
parser.add_argument('-c',
                    type=float, default=0.0,
                    help='bias strength')
parser.add_argument('-b',
                    type=float, default=3.0e40,
                    help='radius of boundary')
parser.add_argument('-s',
                    type=int, default=42,
                    help='random number seed (default: 42)')
parser.add_argument('--rangen',
                    type=int, default=-1,
                    help='generate this many random numbers before starting to build the trail (default: repo size; ignored if less than repo size)')
parser.add_argument('-P0',
                    type=float, default=0.0, dest='P0',
                    help='initial value (for 1D)')
parser.add_argument('-P0x',
                    type=float, default=0.0, dest='P0x',
                    help='initial value of x component (for 3D)')
parser.add_argument('-P0y',
                    type=float, default=0.0, dest='P0y',
                    help='initial value of y component (for 3D)')
parser.add_argument('-P0z',
                    type=float, default=0.0, dest='P0z',
                    help='initial value of z component (for 3D)')
parser.add_argument('-o','--output-file',
                    dest='outfile',
                    type=argparse.FileType('w'),
                    default=sys.stdout,
                    help='output filename (if not given, use stdout)')
parser.add_argument('--numpy', dest='outputformat', action='store_const',
                   const='numpy', default='undecided',
                   help='output in NumPy format (default: ASCII for stdout, NumPy for file')
parser.add_argument('--ascii', dest='outputformat', action='store_const',
                   const='ascii', default='undecided',
                   help='output in ASCII format (default: ASCII for stdout, NumPy for file')

args = parser.parse_args()

seed(args.s)
if args.rangen > args.r :
    dummy = random(args.rangen - args.r)

N = args.N
d = args.d
b = args.b
dP = args.i/np.sqrt(N)

def trail_1d(N, r=0, c=0.0, b=3.0e40) :
    P0 = args.P0
    P = P0
    trl = np.empty(N+1)
    trl[0] = P0
    if r>0 and c!=0.0 : # repository initialization
        use_rep = True
        rep_norm = 1.0/np.sqrt(r)
        rep = random(r)-0.5
        q = np.sum(rep)*rep_norm * c
    else :
        use_rep = False
        dummy = random(r)
        q = 0.0
    for i in range(N) :
        X = random() - 0.5
        if X>q : DP = -dP 
        else   : DP = dP  
        if np.fabs(P+DP) > b :
            P -= DP
        else :
            P += DP
        trl[i+1] = P
        if use_rep :
            rep[i%r] = X
            q = -1.0*np.sum(rep)*rep_norm * c
    return trl

def trail_3d(N, r=0, c=0.0, b=3.0e40) :
    def vec_norm2(a) :
        return a[0]*a[0] + a[1]*a[1] + a[2]*a[2]
    b2 = b*b
    P0x, P0y, P0z =  args.P0x, args.P0y, args.P0z
    Px, Py, Pz = P0x, P0y, P0z
    trl  = np.empty((N+1,3))
    trl[0] = P0x, P0y, P0z
    if r>0 and c!=0.0 : # repository initialization
        use_rep = True
        rep_norm = 1.0/np.sqrt(r)
        repx, repy, repz = random(r)-0.5, random(r)-0.5, random(r)-0.5
        qx = np.sum(repx)*rep_norm * c
        qy = np.sum(repy)*rep_norm * c
        qz = np.sum(repz)*rep_norm * c
    else :
        use_rep = False
        dummy = random(r*3)
        qx, qy, qz = 0.0, 0.0, 0.0
    for i in range(N) :
        Xx = random() - 0.5
        Xy = random() - 0.5
        Xz = random() - 0.5
        if Xx>qx : DPx = -dP 
        else     : DPx = dP 
        if Xy>qy : DPy = -dP 
        else     : DPy = dP 
        if Xz>qz : DPz = -dP 
        else     : DPz = dP 
        Ptry = Px+DPx, Py+DPy, Pz+DPz
        if vec_norm2(Ptry) > b2 : # we'll cross bndry if we take this step
            Px -= DPx           # so we take the opposite step
            Py -= DPy
            Pz -= DPz
        else :                  # we are safe
            Px += DPx           # so we take normal step
            Py += DPy
            Pz += DPz
        trl[i+1] = (Px, Py, Pz)
        if use_rep :
            repx[i%r], repy[i%r], repz[i%r]= Xx, Xy, Xz
            qx = -1.0*np.sum(repx)*rep_norm * c
            qy = -1.0*np.sum(repy)*rep_norm * c
            qz = -1.0*np.sum(repz)*rep_norm * c
    return trl

if args.outputformat == 'undecided' :
    if args.outfile == sys.stdout :
        outputformat = 'ascii'
    else :
        outputformat = 'numpy'
else :
    outputformat = args.outputformat

if d==1 :
    trl = trail_1d(N, r=args.r, c=args.c, b=args.b)
    if outputformat == 'ascii' :
        for p in trl :
            print('%e' % (p), file=args.outfile)
    else :
        np.save(args.outfile, trl)
elif d==3 :
    trl = trail_3d(N, r=args.r, c=args.c, b=args.b)
    if outputformat == 'ascii' :
        for p in trl :
            print('%e %e %e' % (p[0], p[1], p[2]), file=args.outfile)
    else :
        np.save(args.outfile, trl)
else :
    print('illegal dimension given: %d' %(d))
