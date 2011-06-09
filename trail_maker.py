#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import numpy as np
from numpy.random import random, seed
import argparse

parser = argparse.ArgumentParser(description='Creates a random walk trail')

parser.add_argument('-t',
                    type=float, default=1.0,
                    help='duration of the motion')
parser.add_argument('-i',
                    type=float, default=1.0,
                    help='size of increments (will get normalized)')
parser.add_argument('-d',
                    type=int, default=3,
                    help='dimensionality of the trail')
parser.add_argument('-N',
                    type=int, default=1000,
                    help='number of points generated')
parser.add_argument('-r',
                    type=int, default=0,
                    help='repository size')
parser.add_argument('-c',
                    type=int, default=0.0,
                    help='bias strength')
parser.add_argument('-b',
                    type=float, default=3.0e40,
                    help='radius of boundary')
parser.add_argument('-s',
                    type=int, default=42,
                    help='random number seed')

args = parser.parse_args()

# Initial value P_0 is always 1 or (1,1,1)/sqrt(3.0).
# Increment \delta P is (given increment)/sqrt(N).

seed(args.s)

Dt = args.t
N = args.N
d = args.d
b = args.b
dP = args.i/np.sqrt(N)

def trail_1d(N, r=0, c=0.0) :
    P0 = 1.0
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
        q = 0.0
    for i in range(N) :
        X = random() - 0.5
        if X>q : P -= dP # this choice brings positive correlation
        else   : P += dP # for positive c
        trl[i+1] = P
        if use_rep :
            rep[i%r] = X
            q = np.sum(rep)*rep_norm * c
    return trl

def trail_3d(N, r=0, c=0.0, b=3.0e40) :
    def vec_norm2(a) :
        return a[0]*a[0] + a[1]*a[1] + a[2]*a[2]
    b2 = b*b
    P0x, P0y, P0z = 1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0)
    Px, Py, Pz = P0x, P0y, P0z
    trlx, trly, trlz  = np.empty(N+1), np.empty(N+1), np.empty(N+1),
    trlx[0], trly[0], trlz[0] = P0x, P0y, P0z
    if r>0 and c!=0.0 : # repository initialization
        use_rep = True
        rep_norm = 1.0/np.sqrt(r)
        repx, repy, repz = random(r)-0.5, random(r)-0.5, random(r)-0.5
        qx = np.sum(repx)*rep_norm * c
        qy = np.sum(repy)*rep_norm * c
        qz = np.sum(repz)*rep_norm * c
    else :
        use_rep = False
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
            print('# XX', i, Px, Py, Pz)
            Px -= DPx           # so we take the opposite step
            Py -= DPy
            Pz -= DPz
        else :                  # we are safe
            Px += DPx           # so we take normal step
            Py += DPy
            Pz += DPz
        trlx[i+1], trly[i+1], trlz[i+1] = Px, Py, Pz
        if use_rep :
            repx[i%r], repy[i%r], repz[i%r]= Xx, Xy, Xz
            qx = np.sum(repx)*rep_norm * c
            qy = np.sum(repy)*rep_norm * c
            qz = np.sum(repz)*rep_norm * c
    return trlx, trly, trlz

ts = np.linspace(0.0, Dt, N+1)
if d==1 :
    trl = trail_1d(N)
    for i, t in enumerate(ts) :
        print('%e %e' % (t, trl[i]))
elif d==3 :
    trlx, trly, trlz = trail_3d(N,0,0.0,b)
    for i, t in enumerate(ts) :
        print('%e %e %e %e' % (t, trlx[i], trly[i], trlz[i]))

