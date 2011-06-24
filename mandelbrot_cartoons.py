#!/usr/bin/python
# -*- coding: utf_8 -*-

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

# this code produces a diagonal path and the
# increments diagram below
# it also prints the path as a time series

from __future__ import division, print_function

from math import sqrt, fabs
from pyx import *
import numpy as np
from numpy.random import random, seed, randint, shuffle
from copy import copy, deepcopy
import scipy.optimize as so

import argparse, sys

class GenSelAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
#        print('%r %r %r' % (namespace, values, option_string))
        setattr(namespace, 'genchoice', option_string[-2:])

parser = argparse.ArgumentParser(description='Creates fractal cartoons (a la Mandelbrot)')

parser.add_argument('--A4', metavar='Hurst parameter',
                    action=GenSelAction,
                    default = 0.5, nargs='?',
                    help='choose Mandelbrot\'s A4 cartoon with Hurst parameter as generator (this is the default choice for generator; if no value for Hurst parameter is given, default is 1/2, i.e, Brownian motion)')
parser.add_argument('--3b', metavar='x',
                    action=GenSelAction,
                    default = 2.0/3, nargs='?',
                    help='choose Gürkan\'s 3-segment Brownian motion with height of the first section as a parameter (Mandelbrot\'s  A4 generator is the default; if no value for x is given, default is 2/3, i.e., first and third section have same height)')
parser.add_argument('--4a',
                    action=GenSelAction,
                    nargs=0,
                    help='choose Gürkan\'s 4-segment anti-correlated Brownian motion generator with exponent 1/3, ie., H=3 (Mandelbrot\'s  A4 generator is the default)')
parser.add_argument('--4b',
                    action=GenSelAction,
                    nargs=0,
                    help='choose Gürkan\'s 4-segment Brownian motion generator, ie., H=2 (Mandelbrot\'s  A4 generator is the default)')
parser.add_argument('--4c',
                    action=GenSelAction,
                    nargs=0,
                    help='choose Gürkan\'s 4-segment correlated Brownian motion generator with exponent 2/3, ie., H=1.5 (Mandelbrot\'s  A4 generator is the default)')
parser.add_argument('--4j',
                    action=GenSelAction,
                    nargs=0,
                    help='choose Gürkan\'s 4-segment "Brownian motion with jumps" generator (Mandelbrot\'s  A4 generator is the default)')
parser.add_argument('--randomize', action='store_true',
                    dest='randomize', default=False,
                   help='Shuffle the intervals randomly (default: do not shuffle')
parser.add_argument('-s', '--seed', dest='s',
                    type=int, default=42,
                    help='random number seed (default: 42)')
parser.add_argument('-i', '--iterations', dest='i',
                    type=int, default=6,
                    help='number of iterations (default: 6)')
parser.add_argument('-c', '--cartoon', dest='cartoon',
                    metavar='cartoon file', default=None,
                    help='file to put a plot in (default: do not produce a plot)')
args = parser.parse_args()

def make_increment(init, gene) :
    '''  takes a vector (similar to an initiator, but not necessarily
        one; any increment vector will do) and a generator, and
        returns a list of increments.
         this is no more than a multiplication, that is, the scaling
        of the generator list with the initiator vector.'''
    return [[x[0]*init[0], x[1]*init[1]] for x in gene]

def iterate_increments(incr_list_in, gene, lev=1) :
    if lev==0 :
        return incr_list_in
    else :
        lev -= 1
        incr_list_out = []
        for incr in incr_list_in :
            if args.randomize == True :
                shuffle(gene) # XXX simpler randomization
            incr_list_out.extend(make_increment(incr, gene))
        return iterate_increments(incr_list_out, gene, lev)

def refine1_incrs(incr_in) :
    ''' this procedure accepts an array of increments,
    adds the increments with same slope to each other,
    and returns the resulting array. "same slope" is
    defined with a tolerance of 1e-4 (builtin)'''
    incr_out = []
    last_slope = 0
    slope_diff_tol = 1e-4
    for inc in incr_in :
        if inc[0] == 0 :
            slope = np.inf
        else :
            slope = inc[1]/inc[0]
        # XXX two tests are needed below since slope may be inf
        if slope == last_slope or fabs(slope - last_slope)<slope_diff_tol :
            incr_out[-1][0] += inc[0]
            incr_out[-1][1] += inc[1]
        else :
            incr_out.append(copy(inc))
        last_slope = slope
#        if incr_out[-1][0] == 0 :
#            last_slope = np.inf
#        else :
#            last_slope = incr_out[-1][1]/incr_out[-1][0]
    return incr_out


seed(args.s)
origin = [0,0]

if args.genchoice=='A4' :
    # Mandelbrot's A4 cartoon
    def make_func(H) : # note H=1/2 is Brownian motion
        return lambda x: (pow(2*x-1.0, 1.0/H) - 1.0 + 2*pow(x, 1.0/H))
    func = make_func(args.A4)
    x =  so.bisect(func, 0.5, 1.0)
    incpow = 1.0/args.A4
    initiator = [1,1]
    generator = [[pow(x, 1.0/args.A4), x],
                 [1-2*pow(x, 1.0/args.A4), 1.0-2*x],
                 [pow(x, 1.0/args.A4), x]]
elif args.genchoice=='3b' :
    # my own general 3-step brownian motion thingy
    incpow = 1.0/2
    x = args['3b']
    y = (-x+1+sqrt((x-1)**2 - 4*x*(x-1)))/2.0
    z = x + y - 1.0
    initiator = [4,2]
    generator = [[x*x, x], [1-(x*x+y*y), 1-(x+y)], [y*y, y]]        
elif args.genchoice=='4a' :
    # fractional gaussian (1/3)
    incpow = 1.0/3
    initiator = [8,2]
    initiator = [1,1]
    generator = [[1/8, 1/2], [1/8, -1/2], [1/8, 1/2], [1/8, 1/2],
                 [1/8,-1/2], [1/8,  1/2], [1/8,-1/2], [1/8, 1/2]] 
elif args.genchoice=='4b' :
    # gaussian (1/2)
    incpow = 1.0/2
    initiator = [4,2]
    initiator = [1,1]
    generator = [[1/4, 1/2], [1/4, -1/2], [1/4, 1/2], [1/4, 1/2]] 
elif args.genchoice=='4c' :
    # fractional gaussian (2/3)
    incpow = 2.0/3
    initiator = [8,4]
    initiator = [1,1]
    generator = [[1/8, 1/4], [1/8,  1/4], [1/8,-1/4], [1/8, 1/4],
             [1/8, 1/4], [1/8, -1/4], [1/8, 1/4], [1/8, 1/4]] 
elif args.genchoice=='4j' :
    # with jumps
    incpow = 1.0/2
    initiator = [1,1]
    generator = [[1/3, 1/2], [1/3, -1/2], [0, 1/2], [1/3, 1/2]]


stg0 = [initiator]
stgf = iterate_increments(stg0, generator, args.i)

refd_incs = np.asarray(refine1_incrs(stgf))

# TIME SERIES OUTPUT
t = 0
d = 0
print("%22.14e %22.14e" % (t, d))
for inc in refd_incs :
    t += inc[0]
    d += inc[1]
    print("%22.14e %22.14e" % (t, d))

if args.cartoon==None : sys.exit(0)

# CARTOON OUTPUT (UPPER PART)
unit.set(uscale=5)
def make_path(origin, increments) :
    '''  takes an origin point and an increment list and returns
        a path with these increments.'''
    pp = [origin]
    ppp = path.path(path.moveto(origin[0], origin[1]))
    for incr in increments :
        pp.append((pp[-1][0] + incr[0], pp[-1][1] + incr[1]))
        ppp.extend(path.path(path.lineto(pp[-1][0],pp[-1][1])))
    return ppp

p0 = make_path(origin, stg0)
pf = make_path(origin, refine1_incrs(stgf))
#pf_alt = make_path(origin, stgf)

#print(stgf)
#print(refine1_incrs(stgf))
#sys.exit(0)

c = canvas.canvas()
c.stroke(p0, [color.cmyk.RoyalBlue,style.linewidth.THick])
c.stroke(pf, [color.cmyk.PineGreen,style.linewidth.Thin])
#c.stroke(pf_alt,[color.cmyk.Apricot,style.linewidth.THIN])


# CARTOON OUTPUT (LOWER PART)

def increment_path(inc_list, x0, y) :
    x = x0
    inc_p = path.path(path.moveto(x, y))
    for inc in inc_list :
        x += inc[0]/2
        inc_p.extend(path.path(path.moveto(x, y)))
        # XXX the "normalization" below is not ideal, i'll work on it
        inc_p.extend(path.path(path.lineto(x, y+inc[1]*((args.i-1)**incpow))))
        x += inc[0]/2
    return inc_p

# increments w/o refinements
#inc_axis = path.path(path.moveto(0, -1))
#inc_axis.extend(path.path(path.lineto(initiator[0], -1)))

#inc_p = increment_path(stgf, 0, -1)
#c.stroke(inc_p, [color.cmyk.Black,style.linewidth.Thin])
#c.stroke(inc_axis, [color.cmyk.Gray,style.linewidth.Thin])


# increments w/ refinements
inc_axis = path.path(path.moveto(0, -1))
inc_axis.extend(path.path(path.lineto(initiator[0], -1)))

inc_p = increment_path(refd_incs, 0, -1)
c.stroke(inc_p, [color.cmyk.Orange,style.linewidth.THIN])
c.stroke(inc_axis, [color.cmyk.Gray,style.linewidth.Thin])

c.writeEPSfile(args.cartoon)
#c.writeEPSfile("brown")
#c.writePDFfile("brown")
