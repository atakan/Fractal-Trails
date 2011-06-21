#!/usr/bin/python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
setup(
  name = "CythonTrailLength",
  ext_modules=[ 
    Extension("trail_length_calc", ["trail_length_calc.pyx"], libraries = ["gsl", "gslcblas", "m"])
    ],
  cmdclass = {'build_ext': build_ext}
)
