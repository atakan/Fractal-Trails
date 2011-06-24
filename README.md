This is a collection of programs to create and analyze trails of random
walk (in 1D or 3D, optionally with bias and/or bounds), presented as a
time series. The analysis is done by estimating the fractal dimension
of the trails and yields the timescales of correlations as well as
their strength.

I used the techniques for creating and analyzing random walk trails
in a paper that is available in the `doc/` directory.
[This paper](http://arxiv.org/abs/1011.5226) has the formulae for
estimating the fractal dimension. The code in this repository is a
clean-up/rewrite of the code used for that paper.

There is also a little code, called `mandelbrot_cartoons.py`, to create
some fractal trails that Mandelbrot calls cartoons (see, e.g., Chapter
1 in "Multifractals and 1/f Noise: Wild Self-affinity in Physics." New
York NY: Springer, 1999). This is a work in progress, but already
has some material to build various generalizations of Brownian motion.

**trail_maker.py** uses [NumPy](http://numpy.scipy.org/) and
[argparse](http://code.google.com/p/argparse/).

**trail_analyzer** uses [NumPy](http://numpy.scipy.org/),
[argparse](http://code.google.com/p/argparse/),
[Cython](http://cython.org/) and [GSL](http://www.gnu.org/software/gsl/).
To use it, first create `trail_length_calc.so` by running

`python setup.py build_ext --inplace`

**mandelbrot_cartoons.py** uses [NumPy](http://numpy.scipy.org/),
[SciPy](http://www.scipy.org/) optimization routines, 
[argparse](http://code.google.com/p/argparse/) and
[PyX](http://pyx.sourceforge.net/)

Keywords: Fractals -- fractal dimension -- Brownian motion --
random walk

Atakan Gürkan <ato.gurkan@gmail.com>
