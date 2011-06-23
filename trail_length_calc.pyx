from __future__ import print_function
import numpy

cdef extern from "stdlib.h":
    ctypedef long size_t

cdef extern from "gsl/gsl_interp.h":
    ctypedef struct gsl_interp_accel:
        size_t  cache
        size_t  miss_count
        size_t  hit_count

    ctypedef struct gsl_interp_type:
        char * name
        unsigned int min_size
        void *  (*alloc) (size_t size)
        int     (*init)  (void *, double xa[], double ya[], size_t size)
        int     (*eval)  (void *, double xa[], double ya[], size_t size, double x, gsl_interp_accel *, double * y)
        int     (*eval_deriv)  (void *, double xa[], double ya[], size_t size, double x, gsl_interp_accel *, double * y_p)
        int     (*eval_deriv2) (void *, double xa[], double ya[], size_t size, double x, gsl_interp_accel *, double * y_pp)
        int     (*eval_integ)  (void *, double xa[], double ya[], size_t size, gsl_interp_accel *, double a, double b, double * result)
        void    (*free)  (void *)

    ctypedef struct gsl_interp:
        gsl_interp_type * type
        double  xmin
        double  xmax
        size_t  size
        void * state

    gsl_interp_accel * gsl_interp_accel_alloc()

    gsl_interp *gsl_interp_alloc(gsl_interp_type * T, size_t n)
    
    int gsl_interp_init(gsl_interp * obj, double xa[], double ya[], size_t size)

    double gsl_interp_eval(gsl_interp * obj,
                double xa[], double ya[], double x,
                gsl_interp_accel * a)

    void gsl_interp_free(gsl_interp * interp)
    void gsl_interp_accel_free(gsl_interp_accel * a)

    extern gsl_interp_type * gsl_interp_linear
    extern gsl_interp_type * gsl_interp_polynomial
    extern gsl_interp_type * gsl_interp_cspline
    extern gsl_interp_type * gsl_interp_cspline_periodic
    extern gsl_interp_type * gsl_interp_akima
    extern gsl_interp_type * gsl_interp_akima_periodic

cdef extern from "gsl/gsl_spline.h":

    ctypedef struct gsl_spline:
        gsl_interp    *interp
        double        *x
        double        *y
        size_t        size

    gsl_spline *gsl_spline_alloc(gsl_interp_type *T, size_t size)

    int gsl_spline_init(gsl_spline *spline, double xa[], double ya[], size_t size)

    double gsl_spline_eval(gsl_spline *spline, double x, gsl_interp_accel *a)

    void gsl_spline_free(gsl_spline *spline)

cdef extern from "numpy/arrayobject.h":

    ctypedef int intp 

    ctypedef extern class numpy.dtype [object PyArray_Descr]:
        cdef int type_num, elsize, alignment
        cdef char type, kind, byteorder, hasobject
        cdef object fields, typeobj

    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef object base
        cdef dtype descr
        cdef int flags

    void import_array()

import_array()

def get_distance(x, y):
    sep = x-y
    return numpy.sqrt(numpy.sum(sep*sep))

def get_length(x) :
    ''' get the "length" of the array x,
        sum(i: 1 to length of x): distance(x[i],x[i-1])''' 
    length = 0.0
    p = x[0]
    for np in x[1:] :
        length += get_distance(np, p)
        p = np
    return length

def trail_length_1d(ndarray rls, ndarray tt, ndarray dd) :
    cdef double *t_in, *x_in
    cdef int npoi
    cdef ndarray T, X
    
    # we create NumPy arrays from time and trail data, and create C pointers
    T = numpy.asfortranarray(tt)
    X = numpy.asfortranarray(dd)
    t_in = <double *>T.data
    x_in = <double *>X.data
    npoi = T.dimensions[0]

    cdef gsl_interp_accel *accx
    cdef gsl_spline *splinex

    # GSL spline for time and trail data
    accx = gsl_interp_accel_alloc()
    splinex = gsl_spline_alloc(gsl_interp_cspline, npoi);
    gsl_spline_init(splinex, t_in, x_in, npoi);
    
    # For each ruler length, we create a new trail using this spline;
    # then calculate the length of this trail.
    tls = [] # trail lengths
#    print("# tt", list(tt))
#    print("# T", list(T))
    for rl in rls :
        #new_t = numpy.arange(T[0], T[-1], rl)
        new_t = numpy.linspace(T[0], T[-1], round((T[-1]-T[0])/rl))
        new_x = [gsl_spline_eval(splinex, t, accx) for t in new_t]
#        print("# T[0], T[-1], rl", T[0], T[-1], rl)
#        print("# t = ", list(new_t))
#        print("# x = ", list(new_x))
        tls.append(get_length((new_x)))
    gsl_spline_free(splinex)
    gsl_interp_accel_free(accx)
    return numpy.array(tls)

def trail_length_3d(ndarray rls, ndarray tt, ndarray dd) :
    cdef double *t_in, *x_in, *y_in, *z_in
    cdef int npoi
    cdef ndarray T, X, Y, Z

    T = numpy.asfortranarray(tt)
    X = numpy.asfortranarray(dd[:,0])
    Y = numpy.asfortranarray(dd[:,1])
    Z = numpy.asfortranarray(dd[:,2])
    t_in = <double *>T.data
    x_in = <double *>X.data
    y_in = <double *>Y.data
    z_in = <double *>Z.data
    npoi = T.dimensions[0]

    cdef gsl_interp_accel *accx, *accy, *accz
    cdef gsl_spline *splinex, *spliney, *splinez

    accx = gsl_interp_accel_alloc()
    splinex = gsl_spline_alloc(gsl_interp_cspline, npoi);
    gsl_spline_init(splinex, t_in, x_in, npoi);
    accy = gsl_interp_accel_alloc()
    spliney = gsl_spline_alloc(gsl_interp_cspline, npoi);
    gsl_spline_init(spliney, t_in, y_in, npoi);
    accz = gsl_interp_accel_alloc()
    splinez = gsl_spline_alloc(gsl_interp_cspline, npoi);
    gsl_spline_init(splinez, t_in, z_in, npoi);
    
    tls = [] # trail lengths
    for rl in rls :
        new_t = numpy.arange(T[0], T[-1], rl)
        new_x = [gsl_spline_eval(splinex, t, accx) for t in new_t]
        new_y = [gsl_spline_eval(spliney, t, accy) for t in new_t]
        new_z = [gsl_spline_eval(splinez, t, accz) for t in new_t]
        tls.append(get_length(numpy.column_stack((new_x, new_y, new_z))))
    gsl_spline_free(splinex)
    gsl_spline_free(spliney)
    gsl_spline_free(splinez)
    gsl_interp_accel_free(accx)
    gsl_interp_accel_free(accy)
    gsl_interp_accel_free(accz)
    return numpy.array(tls)
