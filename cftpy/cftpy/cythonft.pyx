# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

from libc.stdlib cimport malloc, free
from libc.string cimport memset

from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.modular.etaproducts import qexp_eta 


cpdef object cstring_function_su2(int l, int m, int k, int order):
    cdef Py_ssize_t r, s
    cdef Py_ssize_t ind1, ind2
    cdef int sgn
    cdef int *coeff = <int*> malloc((order+1) * sizeof(int))
    if coeff == NULL:
        raise MemoryError()

    memset(coeff, 0, (order+1) * sizeof(int))

    for r in range(order):
        for s in range(order):
            ind1 = ( r*(r+1+l+m) + s*(s+1+l-m) )//2 + r*s*(k+1)
            ind2 = ( r*(r+1-l-m) + s*(s+1-l+m) )//2 + (k+1)*(r+s+1+r*s) - l

            sgn = -1 if ((r + s) & 1) else 1

            if ind1 <= order:
                coeff[ind1] += sgn
            if ind2 <= order:
                coeff[ind2] -= sgn

    cdef list pycoeff = [coeff[i] for i in range(order+1)]
    free(coeff)
    
    cdef object R, series, eta
    R       = PowerSeriesRing(ZZ, 'q', default_prec=order+1)
    series  = R(pycoeff, order+1)
    eta     = qexp_eta(R, order+1)

    return series * eta**(-3)