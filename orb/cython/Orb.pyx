# cython: language_level=3str
# cython: embedsignature = False

ctypedef double Real

include "Orb.pxi"
include "common/cython/numbers/double.pyx"
include "common/cython/core/basic_conversions.pyx"
include "core/basic.pyx"
include "core/orb_triangulation.pyx"
include "core/orbifold.pyx"
