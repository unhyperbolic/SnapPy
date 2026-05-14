/**
 *  @file double_SnapPy.h
 *  @brief Support for building a kernel that uses standard double precision arithmetic.
 *
 *  Typedefs and structs to enable using standard doubles as the kernel's Real type.
 *
 */

#ifndef _DOUBLE_SNAPPY_
#define _DOUBLE_SNAPPY_

#include "kernel_namespace.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/**
 * Use a standard double as SnapPea's Real type.
 */

typedef double Real;
/**
 * This is used to work around a Cython bug which prevents declaring
 * arrays of C++ objects. See SnapPy.pxi.
 */
typedef double Real_struct;

SNAPPEA_NAMESPACE_END_SCOPE

#define Real_from_string(x) (atof((char *)x))

#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON

#endif
