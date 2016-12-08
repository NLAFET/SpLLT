#ifndef COMMON_HPP
#define COMMON_HPP

#ifdef USE_COMPLEX
#include <complex>
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#endif // USE_COMPLEX

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "mkl.h"

#ifdef USE_FLOAT
#define REAL_BLAS(name) s##name
#define REAL_LAPACK(name) s##name##_
#define CMPLX_BLAS(name) c##name
#define CMPLX_LAPACK(name) c##name##_
#else // USE_DOUBLE
#define REAL_BLAS(name) d##name
#define REAL_LAPACK(name) d##name##_
#define CMPLX_BLAS(name) z##name
#define CMPLX_LAPACK(name) z##name##_
#endif // USE_FLOAT

#endif // !COMMON_HPP
