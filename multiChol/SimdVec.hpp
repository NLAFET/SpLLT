#pragma once

#include <cmath>
#include <cstdio>
#include <limits>

#include <immintrin.h>

/** \brief The SimdVec class isolates use of AVX/whatever intrinsics in a
 *  single place for ease of upgrading to future instruction sets.
 *
 *  Support is only added as required, so don't expect all intrinsics to be
 *  wrapped yet! */
template <typename T, //< Underlying type e.g. int or double
#if defined(__AVX2__) || defined(__AVX__)
          int length = 256 //< Underlying number of bits
#else
          int length = 8*sizeof(T) //< Underlying number of bits
#endif
          >
class SimdVec;

template <>
class SimdVec<int,
#if defined(__AVX2__) || defined(__AVX__)
              128
#else
              8*sizeof(int) // 8 bits to a byte
#endif
             > {
public:
#if defined(__AVX2__) || defined(__AVX__)
   /// Length of underlying vector type
   static const int vector_length = 4;
   /// Typedef for underlying vector type containing doubles
   typedef __m128i simd_int_type;
#else
   /// Length of underlying vector type
   static const int vector_length = 1;
   /// Typedef for underlying vector type containing doubles
   typedef int simd_int_type;
#endif

   /*******************************************
    * Constructors
    *******************************************/

   /// Uninitialized value constructor
   SimdVec()
   {}

   /*******************************************
    * Memory load/store
    *******************************************/

   /// Extract value as array
   void store_aligned(int *dest) const {
#if defined(__AVX2__) || defined(__AVX__)
      _mm_store_si128((__m128i*) dest, val);
#else
      dest[0] = val;
#endif
   }

   /*******************************************
    * Named operations
    *******************************************/

   void set_element(const int idx, int v) {
      val = _mm_insert_epi32( val, v, idx );
   }
   /** Returns element idx of vector value.
    *  idx may be runtime variable */
   /* Messy store to support runtime index.
    * If index is known use extract_element()  */
   int get_element(const int idx) const {
      int __attribute__((aligned(32))) val_as_array[vector_length];
      store_aligned(val_as_array);
      return val_as_array[idx];
   }
   /** Return element idx of vector.
    *  IMPORTANT: idx must be compile time constant */
   int extract_element(const int idx) const {
      return _mm_extract_epi32( val, idx );
   }

   /*******************************************
    * Operators
    *******************************************/

   /// Extract indvidual elements of vector
   /// idx MUST be < vector_length.
   int operator[](const int idx) const {
      return get_element(idx);
   }

   /*******************************************
    * Debug functions
    *******************************************/

   /// Prints the vector (inefficient, use for debug only)
   void print() const {
      for(int i=0; i<vector_length; i++) printf(" %d", (*this)[i]);
   }

private:
   template <typename T_, int length_>
   friend class SimdVec;

   /// Underlying vector that this type wraps
   simd_int_type val;
};

template <>
class SimdVec<double,
#if defined(__AVX2__) || defined(__AVX__)
              256
#else
              8*sizeof(double)
#endif
             > {
public:
   /*******************************************
    * Properties of the type
    *******************************************/

#if defined(__AVX2__) || defined(__AVX__)
   /// Length of underlying vector type
   static const int vector_length = 4;
   /// Typedef for underlying vector type containing doubles
   typedef __m256d simd_double_type;
#else
   /// Length of underlying vector type
   static const int vector_length = 1;
   /// Typedef for underlying vector type containing doubles
   typedef double simd_double_type;
#endif
   typedef SimdVec<int, vector_length*8*sizeof(int)> simd_index_type;

   /*******************************************
    * Constructors
    *******************************************/

   /// Uninitialized value constructor
   SimdVec()
   {}
   /// Initialize all entries in vector to given scalar value
   SimdVec(const double initial_value)
   {
#if defined(__AVX2__) || defined(__AVX__)
      val = _mm256_set1_pd(initial_value);
#else
      val = initial_value;
#endif
   }
#if defined(__AVX2__) || defined(__AVX__)
   /// Initialize with underlying vector type (no version for scalar)
   SimdVec(const simd_double_type &initial_value) {
      val = initial_value;
   }
#endif
   /// Initialize with another SimdVec
   SimdVec(const SimdVec<double> &initial_value) {
      val = initial_value.val;
   }
#if defined(__AVX2__) || defined(__AVX__)
   /// Initialize as a vector by specifying all entries (no version for scalar)
   SimdVec(double x1, double x2, double x3, double x4) {
      val = _mm256_set_pd(x4, x3, x2, x1); // Reversed order expected
   }
#endif

   /*******************************************
    * Memory load/store
    *******************************************/

   /// Load from suitably aligned memory
   static
   const SimdVec load_aligned(const double *src) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec( _mm256_load_pd(src) );
#else
      return SimdVec( src[0] );
#endif
   }

   /// Load from unaligned memory
   static
   const SimdVec load_unaligned(const double *src) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec( _mm256_loadu_pd(src) );
#else
      return SimdVec( src[0] );
#endif
   }

   /// Gather
   static
   const SimdVec gather(const double *base, const simd_index_type idx, const int scale) {
#if defined(__AVX2__)
      //return SimdVec( _mm256_i32gather_pd(base, idx.val, scale) );
      // FIXME: Figure out WHY the gather intrinisic is 50% slower than the
      //        below (prefetch prediction better?)
      return SimdVec(
            base[idx[0]*scale],
            base[idx[1]*scale],
            base[idx[2]*scale],
            base[idx[3]*scale]
            );
#elif defined(__AVX__)
      return SimdVec(
            base[idx[0]*scale],
            base[idx[1]*scale],
            base[idx[2]*scale],
            base[idx[3]*scale]
            );
#else
      return SimdVec( base[idx.val * scale] );
#endif
   }

   /// Scatter
   void scatter(double *base, const simd_index_type idx, const int scale) const {
#if defined(__AVX2__) || defined(__AVX__)
      __m128d low0  = _mm256_extractf128_pd(val, 0);
      __m128d high0 = _mm256_extractf128_pd(val, 1);
      __m128d low1 = _mm_shuffle_pd(low0, low0, _MM_SHUFFLE2(0,1));
      __m128d high1 = _mm_shuffle_pd(high0, high0, _MM_SHUFFLE2(0,1));
      base[idx.extract_element(0)*scale] = _mm_cvtsd_f64(low0);
      base[idx.extract_element(1)*scale] = _mm_cvtsd_f64(low1);
      base[idx.extract_element(2)*scale] = _mm_cvtsd_f64(high0);
      base[idx.extract_element(3)*scale] = _mm_cvtsd_f64(high1);
#else
      base[idx.val * scale] = val;
#endif
   }

   /// Extract value as array
   void store_aligned(double *dest) const {
#if defined(__AVX2__) || defined(__AVX__)
      _mm256_store_pd(dest, val);
#else
      dest[0] = val;
#endif
   }

   /// Extract value as array
   void store_unaligned(double *dest) const {
#if defined(__AVX2__) || defined(__AVX__)
      _mm256_storeu_pd(dest, val);
#else
      dest[0] = val;
#endif
   }

   /*******************************************
    * Named operations
    *******************************************/

   /// Return requested element of vector (arbitrary idx version)
   double get_element(size_t idx) const {
      double __attribute__((aligned(32))) val_as_array[vector_length];
      store_aligned(val_as_array);
      return val_as_array[idx];
   }


   /// Blend operation: returns (mask) ? x2 : x1
   friend
   SimdVec blend(const SimdVec &x1, const SimdVec &x2, const SimdVec &mask) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec( _mm256_blendv_pd(x1.val, x2.val, mask.val) );
#else
      return SimdVec( (mask.val) ? x2 : x1 );
#endif
   }

   /// Returns absolute values
   friend
   SimdVec fabs(const SimdVec &x) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec(
            _mm256_andnot_pd(_mm256_set1_pd(-0.0), x)
         );
#else
      return SimdVec( fabs(x.val) );
#endif
   }

   /// Return(a + b * c)
   friend
   SimdVec fmadd(const SimdVec &a, const SimdVec &b, const SimdVec &c) {
#if defined(__AVX2__)
      return SimdVec(
            _mm256_fmadd_pd(b.val, c.val, a.val) // b*c + 1
         );
#else
      return a + b*c;
#endif
   }

   /// Return (a - b * c)
   friend
   SimdVec fmsub(const SimdVec &a, const SimdVec &b, const SimdVec &c) {
#if defined(__AVX2__)
      return SimdVec(
            _mm256_fnmadd_pd(b.val, c.val, a.val) // -(b*c) + a
         );
#else
      return a - b*c;
#endif
   }

   /// Return sqrt(a)
   friend
   SimdVec sqrt(const SimdVec &a) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec(
            _mm256_sqrt_pd(a.val)
         );
#else
      return SimdVec( sqrt(a.val) );
#endif
   }

   /*******************************************
    * Operators
    *******************************************/

   /// Conversion to underlying type
   operator simd_double_type() const {
      return val;
   }

   /// Extract indvidual elements of vector (messy and inefficient)
   /// idx MUST be < vector_length.
   double operator[](size_t idx) const {
      return get_element(idx);
   }

   /// Vector valued LT comparison
   friend
   SimdVec operator<(const SimdVec &lhs, const SimdVec &rhs) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec( _mm256_cmp_pd(lhs.val, rhs.val, _CMP_LT_OQ) );
#else
      return SimdVec( lhs.val > rhs.val );
#endif
   }

   /// Vector valued GT comparison
   friend
   SimdVec operator>(const SimdVec &lhs, const SimdVec &rhs) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec( _mm256_cmp_pd(lhs.val, rhs.val, _CMP_GT_OQ) );
#else
      return SimdVec( lhs.val > rhs.val );
#endif
   }

   /// Bitwise and
   friend
   SimdVec operator&(const SimdVec &lhs, const SimdVec &rhs) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec( _mm256_and_pd(lhs.val, rhs.val) );
#else
      return SimdVec( lhs.val && rhs.val );
#endif
   }

   /// Divide
   friend
   SimdVec operator/(const SimdVec &lhs, const SimdVec &rhs) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec( _mm256_div_pd(lhs.val, rhs.val) );
#else
      return SimdVec( lhs.val / rhs.val );
#endif
   }

   /// Multiply
   friend
   SimdVec operator*(const SimdVec &lhs, const SimdVec &rhs) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec( _mm256_mul_pd(lhs.val, rhs.val) );
#else
      return SimdVec( lhs.val * rhs.val );
#endif
   }

   /// Add
   friend
   SimdVec operator+(const SimdVec &lhs, const SimdVec &rhs) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec( _mm256_add_pd(lhs.val, rhs.val) );
#else
      return SimdVec( lhs.val + rhs.val );
#endif
   }

   /// Subtract
   friend
   SimdVec operator-(const SimdVec &lhs, const SimdVec &rhs) {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec( _mm256_sub_pd(lhs.val, rhs.val) );
#else
      return SimdVec( lhs.val + rhs.val );
#endif
   }

   /// Divide-in-place
   void operator/=(const SimdVec &rhs) {
      *this = *this / rhs;
   }

   /// Multiply-in-place
   void operator*=(const SimdVec &rhs) {
      *this = *this * rhs;
   }

   /// Returns true if any bit is non-zero
   friend
   bool any(const SimdVec &a) {
#if defined(__AVX2__) || defined(__AVX__)
      return _mm256_testz_pd(a.val, a.val);
#else
      return a.val;
#endif
   }

   /*******************************************
    * Factory functions for special cases
    *******************************************/

   /// Returns an instance initialized to zero using custom instructions
   static
   SimdVec zero() {
#if defined(__AVX2__) || defined(__AVX__)
      return SimdVec(_mm256_setzero_pd());
#else
      return SimdVec(0.0);
#endif
   }

   /// Returns a vector with all positions idx or above set to true, otherwise
   /// false.
   static
   SimdVec gt_mask(int idx) {
#if defined(__AVX2__) || defined(__AVX__)
      const double avx_true  = -std::numeric_limits<double>::quiet_NaN();
      const double avx_false = 0.0;
      switch(idx) {
         case 0:  return SimdVec(avx_true,   avx_true,  avx_true,  avx_true);
         case 1:  return SimdVec(avx_false,  avx_true,  avx_true,  avx_true);
         case 2:  return SimdVec(avx_false, avx_false,  avx_true,  avx_true);
         case 3:  return SimdVec(avx_false, avx_false, avx_false,  avx_true);
         default: return SimdVec(avx_false, avx_false, avx_false, avx_false);
      }
#else
      return (idx>0) ? SimdVec(false) : SimdVec(true);
#endif
   }

   /// Extract elements of vector into scalars
   void unpack(double v[vector_length]) const {
#if defined(__AVX2__) || defined(__AVX__)
      __m128d low0  = _mm256_extractf128_pd(val, 0);
      __m128d high0 = _mm256_extractf128_pd(val, 1);
      __m128d low1 = _mm_shuffle_pd(low0, low0, _MM_SHUFFLE2(0,1));
      __m128d high1 = _mm_shuffle_pd(high0, high0, _MM_SHUFFLE2(0,1));
      v[0] = _mm_cvtsd_f64(low0);
      v[1] = _mm_cvtsd_f64(low1);
      v[2] = _mm_cvtsd_f64(high0);
      v[3] = _mm_cvtsd_f64(high1);
#else
      v[0] = val;
#endif
   }

   /*******************************************
    * Debug functions
    *******************************************/

   /// Prints the vector (inefficient, use for debug only)
   void print() const {
      for(int i=0; i<vector_length; i++) printf(" %e", (*this)[i]);
   }

private:
   /// Underlying vector that this type wraps
   simd_double_type val;
};
