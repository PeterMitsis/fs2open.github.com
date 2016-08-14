#ifndef __SSE_H_
#define __SSE_H_
#include <immintrin.h>

/**
 * @brief Generate a mask for SSE operations
 */
#define sse_vec_mask(i, j, k, l) { .ijk = {i, j, k, l} }

/**
 * @brief Initialize an SSE vector
 */
#define sse_vec_init(x, y, z, w) { .xyz = {x, y, z, w} }

extern const vec3d sse_3d_mask;    /* For clearing 4th element in SSE vector */
extern const vec3d sse_abs_mask;   /* For determining absolute values */

extern const vec3d sse_zero;       /* Zero vector */
extern const vec3d sse_x_unit;     /* Unit vector x direction */
extern const vec3d sse_y_unit;     /* Unit vector y direction */
extern const vec3d sse_z_unit;     /* Unit vector z direction */

/**
 * @return absolute values
 */
static inline __m128 sse_abs(__m128 v)
{
	return _mm_and_ps(v, sse_abs_mask.xyzw);
}

/**
 * @return addend1 + addend2
 */
static inline __m128 sse_vec_add(__m128 addend1, __m128 addend2)
{
	return _mm_add_ps(addend1, addend2);
}

/**
 * @return minuend - subtrahend
 */
static inline __m128 sse_vec_sub(__m128 minuend, __m128 subtrahend)
{
	return _mm_sub_ps(minuend, subtrahend);
}

/**
 * @return factor1 * factor2
 */
static inline __m128 sse_vec_mul(__m128 factor1, __m128 factor2)
{
	return _mm_mul_ps(factor1, factor2);
}

/**
 * @return dividend / divisor
 */
static inline __m128 sse_vec_div(__m128 dividend, __m128 divisor)
{
	return _mm_div_ps(dividend, divisor);
}

/**
 * @return Approximate @a dividend / @a divisor
 */
static inline __m128 sse_vec_div_quick(__m128 dividend, __m128 divisor)
{
	return _mm_mul_ps(dividend, _mm_rcp_ps(divisor));
}

/**
 * @return {konstant, konstant, konstant, konstant}
 */
static inline __m128 sse_vec_constant(float konstant)
{
	return _mm_set1_ps(konstant);
}

/**
 * @return factor * konstant
 */
static inline __m128 sse_vec_scale(__m128 factor, float konstant)
{
	return sse_vec_mul(factor, sse_vec_constant(konstant));
}

#if defined(__SSE3__)
/**
 * @return vector of (@a v1 dot @a v2)
 */
static inline __m128 sse_vec_vdot3(__m128 v1, __m128 v2)
{
	__m128 dp;

	dp = sse_vec_mul(_mm_and_ps(sse_3d_mask.xyzw, v1), v2);
	dp = _mm_hadd_ps(dp, dp);
	dp = _mm_hadd_ps(dp, dp);

	return dp;
}
#elif defined(__SSE2__)
#endif

/**
 * @return float of (v dot v)
 */
static inline float sse_vec_dot3(__m128 v1, __m128 v2)
{
	return _mm_cvtss_f32(sse_vec_vdot3(v1, v2));
}

/**
 * @return vector of squared magnitudes of @a v
 */
static inline __m128 sse_vec_vmag3_squared(__m128 v)
{
	return sse_vec_vdot3(v, v);
}

/**
 * @return float of squared magnitude of @a v
 */
static inline float sse_vec_mag3_squared(__m128 v)
{
	return _mm_cvtss_f32(sse_vec_vmag3_squared(v));
}

/**
 * @return vector of magnitude of @a v
 */
static inline __m128 sse_vec_vmag3(__m128 v)
{
	return _mm_sqrt_ps(sse_vec_vmag3_squared(v));
}

/**
 * @return float of magnitude of @a v
 */
static inline float sse_vec_mag3(__m128 v)
{
	return _mm_cvtss_f32(_mm_sqrt_ss(sse_vec_vmag3_squared(v)));
}

/**
 * @return vector of inverse magnitude of @a v
 */
static inline __m128 sse_vec_vmag3_inv_quick(__m128 v)
{
	return _mm_rsqrt_ps(sse_vec_vmag3_squared(v));
}

/**
 * @return float of inverse magnitude of @a v
 */
static inline float sse_vec_mag3_inv_quick(__m128 v)
{
	return _mm_cvtss_f32(_mm_rsqrt_ss(sse_vec_vmag3_squared(v)));
}

/**
 * @return - quick approximation of vector of magnitudes
 * Except for the denormalized values, error does not exceed 0.047%.
 */
static inline __m128 sse_vec_vmag3_quick(__m128 v)
{
	return _mm_rcp_ps(sse_vec_vmag3_inv_quick(v));
}

/**
 * @return - quick approximation of vector of magnitudes
 * Except for the denormalized values, error does not exceed 0.0472%.
 */
static inline float sse_vec_mag3_quick(__m128 v)
{
	__m128  squared = sse_vec_vmag3_squared(v);

	return _mm_cvtss_f32(_mm_rcp_ss(_mm_rsqrt_ss(squared)));
}

/**
 * @return - v1 cross v2  (4th element is garbage)
 */
static inline __m128 sse_vec_cross3(__m128 v1, __m128 v2)
{
	return _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(v1, v1,
												_MM_SHUFFLE(3, 0, 2, 1)),
								 _mm_shuffle_ps(v2, v2,
												_MM_SHUFFLE(3, 1, 0, 2))),
					  _mm_mul_ps(_mm_shuffle_ps(v1, v1,
												_MM_SHUFFLE(3, 1, 0, 2)),
								 _mm_shuffle_ps(v2, v2,
												_MM_SHUFFLE(3, 0, 2, 1))));
}
#endif /* __SSE_H_ */
#endif /* __SSE_H_ */
