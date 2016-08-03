/*
 * Copyright (C) Volition, Inc. 1999.  All rights reserved.
 *
 * All source code herein is the property of Volition, Inc. You may not sell
 * or otherwise commercially exploit the source or things you created based on
 * the source.
 *
 */



#ifndef _VECMAT_H
#define _VECMAT_H

#include "globalincs/pstypes.h"
#include "math/floating.h"


/**
 * @brief Test if a vector has any elements that are NaN
 */
#define vm_is_vec_nan(v) (_isnan((v)->xyz.x) ||    \
						  _isnan((v)->xyz.y) ||    \
						  _isnan((v)->xyz.z))

/**
 * @brief Test if squaring vector results in a null vector
 */
#define IS_VEC_NULL_SQ_SAFE(v)                              \
		(((v)->xyz.x > -1e-16) && ((v)->xyz.x < 1e-16) &&   \
		 ((v)->xyz.y > -1e-16) && ((v)->xyz.y < 1e-16) &&   \
		 ((v)->xyz.z > -1e-16) && ((v)->xyz.z < 1e-16))

/**
 * @brief Test if vector is a null vector
 */
#define IS_VEC_NULL(v) (((v)->xyz.x > -1e-36) && ((v)->xyz.x < 1e-36) && \
						((v)->xyz.y > -1e-36) && ((v)->xyz.y < 1e-36) && \
						((v)->xyz.z > -1e-36) && ((v)->xyz.z < 1e-36) )

/**
 * @brief Test if matrix is a null matrix
 */
#define IS_MAT_NULL(v) (IS_VEC_NULL(&(v)->vec.fvec) &&  \
						IS_VEC_NULL(&(v)->vec.uvec) &&  \
						IS_VEC_NULL(&(v)->vec.rvec))

/**
 * @brief Set vector elements to zero (no return value)
 *
 * @internal Although this could be done with an inline assembly macro, it
 * is expected that is probably better to let the compiler optimize it.
 */
#define vm_vec_zero(v) (v)->xyz.x = (v)->xyz.y = (v)->xyz.z = 0.0f

extern void vm_set_identity(matrix *m);

/**
 * @brief Initialize vector @v with @a _x, @a _y and @a _z
 */
#define vm_vec_make(v,_x,_y,_z)    ((v)->xyz.x=(_x),   \
									(v)->xyz.y=(_y),   \
									(v)->xyz.z=(_z))

/* Global constants */

extern vec3d vmd_zero_vector;
extern vec3d vmd_x_vector;
extern vec3d vmd_y_vector;
extern vec3d vmd_z_vector;
extern matrix vmd_identity_matrix;

/* Handy constants */

#define ZERO_VECTOR { { { 0.0f, 0.0f, 0.0f } } }

/*
 * First set of inside braces is for union,
 * Second set is for inside union, then for a2d[3][3]
 * (some compiler warning messages just suck)
 */
#define IDENTITY_MATRIX { { { { { { 1.0f, 0.0f, 0.0f } } },        \
							  { { { 0.0f, 1.0f, 0.0f } } },        \
							  { { { 0.0f, 0.0f, 1.0f } } } } } }

/**
 * @brief Fill in fields of an angle vector
 */
#define vm_angvec_make(v,_p,_b,_h) (((v)->p = (_p),   \
									 (v)->b = (_b),   \
									 (v)->h = (_h)),  \
									 (v))

/**
 * @brief Negate a vector
 */
#define vm_vec_negate(v)                 \
		do {                             \
			(v)->xyz.x = - (v)->xyz.x;   \
			(v)->xyz.y = - (v)->xyz.y;   \
			(v)->xyz.z = - (v)->xyz.z;   \
		} while (0)

typedef struct plane {
	float  A, B, C, D;
} plane;

/* Functions in library */

/**
 * @brief Add two vectors
 *
 * @a dest = @a src0 + @a src1.
 *
 * Although @a dest may point to either @a src0 or @a src1, it is recommended
 * that vm_vec_add2() be used instead for that case.
 *
 * @return N/A
 */
void vm_vec_add(vec3d *dest, const vec3d *src0, const vec3d *src1);

/**
 * @brief In-place addition of two vectors
 *
 * @a dest += @a src
 *
 * @return N/A
 */
void vm_vec_add2(vec3d *dest, const vec3d *src);

/**
 * @brief In-place subtraction of a scaled vector
 *
 * @a dest -= @a src * @a k
 *
 * @return N/A
 */
void vm_vec_scale_sub2(vec3d *dest, const vec3d *src, float k);

/**
 * @brief Subtract two vectors
 *
 * @a dest = @a src0 - @a src1
 *
 * Although @a dest may point to either @a src0 or @a src1, it is recommended
 * that vm_vec_sub2() be used instead for that case.
 *
 * @return N/A
 */
void vm_vec_sub(vec3d *dest, const vec3d *src0, const vec3d *src1);

/**
 * @brief In-place subtraction of two vectors
 *
 * @a dest -= @a src
 *
 * @return N/A
 */
void vm_vec_sub2(vec3d *dest, const vec3d *src);

/**
 * @brief Calculate the average of @a n vectors
 *
 * @return Average of the @a n vectors (@a dest)
 */
vec3d *vm_vec_avg_n(vec3d *dest, int n, const vec3d src[]);

/**
 * @brief Calculate the average of two vectors
 *
 * @return Average of the two vectors (@a dest)
 */
vec3d *vm_vec_avg(vec3d *dest, const vec3d *src0, const vec3d *src1);

/**
 * @brief Calculate the average of three vectors
 *
 * @return Average of the three vectors (@a dest)
 */
vec3d *vm_vec_avg3(vec3d *dest, const vec3d *src0, const vec3d *src1,
				   const vec3d *src2);

/**
 * @brief Calculate the average of four vectors
 *
 * @return Average of the four vectors (@a dest)
 */
vec3d *vm_vec_avg4(vec3d *dest, const vec3d *src0, const vec3d *src1,
				   const vec3d *src2, const vec3d *src3);

/**
 * @brief In-place Scaling of a vector
 *
 * @a dest *= @a s
 *
 * @return N/A
 */
void vm_vec_scale(vec3d *dest, float s);

/**
 * @brief Scale a vector
 *
 * @a dest = @a src * @a s
 *
 * @return N/A
 */
void vm_vec_copy_scale(vec3d *dest, const vec3d *src, float s);

/**
 * @brief Add a vector and scaled vector
 *
 * @a dest = @a src1 + @a src2 * @a k
 *
 * @return N/A
 */
void vm_vec_scale_add(vec3d *dest, const vec3d *src1, const vec3d *src2,
					  float k);

/**
 * @a brief Subtract a scaled vector from a vector
 *
 * @a dest = @a src1 - @a src2 * @a k
 *
 * @return N/A
 */
void vm_vec_scale_sub(vec3d *dest, const vec3d *src1, const vec3d *src2,
					  float k);

/**
 * @brief In-place addition of a scaled vector
 *
 * @a dest += @a src * @a k
 *
 * @return N/A
 */
void vm_vec_scale_add2(vec3d *dest, const vec3d *src, float k);

/**
 * @brief In-place scaling of a vector
 *
 * @a dest *= @a n / @a d
 *
 * @return N/A
 */
void vm_vec_scale2(vec3d *dest, float n, float d);

/**
 * @brief Check if two vectors are equal
 *
 * @return true if equal, false if not
 */
bool vm_vec_equal(const vec2d &self, const vec2d &other);

/**
 * @brief Check if two vectors are equal
 *
 * @return true if equal, false if not
 */
bool vm_vec_equal(const vec3d &self, const vec3d &other);

/**
 * @brief Check if two vectors are equal
 *
 * @return true if equal, false if not
 */
bool vm_vec_equal(const vec4 &self, const vec4 &other);

/**
 * @brief Check if two matrices are equal
 *
 * @return true if equal, false if not
 */
bool vm_matrix_equal(const matrix &self, const matrix &other);

/**
 * @brief Check if two matrices are equal
 *
 * @return true if equal, false if not
 */
bool vm_matrix_equal(const matrix4 &self, const matrix4 &other);

/**
 * @brief Calculate projection of source vector along a unit vector
 *
 * @param component Projected vector
 * @param src Source vector
 * @param unit_vector Unit vector
 *
 * @return Magnitude of the resulting @a component
 */
float vm_vec_projection_parallel (vec3d *component, const vec3d *src,
								  const vec3d *unit_vector);

/**
 * @brief Calclate the projection of source vector onto a surface
 *
 * @param projection Projected vector
 * @param src Source vector
 * @param normal Unit vector normal to the plane
 *
 * @return N/A
 */
void vm_vec_projection_onto_plane (vec3d *projection, const vec3d *src,
								   const vec3d *normal);

/**
 * @brief Calculate magnitude of a vector
 *
 * @return Magnitude of vector @a v
 */
float vm_vec_mag(const vec3d *v);

/**
 * @brief Calculate the square of a vector magnitude
 *
 * Similar to vm_vec_mag(), this routine calculates the square of the
 * magnitude of the vector @a v. It can be useful if comparing distances
 * as it avoids a costly sqrt() operation.
 *
 * @return Square of the magnitude of vector @a v
 */
float vm_vec_mag_squared(const vec3d* v);

/**
 * @brief Calculate the square of the distance between two points
 *
 * @return Square of the distance
 */
float vm_vec_dist_squared(const vec3d *v0, const vec3d *v1);

/**
 * @brief Calculates the distance between two points
 *
 * @return Distance between two points
 */
float vm_vec_dist(const vec3d *v0, const vec3d *v1);

/**
 * @brief Approximate the magnitude of the vector
 *
 * The current implementation approximates the magnitude of the vector as
 * follows ...
 *     dist = largest + next_largest*3/8 + smallest*3/16
 *
 * @return Approximate vector magnitude
 */
float vm_vec_mag_quick(const vec3d *v);

/**
 * @brief Approximate the distance between two vectors
 *
 * Uses vm_vec_mag_quick()'s approximation method.
 *
 * @return Approximate distance betwen two points
 */
float vm_vec_dist_quick(const vec3d *v0, const vec3d *v1);


/**
 * @brief Copy normalized @a src into @a dest
 *
 * If @a src is a null vector, then use x-unit vector and generate a Warning().
 *
 * @return Magnitude of @a src
 */
float vm_vec_copy_normalize(vec3d *dest, const vec3d *src);

/**
 * @brief In-place normalization of a vector
 *
 * If @a v is a null vector, then use x-unit vector and generate a Warning().
 *
 * @return Magnitude of @a v
 */
float vm_vec_normalize(vec3d *v);

/**
 * @brief In-place normalization of a vector
 *
 * If @a v is a null vector, then use x-unit vector but do NOT generate
 * a Warning().
 *
 * @return Magnitude of @a v
 */
float vm_vec_normalize_safe(vec3d *v);

/**
 * @brief Copy normalized @a src to @a dest
 *
 * @return Magnitude of @a src
 */
float vm_vec_copy_normalize_quick(vec3d *dest, const vec3d *src);

/**
 * @brief In-place normalization of @a v
 *
 * @return Magnitude of @a src
 */
float vm_vec_normalize_quick(vec3d *v);

/**
 * @brief Copy normalized @a src to @a dest
 *
 * Uses vm_vec_mag_quick() in the calculation.
 *
 * @return Quick magnitude of @a src
 */
float vm_vec_copy_normalize_quick_mag(vec3d *dest, const vec3d *src);

/**
 * @brief In-place normalization of @a v
 *
 * Uses vm_vec_mag_quick() in the calculation.
 *
 * @return Quick magnitude of @a src
 */
float vm_vec_normalize_quick_mag(vec3d *v);

/**
 * @brief Calculate normalized direction vector between two points
 *
 * @a dest = normalized(@a end - @a start). Note that the order of parameters
 * matches that of vm_vec_sub().
 *
 * @return Magnitude of direction vector
 */
float vm_vec_normalized_dir(vec3d *dest,const vec3d *end, const vec3d *start);

/**
 * @brief Calculate normalized direction vector between two points
 *
 * Like vm_vec_normalized_dir() but uses vm_vec_normalize_quick_mag().
 *
 * @return Magnitude of direction vector
 */
float vm_vec_normalized_dir_quick_mag(vec3d *dest, const vec3d *end,
									  const vec3d *start);
/**
 * @brief Normalizes the direction vector @a end - @a start
 *
 * This implementation works with reciprocals which allows an x86 compiler to
 * use the RSQRTSS or RCPSS instructions if they are known to be available.
 *
 * @return Approximate magnitude of the direction vector
 */ 
float vm_vec_normalized_dir_quick(vec3d *dest, const vec3d *end,
								  const vec3d *start);

/**
 * @brief Calculate the dot product of two vectors
 *
 * @return Dot product
 */
float vm_vec_dot(const vec3d *v0, const vec3d *v1);

/**
 * @brief Calculate the dot product of two vectors
 *
 * One of the vectors is explicitly specfied by its three (3) components.
 *
 * @return Dot product
 */

float vm_vec_dot3(float x, float y, float z, vec3d *v);

/**
 * @brief Calclate the cross product of two vectors
 *
 * It is important to note that that @a dest must not point to either
 * @a src0 or @a src1 lest the calculation be corrupted.
 *
 * @return Cross product in @a dest
 */
vec3d *vm_vec_cross(vec3d *dest, const vec3d *src0, const vec3d *src1);

/**
 * @brief Determine if two vectors are parallel
 *
 * @return 0 if not parallel, 1 if parallel
 */
int vm_test_parallel(const vec3d *src0, const vec3d *src1);

/**
 * @brief Calclate surface normal from three points
 *
 * It is important to note that that @a dest must not point to any of
 * @a p0 or @a p1 or @a p2 lest the calculation be corrupted.
 *
 * @return Normalized surface normal stored in @a dest
 */
vec3d *vm_vec_normal(vec3d *dest,const vec3d *p0, const vec3d *p1,
					 const vec3d *p2);

/**
 * @brief Calculate vector perpendicular to three points
 *
 * It is important to note that that @a dest must not point to any of
 * @a p0 or @a p1 or @a p2 lest the calculation be corrupted.
 *
 * @return Non-normalized surface normal
 */
vec3d *vm_vec_perp(vec3d *dest, const vec3d *p0, const vec3d *p1,
				   const vec3d *p2);

/**
 * @brief Calculate the angle between two vectors
 *
 * The forward vector @a fvec is used to control the range of the resulting
 * angle. If @a fvec is NULL, then the range of the angle between @a v0 and
 * @a v1 is in the range of [0..PI). If @a fvec is non-NULL, then the
 * absolute value of the angel is returned and it is in the range of [0..PI/2].
 *
 * Note that neither @a v0 nor @a v1 need to be normalized.
 *
 * @return Angle between vectors @a v0 and @a v1
 */
float vm_vec_delta_ang(const vec3d *v0, const vec3d *v1, const vec3d *fvec);

/**
 * @brief Calclate the angle between two normalized vectors
 *
 * The forward vector @a fvec is used to control the range of the resulting
 * angle. If @a fvec is NULL, then the range of the angle between @a v0 and
 * @a v1 is in the range of [0..PI). If @a fvec is non-NULL, then the
 * absolute value of the angel is returned and it is in the range of [0..PI/2].
 *
 * @return Angle between vectors @a v0 and @a v1
 */
float vm_vec_delta_ang_norm(const vec3d *v0, const vec3d *v1,const vec3d *fvec);

/**
 * @brief Generate matrix @a m from a set of three angles
 *
 * @return Pointer to the generated matrix @a m
 */
matrix *vm_angles_2_matrix(matrix *m, const angles *a);

/**
 * @brief Generate matrix @a m from a single angle @a a
 *
 * @param angle_index 0 for pitch, 1 for bank, 2 for heading
 *
 * @return Pointer to the generated matrx @a m
 */
matrix *vm_angle_2_matrix(matrix *m, float a, int angle_index);

/**
 * @brief Generate matrix @a m from a forward vector and an angle
 *
 * @return Pointer to the generated matrix @a m
 */
matrix *vm_vec_ang_2_matrix(matrix *m, const vec3d *v, float a);

/**
 * @brief Generate matrix @a m from one or more vectors
 *
 * The forward vector @a fvec is required while both the up vector @a uvec
 * and the right vector @a rvec are optional. If both the up and right vectors
 * are passed then the up vector is used. If only the forward is passed then
 * a bank of zero is assumed.
 *
 * @return Pointer to the generated matrix @a m
 */
matrix *vm_vector_2_matrix(matrix *m, const vec3d *fvec, const vec3d *uvec,
						   const vec3d *rvec);

/**
 * @brief Generate matrix @a m from one or more vectors
 *
 * Similar to vm_vector_2_matrix(), this routine requires that the vectors be
 * more-or-less normalized and closed to perpendiclar.
 *
 * @return Pointer to the generated matrix @a m
 */
matrix *vm_vector_2_matrix_norm(matrix *m, const vec3d *fvec,
								const vec3d *uvec = NULL,
								const vec3d *rvec = NULL);

/**
 * @brief Rotate vector @a src through the matrix @a m
 *
 * It is important to note that the @a dest vector must not point to @a src.
 *
 * @return The rotated vector @a dest
 */
vec3d *vm_vec_rotate(vec3d *dest, const vec3d *src, const matrix *m);

/**
 * @brief Rotate vector @a src through the transpose of matrix @a m
 *
 * This routine is meant to replace the following code sequence ...
 *
 *     vm_copy_transpose(&tmp_matrix, src_matrix);
 *     vm_vec_rotate(dest_vec, src_vec, &tmp_matrix);
 *
 * ... with ...
 *
 *     vm_vec_unrotate(dest_vec, src_vec, src_matrix);
 *
 * It is important to note that the @a dest vector must not point to @a src.
 * As this routine does not modify the source matrix @a m, if you need it
 * transposed later on then you should use the following technique instead ...
 *
 *    vm_vec_transpose(...)
 *    vm_vec_rotate(...)
 *
 * @return Pointer to rotated vector @a dest
 */
vec3d *vm_vec_unrotate(vec3d *dest, const vec3d *src, const matrix *m);

/**
 * @brief In-place matrix transposition
 *
 * @return Pointer to the transposed matrix
 */
matrix *vm_transpose(matrix *m);

/**
 * @brief Copy and transpose a matrix
 *
 * It is important to note that the @a dest matrix must not point to @a src.
 *
 * @return Pointer to the transposed matrix @a dest
 */
matrix *vm_copy_transpose(matrix *dest, const matrix *src);

/**
 * @brief Multiply two matrices
 *
 * It is important to note that the @a dest matrix must not point to either
 * @a src0 or @a src1.
 *
 * @return Pointer to destination matrix @a dest
 */
matrix *vm_matrix_x_matrix(matrix *dest, const matrix *src0,
						   const matrix *src1);

/**
 * @brief Extract angles from a matrix
 *
 * @return Pointer to extracted angles @a a
 */
angles *vm_extract_angles_matrix(angles *a, const matrix *m);

/**
 * @brief Alternate method to extract angles from a matrix
 *
 * This method for extracting angles seems to be less susceptible to rounding
 * errors. See section 8.7.2 (pages 278-281) of the 3d Math Primer for Graphics
 * and Game Development, 2nd Edition
 *
 * http://books.google.com/books?id=X3hmuhBoFF0C&printsec=frontcover#v=onepage&q&f=false
 *
 * @return Pointer to extracted angles @a a
 */
angles *vm_extract_angles_matrix_alternate(angles *a, const matrix *m);

/**
 * @brief Extract heading and pitch from a vector (assumes bank == 0)
 *
 * @return Pointer to extracted angles @a a
 */
angles *vm_extract_angles_vector(angles *a, const vec3d *v);

/**
 * @brief In-place orthogonalization of matrix @a m_src
 *
 * This routine computes a matrix from one or more vectors found within it.
 * The forward vector is required while the other two are optional. If both
 * the up and right vectors are non-zero then the up vector is used. If only
 * the forward vector is non-zero, then a bank of zero is assumed.
 *
 * @return Pointer to orthogonalized matrix
 */
void vm_orthogonalize_matrix(matrix *m_src);

/**
 * @brief In-place orthogonalization of matrix @a m_src
 *
 * This routine is ike vm_orthogonalize_matrix(), except that zero vectors
 * may exist within the matrix without causing problems. Valid vectors are
 * created where needed.
 *
 * @return N/A
 */
void vm_fix_matrix(matrix *m);

/**
 * @brief In-place rotation and orthogonalization of matrix @a orient
 *
 * Rotates the matrix @a orient by the angles in @a tangles and then
 * makes sure that the matrix is orthogonal.
 *
 * @return N/A
 */
void vm_rotate_matrix_by_angles(matrix *orient, const angles *tangles);

/**
 * @brief Calculate the distance from a point to a plane
 *
 * @param checkp Point to check
 * @param norm   Normalized normal to the plane
 * @param planep Point on the plane
 *
 * @return Distance to the plane (if < 0, then it is behind the plane)
 */
float vm_dist_to_plane(const vec3d *checkp, const vec3d *norm,
					   const vec3d *planep);

/**
 * @brief Generate 3x3 rotation matrix based on mouse movements
 *
 * Given mouse movements in @a idx and @a idy, generate a 3x3 rotation matrix
 * in @a RotMat.
 *
 * Taken from Graphics Gems III, page 51, "The Rolling Ball".
 * Example:
 *    if ((Mouse.dx != 0) || (Mouse.dy != 0)) {
 *        vm_trackball(Mouse.dx, Mouse.dy, &MouseRotMat);
 *        vm_matrix_x_matrix(&tempm, &LargeView.ev_matrix, &MouseRotMat);
 *        LargeView.ev_matrix = tempm;
 *    }
 *
 * @return N/A
 */
void vm_trackball(int idx, int idy, matrix * RotMat);

/**
 * @brief Find nearest point on line
 *
 * This routine finds the nearest point on the line defined by @a p0 and @a p1
 * to the point @a int_pnt and stores it in @a nearest_point. The magnitude of
 * the return value indicates where on the line segment that point resides:
 *    0.0f <= x <= 1.0f means @a p0 <= @a nearest_point < @a p1
 *    x < 0.0f means @a nearest_point is before @a p0
 *    x > 1.0f means @a nearest_point is beyond @a p1
 * Example: If the return value 'x' is 2.0f, then @a nearest_point is beyond
 * @a p1 by 2 times.
 *
 * @return Value indicating whether @a nearest_point is between @a p0 and @a p1
 */
float find_nearest_point_on_line(vec3d *nearest_point, const vec3d *p0,
								 const vec3d *p1, const vec3d *int_pnt);

/**
 * @brief Calculate dot product of @a dir dot (@a p2 - @a p1)
 *
 * It is important to note that ...
 *   1. @a dir must be a unit vector
 *   2. (@a p2 - @a p1) will be normalized and not zero
 *
 * @return dot product
 */
float vm_vec_dot_to_point(const vec3d *dir, const vec3d *p1, const vec3d *p2);

/**
 * @brief Calculate point on plane @a planep closest to point @a p
 *
 * Puts the result of the calculation in @a q.
 *
 * @return N/A
 */
void compute_point_on_plane(vec3d *q, const plane *planep, const vec3d *p);

/**
 * @brief Finds the point on a plane closest to a given point
 *
 * @param new_point    Point on the plane (result)
 * @param point        Point to which to compute closest plane point
 * @param plane_normal Plane normal
 * @param plane_point  Plane point
 *
 * @return N/A
 */
void vm_project_point_onto_plane(vec3d *new_point, const vec3d *point,
								 const vec3d *plane_normal,
								 const vec3d *plane_point);

/**
 * @brief Generate a fairly random vector that is close to normalized
 *
 * @return N/A
 */
void vm_vec_rand_vec_quick(vec3d *rvec);

/**
 * @brief Rotate a point around a line
 *
 * @param out        Result of point rotation
 * @param in         Point to rotate around the line
 * @param angle      Number of radians to rotate
 * @param line_point Point on the line
 * @param line_dir   Line's normalized direction
 *
 * @return N/A
 */
void vm_rot_point_around_line(vec3d *out, const vec3d *in, float angle,
							  const vec3d *line_point, const vec3d *line_dir);

/**
 * @brief Determine if two position vectors are the same
 *
 * The two vectors are considered to be different if the distance between
 * them is greater than 0.005f.
 *
 * @return 0 if the same; 1 if different
 */
int vm_vec_cmp(const vec3d * a, const vec3d * b);

/**
 * @brief Determine if two orientation matrices are the same
 *
 * @return 0 if the same; 1 if different
 */
int vm_matrix_cmp(const matrix * a, const matrix * b);

/**
 * @brief Move @a h towards @a desired_angle
 *
 * Moves 'heading' angle @a h towards @a desired_angle taking the shortest
 * route possible. It will move a maximum of @a step_size radians each call.
 * All angles are in radians.
 *
 * @return Delta between @h and @a desired_angle
 */
float vm_interp_angle(float *h, float desired_angle, float step_size,
					  bool force_front = false);

/**
 * @brief Calculate the difference between two angles
 *
 * Uses the same method as with vm_interp_angle()
 *
 * @return Delta between @a current_angle and @a desired_angle
 */
float vm_delta_from_interp_angle(float current_angle, float desired_angle);

/**
 * @brief Test a matrix for any zero rows or zero columns
 *
 * @return 0 if no zero rows or zero columns; 1 otherwise
 */
int vm_check_matrix_for_zeros(const matrix *m);

/**
 * @brief Test if two vectors are identical
 *
 * @return 1 if vectors are identical; 0 otherwise
 */
int vm_vec_same(const vec3d *v1, const vec3d *v2);

/**
 * @brief Test if two matrices are identical
 *
 * @return 1 if matrices are identical; 0 otherwise
 */
int vm_matrix_same(matrix *m1, matrix *m2);

/**
 * @brief Interpolate from a start matrix towards a goal matrix
 *
 * Moves at maximum rotational acceleration twoard the goal when far and then
 * max decelearation when close to minimize time between orientations. It is
 * subject to constraints on rotational velocity and angular acceleration.
 * The next orientation is valid at time @a delta_t.
 *
 * @return N/A
 */
// Returns next_orientation valid at time delta_t.
void vm_matrix_interpolate(const matrix *goal_orient,
						   const matrix *start_orient,
						   const vec3d *rotvel_in, float delta_t,
						   matrix *next_orient, vec3d *rotvel_out,
						   const vec3d *rotvel_limit, const vec3d *acc_limit,
						   int no_overshoot = 0);

/**
 * @brief Interpolate from a start forward vector twoard a goal forward vector
 *
 * Moves at maximum rotational acceleration toward the goal when far and then
 * max deceleration when close to minimize time between orientations. It is
 * subject to constaints on rotational velocity and angular accleleration.
 * The next forward vector is valid at time @a delta_t.
 *
 * @return N/A
 */
void vm_forward_interpolate(const vec3d *goal_fvec, const matrix *orient,
							const vec3d *rotvel_in,
							float delta_t, float delta_bank,
							matrix *next_orient, vec3d *rotvel_out,
							const vec3d *vel_limit, const vec3d *acc_limit,
							int no_overshoot = 0);

/**
 * @brief Find the bounding sphere for a set of points
 *
 * Places the calculated center and radius in @a center and @a radius
 * respectively.
 *
 * @return N/A
 */
void vm_find_bounding_sphere(const vec3d *pnts, int num_pnts,
							 vec3d *center, float *radius);

/**
 * @brief Version of atan2() that is safe for optimized builds
 *
 * @return arc-tangent of y/x
 */
float atan2_safe(float x, float y);

/**
 * @brief Translate from world coordinates to body coordinates
 *
 * @return Vector in body coordinates
 */
vec3d* vm_rotate_vec_to_body(vec3d *body_vec, const vec3d *world_vec,
							 const matrix *orient);

/**
 * @brief Translate from body coordinates to world coordinates
 *
 * @return Vector in world coordinates
 */
vec3d* vm_rotate_vec_to_world(vec3d *world_vec, const vec3d *body_vec,
							  const matrix *orient);

/**
 * @brief Estimate next orientation matrix as extrapolation of last and current
 *
 * @return N/A
 */
void vm_estimate_next_orientation(const matrix *last_orient,
								  const matrix *current_orient,
								  matrix *next_orient);

/**
 * @brief Test if any elements in @a vec are NaN
 *
 * @return 0 if any elements are NaN; 1 if no elements are NaN
 */
int is_valid_vec(const vec3d *vec);

/**
 * @brief Test if any elements in @a m are NaN
 *
 * @return 0 if any elements are @a m; 1 if no elements are NaN
 */
int is_valid_matrix(const matrix *m);

/**
 * @brief Generate rotation matrix
 *
 * Finds the rotation matrix corresponding to a rotation of @a theta about 
 * axis @a u.
 *
 * @return N/A
 */
void vm_quaternion_rotate(matrix *m, float theta, const vec3d *u);

/**
 * @brief Calculate angle needed to generate matrix @a m
 *
 * @return N/A
 */
void vm_matrix_to_rot_axis_and_angle(const matrix *m, float *theta,
									 vec3d *rot_axis);

/**
 * @brief Interpolate between two vectors
 *
 * It is important to note that 0.0f <= @a t <= @a 1.0f
 *
 * @return N/A
 */
void vm_vec_interp_constant(vec3d *out, const vec3d *v1, const vec3d *v2,
							float t);

/**
 * @brief Randomly perturb a vector
 *
 * Randomly perturbs a vector around a given (normalized vector) or optional
 * orientation matrix.
 *
 * @return N/A
 */
void vm_vec_random_cone(vec3d *out, const vec3d *in, float max_angle,
						const matrix *orient = NULL);

/**
 * @brief Randomly perturb a vector
 *
 * Randomly perturbs a vector around a given (normalized vector) or optional
 * orientation matrix.
 *
 * @return N/A
 */
void vm_vec_random_cone(vec3d *out, const vec3d *in, float min_angle,
						float max_angle, const matrix *orient = NULL);

/**
 * @brief Generate a random point on the plane of a circle
 *
 * Given the start vector @a in, the orientation matrix @a orient and the
 * radius @a radius, generate a point on the plane of the circle. If
 * @a on_edge is 1 then the point is on the very edge of the circle.
 *
 * @return N/A
 */
void vm_vec_random_in_circle(vec3d *out, const vec3d *in,
							 const matrix *orient, float radius, int on_edge);

/**
 * @brief Generate a random point in a spherical volume
 *
 * Given the start vector @a in, the orientation matrix @a orient and the
 * radius @a radius, generate a point in a spherical volume. If @a on_edge is
 * 1 then the point is on the very edge of the sphere.
 *
 * @return N/A
 */
void vm_vec_random_in_sphere(vec3d *out, const vec3d *in, const matrix *orient,
							 float radius, int on_edge);

/**
 * @brief Find the destance to a line
 *
 * Find the nearest point on the line segment defined by @a l0 to @a l1 to
 * point @a p. Puts the result in @a dest if non-NULL.
 *
 * @retval -1 if the point is 'before' the line segment
 * @retval 0 if the point is inside the line segment
 * @retval 1 if the point is 'after' the line segment
 */
int vm_vec_dist_to_line(const vec3d *p, const vec3d *l0, const vec3d *l1,
						vec3d *nearest, float *dist);

/**
 * @brief Find the square of the distance to a line
 *
 * Like vm_vec_dist_to_line() except that it does not check whether the
 * nearest point is on the line segment.
 *
 * @return N/A
 */
void vm_vec_dist_squared_to_line(const vec3d *p, const vec3d *l0,
								 const vec3d *l1, vec3d *nearest,
								 float *dist_squared);
/**
 * @brief 2D vector "box" scaling
 *
 * In-place scaling of the vector so that the longest dimension = @a scale
 *
 * @return N/A
 */
void vm_vec_boxscale(vec2d *vec, float scale);

/**
 * @brief Attempt to invert a 4x4 matrix
 *
 * @param m       Pointer to matrix to invert
 * @param invOut  Inverted matrix, or nullptr if inversion is impossible
 *
 * @return true if matrix is invertible; false otherwise
 */
bool vm_inverse_matrix4(const matrix4 *m, matrix4 *invOut);

#endif
