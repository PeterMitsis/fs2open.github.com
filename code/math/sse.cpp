/*
 * Copyright (C) Volition, Inc. 1999.  All rights reserved.
 *
 * All source code herein is the property of Volition, Inc. You may not sell
 * or otherwise commercially exploit the source or things you created based on
 * the source.
 *
 */

/* sse_3d_mask: used to zero fourth element */
extern const vec3d sse_3d_mask = sse_vec_imask(0xffffffff, 0xffffffff,
											   0xffffffff, 0x00000000);

/* sse_vabs_mask: used to clear sign bit */
extern const vec3d sse_abs_mask = sse_vec_imask(0x7fffffff, 0x7fffffff,
												0x7fffffff, 0x7fffffff);

extern const vec3d sse_zero = sse_vec_init(0.0f, 0.0f, 0.0f, 0.0f);

extern const vec3d sse_x_unit = sse_vec_init(1.0f, 0.0f, 0.0f, 0.0f);
extern const vec3d sse_y_unit = sse_vec_init(0.0f, 1.0f, 0.0f, 0.0f);
extern const vec3d sse_z_unit = sse_vec_init(0.0f, 0.0f, 1.0f, 0.0f);
