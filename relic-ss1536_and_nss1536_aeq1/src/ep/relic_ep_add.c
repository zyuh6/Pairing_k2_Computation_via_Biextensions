/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2009 RELIC Authors
 *
 * This file is part of RELIC. RELIC is legal property of its developers,
 * whose names are not listed here. Please refer to the COPYRIGHT file
 * for contact information.
 *
 * RELIC is free software; you can redistribute it and/or modify it under the
 * terms of the version 2.1 (or later) of the GNU Lesser General Public License
 * as published by the Free Software Foundation; or version 2.0 of the Apache
 * License as published by the Apache Software Foundation. See the LICENSE files
 * for more details.
 *
 * RELIC is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the LICENSE files for more details.
 *
 * You should have received a copy of the GNU Lesser General Public or the
 * Apache License along with RELIC. If not, see <https://www.gnu.org/licenses/>
 * or <https://www.apache.org/licenses/>.
 */

/**
 * @file
 *
 * Implementation of the point addition on prime elliptic curves.
 *
 * @ingroup ep
 */

#include "relic_core.h"
#include "relic_ep_add_tmpl.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#if EP_ADD == BASIC || !defined(STRIP)

/**
 * Adds two points represented in affine coordinates on an ordinary prime
 * elliptic curve.
 *
 * @param[out] r			- the result.
 * @param[out] s			- the slope.
 * @param[in] p				- the first point to add.
 * @param[in] q				- the second point to add.
 */
TMPL_ADD_BASIC_IMP(ep, fp);

#endif /* EP_ADD == BASIC */

#if EP_ADD == PROJC || !defined(STRIP)

/**
 * Adds a point represented in homogeneous coordinates to a point represented in
 * affine coordinates on an ordinary prime elliptic curve.
 *
 * @param[out] r			- the result.
 * @param[in] p				- the projective point.
 * @param[in] q				- the affine point.
 */
TMPL_ADD_PROJC_MIX(ep, fp);

/**
 * Adds two points represented in homogeneous coordinates on an ordinary prime
 * elliptic curve.
 *
 * @param[out] r			- the result.
 * @param[in] p				- the first point to add.
 * @param[in] q				- the second point to add.
 */
TMPL_ADD_PROJC_IMP(ep, fp);

#endif /* EP_ADD == PROJC */

#if EP_ADD == JACOB || !defined(STRIP)

/**
 * Adds a point represented in Jacobian coordinates to a point represented in
 * affine coordinates on an ordinary prime elliptic curve.
 *
 * @param[out] r			- the result.
 * @param[in] p				- the projective point.
 * @param[in] q				- the affine point.
 */
TMPL_ADD_JACOB_MIX(ep, fp);

/**
 * Adds two points represented in Jacobian coordinates on an ordinary prime
 * elliptic curve.
 *
 * @param[out] r			- the result.
 * @param[in] p				- the first point to add.
 * @param[in] q				- the second point to add.
 */
TMPL_ADD_JACOB_IMP(ep, fp);

#endif /* EP_ADD == JACOB */

/*============================================================================*/
	/* Public definitions                                                         */
/*============================================================================*/

#if EP_ADD == BASIC || !defined(STRIP)

void ep_add_basic(ep_t r, const ep_t p, const ep_t q) {
	if (ep_is_infty(p)) {
		ep_copy(r, q);
		return;
	}

	if (ep_is_infty(q)) {
		ep_copy(r, p);
		return;
	}

	ep_add_basic_imp(r, NULL, p, q);
}

void ep_add_slp_basic(ep_t r, fp_t s, const ep_t p, const ep_t q) {
	if (ep_is_infty(p)) {
		ep_copy(r, q);
		return;
	}

	if (ep_is_infty(q)) {
		ep_copy(r, p);
		return;
	}

	ep_add_basic_imp(r, s, p, q);
}

#endif

#if EP_ADD == PROJC || !defined(STRIP)

void ep_add_projc(ep_t r, const ep_t p, const ep_t q) {
	if (ep_is_infty(p)) {
		ep_copy(r, q);
		return;
	}

	if (ep_is_infty(q)) {
		ep_copy(r, p);
		return;
	}

	ep_add_projc_imp(r, p, q);
}

#endif

#if EP_ADD == JACOB || !defined(STRIP)

void ep_add_jacob(ep_t r, const ep_t p, const ep_t q) {
	if (ep_is_infty(p)) {
		ep_copy(r, q);
		return;
	}

	if (ep_is_infty(q)) {
		ep_copy(r, p);
		return;
	}

	ep_add_jacob_imp(r, p, q);
}

#endif

void ep_sub(ep_t r, const ep_t p, const ep_t q) {
	ep_t t;

	ep_null(t);

	if (p == q) {
		ep_set_infty(r);
		return;
	}

	RLC_TRY {
		ep_new(t);
		ep_neg(t, q);
		ep_add(r, p, t);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep_free(t);
	}
}

#if !defined(EP_SUPER) && FP_PRIME == 1536
//==================FOR G2================	//to be edtied and checked below 0119
// y^2 = x^3 + 1/4 x
static void ep_G2_add_projc_imp(ep_t r, const ep_t p, const ep_t q) {
	fp_t t0, t1, t2, t3, t4, t5;
								
	// if (q->coord == BASIC) {	
	// 	ep_add_projc_mix(r, p, q);
	// 	return;					
	// }							
								
	fp_null(t0);				
	fp_null(t1);				
	fp_null(t2);				
	fp_null(t3);				
	fp_null(t4);				
	fp_null(t5);				
								
	RLC_TRY {					
		fp_new(t0);			
		fp_new(t1);			
		fp_new(t2);			
		fp_new(t3);			
		fp_new(t4);			
		fp_new(t5);			
								
		fp_mul(t0, p->x, q->x);
		fp_mul(t1, p->y, q->y);
		fp_mul(t2, p->z, q->z);
		fp_add(t3, p->x, p->y);
		fp_add(t4, q->x, q->y);
		fp_mul(t3, t3, t4);	
		fp_add(t4, t0, t1);	
		fp_sub(t3, t3, t4);	
			/* Cost of 12M + 3m_a + 2_m3b + 23a. */
			fp_add(t4, p->x, p->z);		
			fp_add(t5, q->x, q->z);		
			fp_mul(t4, t4, t5);
			fp_add(t5, t0, t2);
			fp_sub(t4, t4, t5);
			fp_add(t5, p->y, p->z);		
			fp_add(r->x, q->y, q->z);		
			fp_mul(t5, t5, r->x);
			fp_add(r->x, t1, t2);
			fp_sub(t5, t5, r->x);
			// ep_curve_mul_a(r->z, t4);
			fp_hlv(r->z, t4); fp_hlv(r->z, r->z);
			// fp_dbl(r->x, t2);	
			// fp_add(r->x, r->x, t2);		
			// ep_curve_mul_b(r->x, r->x);	
			fp_zero(r->x);
			fp_add(r->z, r->x, r->z);		
			fp_sub(r->x, t1, r->z);		
			fp_add(r->z, t1, r->z);		
			fp_mul(r->y, r->x, r->z);		
			// fp_dbl(t1, t4);	
			// fp_add(t1, t1, t4);
			// ep_curve_mul_b(t4, t1);	
			fp_zero(t4);
			fp_dbl(t1, t0);	
			fp_add(t1, t1, t0);
			// ep_curve_mul_a(t2, t2);
			fp_hlv(t2, t2); fp_hlv(t2, t2);
			fp_add(t1, t1, t2);
			fp_sub(t2, t0, t2);
			// ep_curve_mul_a(t2, t2);	
			fp_hlv(t2, t2); fp_hlv(t2, t2);	
			fp_add(t4, t4, t2);
			fp_mul(t0, t1, t4);
			fp_add(r->y, r->y, t0);		
			fp_mul(t0, t5, t4);
			fp_mul(r->x, t3, r->x);		
			fp_sub(r->x, r->x, t0);		
			fp_mul(t0, t3, t1);
			fp_mul(r->z, t5, r->z);		
			fp_add(r->z, r->z, t0);		
		r->coord = PROJC;		
	} RLC_CATCH_ANY {			
		RLC_THROW(ERR_CAUGHT);	
	} RLC_FINALLY {				
		fp_free(t0);			
		fp_free(t1);			
		fp_free(t2);			
		fp_free(t3);			
		fp_free(t4);			
		fp_free(t5);			
	}							
}								

void ep_G2_add_projc(ep_t r, const ep_t p, const ep_t q) {
	if (ep_is_infty(p)) {
		ep_copy(r, q);
		return;
	}

	if (ep_is_infty(q)) {
		ep_copy(r, p);
		return;
	}

	ep_G2_add_projc_imp(r, p, q);		//to be edited
}

void ep_G2_sub(ep_t r, const ep_t p, const ep_t q) {
	ep_t t;

	ep_null(t);

	if (p == q) {
		ep_set_infty(r);
		return;
	}

	RLC_TRY {
		ep_new(t);
		ep_neg(t, q);
		ep_G2_add(r, p, t);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep_free(t);
	}
}

#endif
