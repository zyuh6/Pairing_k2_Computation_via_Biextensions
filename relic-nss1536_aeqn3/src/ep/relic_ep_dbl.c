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
 * Implementation of the point doubling on prime elliptic curves.
 *
 * @ingroup ep
 */

#include "relic_core.h"
#include "relic_ep_dbl_tmpl.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#if EP_ADD == BASIC || !defined(STRIP)

/**
 * Doubles a point represented in affine coordinates on an ordinary prime
 * elliptic curve.
 *
 * @param[out] r			- the result.
 * @param[out] s			- the slope.
 * @param[in] p				- the point to double.
 */
TMPL_DBL_BASIC_IMP(ep, fp);

#endif /* EP_ADD == BASIC */

#if EP_ADD == PROJC || !defined(STRIP)

/**
 * Doubles a point represented in projective coordinates on an ordinary prime
 * elliptic curve.
 *
 * @param r					- the result.
 * @param p					- the point to double.
 */
TMPL_DBL_PROJC_IMP(ep, fp);

#endif /* EP_ADD == PROJC */

#if EP_ADD == JACOB || !defined(STRIP)

/**
 * Doubles a point represented in Jacobian coordinates on an ordinary prime
 * elliptic curve.
 *
 * @param r					- the result.
 * @param p					- the point to double.
 */
TMPL_DBL_JACOB_IMP(ep, fp);

#endif /* EP_ADD == JACOB */

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

#if EP_ADD == BASIC || !defined(STRIP)

void ep_dbl_basic(ep_t r, const ep_t p) {
	if (ep_is_infty(p)) {
		ep_set_infty(r);
		return;
	}
	ep_dbl_basic_imp(r, NULL, p);
}

void ep_dbl_slp_basic(ep_t r, fp_t s, const ep_t p) {
	if (ep_is_infty(p)) {
		ep_set_infty(r);
		return;
	}

	ep_dbl_basic_imp(r, s, p);
}

#endif

#if EP_ADD == PROJC || !defined(STRIP)

void ep_dbl_projc(ep_t r, const ep_t p) {
	if (ep_is_infty(p)) {
		ep_set_infty(r);
		return;
	}

	ep_dbl_projc_imp(r, p);
}

#endif

#if EP_ADD == JACOB || !defined(STRIP)

void ep_dbl_jacob(ep_t r, const ep_t p) {
	if (ep_is_infty(p)) {
		ep_set_infty(r);
		return;
	}

	ep_dbl_jacob_imp(r, p);
}

#endif

#if !defined(EP_SUPER) && FP_PRIME == 1536
//==================FOR G2================	//to be edtied and checked below 0119
// y^2 = x^3 - 3/4 x
static void ep_G2_dbl_projc_imp(ep_t r, const ep_t p) {
	fp_t t0, t1, t2, t3, t4, t5;
								
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
								
			
			fp_sqr(t0, p->x);	
			fp_sqr(t1, p->y);	
			fp_mul(t3, p->x, p->y);				
			fp_dbl(t3, t3);	
			fp_mul(t4, p->y, p->z);				
											
				/* Common cost of 8M + 3S + 3m_a + 2m_3b + 15a. */	
				// if (p->coord == BASIC) {			
				// 	/* Save 1S + 1m_b + 1m_a if z1 = 1. */			
				// 	fp_dbl(r->y, ep_curve_get_b());
				// 	fp_add(r->y, r->y, ep_curve_get_b());			
				// 	fp_copy(t2, ep_curve_get_a());
				// } else {		
					fp_sqr(t2, p->z);				
					fp_dbl(t5, t2);				
					fp_add(t5, t5, t2);			
					// ep_curve_mul_b(r->y, t5);
					fp_zero(r->y);
					// ep_curve_mul_a(t2, t2);
					fp_dbl(t5, t2); fp_add(t2, t5, t2); fp_neg(t2, t2); fp_hlv(t2, t2); fp_hlv(t2, t2); 
				// }				
				fp_mul(r->z, p->x, p->z);			
				fp_dbl(r->z, r->z);				
				// ep_curve_mul_a(r->x, r->z);
				fp_dbl(r->x, r->z); fp_add(r->x, r->x, r->z); fp_neg(r->x, r->x); fp_hlv(r->x, r->x); fp_hlv(r->x, r->x);
				fp_add(r->y, r->x, r->y);			
				fp_sub(r->x, t1, r->y);			
				fp_add(r->y, t1, r->y);			
				fp_mul(r->y, r->x, r->y);			
				fp_mul(r->x, t3, r->x);			
				// fp_dbl(t5, r->z);
				// fp_add(t5, t5, r->z);				
				// ep_curve_mul_b(r->z, t5);
				fp_zero(r->z);			
				fp_sub(t3, t0, t2);				
				// ep_curve_mul_a(t3, t3);	
				fp_dbl(t5, t3); fp_add(t3, t5, t3); fp_neg(t3, t3); fp_hlv(t3, t3); fp_hlv(t3, t3); 		
				fp_add(t3, t3, r->z);				
				fp_dbl(r->z, t0);
				fp_add(t0, t0, r->z);				
				fp_add(t0, t0, t2);				
			
			/* Common part with renamed variables. */
			fp_mul(t0, t0, t3);
			fp_add(r->y, r->y, t0);				
			fp_dbl(t2, t4);	
			fp_mul(t0, t2, t3);
			fp_sub(r->x, r->x, t0);				
			fp_mul(r->z, t2, t1);
			fp_dbl(r->z, r->z);
			fp_dbl(r->z, r->z);
				
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

void ep_G2_dbl_projc(ep_t r, const ep_t p) {
	if (ep_is_infty(p)) {
		ep_set_infty(r);
		return;
	}

	ep_G2_dbl_projc_imp(r, p);
}

#endif

