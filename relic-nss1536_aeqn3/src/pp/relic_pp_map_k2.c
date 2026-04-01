/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2019 RELIC Authors
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
 * Implementation of pairing computation for curves with embedding degree 2.
 *
 * @ingroup pp
 */

#include "relic_core.h"
#include "relic_pp.h"
#include "relic_util.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Compute the Miller loop for pairings of type G_2 x G_1 over the bits of a
 * given parameter.
 *
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] p				- the first pairing argument in affine coordinates.
 * @param[in] q				- the second pairing argument in affine coordinates.
 * @param[in] m 			- the number of pairings to evaluate.
 * @param[in] a				- the loop parameter.
 */
static void pp_mil_k2(fp2_t r, ep_t *t, ep_t *p, ep_t *q, int m, bn_t a) {
	fp2_t l;
	ep_t *_q = RLC_ALLOCA(ep_t, m);
	int i, j;

	fp2_null(l);

	RLC_TRY {
		if (_q == NULL) {
			RLC_THROW(ERR_NO_MEMORY);
		}
		fp2_new(l);
		for (j = 0; j < m; j++) {
			ep_null(_q[j]);
			ep_new(_q[j]);
			ep_copy(t[j], p[j]);
			ep_neg(_q[j], q[j]);
		}

		fp2_zero(l);
		for (i = bn_bits(a) - 2; i >= 0; i--) {
			fp2_sqr(r, r);
			for (j = 0; j < m; j++) {
				pp_dbl_k2(l, t[j], t[j], _q[j]);
				fp2_mul(r, r, l);
				if (bn_get_bit(a, i)) {
					pp_add_k2(l, t[j], p[j], q[j]);
					fp2_mul(r, r, l);
				}
			}
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp2_free(l);
		for (j = 0; j < m; j++) {
			ep_free(_q[j]);
		}
		RLC_FREE(_q);
	}
}

#if !defined(EP_SUPER) && FP_PRIME == 1536

/**
 * Compute the Miller loop for pairings of type G_2 x G_1 over the bits of a
 * given parameter.
 * 
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] p				- the first pairing argument in affine coordinates.
 * @param[in] q				- the second pairing argument in affine coordinates.
 * @param[in] m 			- the number of pairings to evaluate.
 * @param[in] a				- the loop parameter.
 */
static void pp_mil_k2_NSS(fp2_t r, ep_t *t, ep_t *p, ep_t *q, int m, bn_t a) {
	fp2_t l;
	fp2_t ll[m][bn_bits(a)-1];
	fp2_t lll[m];
	fp2_t T;
	int i, j;

	fp2_null(l);
	fp2_null(T);
	for (i=0;i<m;i++){
		fp2_null(lll[i]);
		for (j=0;j<bn_bits(a)-1;j++)
			fp2_null(ll[i][j]);
	}
		
	RLC_TRY {
		fp2_new(l);
		fp2_new(T);
		for (j = 0; j < m; j++) {
			ep_copy(t[j], p[j]);
			fp2_new(lll[j]);
			for (i=0;i<bn_bits(a)-1;i++) fp2_new(ll[j][i]);
		}

		fp2_zero(l);
		for (i=0;i<m;i++){
			fp2_zero(lll[i]);
			for (j=0;j<bn_bits(a)-1;j++)
				fp2_zero(ll[i][j]);
		}
		
		for (i = bn_bits(a) - 2; i >= 0; i--) {
			fp2_sqr(r, r);
			for (j = 0; j < m; j++) {
				pp_dbl_k2_NSS(l, ll[j][i], t[j], t[j], q[j]);
				fp2_mul(r, r, l);
				if (bn_get_bit(a, i)) {
					pp_add_k2_NSS(l, lll[j], t[j], p[j], q[j]);
					fp2_mul(r, r, l);
				}
			}
		}

		fp2_copy(T,r);
		for (i = bn_bits(a) - 2; i >= 0; i--) {
			fp2_sqr(r, r);
			for (j = 0; j < m; j++) {
				fp2_mul(r, r, ll[j][i]);
				if (bn_get_bit(a, i)) {
					fp2_mul(r, r, lll[j]);
				}
			}
			if (bn_get_bit(a, i)) fp2_mul(r, r, T);
		}

	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp2_free(l);
	}
}


/**
 * Compute the Miller loop for pairings of type G_1 x G_2 over the bits of a
 * given parameter.
 *
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] p				- the first pairing argument in affine coordinates.
 * @param[in] q				- the second pairing argument in affine coordinates.
 * @param[in] a				- the loop parameter.
 */

static void pp_mil_k2_NSS_omega(fp2_t r, ep_t *t, ep_t *_t, ep_t *p, ep_t *q, int m, bn_t a) {
	fp2_t l, _l;
	int i, j;

	fp2_null(l);
	fp2_null(_l);

	RLC_TRY {
		fp2_new(l);
		fp2_new(_l);
		for (j = 0; j < m; j++) {
			ep_copy(t[j], p[j]);
			ep_copy(_t[j], q[j]);

			fp_hlv(_t[j]->z, _t[j]->z);
			fp_neg(_t[j]->z, _t[j]->z);
		}

		for (i = bn_bits(a) - 2; i >= 0; i--) {
			fp2_sqr(r, r);
			for (j = 0; j < m; j++) {
				pp_dbl_k2_projc_lazyr_NSS_omega_num(l, t[j], t[j], q[j]);
				pp_dbl_k2_projc_lazyr_NSS_omega_den(_l, _t[j], _t[j], p[j]);
				fp2_frb(_l, _l, 1);
				fp2_mul(r, r, l);
				fp2_mul(r, r, _l);
				if (bn_get_bit(a, i)) {
					pp_add_k2_projc_lazyr_NSS_omega_num(l, t[j], p[j], q[j]);
					pp_add_k2_projc_lazyr_NSS_omega_den(_l, _t[j], q[j], p[j]);
					fp2_frb(_l, _l, 1);
					fp2_mul(r, r, l);
					fp2_mul(r, r, _l);
				}
				
			}
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp2_free(l);
		fp2_free(_l);
		for (j = 0; j < m; j++) {
			ep_null(_q[j]);
		}
	}
}


#endif

/**
 * Compute the Miller loop for pairings of type G_1 x G_2 over the bits of a
 * given parameter.
 *
 * @param[out] r			- the result.
 * @param[out] t			- the resulting point.
 * @param[in] p				- the first pairing argument in affine coordinates.
 * @param[in] q				- the second pairing argument in affine coordinates.
 * @param[in] a				- the loop parameter.
 */
static void pp_mil_lit_k2(fp2_t r, ep_t *t, ep_t *p, ep_t *q, int m, bn_t a) {
	fp2_t l, _l;
	ep_t *_q = RLC_ALLOCA(ep_t, m);
	int i, j;

	fp2_null(l);
	fp2_null(_l);

	RLC_TRY {
		if (_q == NULL) {
			RLC_THROW(ERR_NO_MEMORY);
		}
		fp2_new(l);
		fp2_new(_l);
		for (j = 0; j < m; j++) {
			ep_null(_q[j]);
			ep_new(_q[j]);
			ep_copy(t[j], p[j]);
			ep_neg(_q[j], q[j]);
		}

		for (i = bn_bits(a) - 2; i >= 0; i--) {
			fp2_sqr(r, r);
			for (j = 0; j < m; j++) {
				pp_dbl_k2(l, t[j], t[j], _q[j]);
				fp_copy(_l[0], l[1]);
				fp_copy(_l[1], l[0]);
				fp2_mul(r, r, _l);
				if (bn_get_bit(a, i)) {
					pp_add_k2(l, t[j], p[j], q[j]);
					fp_copy(_l[0], l[1]);
					fp_copy(_l[1], l[0]);
					fp2_mul(r, r, _l);
				}
			}
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp2_free(l);
		fp2_free(_l);
		for (j = 0; j < m; j++) {
			ep_null(_q[j]);
		}
		RLC_FREE(_q);
	}
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

#if PP_MAP == TATEP || PP_MAP == OATEP || !defined(STRIP)

void pp_map_tatep_k2(fp2_t r, const ep_t p, const ep_t q) {
	ep_t _p[1], _q[1], t[1];
	bn_t n;

	ep_null(_p[0]);
	ep_null(_q[0]);
	ep_null(t[0]);
	bn_null(n);

	RLC_TRY {
		ep_new(t[0]);
		bn_new(n);

		ep_norm(_p[0], p);
		ep_norm(_q[0], q);
		ep_curve_get_ord(n);
		/* Since p has order n, we do not have to perform last iteration. */
		// bn_sub_dig(n, n, 1);
		fp2_set_dig(r, 1);

		if (!ep_is_infty(p) && !ep_is_infty(q)) {
			pp_mil_k2(r, t, _p, _q, 1, n);
			pp_exp_k2(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep_free(_p[0]);
		ep_free(_q[0]);
		ep_free(t[0]);
		bn_free(n);
	}
}

void pp_map_sim_tatep_k2(fp2_t r, const ep_t *p, const ep_t *q, int m) {
	ep_t *_p = RLC_ALLOCA(ep_t, m),
			*_q = RLC_ALLOCA(ep_t, m), *t = RLC_ALLOCA(ep_t, m);
	bn_t n;
	int i, j;

	bn_null(n);

	RLC_TRY {
		bn_new(n);
		if (_p == NULL || _q == NULL || t == NULL) {
			RLC_THROW(ERR_NO_MEMORY);
		}
		for (i = 0; i < m; i++) {
			ep_null(_p[i]);
			ep_null(_q[i]);
			ep_null(t[i]);
			ep_new(_p[i]);
			ep_new(_q[i]);
			ep_new(t[i]);
		}

		j = 0;
		for (i = 0; i < m; i++) {
			if (!ep_is_infty(p[i]) && !ep_is_infty(q[i])) {
				ep_norm(_p[j], p[i]);
				ep_norm(_q[j++], q[i]);
			}
		}

		ep_curve_get_ord(n);
		bn_sub_dig(n, n, 1);
		fp2_set_dig(r, 1);
		if (j > 0) {
			pp_mil_k2(r, t, _p, _q, j, n);
			pp_exp_k2(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		for (i = 0; i < m; i++) {
			ep_free(_p[i]);
			ep_free(_q[i]);
			ep_free(t[i]);
		}
		RLC_FREE(_p);
		RLC_FREE(_q);
		RLC_FREE(t);
	}
}

#if !defined(EP_SUPER) && FP_PRIME == 1536

void pp_map_tatep_k2_NSS(fp2_t r, const ep_t p, const ep_t q) {
	ep_t _p[1], _q[1], t[1];
	bn_t n;

	ep_null(_p[0]);
	ep_null(_q[0]);
	ep_null(t[0]);
	bn_null(n);

	RLC_TRY {
		ep_new(t[0]);
		bn_new(n);

		ep_norm(_p[0], p);
		ep_norm(_q[0], q);
		// ep_curve_get_ord(n);
		// /* Since p has order n, we do not have to perform last iteration. */
		// bn_sub_dig(n, n, 1);
		fp_prime_get_par(n);
		fp2_set_dig(r, 1);

		if (!ep_is_infty(p) && !ep_is_infty(q)) {
			pp_mil_k2_NSS(r, t, _p, _q, 1, n);
			pp_exp_k2(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep_free(_p[0]);
		ep_free(_q[0]);
		ep_free(t[0]);
		bn_free(n);
	}
}

void pp_map_sim_tatep_k2_NSS(fp2_t r, const ep_t *p, const ep_t *q, int m) {
	ep_t *_p = RLC_ALLOCA(ep_t, m),
			*_q = RLC_ALLOCA(ep_t, m), *t = RLC_ALLOCA(ep_t, m);
	bn_t n;
	int i, j;

	bn_null(n);

	RLC_TRY {
		bn_new(n);
		if (_p == NULL || _q == NULL || t == NULL) {
			RLC_THROW(ERR_NO_MEMORY);
		}
		for (i = 0; i < m; i++) {
			ep_null(_p[i]);
			ep_null(_q[i]);
			ep_null(t[i]);
			ep_new(_p[i]);
			ep_new(_q[i]);
			ep_new(t[i]);
		}

		j = 0;
		for (i = 0; i < m; i++) {
			if (!ep_is_infty(p[i]) && !ep_is_infty(q[i])) {
				ep_norm(_p[j], p[i]);
				ep_norm(_q[j++], q[i]);
			}
		}

		// ep_curve_get_ord(n);
		// bn_sub_dig(n, n, 1);
		fp_prime_get_par(n);
		fp2_set_dig(r, 1);
		if (j > 0) {
			pp_mil_k2_NSS(r, t, _p, _q, j, n);
			pp_exp_k2(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		for (i = 0; i < m; i++) {
			ep_free(_p[i]);
			ep_free(_q[i]);
			ep_free(t[i]);
		}
		RLC_FREE(_p);
		RLC_FREE(_q);
		RLC_FREE(t);
	}
}

#endif

#endif

#if PP_MAP == WEILP || !defined(STRIP)

void pp_map_weilp_k2(fp2_t r, const ep_t p, const ep_t q) {
	ep_t _p[1], _q[1], t0[1], t1[1];
	fp2_t r0, r1;
	bn_t n;

	ep_null(_p[0]);
	ep_null(_q[0]);
	ep_null(t0[0]);
	ep_null(t1[0]);
	fp2_null(r0);
	fp2_null(r1);
	bn_null(n);

	RLC_TRY {
		ep_new(_p[0]);
		ep_new(_q[0]);
		ep_new(t0[0]);
		ep_new(t1[0]);
		fp2_new(r0);
		fp2_new(r1);
		bn_new(n);

		ep_norm(_p[0], p);
		ep_norm(_q[0], q);
		ep_curve_get_ord(n);
		/* Since p has order n, we do not have to perform last iteration. */
		bn_sub_dig(n, n, 1);
		fp2_set_dig(r0, 1);
		fp2_set_dig(r1, 1);

		if (!ep_is_infty(_p[0]) && !ep_is_infty(_q[0])) {
			pp_mil_lit_k2(r0, t0, _p, _q, 1, n);
			pp_mil_k2(r1, t1, _q, _p, 1, n);
			fp2_inv(r1, r1);
			fp2_mul(r0, r0, r1);
			fp2_inv(r1, r0);
			fp2_inv_cyc(r0, r0);
		}
		fp2_mul(r, r0, r1);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep_free(_p[0]);
		ep_free(_q[0]);
		ep_free(t0[0]);
		ep_free(t1[0]);
		fp2_free(r0);
		fp2_free(r1);
		bn_free(n);
	}
}

void pp_map_sim_weilp_k2(fp2_t r, const ep_t *p, const ep_t *q, int m) {
	ep_t *_p = RLC_ALLOCA(ep_t, m),
			*_q = RLC_ALLOCA(ep_t, m),
			*t0 = RLC_ALLOCA(ep_t, m), *t1 = RLC_ALLOCA(ep_t, m);
	fp2_t r0, r1;
	bn_t n;
	int i, j;

	fp2_null(r0);
	fp2_null(r1);
	bn_null(r);

	RLC_TRY {
		fp2_new(r0);
		fp2_new(r1);
		bn_new(n);
		if (_p == NULL || _q == NULL || t0 == NULL || t1 == NULL) {
			RLC_THROW(ERR_NO_MEMORY);
		}
		for (i = 0; i < m; i++) {
			ep_null(_p[i]);
			ep_null(_q[i]);
			ep_null(t0[i]);
			ep_null(t1[i]);
			ep_new(_p[i]);
			ep_new(_q[i]);
			ep_new(t0[i]);
			ep_new(t1[i]);
		}

		j = 0;
		for (i = 0; i < m; i++) {
			if (!ep_is_infty(p[i]) && !ep_is_infty(q[i])) {
				ep_norm(_p[j], p[i]);
				ep_norm(_q[j++], q[i]);
			}
		}

		ep_curve_get_ord(n);
		bn_sub_dig(n, n, 1);
		fp2_set_dig(r0, 1);
		fp2_set_dig(r1, 1);

		if (j > 0) {
			pp_mil_lit_k2(r0, t0, _p, _q, j, n);
			pp_mil_k2(r1, t1, _q, _p, j, n);
			fp2_inv(r1, r1);
			fp2_mul(r0, r0, r1);
			fp2_inv(r1, r0);
			fp2_inv_cyc(r0, r0);
		}
		fp2_mul(r, r0, r1);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp2_free(r0);
		fp2_free(r1);
		bn_free(n);
		for (i = 0; i < m; i++) {
			ep_free(_p[i]);
			ep_free(_q[i]);
			ep_free(t0[i]);
			ep_free(t1[i]);
		}
		RLC_FREE(_p);
		RLC_FREE(_q);
		RLC_FREE(t0);
		RLC_FREE(t1);
	}
}

#endif

#if !defined(EP_SUPER) && FP_PRIME == 1536

void pp_map_omegap_k2_NSS(fp2_t r, const ep_t p, const ep_t q) {
	ep_t _p[1], _q[1], t[1], _t[1];
	bn_t n;

	ep_null(_p[0]);
	ep_null(_q[0]);
	ep_null(t[0]);
	ep_null(_t[0]);
	bn_null(n);

	RLC_TRY {
		ep_new(t[0]);
		ep_new(_t[0]);
		bn_new(n);

		ep_norm(_p[0], p);
		ep_norm(_q[0], q);
		// ep_curve_get_ord(n);
		// /* Since p has order n, we do not have to perform last iteration. */
		// bn_sub_dig(n, n, 1);	
		fp_prime_get_par(n);
		fp2_set_dig(r, 1);

		if (!ep_is_infty(p) && !ep_is_infty(q)) {
			pp_mil_k2_NSS_omega(r, t, _t, _p, _q, 1, n);
			fp2_conv_cyc(r,r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep_free(_p[0]);
		ep_free(_q[0]);
		ep_free(t[0]);
		ep_free(_t[0]);
		bn_free(n);
	}
}

void pp_map_sim_omegap_k2_NSS(fp2_t r, const ep_t *p, const ep_t *q, int m) {
	ep_t *_p = RLC_ALLOCA(ep_t, m),
			*_q = RLC_ALLOCA(ep_t, m), *t = RLC_ALLOCA(ep_t, m),  *_t = RLC_ALLOCA(ep_t, m);
	bn_t n;
	int i, j;

	bn_null(n);

	RLC_TRY {
		bn_new(n);
		if (_p == NULL || _q == NULL || t == NULL) {
			RLC_THROW(ERR_NO_MEMORY);
		}
		for (i = 0; i < m; i++) {
			ep_null(_p[i]);
			ep_null(_q[i]);
			ep_null(t[i]);
			ep_null(_t[i]);
			ep_new(_p[i]);
			ep_new(_q[i]);
			ep_new(t[i]);
			ep_new(_t[i]);
		}

		j = 0;
		for (i = 0; i < m; i++) {
			if (!ep_is_infty(p[i]) && !ep_is_infty(q[i])) {
				ep_norm(_p[j], p[i]);
				ep_norm(_q[j++], q[i]);
			}
		}
		fp_prime_get_par(n);
		fp2_set_dig(r, 1);
		if (j > 0) {
			pp_mil_k2_NSS_omega(r, t, _t, _p, _q, j, n);
			fp2_conv_cyc(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		for (i = 0; i < m; i++) {
			ep_free(_p[i]);
			ep_free(_q[i]);
			ep_free(t[i]);
			ep_free(_t[i]);
		}
		RLC_FREE(_p);
		RLC_FREE(_q);
		RLC_FREE(t);
		RLC_FREE(_t);
	}
}

#endif
