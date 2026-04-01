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
 * Tests for pairings defined over prime elliptic curves.
 *
 * @ingroup test
 */

#include <stdio.h>

#include "relic.h"
#include "relic_test.h"
#include "relic_bench.h"

static int pairing_tate(void) {
	int j, code = RLC_ERR;
	bn_t k, n;
	ep_t p[2], q[2], r;
	fp2_t e1, e2;

	bn_null(k);
	bn_null(n);
	ep_null(r);
	fp2_null(e1);
	fp2_null(e2);

	RLC_TRY {
		bn_new(n);
		bn_new(k);
		ep_new(r);
		fp2_new(e1);
		fp2_new(e2);

		for (j = 0; j < 2; j++) {
			ep_null(p[j]);
			ep_null(q[j]);
			ep_new(p[j]);
			ep_new(q[j]);
		}

		ep_curve_get_ord(n);

		TEST_CASE("Tate pairing non-degeneracy is correct (NSS-1536(a=-3))") {
			ep_rand(p[0]);
			// ep_rand(q[0]);
			ep_G2_rand(q[0]);
			pp_map_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) != RLC_EQ, end);
			ep_set_infty(p[0]);
			pp_map_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) == RLC_EQ, end);
			ep_rand(p[0]);
			ep_set_infty(q[0]);
			pp_map_k2(e1, p[0], q[0]);
			TEST_ASSERT(fp2_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("Tate pairing is bilinear (NSS-1536(a=-3))") {
			ep_rand(p[0]);
			ep_G2_rand(q[0]);
			bn_rand_mod(k, n);
			ep_G2_mul(r, q[0], k);
			pp_map_k2(e1, p[0], r);
			pp_map_k2(e2, p[0], q[0]);
			fp2_exp(e2, e2, k);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_mul(p[0], p[0], k);
			pp_map_k2(e2, p[0], q[0]);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			ep_dbl(p[0], p[0]);
			pp_map_k2(e2, p[0], q[0]);
			fp2_sqr(e1, e1);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
			// ep_dbl(q[0], q[0]);
			ep_G2_dbl(q[0], q[0]);
			pp_map_k2(e2, p[0], q[0]);
			fp2_sqr(e1, e1);
			TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

        TEST_CASE("Tate multi-pairing is correct (NSS-1536(a=-3))") {
            ep_rand(p[i % 2]);
            ep_rand(q[i % 2]);
            pp_map_k2(e1, p[i % 2], q[i % 2]);
            ep_rand(p[1 - (i % 2)]);
            ep_set_infty(q[1 - (i % 2)]);
            pp_map_sim_k2(e2, p, q, 2);
            TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
            ep_set_infty(p[1 - (i % 2)]);
            ep_rand(q[1 - (i % 2)]);
            pp_map_sim_k2(e2, p, q, 2);
            TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
            ep_set_infty(q[i % 2]);
            pp_map_sim_k2(e2, p, q, 2);
            TEST_ASSERT(fp2_cmp_dig(e2, 1) == RLC_EQ, end);
            ep_rand(p[0]);
            ep_rand(q[0]);
            pp_map_k2(e1, p[0], q[0]);
            ep_rand(p[1]);
            ep_rand(q[1]);
            pp_map_k2(e2, p[1], q[1]);
            fp2_mul(e1, e1, e2);
            pp_map_sim_k2(e2, p, q, 2);
            TEST_ASSERT(fp2_cmp(e1, e2) == RLC_EQ, end);
        } TEST_END;

	}
	RLC_CATCH_ANY {
		util_print("FATAL ERROR!\n");
		RLC_ERROR(end);
	}
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	ep_free(r);
	fp2_free(e1);
	fp2_free(e2);

	for (j = 0; j < 2; j++) {
		ep_free(p[j]);
		ep2_free(q[j]);
	}

    return code;
}

static int pairing_omega(void) {
	int j, code = RLC_ERR;
	g1_t p[2];
	g2_t q[2];
	gt_t e1, e2;
	bn_t k, n;

	gt_null(e1);
	gt_null(e2);
	bn_null(k);
	bn_null(n);

	RLC_TRY {
		gt_new(e1);
		gt_new(e2);
		bn_new(k);
		bn_new(n);

		for (j = 0; j < 2; j++) {
			g1_null(p[j]);
			g2_null(q[j]);
			g1_new(p[j]);
			g2_new(q[j]);
		}

		pc_get_ord(n);

		TEST_CASE("Omega pairing non-degeneracy is correct (NSS-1536(a=-3))") {
			g1_rand(p[0]);
			g2_rand(q[0]);
			pp_map_omegap_k2_NSS(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) != RLC_EQ, end);
			g1_set_infty(p[0]);
			pp_map_omegap_k2_NSS(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) == RLC_EQ, end);
			g1_rand(p[0]);
			g2_set_infty(q[0]);
			pp_map_omegap_k2_NSS(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("Omega pairing is bilinear (NSS-1536(a=-3))") {
			g1_rand(p[0]);
			g2_rand(q[0]);
			g2_mul(q[1], q[0], k);
			pp_map_omegap_k2_NSS(e1, p[0], q[1]);
			pp_map_omegap_k2_NSS(e2, p[0], q[0]);
			gt_exp(e2, e2, k);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
			g1_mul(p[0], p[0], k);
			pp_map_omegap_k2_NSS(e2, p[0], q[0]);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
			g1_dbl(p[0], p[0]);
			pp_map_omegap_k2_NSS(e2, p[0], q[0]);
			gt_sqr(e1, e1);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
			g2_dbl(q[0], q[0]);
			pp_map_omegap_k2_NSS(e2, p[0], q[0]);
			gt_sqr(e1, e1);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("Omega multi-pairing is correct (NSS-1536(a=-3))") {
			g1_rand(p[i % 2]);
			g2_rand(q[i % 2]);
			pp_map_omegap_k2_NSS(e1, p[i % 2], q[i % 2]);
			g1_rand(p[1 - (i % 2)]);
			g2_set_infty(q[1 - (i % 2)]);
			pp_map_sim_omegap_k2_NSS(e2, p, q, 2);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
			g1_set_infty(p[1 - (i % 2)]);
			g2_rand(q[1 - (i % 2)]);
			pp_map_sim_omegap_k2_NSS(e2, p, q, 2);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
			g2_set_infty(q[i % 2]);
			pp_map_sim_omegap_k2_NSS(e2, p, q, 2);
			TEST_ASSERT(gt_cmp_dig(e2, 1) == RLC_EQ, end);
			g1_rand(p[0]);
			g2_rand(q[0]);
			pp_map_omegap_k2_NSS(e1, p[0], q[0]);
			g1_rand(p[1]);
			g2_rand(q[1]);
			pp_map_omegap_k2_NSS(e2, p[1], q[1]);
			gt_mul(e1, e1, e2);
			pp_map_sim_omegap_k2_NSS(e2, p, q, 2);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
			g1_neg(p[1], p[0]);
			g2_copy(q[1], q[0]);
			pp_map_sim_omegap_k2_NSS(e1, p, q, 2);
			TEST_ASSERT(gt_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;
	}
	RLC_CATCH_ANY {
		util_print("FATAL ERROR!\n");
		RLC_ERROR(end);
	}
	code = RLC_OK;
  end:
	gt_free(e1);
	gt_free(e2);
	bn_free(k);
	bn_free(n);
	for (j = 0; j < 2; j++) {
		g1_free(p[j]);
		g2_free(q[j]);
	}
	return code;
}

int main(void) {
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}

	util_banner("Tests for the PP module", 0);

	if (ep_param_set_any_pairf() == RLC_ERR) {
		RLC_THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}

	ep_param_print();

	util_banner("Pairing", 1);

	if (ep_curve_embed() == 2) {

		if (pairing_tate() != RLC_OK) {
			core_clean();
			return 1;
		}
		if (pairing_omega() != RLC_OK) {
			core_clean();
			return 1;
		}
	}

	util_banner("All tests have passed.\n", 0);

	core_clean();
	return 0;
}

