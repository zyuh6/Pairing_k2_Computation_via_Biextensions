/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2010 RELIC Authors
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
 * Tests for ss1536.
 *
 * @ingroup test
 */

#include <stdio.h>

#include "relic.h"
#include "relic_test.h"

#if !defined(EP_SUPER) && FP_PRIME == 1536
static int pairing_nss(void) {
	int j, code = RLC_ERR;
	g1_t p[2];
	g2_t q[2];
	gt_t e1, e2;
	bn_t k, n;
	fp2_t t1, t2;

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

		TEST_CASE("pairing non-degeneracy is correct (NSS1536 Tate)") {
			g1_rand(p[0]);
			g2_rand(q[0]);
			pp_map_tatep_biext_nss(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) != RLC_EQ, end);
			g1_set_infty(p[0]);
			pp_map_tatep_biext_nss(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) == RLC_EQ, end);
			g1_rand(p[0]);
			g2_set_infty(q[0]);
			pp_map_tatep_biext_nss(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("cubical-based pairing is bilinear (NSS1536 Tate)") {
			g1_rand(p[0]);
			g2_rand(q[0]);
			bn_rand_mod(k, n);
			g1_mul(p[1], p[0], k);
			g2_mul(q[1], q[0], k);

			pp_map_tatep_biext_nss(e2, p[0], q[0]);
			pp_map_tatep_biext_nss(e1, p[0], q[1]);
			gt_exp(e2, e2, k);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
			pp_map_tatep_biext_nss(e1, p[1], q[0]);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
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

static int pairing_nss_omega(void) {
	int j, code = RLC_ERR;
	g1_t p[2];
	g2_t q[2];
	gt_t e1, e2;
	bn_t k, n;
	fp2_t t1, t2;

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

		TEST_CASE("pairing non-degeneracy is correct (NSS1536 Omega)") {
			g1_rand(p[0]);
			g2_rand(q[0]);
			pp_map_omegap_biext_nss(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) != RLC_EQ, end);
			g1_set_infty(p[0]);
			pp_map_omegap_biext_nss(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) == RLC_EQ, end);
			g1_rand(p[0]);
			g2_set_infty(q[0]);
			pp_map_omegap_biext_nss(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("cubical-based pairing is bilinear (NSS1536 Omega)") {
			g1_rand(p[0]);
			g2_rand(q[0]);
			bn_rand_mod(k, n);
			g1_mul(p[1], p[0], k);
			g2_mul(q[1], q[0], k);

			pp_map_omegap_biext_nss(e2, p[0], q[0]);
			pp_map_omegap_biext_nss(e1, p[0], q[1]);
			gt_exp(e2, e2, k);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
			pp_map_omegap_biext_nss(e1, p[1], q[0]);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
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
#elif defined(EP_SUPER) && FP_PRIME == 1536
static int pairing_ssmont(void) {
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

		TEST_CASE("pairing non-degeneracy is correct (SSmont1536)") {
			g1_rand(p[0]);
			g2_rand(q[0]);
			pp_map_tatep_biext_ssmont1536(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) != RLC_EQ, end);
			g1_set_infty(p[0]);
			pp_map_tatep_biext_ssmont1536(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) == RLC_EQ, end);
			g1_rand(p[0]);
			g2_set_infty(q[0]);
			pp_map_tatep_biext_ssmont1536(e1, p[0], q[0]);
			TEST_ASSERT(gt_cmp_dig(e1, 1) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("cubical-based pairing is bilinear (SSmont1536)") {
			g1_rand(p[0]);
			g2_rand(q[0]);
			bn_rand_mod(k, n);
			g1_mul(p[1], p[0], k);
			g2_mul(q[1], q[0], k);

			ep_ss_to_ssmont(p[0], p[0]);
			ep_ss_to_ssmont(q[0], q[0]);
			ep_ss_to_ssmont(p[1], p[1]);
			ep_ss_to_ssmont(q[1], q[1]);

			pp_map_tatep_biext_ssmont1536(e2, p[0], q[0]);
			pp_map_tatep_biext_ssmont1536(e1, p[0], q[1]);
			gt_exp(e2, e2, k);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
			pp_map_tatep_biext_ssmont1536(e2, p[1], q[0]);
			TEST_ASSERT(gt_cmp(e1, e2) == RLC_EQ, end);
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
#endif

int main(void) {
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}

	util_banner("Tests for the PC module on SSmont1536 and NSS1536:", 0);

	if (pc_param_set_any() != RLC_OK) {
		RLC_THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}

	pc_param_print();
    ep_param_print();
#if !defined(EP_SUPER) && FP_PRIME == 1536
	pairing_nss();
	pairing_nss_omega();
#elif defined(EP_SUPER) && FP_PRIME == 1536
	pairing_ssmont();
#endif

	util_banner("All tests have passed.\n", 0);

	core_clean();
	return 0;

}