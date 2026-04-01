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
 * Benchmarks for Pairing-Based Cryptography, especially for NSS-1536(a=-3).
 *
 * @ingroup bench
 */

#include <stdio.h>

#include "relic.h"
#include "relic_bench.h"

#if !defined(EP_SUPER) && FP_PRIME == 1536
static void pairing_nss(void) {
	g1_t p[2];
	g2_t q[2];
	gt_t r;

	g1_new(p[0]);
	g2_new(q[0]);
	g1_new(p[1]);
	g2_new(q[1]);
	gt_new(r);

	BENCH_RUN("The Variant of the Tate Pairing on NSS1536(a=-3) via Miller's Algorithm") {
		g1_rand(p[0]);
		g2_rand(q[0]);
		BENCH_ADD(pp_map_tatep_k2_NSS(r, p[0], q[0]));
	}
	BENCH_END;

	BENCH_RUN("The Final Exponentiation of the Variant of the Tate Pairing on NSS1536(a=-3)") {
		gt_rand(r);
		BENCH_ADD(pp_exp_k2(r, r));
	}
	BENCH_END;

	BENCH_RUN("The Omega Pairing on NSS1536(a=-3) via Miller's Algorithm") {
		g1_rand(p[0]);
		g2_rand(q[0]);
		BENCH_ADD(pp_map_omegap_k2_NSS(r, p[0], q[0]));
	}
	BENCH_END;

	BENCH_RUN("The Final Exponentiation of the Omega Pairing on NSS1536(a=-3)") {
		gt_rand(r);
		BENCH_ADD(fp2_conv_cyc(r, r));
	}
	BENCH_END;

	g1_free(p[0]);
	g2_free(q[0]);
	g1_free(p[1]);
	g2_free(q[1]);
	gt_free(r);
}
#endif


int main(void) {
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}

	conf_print();
	util_banner("Benchmarks for the PC module:", 0);

	if (pc_param_set_any() != RLC_OK) {
		RLC_THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}

	pc_param_print();

#if !defined(EP_SUPER) && FP_PRIME == 1536
	util_banner("Pairing on Non-Supersingular Curve NSS1536(a=-3):", 0);
	pairing_nss();
#endif
	core_clean();
	return 0;
}
