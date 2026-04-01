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
#include "relic_fp_low.h"
#include "relic_fpx_low.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

//P=[px,py]  twist(Q)=[-2*qx,-2*u*qy],  phihat(twist(Q))=[2*qx, -2*u*qy/beta]   //compute P+twist(Q) and P+phihat(twist(Q)) on E(Fp2)  //a=1
//5m + 2s on Fp
void ep2_add_PQ_PphihatQ_nss(fp2_t XPQ, fp_t ZPQ, fp2_t XPphihatQ, fp_t ZPphihatQ, const ep_t P, const ep_t Q){
    fp_t t0, t1, t2, t3;
    fp_null(t0);
    fp_null(t1);
    fp_null(t2);
    fp_null(t3);
    
    fp_new(t0);
    fp_new(t1);
    fp_new(t2);
    fp_new(t3);
    
    fp_dbl(t3, Q->x);
    fp_sub(t0, t3, P->x);
    fp_mul(t1, P->x, t3);
    fp_sub_dig(t2, t1, 1);
    // fp_mul(XPphihatQ[0], t0, t2);
    fp_mul(XPQ[0], t0, t2);
    // fp_sqr(ZPQ, t0);
    fp_sqr(ZPphihatQ, t0);

    fp_add(t0, t3, P->x);
    fp_add_dig(t2, t1, 1);
    // fp_mul(XPQ[0], t0, t2);
    fp_mul(XPphihatQ[0], t0, t2);
    // fp_sqr(ZPphihatQ, t0);
    fp_sqr(ZPQ, t0);

    fp_mul(t1, P->y, Q->y);
    // fp_dbl(XPphihatQ[1], t1);
    // fp_dbl(XPphihatQ[1], XPphihatQ[1]);
    fp_dbl(XPQ[1], t1);
    fp_dbl(XPQ[1], XPQ[1]);

    fp_neg(XPphihatQ[1], XPQ[1]);
    fp_mul(XPphihatQ[1], XPphihatQ[1], ep_curve_get_beta());

    fp_free(t0);
    fp_free(t1);
    fp_free(t2);
    fp_free(t3);
}

//a=1 all below
void ep_biex_doubleadd_saved_nss(fp_t r0, fp_t r1, fp_t drx, fp_t drz, fp2_t arx, fp2_t arz, const fp_t p1x, const fp_t p1z, const fp2_t p2x, const fp2_t p2z, const fp_t qx){
    fp_t t0, x4, z4, c24z4;
    fp2_t t2, t3, x5, z5, tt0, tt1;
    RLC_TRY {
        fp_add(t0, p1x, p1z);
        fp_copy(r0, t0);
        fp_sub(r1, p1x, p1z);
        fp_sqr(x4, t0);
        fp2_sub(t2, p2x, p2z);
        fp2_add(x5, p2x, p2z);

        for(int i=0;i<2;i++) fp_mul(t2[i], t0, t2[i]);
        fp_sqr(z4, r1);
        fp_dbl(c24z4, z4);
        for(int i=0;i<2;i++) fp_mul(t3[i], r1, x5[i]);
        fp_sub(t0, x4, z4);
        fp_mul(x4, x4, c24z4);
        
        fp2_sub(z5, t2, t3);
        fp_add(z4, t0, c24z4);
        fp2_add(x5, t2, t3);
        fp_mul(z4, z4, t0);

        fp2_sqr(z5, z5);
        fp2_sqr(x5, x5);

        // fp_neg(qx, qx);
        for(int i=0;i<2;i++) fp_mul(z5[i], qx, z5[i]);
        fp2_copy(arz, z5);
        fp2_copy(arx, x5);
        fp_copy(drz, z4);
        fp_copy(drx, x4);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
	}
}

//just add_dig and sub_dig while use saved
void ep_biex_doubleadd_nss_first(/*fp_t r0, fp_t r1,*/ fp_t drx, fp_t drz, fp2_t arx, fp2_t arz, const fp_t p1x, /*const fp_t p1z,*/ const fp2_t p2x, const fp_t p2z, const fp_t qx){
    fp_t t0, x4, z4, c24z4;
    fp2_t t2, t3, x5, z5, tt0, tt1;
    RLC_TRY {
        // fp_add(t0, p1x, p1z);
        // fp_sub(t1, p1x, p1z);
        // fp_add_dig(t0, p1x, 1);
        // fp_sub_dig(t1, p1x, 1);
        // fp_sqr(x4, t0);
        fp_sqr(z4, p1x);
        fp_add_dig(z4, z4, 1);
        fp_dbl(t0, p1x);
        fp_add(x4, z4, t0);
        // fp2_sub(t2, p2x, p2z);
        // fp2_add(x5, p2x, p2z);

        // for(int i=0;i<2;i++) fp_mul(t2[i], t0, t2[i]);
        // fp_sqr(z4, t1);
        fp_sub(z4, z4, t0);

        fp_dbl(c24z4, z4);
        fp_sub(t0, x4, z4);
        fp_mul(x4, x4, c24z4);
        
        fp_mul(z5[0], p1x, p2z);
        fp_sub(z5[0], p2x[0], z5[0]);
        fp_copy(z5[1], p2x[1]);

        fp_add(z4, t0, c24z4);
        for(int i=0;i<2;i++) fp_mul(x5[i], p1x, p2x[i]);
        fp_sub(x5[0], x5[0], p2z);

        fp_mul(z4, z4, t0);

        fp2_sqr(z5, z5);
        fp2_sqr(x5, x5);

        for(int i=0;i<2;i++) fp_mul(z5[i], qx, z5[i]);

        fp2_copy(arz, z5);
        fp2_copy(arx, x5);
        fp_copy(drz, z4);
        fp_copy(drx, x4);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
	}
}

void ep_biex_doubleadd_saved_nss_mid(fp2_t mid, fp_t r0, fp_t r1, fp_t drx, fp_t drz, fp2_t arx, fp2_t arz, const fp_t p1x, const fp_t p1z, const fp2_t p2x, const fp2_t p2z, const fp_t qx){
    fp_t t0, x4, z4, c24z4;
    fp2_t t2, t3, x5, z5, tt0, tt1;
    RLC_TRY {
        fp_add(t0, p1x, p1z);
        fp_copy(r0, t0);
        fp_sub(r1, p1x, p1z);
        fp_sqr(x4, t0);
        fp2_sub(t2, p2x, p2z);
        fp2_add(x5, p2x, p2z);

        for(int i=0;i<2;i++) fp_mul(t2[i], t0, t2[i]);
        fp_sqr(z4, r1);
        fp_dbl(c24z4, z4);
        for(int i=0;i<2;i++) fp_mul(t3[i], r1, x5[i]);
        fp_sub(t0, x4, z4);
        fp_mul(x4, x4, c24z4);
        
        fp2_sub(z5, t2, t3);
        fp_add(z4, t0, c24z4);
        fp2_add(x5, t2, t3);
        fp_mul(z4, z4, t0);

        fp2_copy(mid, z5);//new for test
        fp2_sqr(z5, z5);
        fp2_sqr(x5, x5);

        // fp_neg(qx, qx);
        for(int i=0;i<2;i++) fp_mul(z5[i], qx, z5[i]);
        fp2_copy(arz, z5);
        // fp2_copy(mid, arz);//new for test
        fp2_copy(arx, x5);
        fp_copy(drz, z4);
        fp_copy(drx, x4);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
	}
}

void ep_biex_doubleadd_saved_nss_last(fp_t r0, fp_t r1, fp_t drx, fp_t drz, fp2_t arzlast, const fp_t p1x, const fp_t p1z, const fp2_t p2x, const fp2_t p2z, const fp_t qx){
    fp_t t0, x4, z4, c24z4;
    fp2_t t2, t3, x5, z5, tt0, tt1;
    RLC_TRY {
        fp_add(t0, p1x, p1z);
        fp_copy(r0, t0);
        fp_sub(r1, p1x, p1z);
        fp_sqr(x4, t0);
        fp2_sub(t2, p2x, p2z);
        fp2_add(x5, p2x, p2z);

        for(int i=0;i<2;i++) fp_mul(t2[i], t0, t2[i]);
        fp_sqr(z4, r1);
        fp_dbl(c24z4, z4);
        for(int i=0;i<2;i++) fp_mul(t3[i], r1, x5[i]);
        fp_sub(t0, x4, z4);
        fp_mul(x4, x4, c24z4);
        
        fp2_sub(z5, t2, t3);    //lzyr
        fp_add(z4, t0, c24z4);
        // fp2_add(x5, t2, t3);
        fp_mul(z4, z4, t0);

        // fp2_sqr(z5, z5);
        // fp2_sqr(x5, x5);

        // fp_neg(qx, qx);
        // for(int i=0;i<2;i++) fp_mul(z5[i], qx, z5[i]);
        // fp2_copy(arz, z5);
        fp2_copy(arzlast, z5);
        // fp2_copy(arx, x5);
        fp_copy(drz, z4);
        fp_copy(drx, x4);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
	}
}

//a=1 all below
//6m+2S
void ep_biex_diff_add_use_saved_nss(fp2_t rx, fp2_t rz, const fp_t r0, const fp_t r1, const fp2_t p2x, const fp2_t p2z, const fp_t XphihatQ){
    fp2_t t0, t1, t2, X5, Z5;

    RLC_TRY {
        fp2_sub(t2, p2x, p2z);
        fp2_add(X5, p2x, p2z);
        for(int i=0;i<2;i++) fp_mul(t0[i], r0, t2[i]);
        for(int i=0;i<2;i++) fp_mul(t1[i], r1, X5[i]);

        fp2_sub(Z5, t0, t1);
        fp2_add(X5, t0, t1);

        fp2_sqr(Z5, Z5);
        for(int i=0;i<2;i++) fp_mul(rz[i], XphihatQ, Z5[i]);
        fp2_sqr(rx, X5);
    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
	}
}
//6m+2S
void ep_biex_diff_add_nss_first(fp2_t rx, fp2_t rz, const fp_t p1x, /*const fp_t p1z,*/ const fp2_t p2x, const fp2_t p2z, const fp_t XphihatQ){
    fp2_t x5, z5;
    RLC_TRY {
        for(int i=0;i<2;i++) fp_mul(z5[i], p2z[i], p1x);
        fp2_sub(z5, p2x, z5);

        for(int i=0;i<2;i++) fp_mul(x5[i], p2x[i], p1x);
        fp2_sub(x5, x5, p2z);

        fp2_sqr(z5, z5);
        for(int i=0;i<2;i++) fp_mul(rz[i], XphihatQ, z5[i]);
        fp2_sqr(rx, x5);
    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
	}
}

// mid and last
void ep_biex_diff_add_use_saved_nss_mid(fp2_t mid, fp2_t rx, fp2_t rz, const fp_t r0, const fp_t r1, const fp2_t p2x, const fp2_t p2z, const fp_t XphihatQ){
    fp2_t t0, t1, t2, X5, Z5;

    RLC_TRY {
        fp2_sub(t2, p2x, p2z);
        fp2_add(X5, p2x, p2z);
        for(int i=0;i<2;i++) fp_mul(t0[i], r0, t2[i]);
        for(int i=0;i<2;i++) fp_mul(t1[i], r1, X5[i]);

        fp2_sub(Z5, t0, t1);
        fp2_add(X5, t0, t1);

        fp2_copy(mid, Z5);

        fp2_sqr(Z5, Z5);
        for(int i=0;i<2;i++) fp_mul(rz[i], XphihatQ, Z5[i]);
        fp2_sqr(rx, X5);
    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
	}
}
//4m 
void ep_biex_diff_add_use_saved_nss_last(/*fp2_t rx,*/ fp2_t rz, const fp_t r0, const fp_t r1, const fp2_t p2x, const fp2_t p2z, const fp_t XphihatQ){
    fp2_t t0, t1, t2, X5, Z5;

    RLC_TRY {
        fp2_sub(t2, p2x, p2z);
        fp2_add(X5, p2x, p2z);
        for(int i=0;i<2;i++) fp_mul(t0[i], r0, t2[i]);
        for(int i=0;i<2;i++) fp_mul(t1[i], r1, X5[i]);

        fp2_sub(Z5, t0, t1);    //lzyr
        fp2_copy(rz, Z5);
    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
	}
}
//line functions two
//compute line_{[2^41]P, -P}(Q)     //a=1
//15m+1s
void ep2_line_function_tate_nss(fp2_t resultQ, fp2_t resultphihatQ, const fp_t Xn2P, const fp_t Zn2P, const fp_t Xn1P, const fp_t Zn1P, const fp_t px, const fp_t py, const fp_t qx, const fp_t qy){
    fp_t k0, k1, k2, k3, k4, k5, k6, rQr, rQi, rphihatQr, rphihatQi;
    RLC_TRY {
        fp_mul(k0, px, Zn2P);
        fp_add(k1, Xn2P, k0);
        fp_mul(k2, ep_curve_get_beta(), py);
        fp_mul(k3, k2, Zn2P);
        fp_mul(k4, Zn1P, k3);
        fp_dbl(k4, k4);
        fp_mul(k3, k3, k4);
        fp_mul(k5, k1, k4);
        fp_neg(k5, k5);
        fp_sqr(k1, k1);
        fp_mul(k1, k1, Xn1P);
        fp_sub(k0, Xn2P, k0);
        fp_mul(k4, px, Xn2P);
        fp_sub(k4, Zn2P, k4);
        fp_mul(k6, k4, k0);
        fp_mul(k6, k6, Zn1P);
        fp_sub(k6, k6, k1);
        fp_sub(k6, k6, k3);
        fp_mul(rQr, k2, k5);
        fp_mul(k5, k5, qy);
        fp_dbl(k5, k5);
        fp_neg(rQi, k5);
        fp_mul(rphihatQi, k5, ep_curve_get_beta());
        fp_dbl(k1, qx);
        fp_neg(k1, k1);
        fp_add(k3, k1, px);
        fp_sub(k4, k1, px);
        fp_mul(k3, k3, k6);
        fp_mul(k4, k4, k6);
        fp_add(rphihatQr, rQr, k4);
        fp_sub(rQr, rQr, k3);

        fp_copy(resultQ[0], rQr);
        fp_copy(resultQ[1], rQi);

        fp_copy(resultphihatQ[0], rphihatQr);
        fp_copy(resultphihatQ[1], rphihatQi);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
	}
}

void ep2_doubleadd_ladder_tate_nss(fp2_t result, const bn_t n, const ep_t P, const ep_t Q, const fp2_t XQP, const fp_t ZQP, const fp2_t XphihatQP, const fp_t ZphihatQP){
    fp_t Xn1P, Zn1P, Xn2P, Zn2P, XQ, XphihatQ, XP;
    fp2_t XQn1P, ZQn1P, ZQn2P, lineQ, linephihatQ;
    size_t bits;
    bits = bn_bits(n);
    fp_t R0[bits-2], R1[bits-2];
    int i;

    for (i=0;i<bits-2;i++){
        fp_null(R0[i]);
        fp_null(R1[i]);
    }

    RLC_TRY {
        fp_dbl(XphihatQ, Q->x);
        fp_neg(XQ, XphihatQ);
        fp_copy(XP, P->x);
        
        
        //first
        ep_biex_doubleadd_nss_first(Xn1P, Zn1P, XQn1P, ZQn1P, XP, XQP, ZQP, XQ);

        for (i = bits - 3; i >= bits - 23 + 1; i--) {       //n2
            ep_biex_doubleadd_saved_nss(R0[bits-2-1-i], R1[bits-2-1-i], Xn1P, Zn1P, XQn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
        }
        ep_biex_doubleadd_saved_nss_mid(result, R0[bits-2-1-i], R1[bits-2-1-i], Xn1P, Zn1P, XQn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
        // fp2_copy(result, ZQn1P);
        i--;
        fp_copy(Zn2P, Zn1P);
        fp_copy(Xn2P, Xn1P);
        for(;i>=0 + 1;i--){     //n1
            ep_biex_doubleadd_saved_nss(R0[bits-2-1-i], R1[bits-2-1-i], Xn1P, Zn1P, XQn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
        }
        ep_biex_doubleadd_saved_nss_last(R0[bits-2-1-i], R1[bits-2-1-i], Xn1P, Zn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
        
        ep2_line_function_tate_nss(lineQ, linephihatQ, Xn1P, Zn1P, Xn2P, Zn2P, P->x, P->y, Q->x, Q->y);
        fp2_mul(result, result, ZQn1P);
        fp2_mul(result, result, lineQ);  //millerQ
        fp2_sqr(result, result);

        fp2_mul(XQn1P, XphihatQP, result);
        for(i=0; i<2; i++) fp_mul(ZQn1P[i], ZphihatQP, result[i]);

        ep_biex_diff_add_nss_first(XQn1P, ZQn1P, XP, XQn1P, ZQn1P, XphihatQ);
        for (i = bits - 3; i >= bits - 23 + 1; i--) {       //n2
            ep_biex_diff_add_use_saved_nss(XQn1P, ZQn1P, R0[bits-2-1-i], R1[bits-2-1-i], XQn1P, ZQn1P, XphihatQ);
        }
        ep_biex_diff_add_use_saved_nss_mid(lineQ, XQn1P, ZQn1P, R0[bits-2-1-i], R1[bits-2-1-i], XQn1P, ZQn1P, XphihatQ);    //borrow |lineQ| for the midresult
        i--;

        for(;i>=0 + 1;i--){     //n1
            ep_biex_diff_add_use_saved_nss(XQn1P, ZQn1P, R0[bits-2-1-i], R1[bits-2-1-i], XQn1P, ZQn1P, XphihatQ);
        }
        ep_biex_diff_add_use_saved_nss_last(ZQn1P, R0[bits-2-1-i], R1[bits-2-1-i], XQn1P, ZQn1P, XphihatQ);

        fp2_mul(result, ZQn1P, lineQ);
        fp2_mul(result, result, linephihatQ);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
	}
}

void pp_map_tatep_biext_nss(fp2_t result, ep_t P, ep_t Q){
    fp2_t XPQ, XPphihatQ;
    fp_t ZPQ, ZPphihatQ;
    bn_t n;
    fp2_null(XPQ);
    fp2_null(XPphihatQ);
    fp_null(ZPQ);
    fp_null(ZPphihatQ);
    RLC_TRY {
        fp2_new(XPQ);
        fp2_new(XPphihatQ);
        fp_new(ZPQ);
        fp_new(ZPphihatQ);

        fp2_set_dig(result,1);

        if (!ep_is_infty(P) && !ep_is_infty(Q)) {
            fp_prime_get_par(n);
            ep2_add_PQ_PphihatQ_nss(XPQ, ZPQ, XPphihatQ, ZPphihatQ, P, Q);
            ep2_doubleadd_ladder_tate_nss(result, n, P, Q, XPQ, ZPQ, XPphihatQ, ZPphihatQ);
            pp_exp_k2(result, result);
        }
    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
        fp2_free(XPQ);
        fp2_free(XPphihatQ);
        fp_free(ZPQ);
        fp_free(ZPphihatQ);
	}
}


//=========================================
//          Omega Pairing
//=========================================

void ep2_add_PQ_nss(fp2_t XPQ, fp_t ZPQ, const ep_t P, const ep_t Q){
    fp_t t0, t1, t2, t3;
    fp_null(t0);
    fp_null(t1);
    fp_null(t2);
    fp_null(t3);
    
    fp_new(t0);
    fp_new(t1);
    fp_new(t2);
    fp_new(t3);
    
    fp_dbl(t3, Q->x);
    fp_sub(t0, t3, P->x);
    fp_mul(t1, P->x, t3);
    fp_sub_dig(t2, t1, 1);
    fp_mul(XPQ[0], t0, t2);

    fp_add(t0, t3, P->x);
    fp_add_dig(t2, t1, 1);
    fp_sqr(ZPQ, t0);

    fp_mul(t1, P->y, Q->y);
    fp_dbl(XPQ[1], t1);
    fp_dbl(XPQ[1], XPQ[1]);

    fp_free(t0);
    fp_free(t1);
    fp_free(t2);
    fp_free(t3);
}

//a=1 all below
void ep_biex_doubleadd_nss(fp_t drx, fp_t drz, fp2_t arx, fp2_t arz, const fp_t p1x, const fp_t p1z, const fp2_t p2x, const fp2_t p2z, const fp_t qx){
    fp_t t0, t1, x4, z4, c24z4;
    fp2_t t2, t3, x5, z5;
    fp_null(t0);
    fp_null(t1);
    fp_null(x4);
    fp_null(z4);
    fp_null(c24z4);
    fp2_null(t2);
    fp2_null(t3);
    fp2_null(x5);
    fp2_null(z5);
    RLC_TRY {
        fp_new(t0);
        fp_new(t1);
        fp_new(x4);
        fp_new(z4);
        fp_new(c24z4);
        fp2_new(t2);
        fp2_new(t3);
        fp2_new(x5);
        fp2_new(z5);
        fp_add(t0, p1x, p1z);
        fp_sub(t1, p1x, p1z);
        fp_sqr(x4, t0);
        fp2_sub(t2, p2x, p2z);
        fp2_add(x5, p2x, p2z);

        for(int i=0;i<2;i++) fp_mul(t2[i], t0, t2[i]);
        fp_sqr(z4, t1);
        fp_dbl(c24z4, z4);
        for(int i=0;i<2;i++) fp_mul(t3[i], t1, x5[i]);
        fp_sub(t0, x4, z4);
        fp_mul(x4, x4, c24z4);
        
        fp2_sub(z5, t2, t3);
        fp_add(z4, t0, c24z4);
        fp2_add(x5, t2, t3);
        fp_mul(z4, z4, t0);

        fp2_sqr(z5, z5);
        fp2_sqr(x5, x5);

        for(int i=0;i<2;i++) fp_mul(z5[i], qx, z5[i]);
        fp2_copy(arz, z5);
        fp2_copy(arx, x5);
        fp_copy(drz, z4);
        fp_copy(drx, x4);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
        fp_free(t0);
        fp_free(t1);
        fp_free(x4);
        fp_free(z4);
        fp_free(c24z4);
        fp2_free(t2);
        fp2_free(t3);
        fp2_free(x5);
        fp2_free(z5);
	}
}

// void ep_biex_doubleadd_nss_first(/*fp_t r0, fp_t r1,*/ fp_t drx, fp_t drz, fp2_t arx, fp2_t arz, const fp_t p1x, /*const fp_t p1z,*/ const fp2_t p2x, const fp_t p2z, const fp_t qx)

void ep_biex_doubleadd_nss_mid(fp2_t mid, fp_t drx, fp_t drz, fp2_t arx, fp2_t arz, const fp_t p1x, const fp_t p1z, const fp2_t p2x, const fp2_t p2z, const fp_t qx){
    fp_t t0, t1, x4, z4, c24z4;
    fp2_t t2, t3, x5, z5;
    fp_null(t0);
    fp_null(t1);
    fp_null(x4);
    fp_null(z4);
    fp_null(c24z4);
    fp2_null(t2);
    fp2_null(t3);
    fp2_null(x5);
    fp2_null(z5);
    RLC_TRY {
        fp_new(t0);
        fp_new(t1);
        fp_new(x4);
        fp_new(z4);
        fp_new(c24z4);
        fp2_new(t2);
        fp2_new(t3);
        fp2_new(x5);
        fp2_new(z5);
        fp_add(t0, p1x, p1z);
        fp_sub(t1, p1x, p1z);
        fp_sqr(x4, t0);
        fp2_sub(t2, p2x, p2z);
        fp2_add(x5, p2x, p2z);

        for(int i=0;i<2;i++) fp_mul(t2[i], t0, t2[i]);
        fp_sqr(z4, t1);
        fp_dbl(c24z4, z4);
        for(int i=0;i<2;i++) fp_mul(t3[i], t1, x5[i]);
        fp_sub(t0, x4, z4);
        fp_mul(x4, x4, c24z4);
        
        fp2_sub(z5, t2, t3);
        fp_add(z4, t0, c24z4);
        fp2_add(x5, t2, t3);
        fp_mul(z4, z4, t0);

        fp2_copy(mid, z5);//new for test
        fp2_sqr(z5, z5);
        fp2_sqr(x5, x5);

        // fp_neg(qx, qx);
        for(int i=0;i<2;i++) fp_mul(z5[i], qx, z5[i]);
        fp2_copy(arz, z5);
        // fp2_copy(mid, arz);//new for test
        fp2_copy(arx, x5);
        fp_copy(drz, z4);
        fp_copy(drx, x4);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
        fp_free(t0);
        fp_free(t1);
        fp_free(x4);
        fp_free(z4);
        fp_free(c24z4);
        fp2_free(t2);
        fp2_free(t3);
        fp2_free(x5);
        fp2_free(z5);
	}
}

//no result point for QnP, just the Z not squared   //but nP needed as line function need it
void ep_biex_doubleadd_nss_last(fp_t drx, fp_t drz, fp2_t arzlast, const fp_t p1x, const fp_t p1z, const fp2_t p2x, const fp2_t p2z, const fp_t qx){
    fp_t t0, t1, x4, z4, c24z4;
    fp2_t t2, t3, x5, z5;
    fp_null(t0);
    fp_null(t1);
    fp_null(x4);
    fp_null(z4);
    fp_null(c24z4);
    fp2_null(t2);
    fp2_null(t3);
    fp2_null(x5);
    fp2_null(z5);
    RLC_TRY {
        fp_new(t0);
        fp_new(t1);
        fp_new(x4);
        fp_new(z4);
        fp_new(c24z4);
        fp2_new(t2);
        fp2_new(t3);
        fp2_new(x5);
        fp2_new(z5);
        fp_add(t0, p1x, p1z);
        fp_sub(t1, p1x, p1z);
        fp_sqr(x4, t0);
        fp2_sub(t2, p2x, p2z);
        fp2_add(x5, p2x, p2z);

        for(int i=0;i<2;i++) fp_mul(t2[i], t0, t2[i]);
        fp_sqr(z4, t1);
        fp_dbl(c24z4, z4);
        for(int i=0;i<2;i++) fp_mul(t3[i], t1, x5[i]);
        fp_sub(t0, x4, z4);
        fp_mul(x4, x4, c24z4);
        
        fp2_sub(z5, t2, t3);    //lzyr
        fp_add(z4, t0, c24z4);
        fp_mul(z4, z4, t0);

        fp2_copy(arzlast, z5);
        fp_copy(drz, z4);
        fp_copy(drx, x4);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
        fp_free(t0);
        fp_free(t1);
        fp_free(x4);
        fp_free(z4);
        fp_free(c24z4);
        fp2_free(t2);
        fp2_free(t3);
        fp2_free(x5);
        fp2_free(z5);
	}
}

void ep2_line_function_omega_nss(fp2_t resultQ, fp2_t resultP, const fp_t Xn1P, const fp_t Zn1P, const fp_t Xn2P, const fp_t Zn2P, const fp_t Xn1Q, const fp_t Zn1Q, const fp_t Xn2Q, const fp_t Zn2Q, const fp_t px, const fp_t py, const fp_t qx, const fp_t qy){
    fp_t k00, k0, k1, k2, k3, k4, k5, k6, rQr, rQi, rPr, rPi;
    fp_null(k00);
    fp_null(k0);
    fp_null(k1);
    fp_null(k2);
    fp_null(k3);
    fp_null(k4);
    fp_null(k5);
    fp_null(k6);
    fp_null(rQr);
    fp_null(rQi);
    fp_null(rPr);
    fp_null(rPi);
    RLC_TRY {
        fp_new(k00);
        fp_new(k0);
        fp_new(k1);
        fp_new(k2);
        fp_new(k3);
        fp_new(k4);
        fp_new(k5);
        fp_new(k6);
        fp_new(rQr);
        fp_new(rQi);
        fp_new(rPr);
        fp_new(rPi);
        fp_mul(k0, px, Zn1P);
        fp_add(k1, Xn1P, k0);
        fp_mul(k2, ep_curve_get_beta(), py);
        fp_mul(k3, k2, Zn1P);
        fp_mul(k4, Zn2P, k3);
        fp_dbl(k4, k4);
        fp_mul(k3, k3, k4);
        fp_mul(k5, k1, k4);
        fp_neg(k5, k5);
        fp_sqr(k1, k1);
        fp_mul(k1, k1, Xn2P);
        fp_sub(k0, Xn1P, k0);
        fp_mul(k4, px, Xn1P);
        fp_sub(k4, Zn1P, k4);
        fp_mul(k6, k4, k0);
        fp_mul(k6, k6, Zn2P);
        fp_sub(k6, k6, k1);
        fp_sub(k6, k6, k3);
        fp_mul(rQr, k2, k5);
        fp_mul(k5, k5, qy);
        fp_dbl(k5, k5);
        fp_neg(rQi, k5);
        fp_dbl(k1, qx);
        fp_neg(k1, k1);
        fp_add(k3, k1, px);
        fp_copy(rPi, k3);   //px-2qx
        fp_sub(k4, k1, px);
        fp_mul(k3, k3, k6);
        fp_mul(k4, k4, k6);
        fp_sub(rQr, rQr, k3);
        fp_copy(resultQ[0], rQr);
        fp_copy(resultQ[1], rQi);

        fp_dbl(k00, qx);
        fp_neg(k00, k00);
        fp_mul(k0, k00, Zn1Q);
        fp_add(k1, Xn1Q, k0);
        fp_sub(k2, Xn1Q, k0);
        fp_mul(k2, k2, Zn2Q);
        fp_mul(k3, ep_curve_get_beta(), qy);
        fp_neg(k3, k3);
        fp_dbl(k3, k3);
        fp_mul(k4, k3, Zn1Q);
        fp_mul(k5, Zn2Q, k4);
        fp_dbl(k5, k5);
        fp_mul(k4, k4, k5);
        fp_mul(k5, k1, k5);
        fp_mul(k3, k3, k5);
        fp_mul(k6, k00, Xn1Q);
        fp_sub(k6, Zn1Q, k6);
        fp_mul(k2, k2, k6);
        fp_sqr(k1, k1);
        fp_mul(k1, Xn2Q, k1);

        fp_mul(rPr, k5, py);
        fp_dbl(rPr, rPr);
        fp_neg(rPr, rPr);
        fp_dbl(k4, k4);
        fp_add(k4, k4, k2);
        fp_sub(k4, k4, k1);
        fp_mul(rPi, rPi, k4);
        fp_dbl(k3, k3);
        fp_sub(rPi, k3, rPi);

        fp_copy(resultP[0], rPr);
        fp_copy(resultP[1], rPi);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
        fp_free(k00);
        fp_free(k0);
        fp_free(k1);
        fp_free(k2);
        fp_free(k3);
        fp_free(k4);
        fp_free(k5);
        fp_free(k6);
        fp_free(rQr);
        fp_free(rQi);
        fp_free(rPr);
        fp_free(rPi);
	}
}

void ep2_doubleadd_ladder_omega_nss(fp2_t result, const bn_t n, const ep_t P, const ep_t Q, const fp2_t XQP, const fp_t ZQP){
    fp_t Xn1P, Zn1P, Xn2P, Zn2P, Xn1Q, Zn1Q, Xn2Q, Zn2Q, XQ, XP;
    fp2_t XQn1P, ZQn1P, XPn1Q, ZPn1Q, lineQ, lineP, resultQ, resultP;
    size_t bits;
    bits = bn_bits(n);
    int i;
    fp_null(Xn1P);
    fp_null(Zn1P);
    fp_null(Xn2P);
    fp_null(Zn2P);
    fp_null(Xn1Q);
    fp_null(Zn1Q);
    fp_null(Xn2Q);
    fp_null(Zn2Q);
    fp_null(XQ);
    fp_null(XP);
    fp2_null(XQn1P);
    fp2_null(ZQn1P);
    fp2_null(XPn1Q);
    fp2_null(ZPn1Q);
    fp2_null(lineQ);
    fp2_null(lineP);
    fp2_null(resultQ);
    fp2_null(resultP);

    RLC_TRY {
        fp_new(Xn1P);
        fp_new(Zn1P);
        fp_new(Xn2P);
        fp_new(Zn2P);
        fp_new(Xn1Q);
        fp_new(Zn1Q);
        fp_new(Xn2Q);
        fp_new(Zn2Q);
        fp_new(XQ);
        fp_new(XP);
        fp2_new(XQn1P);
        fp2_new(ZQn1P);
        fp2_new(XPn1Q);
        fp2_new(ZPn1Q);
        fp2_new(lineQ);
        fp2_new(lineP);
        fp2_new(resultQ);
        fp2_new(resultP);

        fp_dbl(XQ, Q->x);
        fp_neg(XQ, XQ);
        fp_copy(XP, P->x);
        
        ep_biex_doubleadd_nss_first(Xn1P, Zn1P, XQn1P, ZQn1P, XP, XQP, ZQP, XQ);
        ep_biex_doubleadd_nss_first(Xn1Q, Zn1Q, XPn1Q, ZPn1Q, XQ, XQP, ZQP, XP);

        for (i = bits - 3; i >= bits - 23 + 1; i--) {       //n2
            ep_biex_doubleadd_nss(Xn1P, Zn1P, XQn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
            ep_biex_doubleadd_nss(Xn1Q, Zn1Q, XPn1Q, ZPn1Q, Xn1Q, Zn1Q, XPn1Q, ZPn1Q, XP);
        }
        ep_biex_doubleadd_nss_mid(resultQ, Xn1P, Zn1P, XQn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
        ep_biex_doubleadd_nss_mid(resultP, Xn1Q, Zn1Q, XPn1Q, ZPn1Q, Xn1Q, Zn1Q, XPn1Q, ZPn1Q, XP);
        i--;
        fp_copy(Zn2P, Zn1P);
        fp_copy(Xn2P, Xn1P);
        fp_copy(Zn2Q, Zn1Q);
        fp_copy(Xn2Q, Xn1Q);
        
        for(;i>=0 + 1;i--){     //n1
            ep_biex_doubleadd_nss(Xn1P, Zn1P, XQn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
            ep_biex_doubleadd_nss(Xn1Q, Zn1Q, XPn1Q, ZPn1Q, Xn1Q, Zn1Q, XPn1Q, ZPn1Q, XP);
        }
        ep_biex_doubleadd_nss_last(Xn1P, Zn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
        ep_biex_doubleadd_nss_last(Xn1Q, Zn1Q, ZPn1Q, Xn1Q, Zn1Q, XPn1Q, ZPn1Q, XP);

        ep2_line_function_omega_nss(lineQ, lineP, Xn1P, Zn1P, Xn2P, Zn2P, Xn1Q, Zn1Q, Xn2Q, Zn2Q, P->x, P->y, Q->x, Q->y);
        fp2_mul(resultQ, resultQ, ZQn1P);
        fp2_mul(resultQ, resultQ, lineQ);
        fp2_mul(resultP, resultP, ZPn1Q);
        fp2_mul(resultP, resultP, lineP);
        fp2_inv_cyc(resultP, resultP);
        fp2_mul(result, resultQ, resultP);
    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
        fp_free(Xn1P);
        fp_free(Zn1P);
        fp_free(Xn2P);
        fp_free(Zn2P);
        fp_free(Xn1Q);
        fp_free(Zn1Q);
        fp_free(Xn2Q);
        fp_free(Zn2Q);
        fp_free(XQ);
        fp_free(XP);
        fp2_free(XQn1P);
        fp2_free(ZQn1P);
        fp2_free(XPn1Q);
        fp2_free(ZPn1Q);
        fp2_free(lineQ);
        fp2_free(lineP);
        fp2_free(resultQ);
        fp2_free(resultP);
	}
}

void pp_map_omegap_biext_nss(fp2_t result, ep_t P, ep_t Q){
    fp2_t XPQ;
    fp_t ZPQ;
    bn_t n;
    fp2_null(XPQ);
    fp_null(ZPQ);
    RLC_TRY {
        fp2_new(XPQ);
        fp_new(ZPQ);

        fp2_set_dig(result,1);

        if (!ep_is_infty(P) && !ep_is_infty(Q)) {
            fp_prime_get_par(n);
            ep2_add_PQ_nss(XPQ, ZPQ, P, Q);
            ep2_doubleadd_ladder_omega_nss(result, n, P, Q, XPQ, ZPQ);
            fp2_conv_cyc(result, result);
        }
    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
        fp2_free(XPQ);
        fp_free(ZPQ);
	}
}
