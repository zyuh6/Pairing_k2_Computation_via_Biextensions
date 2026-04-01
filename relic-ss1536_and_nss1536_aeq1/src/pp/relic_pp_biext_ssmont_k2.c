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

//convert from y^2=x^3-x to y^2=x^3+x through a 2-isogeny and an isomorphism
//this can be precomputed as the generator of the curve is fixed
void ep_ss_to_ssmont(ep_t Pmont, const ep_t P){
    fp_t a1, a2;
    ep_t tmp;
    fp_null(a1);
    fp_null(a2);
    ep_null(tmp);

    RLC_TRY {
        fp_new(a1);
        fp_new(a2);
        ep_new(tmp);

        if(P->coord != BASIC) ep_norm(P, P);

        fp_set_dig(a1, 12);
        fp_inv(a1, a1);
        fp_srt(a1, a1);
        fp_srt(a2, a1);
        fp_exp_dig(a2, a2, 3);

        fp_inv(tmp->x, P->x);
        fp_sqr(tmp->y, tmp->x);
        fp_mul_dig(tmp->y, tmp->y, 3);
        fp_add_dig(tmp->y, tmp->y, 1);
        fp_mul(tmp->y, tmp->y, P->y);
        fp_mul(tmp->x, tmp->x, P->y);
        fp_sqr(tmp->x, tmp->x);

        fp_mul(tmp->x, tmp->x, a1);
        fp_mul(tmp->y, tmp->y, a2);
        fp_copy(Pmont->x, tmp->x);
        fp_copy(Pmont->y, tmp->y);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
	}

}

//P=[px,py],  twist(Q)=[-qx,qy*u]   //compute P+twist(Q) on E(Fp2)  //a=1
//3m + 1s on Fp
void ep2_add_PQ_ssmont(fp2_t XPQ, fp_t ZPQ, const ep_t P, const ep_t Q){
    fp_t t0, t1;
    fp_null(t0);
    fp_null(t1);
    
    fp_new(t0);
    fp_new(t1);
    
    fp_sub(t0, Q->x, P->x);
    fp_mul(t1, P->x, Q->x);
    fp_sub_dig(t1, t1, 1);
    fp_mul(XPQ[0], t0, t1);

    fp_mul(t1, P->y, Q->y);
    fp_neg(t1, t1);
    fp_dbl(XPQ[1], t1);

    fp_add(t0, Q->x, P->x);
    fp_sqr(ZPQ, t0);

    fp_free(t0);
    fp_free(t1);
}

//8m+2s+2S
void ep_biex_doubleadd_ssmont(fp_t drx, fp_t drz, fp2_t arx, fp2_t arz, const fp_t p1x, const fp_t p1z, const fp2_t p2x, const fp2_t p2z, const fp_t qx){
    fp_t t0, t1, x4, z4, c24z4;
    fp2_t t2, t3, x5, z5;
    fp_null(t0);
    fp_null(x4);
    fp_null(z4);
    fp_null(c24z4);
    fp2_null(t2);
    fp2_null(t3);
    fp2_null(x5);
    fp2_null(z5);

    RLC_TRY {
        fp_new(t0);
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
        fp_free(t0);
        fp_free(x4);
        fp_free(z4);
        fp_free(c24z4);
        fp2_free(t2);
        fp2_free(t3);
        fp2_free(x5);
        fp2_free(z5);
	}
}

//7m+1s+2S
void ep_biex_doubleadd_ssmont_first(fp_t drx, fp_t drz, fp2_t arx, fp2_t arz, const fp_t p1x, /*const fp_t p1z,*/ const fp2_t p2x, const fp_t p2z, const fp_t qx){
    fp_t t0, x4, z4, c24z4;
    fp2_t t2, t3, x5, z5;
    RLC_TRY {
        fp_sqr(z4, p1x);
        fp_add_dig(z4, z4, 1);
        fp_dbl(t0, p1x);
        fp_add(x4, z4, t0);
        fp_sub(z4, z4, t0);

        fp_dbl(c24z4, z4);
        fp_sub(t0, x4, z4);
        fp_mul(x4, x4, c24z4);
        
        fp_copy(z5[1], p2x[1]);
        fp_mul(z5[0], p1x, p2z);
        fp_sub(z5[0], p2x[0], z5[0]);

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

//8m+2s+2S
void ep_biex_doubleadd_ssmont_mid(fp2_t mid, fp_t drx, fp_t drz, fp2_t arx, fp2_t arz, const fp_t p1x, const fp_t p1z, const fp2_t p2x, const fp2_t p2z, const fp_t qx){
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

        fp2_copy(mid, z5);
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

//6m+2s     //to be edited with lzyr
void ep_biex_doubleadd_ssmont_last(fp_t drx, fp_t drz, fp2_t arzlast, const fp_t p1x, const fp_t p1z, const fp2_t p2x, const fp2_t p2z, const fp_t qx){
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

//compute line_{[2^41]P, -P}(Q)     //a=1
void ep2_line_function_ssmont(fp2_t result, const fp_t Xn1P, const fp_t Zn1P, const fp_t Xn2P, const fp_t Zn2P, const fp_t px, const fp_t py, const fp_t qx, const fp_t qy){
    fp_t t0, t1, t2, t3, t4, t5, rr, ri;
    fp_null(t0);
    fp_null(t1);
    fp_null(t2);
    fp_null(t3);
    fp_null(t4);
    fp_null(t5);
    fp_null(rr);
    fp_null(ri);
    RLC_TRY {
        fp_new(t0);
        fp_new(t1);
        fp_new(t2);
        fp_new(t3);
        fp_new(t4);
        fp_new(t5);
        fp_new(rr);
        fp_new(ri);

        fp_add(t0, px, qx);
        fp_mul(t1, Zn2P, px);
        fp_sub(t2, Xn2P, t1);
        fp_mul(t3, Zn1P, Zn2P);
        fp_sqr(t4, py);
        fp_mul(t4, t4, t3);
        fp_dbl(t4, t4);

        fp_mul(ri, py, qy);
        fp_mul(ri, ri, t3);
        fp_mul(ri, ri, t2);
        fp_dbl(ri, ri);

        fp_add(t1, t1, Xn2P);
        fp_mul(t1, t1, Zn1P);
        fp_mul(t5, px, Xn2P);
        fp_add(t5, t5, Zn2P);
        fp_mul(t5, t5, t1);

        fp_sqr(t3, t2);
        fp_mul(t3, t3, Xn1P);
        fp_mul(t1, Zn2P, t4);
        fp_add(t1, t3, t1);

        fp_mul(t4, t4, t2);

        fp_sub(t1, t5, t1);
        fp_mul(t1, t1, t0);
        fp_sub(rr, t1, t4);

        fp_copy(result[0], rr);
        fp_copy(result[1], ri);

    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
        fp_free(t0);
        fp_free(t1);
        fp_free(t2);
        fp_free(t3);
        fp_free(t4);
        fp_free(t5);
        fp_free(rr);
        fp_free(ri);
	}
}

//n = 2^255 + 2^41 + 1 = n1 + n2 + 1
void ep2_doubleadd_ladder_ssmont(fp2_t result, fp_t Xn1P, fp_t Zn1P, fp_t Xn2P, fp_t Zn2P, /*fp2_t ZQn1P, fp2_t ZQn2P,*/ const bn_t n, const fp_t XQ, const fp_t XP, const fp2_t XQP, const fp_t ZQP){
    fp2_t XQn1P, XQn2P, ZQn1P, ZQn2P;
    size_t bits;
    int i;
    fp2_null(XQn1P);
    fp2_null(XQn2P);
    fp2_null(ZQn1P);
    fp2_null(ZQn2P);

    RLC_TRY {
        fp2_new(XQn1P);
        fp2_new(XQn2P);
        fp2_new(ZQn1P);
        fp2_new(ZQn2P);
        bits = bn_bits(n);
        ep_biex_doubleadd_ssmont_first(Xn1P, Zn1P, XQn1P, ZQn1P, XP, XQP, ZQP, XQ);

        for (i = bits - 3; i >= bits - 42 + 1; i--) {       //n2
            ep_biex_doubleadd_ssmont(Xn1P, Zn1P, XQn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
        }
        ep_biex_doubleadd_ssmont_mid(result, Xn1P, Zn1P, XQn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
        i--;
        fp_copy(Xn2P, Xn1P);
        fp_copy(Zn2P, Zn1P);
        for(;i>=0 + 1;i--){     //n1
            ep_biex_doubleadd_ssmont(Xn1P, Zn1P, XQn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
        }
        ep_biex_doubleadd_ssmont_last(Xn1P, Zn1P, ZQn1P, Xn1P, Zn1P, XQn1P, ZQn1P, XQ);
        fp2_mul(result, result, ZQn1P);
    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
        fp2_free(XQn1P);
        fp2_free(XQn2P);
        fp2_free(ZQn1P);
        fp2_free(ZQn2P);
	}
}

void pp_map_tatep_biext_ssmont1536(fp2_t result, ep_t P, ep_t Q){
    fp2_t XPQ, millerresult;
    fp_t ZPQ, XQ, Xn1P, Zn1P, Xn2P, Zn2P;
    bn_t n;
    fp2_null(XPQ);
    fp2_null(millerresult);
    fp_null(ZPQ);
    fp_null(XQ);
    fp_null(Xn1P);
    fp_null(Zn1P);
    fp_null(Xn2P);
    fp_null(Zn2P);
    bn_null(n);

    RLC_TRY {
        fp2_new(XPQ);
        fp2_new(millerresult);
        fp_new(ZPQ);
        fp_new(XQ);
        fp_new(Xn1P);
        fp_new(Zn1P);
        fp_new(Xn2P);
        fp_new(Zn2P);
        bn_new(n);

        fp2_set_dig(result,1);
        if (!ep_is_infty(P) && !ep_is_infty(Q)) {
            fp_prime_get_par(n);
            fp_neg(XQ, Q->x);
            ep2_add_PQ_ssmont(XPQ, ZPQ, P, Q);
            ep2_doubleadd_ladder_ssmont(millerresult, Xn1P, Zn1P, Xn2P, Zn2P, n, XQ, P->x, XPQ, ZPQ);
            ep2_line_function_ssmont(result, Xn1P, Zn1P, Xn2P, Zn2P, P->x, P->y, Q->x, Q->y);
            fp2_mul(result, result, millerresult);
            pp_exp_k2(result, result);
        }
    }
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
        fp2_free(XPQ);
        fp2_free(millerresult);
        fp_free(ZPQ);
        fp_free(XQ);
        fp_free(Xn1P);
        fp_free(Zn1P);
        fp_free(Xn2P);
        fp_free(Zn2P);
        bn_free(n);
	}
}
