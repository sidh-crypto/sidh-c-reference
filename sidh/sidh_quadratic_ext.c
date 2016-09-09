/*
 * Copyright (C) 2016 Javad Doliskani, javad.doliskani@uwaterloo.ca
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "sidh_quadratic_ext.h"
#include "sidh_util.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void sidh_fp_init_chararacteristic_ui(long p) {
    mpz_init_set_ui(characteristic, p);
}

void sidh_fp_init_chararacteristic_str(const char *value) {
    mpz_init_set_str(characteristic, value, 10);
}

void sidh_fp_init_chararacteristic(const mpz_t p) {
    mpz_init_set(characteristic, p);
}

void sidh_fp_set(mpz_t x,
                 const mpz_t a) {
    mpz_mod(x, a, characteristic);
}

void sidh_fp_add(mpz_t x,
                 const mpz_t a,
                 const mpz_t b) {
    mpz_add(x, a, b);
    mpz_mod(x, x, characteristic);
}

void sidh_fp_add_ui(mpz_t x,
                    const mpz_t a,
                    unsigned long b) {
    mpz_add_ui(x, a, b);
    mpz_mod(x, x, characteristic);
}

void sidh_fp_sub(mpz_t x,
                 const mpz_t a,
                 const mpz_t b) {
    mpz_sub(x, a, b);
    mpz_mod(x, x, characteristic);
}

void sidh_fp_sub_ui(mpz_t x,
                    const mpz_t a,
                    unsigned long b) {
    mpz_sub_ui(x, a, b);
    mpz_mod(x, x, characteristic);
}

void sidh_fp_mul(mpz_t x,
                 const mpz_t a,
                 const mpz_t b) {
    mpz_mul(x, a, b);
    mpz_mod(x, x, characteristic);
}

void sidh_fp_mul_si(mpz_t x,
                    const mpz_t a,
                    long b) {
    mpz_mul_si(x, a, b);
    mpz_mod(x, x, characteristic);
}

void sidh_fp_inv(mpz_t x,
                 const mpz_t a) {
    mpz_invert(x, a, characteristic);
}

void sidh_fp_div(mpz_t x,
                 const mpz_t a,
                 const mpz_t b) {
    sidh_fp_inv(x, b);
    sidh_fp_mul(x, a, x);
}

void sidh_fp_neg(mpz_t x,
                 const mpz_t a) {
    sidh_fp_sub(x, characteristic, a);
}

void sidh_fp_sqrt(mpz_t x,
                  const mpz_t a) {
    mpz_t exponent;
    mpz_init(exponent);

    // compute (p + 1) / 4
    mpz_add_ui(exponent, characteristic, 1);
    mpz_divexact_ui(exponent, exponent, 4);

    mpz_powm(x, a, exponent, characteristic);
    mpz_clear(exponent);
}

//////////////// fp2 methods //////////////////////////

void sidh_fp2_init(fp2_element_t x) {
    mpz_inits(x->a, x->b, NULL);
}

void sidh_fp2_init_set_si(fp2_element_t x,
                          long a,
                          long b) {
    mpz_init_set_si(x->a, a);
    mpz_init_set_si(x->b, b);
}

void sidh_fp2_init_set_str(fp2_element_t x,
                           const char *a,
                           const char *b) {
    mpz_init_set_str(x->a, a, 10);
    mpz_init_set_str(x->b, b, 10);
}

void sidh_fp2_init_set(fp2_element_t x,
                       const fp2_element_t a) {
    mpz_init_set(x->a, a->a);
    mpz_init_set(x->b, a->b);
}

void sidh_fp2_clear(fp2_element_t x) {
    mpz_clears(x->a, x->b, NULL);
}

void sidh_fp2_set(fp2_element_t x,
                  const fp2_element_t b) {
    mpz_set(x->a, b->a);
    mpz_set(x->b, b->b);
}

void sidh_fp2_zero(fp2_element_t x) {
    mpz_set_si(x->a, 0);
    mpz_set_si(x->b, 0);
}

void sidh_fp2_one(fp2_element_t x) {
    mpz_set_si(x->a, 0);
    mpz_set_si(x->b, 1);
}

char *sidh_fp2_get_str(const fp2_element_t a) {

    if (mpz_cmp_si(a->a, 0) == 0 && mpz_cmp_si(a->b, 0) == 0) {
        return "0";
    }

    if (mpz_cmp_si(a->a, 0) == 0) {
        return mpz_get_str(NULL, 10, a->b);
    }

    char *result = "";

    if (mpz_cmp_si(a->b, 0) == 0) {
        result = sidh_concat(result, mpz_get_str(NULL, 10, a->a));
        result = sidh_concat(result, " * i");
        return result;
    }

    result = sidh_concat(result, mpz_get_str(NULL, 10, a->a));
    result = sidh_concat(result, " * i + ");
    result = sidh_concat(result, mpz_get_str(NULL, 10, a->b));

    return result;
}

void sidh_fp2_add(fp2_element_t x,
                  const fp2_element_t a,
                  const fp2_element_t b) {
    sidh_fp_add(x->a, a->a, b->a);
    sidh_fp_add(x->b, a->b, b->b);
}

void sidh_fp2_add_ui(fp2_element_t x,
                     const fp2_element_t a,
                     unsigned long b) {
    sidh_fp_add_ui(x->b, a->b, b);
    sidh_fp_set(x->a, a->a);
}

void sidh_fp2_sub(fp2_element_t x,
                  const fp2_element_t a,
                  const fp2_element_t b) {
    sidh_fp_sub(x->a, a->a, b->a);
    sidh_fp_sub(x->b, a->b, b->b);
}

void sidh_fp2_sub_ui(fp2_element_t x,
                     const fp2_element_t a,
                     unsigned long b) {
    sidh_fp_sub_ui(x->b, a->b, b);
    sidh_fp_set(x->a, a->a);
}

void sidh_fp2_mul(fp2_element_t x,
                  const fp2_element_t a,
                  const fp2_element_t b) {
    mpz_t temp1;
    mpz_t temp2;

    mpz_init(temp1);
    mpz_init(temp2);

    fp2_element_t result;
    sidh_fp2_init(result);

    // (a + b) * (c + d)
    sidh_fp_add(temp1, a->a, a->b);
    sidh_fp_add(temp2, b->a, b->b);
    sidh_fp_mul(result->a, temp1, temp2);

    // a * c
    sidh_fp_mul(temp1, a->a, b->a);
    // b * d
    sidh_fp_mul(temp2, a->b, b->b);

    sidh_fp_sub(result->a, result->a, temp1);
    sidh_fp_sub(result->a, result->a, temp2);
    sidh_fp_sub(result->b, temp2, temp1);
    sidh_fp2_set(x, result);

    mpz_clear(temp1);
    mpz_clear(temp2);
    sidh_fp2_clear(result);
}

void sidh_fp2_square(fp2_element_t x,
                     const fp2_element_t a) {
    mpz_t temp1;
    mpz_t temp2;

    mpz_init(temp1);
    mpz_init(temp2);

    fp2_element_t result;
    sidh_fp2_init(result);

    // (b + a) * (b - a)
    sidh_fp_add(temp1, a->a, a->b);
    sidh_fp_sub(temp2, a->b, a->a);
    sidh_fp_mul(result->b, temp1, temp2);

    // 2 * a * b
    sidh_fp_mul(result->a, a->a, a->b);
    sidh_fp_mul_si(result->a, result->a, 2);

    sidh_fp2_set(x, result);

    mpz_clear(temp1);
    mpz_clear(temp2);
    sidh_fp2_clear(result);
}

void sidh_fp2_pow_ui(fp2_element_t x,
                     const fp2_element_t a,
                     unsigned long n) {
    mpz_t temp_n;
    mpz_init_set_ui(temp_n, n);
    sidh_fp2_pow(x, a, temp_n);
    mpz_clear(temp_n);
}

void sidh_fp2_pow(fp2_element_t x,
                  const fp2_element_t a,
                  const mpz_t n) {
    if (mpz_cmp_ui(n, 0) == 0) {
        sidh_fp2_one(x);
        return;
    }

    fp2_element_t temp1;
    fp2_element_t temp2;
    sidh_fp2_init_set_si(temp1, 0, 1);
    sidh_fp2_init_set(temp2, a);

    long num_bits = mpz_sizeinbase(n, 2);
    for (long i = 0; i < num_bits; i++) {
        if (mpz_tstbit(n, i) == 1)
            sidh_fp2_mul(temp1, temp1, temp2);
        sidh_fp2_square(temp2, temp2);
    }

    sidh_fp2_set(x, temp1);

    sidh_fp2_clear(temp1);
    sidh_fp2_clear(temp2);
}

void sidh_fp2_conjugate(fp2_element_t x,
                        const fp2_element_t a) {
    sidh_fp2_set(x, a);
    sidh_fp_neg(x->a, x->a);
}

void sidh_fp2_negate(fp2_element_t x,
                     const fp2_element_t a) {
    sidh_fp2_set(x, a);
    sidh_fp_neg(x->a, x->a);
    sidh_fp_neg(x->b, x->b);
}

void sidh_fp2_mul_scaler(fp2_element_t x,
                         const fp2_element_t a,
                         const mpz_t scaler) {
    sidh_fp_mul(x->a, a->a, scaler);
    sidh_fp_mul(x->b, a->b, scaler);
}

void sidh_fp2_mul_scaler_si(fp2_element_t x,
                            const fp2_element_t a,
                            long scaler) {
    sidh_fp_mul_si(x->a, a->a, scaler);
    sidh_fp_mul_si(x->b, a->b, scaler);
}

void sidh_fp2_inv(fp2_element_t x,
                  const fp2_element_t a) {
    mpz_t temp;
    fp2_element_t result;

    mpz_init(temp);
    sidh_fp2_init(result);

    sidh_fp2_conjugate(result, a);
    sidh_fp2_norm(temp, a);
    sidh_fp_inv(temp, temp);
    sidh_fp2_mul_scaler(result, result, temp);
    sidh_fp2_set(x, result);

    mpz_clear(temp);
    sidh_fp2_clear(result);
}

void sidh_fp2_div(fp2_element_t x,
                  const fp2_element_t a,
                  const fp2_element_t b) {
    fp2_element_t result;
    sidh_fp2_init(result);

    sidh_fp2_inv(result, b);
    sidh_fp2_mul(result, a, result);
    sidh_fp2_set(x, result);

    sidh_fp2_clear(result);
}

int sidh_fp2_is_zero(const fp2_element_t a) {
    return !mpz_cmp_si(a->a, 0) && !mpz_cmp_si(a->b, 0);
}

int sidh_fp2_is_one(const fp2_element_t a) {
    return !mpz_cmp_si(a->a, 0) && !mpz_cmp_si(a->b, 1);
}

int sidh_fp2_equals(const fp2_element_t a,
                    const fp2_element_t b) {
    return (mpz_cmp(a->a, b->a) == 0) && (mpz_cmp(a->b, b->b) == 0);
}

void sidh_fp2_random(fp2_element_t x,
                     gmp_randstate_t randstate) {
    mpz_urandomm(x->a, randstate, characteristic);
    mpz_urandomm(x->b, randstate, characteristic);
}

void sidh_fp2_sqrt(fp2_element_t x,
                   const fp2_element_t a) {
    mpz_t exponent;
    fp2_element_t temp_a;
    fp2_element_t b;
    fp2_element_t c;
    fp2_element_t beta;
    mpz_t base_root;
    gmp_randstate_t randstate;

    mpz_init(exponent);
    sidh_fp2_init(temp_a);
    sidh_fp2_init(b);
    sidh_fp2_init(c);
    sidh_fp2_init(beta);
    mpz_init(base_root);
    gmp_randinit_default(randstate);

    // compute (p - 1) / 2
    mpz_sub_ui(exponent, characteristic, 1);
    mpz_divexact_ui(exponent, exponent, 2);

    while (sidh_fp2_is_zero(b)) {
        sidh_fp2_random(c, randstate);
        sidh_fp2_square(temp_a, c);
        sidh_fp2_mul(temp_a, temp_a, a);

        // compute 1 + temp_a^((p - 1) / 2)
        sidh_fp2_pow(b, temp_a, exponent);
        sidh_fp2_add_ui(b, b, 1);
    }

    // compute temp_a * b^2
    sidh_fp2_square(beta, b);
    sidh_fp2_mul(beta, beta, temp_a);

    // beta is now in the prime field
    sidh_fp_sqrt(base_root, beta->b);
    sidh_fp2_inv(b, b);
    sidh_fp2_mul_scaler(b, b, base_root);
    sidh_fp2_div(x, b, c);


    mpz_clear(exponent);
    sidh_fp2_clear(temp_a);
    sidh_fp2_clear(b);
    sidh_fp2_clear(c);
    sidh_fp2_clear(beta);
    mpz_clear(base_root);
    gmp_randclear(randstate);
}

int sidh_fp2_is_square(const fp2_element_t a) {
    mpz_t exponent;
    mpz_t norm;
    fp2_element_t temp;

    mpz_init(exponent);
    mpz_init(norm);
    sidh_fp2_init(temp);

    // a^((p - 1) / 2)
    mpz_sub_ui(exponent, characteristic, 1);
    mpz_divexact_ui(exponent, exponent, 2);
    sidh_fp2_pow(temp, a, exponent);

    sidh_fp2_norm(norm, temp);
    int result = (mpz_cmp_si(norm, 1) == 0);

    mpz_clear(exponent);
    mpz_clear(norm);
    sidh_fp2_clear(temp);

    return result;
}

void sidh_fp2_norm(mpz_t x,
                   const fp2_element_t a) {
    mpz_t temp1;
    mpz_t temp2;
    mpz_inits(temp1, temp2, NULL);

    sidh_fp_mul(temp1, a->a, a->a);
    sidh_fp_mul(temp2, a->b, a->b);
    sidh_fp_add(temp1, temp1, temp2);

    mpz_set(x, temp1);
    mpz_clears(temp1, temp2, NULL);
}