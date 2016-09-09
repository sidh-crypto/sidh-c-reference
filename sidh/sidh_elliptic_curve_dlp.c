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

#include "sidh_elliptic_curve_dlp.h"
#include <stdio.h>

void sidh_elliptic_curve_prime_power_dlp(mpz_t x,
                                         const point_t P,
                                         const point_t Q,
                                         const elliptic_curve_t E,
                                         long l,
                                         long e) {
    mpz_t exponent1;
    mpz_t exponent2;
    point_t temp_P;
    point_t temp_Q;
    point_t temp_R;
    point_t PP;

    mpz_init(exponent1);
    mpz_init(exponent2);
    sidh_point_init(temp_P);
    sidh_point_init(temp_Q);
    sidh_point_init(temp_R);
    sidh_point_init(PP);

    int ladic_rep[e];
    mpz_ui_pow_ui(exponent1, l, e - 1);

    // PP = l^(e - 1) * P once and for all
    sidh_point_mul_scaler(PP, P, exponent1, E);

    // compute the first ladic coefficient
    sidh_point_mul_scaler(temp_Q, Q, exponent1, E);
    long ladic_coeff = sidh_elliptic_curve_prime_dlp(PP, temp_Q, E, l);

    for (int j = 1; j < e; j++) {
        if (ladic_coeff >= 0) {
            ladic_rep[j - 1] = ladic_coeff;
        } else {
            break;
        }

        mpz_ui_pow_ui(exponent2, l, j - 1);
        mpz_mul_ui(exponent2, exponent2, ladic_rep[j - 1]);
        mpz_divexact_ui(exponent1, exponent1, l);
        sidh_point_mul_scaler(temp_P, P, exponent2, E);
        sidh_point_add(temp_R, temp_R, temp_P, E);
        sidh_point_sub(temp_Q, Q, temp_R, E);
        sidh_point_mul_scaler(temp_Q, temp_Q, exponent1, E);
        ladic_coeff = sidh_elliptic_curve_prime_dlp(PP, temp_Q, E, l);
    }

    if (ladic_coeff >= 0) {
        ladic_rep[e - 1] = ladic_coeff;

        // set x = l_{e - 1}l^{e - 1} + ... + l_1l + l_0
        mpz_set_ui(x, ladic_rep[e - 1]);
        for (long i = e - 2; i >= 0; i--) {
            mpz_mul_ui(x, x, l);
            mpz_add_ui(x, x, ladic_rep[i]);
        }
    } else {
        mpz_set_si(x, -1);
    }

    mpz_clear(exponent1);
    mpz_clear(exponent2);
    sidh_point_clear(temp_P);
    sidh_point_clear(temp_Q);
    sidh_point_clear(temp_R);
    sidh_point_clear(PP);
}

long sidh_elliptic_curve_prime_dlp(const point_t P,
                                   const point_t Q,
                                   const elliptic_curve_t E,
                                   long l) {
    if (sidh_point_is_zero(Q))
        return 0;

    if (sidh_point_equals(P, Q))
        return 1;

    point_t temp;
    sidh_point_init(temp);
    sidh_point_set(temp, P);

    long result = -1;
    for (long i = 2; i < l; i++) {
        sidh_point_add(temp, temp, P, E);
        if (sidh_point_equals(temp, Q)) {
            result = i;
            break;
        }
    }

    sidh_point_clear(temp);
    return result;
}