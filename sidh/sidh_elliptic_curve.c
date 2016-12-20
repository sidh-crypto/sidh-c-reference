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


#include <stdlib.h>

#include "sidh_elliptic_curve.h"
#include "sidh_util.h"
#include <string.h>

void sidh_elliptic_curve_init(elliptic_curve_t E) {
    sidh_fp2_init_set_si(E->a, 0, 1);
    sidh_fp2_init_set_si(E->b, 0, 1);
}

void sidh_elliptic_curve_set(elliptic_curve_t E,
                             const elliptic_curve_t T) {
    sidh_fp2_set(E->a, T->a);
    sidh_fp2_set(E->b, T->b);
}

void sidh_elliptic_curve_set_coeffs(elliptic_curve_t E,
                                    const fp2_element_t a,
                                    const fp2_element_t b) {
    sidh_fp2_set(E->a, a);
    sidh_fp2_set(E->b, b);
}

void sidh_point_init(point_t P) {
    sidh_fp2_init(P->x);
    sidh_fp2_init(P->y);
    sidh_point_zero(P);
}

void sidh_point_set_coordinates(point_t P,
                                const fp2_element_t x,
                                const fp2_element_t y,
                                int z) {
    sidh_fp2_set(P->x, x);
    sidh_fp2_set(P->y, y);
    P->z = z;
}

void sidh_point_set(point_t P,
                    const point_t Q) {
    sidh_point_set_coordinates(P, Q->x, Q->y, Q->z);
}

void sidh_point_zero(point_t P) {
    sidh_fp2_zero(P->x);
    sidh_fp2_one(P->y);
    P->z = 0;
}

int sidh_point_is_zero(const point_t P) {
    return P->z == 0;
}

void sidh_point_negate(point_t P,
                       const point_t Q) {
    sidh_point_set(P, Q);
    sidh_fp2_negate(P->y, P->y);
}

int sidh_point_has_order_2(const point_t P) {
    return sidh_fp2_is_zero(P->y);
}

void sidh_elliptic_curve_clear(elliptic_curve_t E) {
    sidh_fp2_clear(E->a);
    sidh_fp2_clear(E->b);
}

void sidh_point_clear(point_t P) {
    sidh_fp2_clear(P->x);
    sidh_fp2_clear(P->y);
}

int sidh_point_equals(const point_t P,
                      const point_t Q) {
    return sidh_fp2_equals(P->x, Q->x) &&
           sidh_fp2_equals(P->y, Q->y) &&
           (P->z == Q->z);
}

char *sidh_elliptic_curve_get_str(const elliptic_curve_t E) {
    char *result = "";
    result = sidh_concat(result, "y^2 = x^3");
    if (!sidh_fp2_is_zero(E->a)) {
        result = sidh_concat(result, " + (");
        result = sidh_concat(result, sidh_fp2_get_str(E->a));
        result = sidh_concat(result, ")");
        result = sidh_concat(result, " * x");
    }

    if (!sidh_fp2_is_zero(E->b)) {
        result = sidh_concat(result, " + (");
        result = sidh_concat(result, sidh_fp2_get_str(E->b));
        result = sidh_concat(result, ")");
    }

    return result;
}

char *sidh_point_get_str(const point_t P) {
    char *result = "";
    result = sidh_concat(result, "(");
    result = sidh_concat(result, sidh_fp2_get_str(P->x));
    result = sidh_concat(result, " : ");
    result = sidh_concat(result, sidh_fp2_get_str(P->y));
    result = sidh_concat(result, " : ");
    result = sidh_concat(result, (P->z == 1 ? "1" : "0"));
    result = sidh_concat(result, ")");

    return result;
}

void sidh_point_add_with_lambda(point_t R,
                                const point_t P,
                                const point_t Q,
                                const fp2_element_t lambda) {
    point_t result;
    sidh_point_init(result);
    result->z = 1;

    // x_R = lambda^2 - x_P - x_Q
    sidh_fp2_square(result->x, lambda);
    sidh_fp2_sub(result->x, result->x, P->x);
    sidh_fp2_sub(result->x, result->x, Q->x);

    // y_R = lambda * (x_P - x_R) - y_P
    sidh_fp2_sub(result->y, P->x, result->x);
    sidh_fp2_mul(result->y, result->y, lambda);
    sidh_fp2_sub(result->y, result->y, P->y);
    sidh_point_set(R, result);

    sidh_point_clear(result);
}

void sidh_point_double(point_t R,
                       const point_t P,
                       const elliptic_curve_t E) {
    if (sidh_point_is_zero(P)) {
        sidh_point_zero(R);
        return;
    }

    // check if the point is of order 2
    if (sidh_point_has_order_2(P)) {
        sidh_point_zero(R);
        return;
    }

    fp2_element_t temp;
    fp2_element_t lambda;

    sidh_fp2_init(temp);
    sidh_fp2_init(lambda);

    // lambda = (3(x_P)^2 + a) / (2y_p)
    sidh_fp2_square(lambda, P->x);
    sidh_fp2_mul_scaler_si(lambda, lambda, 3);
    sidh_fp2_add(lambda, lambda, E->a);
    sidh_fp2_mul_scaler_si(temp, P->y, 2);
    sidh_fp2_div(lambda, lambda, temp);

    sidh_point_add_with_lambda(R, P, P, lambda);

    sidh_fp2_clear(temp);
    sidh_fp2_clear(lambda);
}

void sidh_point_add(point_t R,
                    const point_t P,
                    const point_t Q,
                    const elliptic_curve_t E) {
    if (sidh_point_is_zero(P)) {
        sidh_point_set(R, Q);
        return;
    }

    if (sidh_point_is_zero(Q)) {
        sidh_point_set(R, P);
        return;
    }

    if (sidh_fp2_equals(P->x, Q->x)) {
        if (sidh_fp2_equals(P->y, Q->y)) {
            sidh_point_double(R, P, E);
            return;
        }

        sidh_point_zero(R);
        return;
    }

    fp2_element_t temp;
    fp2_element_t lambda;

    sidh_fp2_init(temp);
    sidh_fp2_init(lambda);

    // lambda = (y_Q - y_P) / (x_Q - x_P)
    sidh_fp2_sub(lambda, Q->y, P->y);
    sidh_fp2_sub(temp, Q->x, P->x);
    sidh_fp2_div(lambda, lambda, temp);

    sidh_point_add_with_lambda(R, P, Q, lambda);

    sidh_fp2_clear(temp);
    sidh_fp2_clear(lambda);
}

void sidh_point_sub(point_t R,
                    const point_t P,
                    const point_t Q,
                    const elliptic_curve_t E) {
    point_t temp;
    sidh_point_init(temp);
    sidh_point_negate(temp, Q);
    sidh_point_add(R, P, temp, E);
    sidh_point_clear(temp);
}

void sidh_point_mul_scaler(point_t R,
                           const point_t P,
                           const mpz_t scaler,
                           const elliptic_curve_t E) {
    if (mpz_cmp_ui(scaler, 0) == 0) {
        sidh_point_zero(R);
        return;
    }

    if (mpz_cmp_ui(scaler, 1) == 0) {
        sidh_point_set(R, P);
        return;
    }

    point_t R0;
    point_t R1;

    sidh_point_init(R0);
    sidh_point_init(R1);
    sidh_point_set(R1, P);

    long num_bits = mpz_sizeinbase(scaler, 2);
    for (long i = 0; i < num_bits; i++) {
        if (mpz_tstbit(scaler, i) == 1)
            sidh_point_add(R0, R0, R1, E);
        sidh_point_double(R1, R1, E);
    }

    if (mpz_sgn(scaler) < 0)
        sidh_point_negate(R0, R0);

    sidh_point_set(R, R0);
    sidh_point_clear(R0);
    sidh_point_clear(R1);
}

void sidh_point_mul_scaler_si(point_t R,
                              const point_t P,
                              long scaler,
                              const elliptic_curve_t E) {
    mpz_t temp;
    mpz_init_set_si(temp, scaler);
    sidh_point_mul_scaler(R, P, temp, E);
    mpz_clear(temp);
}

void sidh_elliptic_curve_compute_j_inv(fp2_element_t j_inv,
                                       const elliptic_curve_t E) {
    fp2_element_t result;
    fp2_element_t temp;
    sidh_fp2_init(result);
    sidh_fp2_init(temp);

    sidh_fp2_pow_ui(temp, E->a, 3);
    sidh_fp2_mul_scaler_si(temp, temp, 4);
    sidh_fp2_square(result, E->b);
    sidh_fp2_mul_scaler_si(result, result, 27);
    sidh_fp2_add(result, result, temp);
    sidh_fp2_inv(result, result);
    sidh_fp2_mul(result, result, temp);
    sidh_fp2_mul_scaler_si(result, result, 1728);
    sidh_fp2_set(j_inv, result);

    sidh_fp2_clear(result);
    sidh_fp2_clear(temp);
}

int sidh_point_is_on_curve(const point_t P,
                           const elliptic_curve_t E) {

    if (sidh_point_is_zero(P))
        return 1;

    fp2_element_t temp_x;
    sidh_fp2_init(temp_x);

    // compute x^3 + a * x + b = x * (x^2 + a) + b
    sidh_fp2_square(temp_x, P->x);
    sidh_fp2_add(temp_x, temp_x, E->a);
    sidh_fp2_mul(temp_x, temp_x, P->x);
    sidh_fp2_add(temp_x, temp_x, E->b);

    fp2_element_t temp_y;
    sidh_fp2_init(temp_y);
    sidh_fp2_square(temp_y, P->y);

    int result = sidh_fp2_equals(temp_y, temp_x);

    sidh_fp2_clear(temp_x);
    sidh_fp2_clear(temp_y);

    return result;
}

void sidh_elliptic_curve_random_point(point_t P,
                                      const elliptic_curve_t E) {
    point_t result;
    sidh_point_init(result);
    result->z = 1;

    fp2_element_t temp_x;
    sidh_fp2_init(temp_x);

    fp2_element_t temp_y;
    sidh_fp2_init(temp_y);

    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);

    while (1) {
        sidh_fp2_random(result->x, randstate);

        // compute x^3 + a * x + b = x * (x^2 + a) + b
        sidh_fp2_square(temp_x, result->x);
        sidh_fp2_add(temp_x, temp_x, E->a);
        sidh_fp2_mul(temp_x, temp_x, result->x);
        sidh_fp2_add(temp_x, temp_x, E->b);

        if (sidh_fp2_is_square(temp_x)) {
            sidh_fp2_sqrt(result->y, temp_x);
            break;
        }
    }

    sidh_point_set(P, result);

    sidh_point_clear(result);
    sidh_fp2_clear(temp_x);
    sidh_fp2_clear(temp_y);
    gmp_randclear(randstate);
}
