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


#include "sidh_public_key_validation.h"
#include "sidh_elliptic_curve_dlp.h"
#include <stdio.h>

int sidh_public_key_is_valid(const public_key_t public_key,
                             const public_params_t params) {
    if (!sidh_public_key_check_order(public_key->P, public_key->E, params))
        return 0;

    if (!sidh_public_key_check_order(public_key->Q, public_key->E, params))
        return 0;

    if (!sidh_public_key_check_dependency(public_key, params))
        return 0;

    if (!sidh_public_key_check_curve(public_key->E))
        return 0;

    return 1;
}

int sidh_public_key_check_order(const point_t P,
                                const elliptic_curve_t E,
                                const public_params_t params) {
    mpz_t order;
    point_t temp;

    mpz_init_set(order, params->le);
    sidh_point_init(temp);

    int result = 0;
    mpz_divexact_ui(order, order, params->l);
    sidh_point_mul_scaler(temp, P, order, E);
    if (!sidh_point_is_zero(temp)) {
        sidh_point_mul_scaler_si(temp, temp, params->l, E);
        if (sidh_point_is_zero(temp))
            result = 1;
    }

    mpz_clear(order);
    sidh_point_clear(temp);
    return result;
}

int sidh_public_key_check_dependency(const public_key_t public_key,
                                     const public_params_t params) {
    mpz_t x;
    mpz_init(x);

    int result = 0;
    sidh_elliptic_curve_prime_power_dlp(x,
                                        public_key->P,
                                        public_key->Q,
                                        public_key->E,
                                        params->l,
                                        params->e);

    if (mpz_cmp_si(x, -1) == 0) {
        sidh_elliptic_curve_prime_power_dlp(x,
                                            public_key->Q,
                                            public_key->P,
                                            public_key->E,
                                            params->l,
                                            params->e);
        if (mpz_cmp_si(x, -1) == 0)
            result = 1;
    }

    mpz_clear(x);
    return result;
}

int sidh_public_key_check_curve(const elliptic_curve_t E) {
    point_t temp;
    mpz_t exponent;

    sidh_point_init(temp);
    mpz_init_set(exponent, characteristic);
    mpz_add_ui(exponent, exponent, 1);

    sidh_elliptic_curve_random_point(temp, E);
    sidh_point_mul_scaler(temp, temp, exponent, E);
    int result = sidh_point_is_zero(temp);

    sidh_point_clear(temp);
    mpz_clear(exponent);

    return result;
}