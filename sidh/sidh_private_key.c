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

#include "sidh_private_key.h"
#include "sidh_util.h"
#include "sidh_public_param.h"
#include <stdio.h>

void sidh_private_key_init(private_key_t private_key) {
    mpz_inits(private_key->m, private_key->n, NULL);
}

void sidh_private_key_clear(private_key_t private_key) {
    mpz_clears(private_key->m, private_key->n, NULL);
}

void sidh_private_key_generate(private_key_t private_key,
                               const public_params_t params) {
    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);
    mpz_t seed;
    mpz_init(seed);
    sidh_get_random_mpz(seed);
    gmp_randseed(randstate, seed);

    while (1) {
        mpz_urandomm(private_key->m, randstate, params->le);
        mpz_urandomm(private_key->n, randstate, params->le);

        if (!mpz_divisible_ui_p(private_key->m, params->l))
            break;

        if (!mpz_divisible_ui_p(private_key->n, params->l)) {
            mpz_swap(private_key->m, private_key->n);
            break;
        }
    }

    gmp_randclear(randstate);
    mpz_clear(seed);
}

void sidh_private_key_compute_kernel_gen(point_t gen,
                                         const private_key_t private_key,
                                         const point_t P,
                                         const point_t Q,
                                         const mpz_t le,
                                         const elliptic_curve_t E) {
    mpz_t temp_m;
    mpz_t temp_n;
    mpz_init_set(temp_m, private_key->m);
    mpz_init_set(temp_n, private_key->n);

    point_t result;
    sidh_point_init(result);

    mpz_invert(temp_m, temp_m, le);
    mpz_mul(temp_n, temp_m, temp_n);
    mpz_mod(temp_n, temp_n, le);

    sidh_point_mul_scaler(result, Q, temp_n, E);
    sidh_point_add(result, result, P, E);
    sidh_point_set(gen, result);

    mpz_clears(temp_m, temp_n, NULL);
    sidh_point_clear(result);
}

void sidh_private_key_print(const private_key_t private_key) {
    printf("m: %s\n", mpz_get_str(NULL, 10, private_key->m));
    printf("n: %s\n", mpz_get_str(NULL, 10, private_key->n));
}

