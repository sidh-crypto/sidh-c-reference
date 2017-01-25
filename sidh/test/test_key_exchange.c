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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sidh_elliptic_curve.h"
#include "sidh_public_param.h"
#include "sidh_isogeny.h"
#include "sidh_private_key.h"
#include "sidh_public_key.h"
#include "sidh_shared_key.h"
#include "sidh_util.h"
#include "sidh_public_key_validation.h"
#include <math.h>
#include <time.h>
#include <stdint.h>

void test_key_exchange(const public_params_t paramsA,
                       const public_params_t paramsB) {
    sidh_fp_init_chararacteristic(paramsA->characteristic);

    private_key_t Alice_private_key;
    sidh_private_key_init(Alice_private_key);
//    sidh_private_key_generate(Alice_private_key, paramsA);
    mpz_set_ui(Alice_private_key->m, 10);
    mpz_set_ui(Alice_private_key->n, 11);
    
    printf("Alice's private key: \n");
    sidh_private_key_print(Alice_private_key);
    printf("\n--------------------------------------------\n\n");

    public_key_t Alice_public_key;
    sidh_public_key_init(Alice_public_key);
    point_t kernel_gen;
    sidh_point_init(kernel_gen);
    sidh_private_key_compute_kernel_gen(kernel_gen,
                                        Alice_private_key,
                                        paramsA->P,
                                        paramsA->Q,
                                        paramsA->le,
                                        paramsA->E);
    sidh_public_key_generate(Alice_public_key, kernel_gen, paramsA, paramsB);
    printf("Alice's public key: \n");
    sidh_public_key_print(Alice_public_key);
    printf("\n--------------------------------------------\n\n");

    private_key_t Bob_private_key;
    sidh_private_key_init(Bob_private_key);
//    sidh_private_key_generate(Bob_private_key, paramsB);
    mpz_set_ui(Bob_private_key->m, 105);
    mpz_set_ui(Bob_private_key->n, 11945);
    printf("Bob's private key: \n");
    sidh_private_key_print(Bob_private_key);
    printf("\n--------------------------------------------\n\n");

    public_key_t Bob_public_key;
    sidh_public_key_init(Bob_public_key);
    sidh_private_key_compute_kernel_gen(kernel_gen,
                                        Bob_private_key,
                                        paramsB->P,
                                        paramsB->Q,
                                        paramsB->le,
                                        paramsB->E);
    sidh_public_key_generate(Bob_public_key, kernel_gen, paramsB, paramsA);
    printf("Bob's public key: \n");
    sidh_public_key_print(Bob_public_key);
    printf("\n--------------------------------------------\n\n");

    fp2_element_t Alice_shared_key;
    sidh_fp2_init(Alice_shared_key);
    sidh_shared_key_generate(Alice_shared_key,
                             Bob_public_key,
                             Alice_private_key,
                             paramsA);
    printf("Alice's shared key: \n");
    printf("%s\n", sidh_fp2_get_str(Alice_shared_key));
    printf("\n--------------------------------------------\n\n");

    fp2_element_t Bob_shared_key;
    sidh_fp2_init(Bob_shared_key);
    sidh_shared_key_generate(Bob_shared_key,
                             Alice_public_key,
                             Bob_private_key,
                             paramsB);
    printf("Bob's shared key: \n");
    printf("%s\n", sidh_fp2_get_str(Bob_shared_key));
    printf("\n--------------------------------------------\n\n");

    long prime_size = (mpz_sizeinbase(characteristic, 2) + 7) / 8;
    long key_size = 12 * prime_size;
    uint8_t keyA[key_size];
    uint8_t keyB[key_size];
    sidh_public_key_to_bytes(keyA, Alice_public_key, prime_size);
    sidh_public_key_to_bytes(keyB, Bob_public_key, prime_size);
    
    for (long i = 0; i < key_size; i++)
        printf("%02X", keyA[i]);
    printf("\n");
    for (long i = 0; i < key_size; i++)
        printf("%02X", keyB[i]);
    printf("\n\n");
    
    public_key_t temp;
    sidh_public_key_init(temp);
    sidh_bytes_to_public_key(temp, keyA, prime_size);
    sidh_public_key_print(temp);
        
    sidh_private_key_clear(Alice_private_key);
    sidh_public_key_clear(Alice_public_key);
    sidh_private_key_clear(Bob_private_key);
    sidh_public_key_clear(Bob_public_key);
    sidh_point_clear(kernel_gen);
    sidh_fp2_clear(Alice_shared_key);
    sidh_fp2_clear(Bob_shared_key);
}

int main(int argc, char** argv) {
    char *input_params = "public_params_521";
    if (argc > 1) {
        input_params = argv[1];
    }

    public_params_t paramsA;
    public_params_t paramsB;
    sidh_public_params_init(paramsA);
    sidh_public_params_init(paramsB);

    if (!sidh_public_params_read(paramsA, paramsB, input_params))
        return (EXIT_FAILURE);

    printf("Public parameters:\n");
    sidh_public_params_print(paramsA, 0);
    sidh_public_params_print(paramsB, 1);
    printf("\n--------------------------------------------\n\n");

    test_key_exchange(paramsA, paramsB);

    sidh_public_params_clear(paramsA);
    sidh_public_params_clear(paramsB);

    return (EXIT_SUCCESS);
}

