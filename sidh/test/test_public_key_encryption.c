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
#include <time.h>
#include "sidh_public_key_encryption.h"
#include "sidh_util.h"
#include "sidh_elliptic_curve_dlp.h"
#include "sidh_public_key_validation.h"

void test_public_key_encryption(const public_params_t paramsA,
                                const public_params_t paramsB) {
    sidh_fp_init_chararacteristic(paramsA->characteristic);
    private_key_t Alice_private_key;
    sidh_private_key_init(Alice_private_key);
    sidh_private_key_generate(Alice_private_key, paramsA);

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

    plaintext_t plaintext;
    sidh_public_key_plaintext_init(plaintext);
    plaintext->content = "The quick brown fox jumps over the lazy dog.";
    plaintext->size = strlen(plaintext->content);
    printf("plaintext:\n%s", plaintext->content);
    printf("\n--------------------------------------------\n\n");

    ciphertext_t ciphertext;
    sidh_public_key_ciphertext_init(ciphertext);
    sidh_public_key_pad_plaintext(plaintext, plaintext);
    sidh_public_key_encrypt(ciphertext,
                            plaintext,
                            Alice_public_key,
                            paramsA,
                            paramsB);
    printf("ciphertext using Alice's public-key: \n");
    for (long i = 0; i < ciphertext->size; i++)
        printf("%x", 0xff & ciphertext->content[i]);
    printf("\n--------------------------------------------\n\n");

    plaintext_t decrypted_text;
    sidh_public_key_plaintext_init(decrypted_text);
    sidh_public_key_decrypt(decrypted_text,
                            ciphertext,
                            Alice_private_key,
                            paramsA);
    printf("decrypted text using Alice's private-key:\n");
    for (long i = 0; i < decrypted_text->size; i++)
        printf("%c", decrypted_text->content[i]);
    printf("\n--------------------------------------------\n\n");

    sidh_public_key_ciphertext_clear(ciphertext);
    sidh_private_key_clear(Alice_private_key);
    sidh_public_key_clear(Alice_public_key);
    sidh_point_clear(kernel_gen);
    sidh_public_key_plaintext_clear(plaintext);
    sidh_public_key_plaintext_clear(decrypted_text);
}

int main(int argc, char** argv) {
    char *input_params = "public_params_263";
    if (argc > 1) {
        input_params = argv[1];
    }

    public_params_t paramsA;
    public_params_t paramsB;
    sidh_public_params_init(paramsA);
    sidh_public_params_init(paramsB);

    if (!sidh_public_params_read(paramsA, paramsB, input_params))
        return (EXIT_FAILURE);

    test_public_key_encryption(paramsA, paramsB);

    sidh_public_params_clear(paramsA);
    sidh_public_params_clear(paramsB);

    return (EXIT_SUCCESS);
}

