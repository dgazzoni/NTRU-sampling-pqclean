#ifndef PQCLEAN_NTRUHPS2048677SHUFFLING_CLEAN_API_H
#define PQCLEAN_NTRUHPS2048677SHUFFLING_CLEAN_API_H

#include <stdint.h>

#define PQCLEAN_NTRUHPS2048677SHUFFLING_CLEAN_CRYPTO_SECRETKEYBYTES 1234
#define PQCLEAN_NTRUHPS2048677SHUFFLING_CLEAN_CRYPTO_PUBLICKEYBYTES 930
#define PQCLEAN_NTRUHPS2048677SHUFFLING_CLEAN_CRYPTO_CIPHERTEXTBYTES 930
#define PQCLEAN_NTRUHPS2048677SHUFFLING_CLEAN_CRYPTO_BYTES 32

#define PQCLEAN_NTRUHPS2048677SHUFFLING_CLEAN_CRYPTO_ALGNAME "ntruhps2048677"

int PQCLEAN_NTRUHPS2048677SHUFFLING_CLEAN_crypto_kem_keypair(uint8_t *pk,
                                                             uint8_t *sk);

int PQCLEAN_NTRUHPS2048677SHUFFLING_CLEAN_crypto_kem_enc(uint8_t *c, uint8_t *k,
                                                         const uint8_t *pk);

int PQCLEAN_NTRUHPS2048677SHUFFLING_CLEAN_crypto_kem_dec(uint8_t *k,
                                                         const uint8_t *c,
                                                         const uint8_t *sk);

#endif
