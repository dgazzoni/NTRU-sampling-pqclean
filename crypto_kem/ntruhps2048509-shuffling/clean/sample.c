// clang-format off

#include "sample.h"

void PQCLEAN_NTRUHPS2048509SHUFFLING_CLEAN_sample_fg(
    poly *f, poly *g, const unsigned char uniformbytes[NTRU_SAMPLE_FG_BYTES]) {

  PQCLEAN_NTRUHPS2048509SHUFFLING_CLEAN_sample_iid(f, uniformbytes);
  PQCLEAN_NTRUHPS2048509SHUFFLING_CLEAN_sample_fixed_type(
      g, uniformbytes + NTRU_SAMPLE_IID_BYTES);
}

void PQCLEAN_NTRUHPS2048509SHUFFLING_CLEAN_sample_rm(
    poly *r, poly *m, const unsigned char uniformbytes[NTRU_SAMPLE_RM_BYTES]) {

  PQCLEAN_NTRUHPS2048509SHUFFLING_CLEAN_sample_iid(r, uniformbytes);
  PQCLEAN_NTRUHPS2048509SHUFFLING_CLEAN_sample_fixed_type(
      m, uniformbytes + NTRU_SAMPLE_IID_BYTES);
}

// Unlike Algorithm 5 in the paper, we use the batch random number generation idea, discussed in Section 4. Moreover,
// although not required in non-SIMD versions, we use disjoint ranges in the random number array for the initial
// sampling and the rejection sampling fixup procedure, as discussed in "SIMD implementation of Algorithm 5" in
// Section 4, to ensure interoperability of KATs between the scalar and SIMD versions.
void PQCLEAN_NTRUHPS2048509SHUFFLING_CLEAN_rejsamplingmod(
  int16_t shuffle_indices[], const uint16_t u[NTRU_SAMPLE_FT_BYTES / 2]) {
  int i, j = NTRU_N - 1;

  for (i = 0; i < NTRU_N - 1; i++)
  {
    uint32_t m;
    uint16_t s, t, l;

    s = NTRU_N - 1 - i;
    t = 65536 % s;

    m = (uint32_t)u[i] * s;
    l = m;

    while (l < t)
    {
      m = (uint32_t)u[j++] * s;
      l = m;
    }

    shuffle_indices[i] = m >> 16;
  }
}

void PQCLEAN_NTRUHPS2048509SHUFFLING_CLEAN_sample_fixed_type(
  poly *r, const unsigned char u[NTRU_SAMPLE_FT_BYTES]) {
  int16_t shuffle_indices[NTRU_N - 1];

  PQCLEAN_NTRUHPS2048509SHUFFLING_CLEAN_rejsamplingmod(shuffle_indices, (const uint16_t*)u);

  int c0 = NTRU_N - 1 - NTRU_WEIGHT, c1 = NTRU_WEIGHT / 2;

  for (int i = 0; i < NTRU_N - 1; i++)
  {
    int16_t p = shuffle_indices[i];
    if (p < c0)
    {
      r->coeffs[i] = 0;
      c0--;
    }
    else if (p < c0 + c1)
    {
      r->coeffs[i] = 1;
      c1--;
    }
    else
    {
      r->coeffs[i] = 2;
    }
  }

  r->coeffs[NTRU_N-1] = 0;
}

// clang-format on
