#ifndef RSKETCH_H
#define RSKETCH_H

#include <string.h> //for memset
#include "rutils.h"

//A large float value can be used for masking if needed
#define RI_MASK_SIGNAL 3.402823466e+32F

#define LAST_SIG_DIFF 0.3F

//To store sketches in vectors
#define RI_HASH_SHIFT 6
#define RI_ID_SHIFT 32
#define RI_POS_SHIFT 1

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Generate sketches (hash values) by quantizing and concatanating normalized event values
 *
 * @param km       thread-local memory pool; using NULL falls back to malloc()
 * @param s_values normalized event values 
 * @param id       ID of the signals (e.g., ref ID or read ID); will be copied to the output $p array
 * @param strand   strand of signal values. Either 1 (forw) or 0 (reverse); will be copied to the output $p array
 * @param len      length of $s_values
 * @param w        if w>0, find a minimizer for every $w consecutive k-mers
 * @param e        number of packed events in a single 32/64-bit hash value
 * @param n        [Currently disabled] if n>0 Enables BLEND. Number of neighbors hash values to use with SimHash.
 * 				   Generates a single hash value from $n hash values; will be stored as a hash value if enabled.
 * 				   Please see the BLEND paper for details of this hash value calculation.
 * @param q        most significant Q bits of normalized event values to use in the quantization step
 * @param lq       least significant lq bits within the most significant Q bits extracted from the normalized event values.
 * 				   If you want to use n as described in the RawHash manuscript, calculate lq as follows: $lq = $Q - 2 - n
 * @param k        k-mer size that a single event represents (usually given in the k-mer model file)
 * @param p        vector of seeds/sketches
 *                 p->a[i].x = hash value<<RI_HASH_SHIFT | span
 *                 p->a[i].y = $id<<RI_ID_SHIFT | lastPos<<RI_POS_SHIFT | $strand
 *                 where lastPos is the position of the i-th signal,
 * 				   span indicates the region that the hash value covers ($k+$e-1),
 *                 and strand indicates whether the seed comes from the top or the bottom $strand.
 *                 Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void ri_sketch(void *km, const float* s_values, uint32_t id, int strand, int len, int w, int e, int n, int q, int lq, int k, mm128_v *p);

#ifdef __cplusplus
}
#endif
#endif //RSKETCH_H