#ifndef RSEED_H
#define RSEED_H

#include <string.h> //for memset
#include "rindex.h"
#include "rutils.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RI_SEED_TANDEM (1ULL<<38) //32+6
#define RI_SEED_SELF (1ULL<<39)
#define RI_SEED_SEG_SHIFT 40
#define RI_SEED_SEG_MASK   (0xffULL<<(RI_SEED_SEG_SHIFT))

typedef struct {
	uint32_t n;
	uint32_t q_pos;
	uint32_t q_span:31, flt:1;
	uint32_t seg_id:31, is_tandem:1;
	const uint64_t *cr;
} ri_seed_t;

void ri_seed_mz_flt(void *km, mm128_v *riv, int32_t q_occ_max, float q_occ_frac);

/**
 * @brief Collects all seeds that match the given query vector and index.
 *
 * @param km Pointer to memory pool.
 * @param _n_m Pointer to the number of collected seeds.
//  * @param qlen Length of the query vector.
//  * @param max_occ Maximum number of occurrences.
//  * @param max_max_occ Maximum number of maximum occurrences.
//  * @param dist Distance between seeds.
 * @param ri Pointer to the index.
 * @param riv Pointer to the query vector.
 * @param n_a Pointer to the number of collected seeds.
 * @param rep_len Pointer to the length of the repetition.
 * @param n_seed_pos Pointer to the number of seed positions.
 * @param seed_pos Pointer to the seed positions.
 *
 * @return Pointer to the collected seeds.
 */
ri_seed_t *ri_collect_matches(void *km,
                              int *_n_seed_m,
                              uint32_t qlen,
                              int max_occ,
                              int max_max_occ,
                              int dist,
                              const ri_idx_t *ri,
                              const mm128_v *riv,
                              int64_t *n_seed_pos,
                              int *rep_len);
                            //   int *n_seed_mini,
                            //   uint64_t **seed_mini)

#ifdef __cplusplus
}
#endif
#endif //RSEED_H