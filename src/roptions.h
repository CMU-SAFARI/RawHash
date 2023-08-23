#include <stdint.h>
#include <string.h> //for memset
// #include "rindex.h"

#ifndef ROPTIONS_H
#define ROPTIONS_H

#define RI_I_NAIVE		0x1
#define RI_I_MIN		0x2
#define RI_I_BLEND		0x4
#define RI_I_SYNCMER	0x8

#define RI_M_SEQUENCEUNTIL	0x1
#define RI_M_RMQ			0x2
#define RI_M_HARD_MLEVEL	0x4

#ifdef __cplusplus
extern "C" {
#endif

// indexing and mapping options
typedef struct ri_idxopt_s{
	short b, w, e, n, q, lq, k, flag, lev_col;
	int64_t mini_batch_size;
	uint64_t batch_size;
} ri_idxopt_t;

typedef struct ri_mapopt_s{

	//ONT Device specific parameters
	uint32_t bp_per_sec;
	uint32_t sample_rate;
	uint32_t chunk_size;
	float sample_per_base;

	//Seeding parameters
	float mid_occ_frac;
	int32_t min_mid_occ, max_mid_occ;
	// int32_t q_mid_occ;
	float q_occ_frac;

	int32_t mid_occ;     // ignore seeds with occurrences above this threshold
	int32_t max_occ, max_max_occ, occ_dist;

	//Chaining parameters
	int min_events;
	int bw;
	int bw_long;
	int max_target_gap_length;
	int max_query_gap_length;
	int max_chain_iter;
	int rmq_inner_dist;
	int rmq_size_cap;
	int max_num_skips;
	int min_num_anchors;
	// uint32_t num_best_chains;
	int min_chaining_score;
	float chain_gap_scale;
	float chain_skip_scale;

	float mask_level;
	int mask_len;
	float pri_ratio;
	int best_n;      // top best_n chains are subjected to DP alignment

	float alt_drop;

	int a, b;
	// int q, e, q2, e2; // matching score, mismatch, gap-open and gap-ext penalties
	// int sc_ambi; // score when one or both bases are "N"
	// int noncan;      // cost of non-canonical splicing sites
	// int junc_bonus;
	// int zdrop, zdrop_inv;   // break alignment if alignment score drops too fast along the diagonal
	// int end_bonus;
	// int min_dp_max;  // drop an alignment if the score of the max scoring segment is below this threshold
	// int min_ksw_len;
	// int anchor_ext_len, anchor_ext_shift;
	// float max_clip_ratio; // drop an alignment if BOTH ends are clipped above this ratio

	// int rank_min_len;
	// float rank_frac;

	// int pe_ori, pe_bonus;

	// int64_t max_sw_mat;
	// int64_t cap_kalloc;

	//Mapping parameters
	uint32_t step_size;
	uint32_t max_num_chunk;
	// uint32_t min_chain_anchor;

	int min_mapq, min_bestmapq;
	float min_bestmapq_ratio, min_meanmapq_ratio;

	float min_bestchain_ratio, min_meanchain_ratio;

	float t_threshold;
	uint32_t tn_samples;
	uint32_t ttest_freq;
	uint32_t tmin_reads;
	
	int64_t flag;    // see ri_F_* macros
	int64_t mini_batch_size; // size of a batch of query bases to process in parallel

	//Event detector options
	uint32_t window_length1;
	uint32_t window_length2;
	float threshold1;
	float threshold2;
	float peak_height;
} ri_mapopt_t;

/**
 * Initializes the default indexing options
 *
 * @param opt	pointer to the indexing options
 * 
 */
void ri_idxopt_init(ri_idxopt_t *opt);

/**
 * Initializes the default mapping options
 *
 * @param opt	pointer to the mapping options
 * 
 */
void ri_mapopt_init(ri_mapopt_t *opt);

#ifdef __cplusplus
}
#endif
#endif //ROPTIONS_H