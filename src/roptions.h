#include <stdint.h>
#include <string.h> //for memset
// #include "rindex.h"

#ifndef ROPTIONS_H
#define ROPTIONS_H

#define RI_I_NAIVE		0x1
#define RI_I_MIN		0x2
#define RI_I_BLEND		0x4
#define RI_I_SYNCMER	0x8
#define RI_I_STORE_SIG	0x10
#define RI_I_SIG_TARGET	0x20
#define RI_I_REV_QUERY	0x40

#define RI_M_SEQUENCEUNTIL	0x1
#define RI_M_RMQ			0x2
#define RI_M_HARD_MLEVEL	0x4
#define RI_M_NO_SPAN		0x8
#define RI_M_ALIGN			0x10
#define RI_M_NO_ADAPTIVE	0x20
//DTW related
#define RI_M_DTW_EVALUATE_CHAINS 0x40
#define RI_M_DTW_OUTPUT_CIGAR	0x80
#define RI_M_DTW_LOG_SCORES		0x100
#define RI_M_DISABLE_CHAININGSCORE_FILTERING 0x200
#define RI_M_OUTPUT_CHAINS		0x400
#define RI_M_LOG_ANCHORS		0x800
#define RI_M_LOG_NUM_ANCHORS	0x1000
//Overlapping related
#define RI_M_ALL_CHAINS			0x2000

//Characterization related
#define RI_M_OUT_ALL_CHAINS 0x4000

//DTW related
#define RI_M_DTW_BORDER_CONSTRAINT_GLOBAL	0
#define RI_M_DTW_BORDER_CONSTRAINT_SPARSE	1
#define RI_M_DTW_BORDER_CONSTRAINT_LOCAL	2
#define RI_M_DTW_FILL_METHOD_FULL		0
#define RI_M_DTW_FILL_METHOD_BANDED		1

#ifdef __cplusplus
extern "C" {
#endif

// indexing and mapping options
typedef struct ri_idxopt_s{
	short b, w, e, n, q, lq, k, flag, lev_col;
	int64_t mini_batch_size;
	uint64_t batch_size;

	float diff;

	uint32_t window_length1;
	uint32_t window_length2;
	float threshold1;
	float threshold2;
	float peak_height;
	float sample_per_base;
	uint32_t bp_per_sec;
	uint32_t sample_rate;

} ri_idxopt_t;

// mapping options for a reference
typedef struct ri_mapopt_s{

	//ONT Device specific parameters
	uint32_t bp_per_sec;
	uint32_t sample_rate;
	uint32_t chunk_size;
	float sample_per_base;

	//Seeding parameters
	float mid_occ_frac;
	int32_t min_mid_occ, max_mid_occ;
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
	int min_chaining_score;
	float chain_gap_scale;
	float chain_skip_scale;

	float w_bestq, w_besta, w_bestma, w_bestmq, w_bestmc, w_threshold;

	float mask_level;
	int mask_len;
	float pri_ratio;
	int best_n;

	int top_n_mean;

	float alt_drop;

	//Mapping parameters
	uint32_t step_size;
	uint32_t max_num_chunk;
	// uint32_t min_chain_anchor;

	int min_mapq;

	uint32_t dtw_border_constraint;
	uint32_t dtw_fill_method;
	float dtw_band_radius_frac;
	float dtw_match_bonus;
	float dtw_min_score;

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

	char* model_map_path;
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