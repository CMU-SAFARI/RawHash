#include "roptions.h"
#include <limits.h>

void ri_idxopt_init(ri_idxopt_t *opt)
{
	memset(opt, 0, sizeof(ri_idxopt_t));
	opt->e = 6; opt->w = 0; opt->q = 9; opt->lq = 3; opt->n = 0; opt->k = 6, opt->lev_col = 1;
	opt->b = 14;
	opt->mini_batch_size = 50000000;
	opt->batch_size = 4000000000ULL;
}

void ri_mapopt_init(ri_mapopt_t *opt)
{
	memset(opt, 0, sizeof(ri_mapopt_t));

	opt->bp_per_sec = 450;
	opt->sample_rate = 4000;
	opt->chunk_size = 4000;
	opt->sample_per_base = (float)opt->sample_rate / opt->bp_per_sec;

	//seeding
	opt->mid_occ_frac = 1e-2f;
	opt->min_mid_occ = 50;
	opt->max_mid_occ = 500000;
	// opt->q_mid_occ = 10;
	opt->q_occ_frac = 0.01f;

	opt->max_max_occ = 32767;
	opt->occ_dist = 500;

	//chaining
	opt->bw = 4000;
	opt->bw_long = 20000;
	opt->max_target_gap_length = 5000;
	opt->max_query_gap_length = 5000;
	opt->max_chain_iter = 5000;
	opt->max_num_skips = 25;
	opt->min_num_anchors = 2;
	opt->num_best_chains = 3;
	opt->min_chaining_score = 10;
	opt->rmq_inner_dist = 1000;
	opt->rmq_size_cap = 100000;
	opt->chain_gap_scale = 0.8f;
	opt->chain_skip_scale = 0.0f;

	opt->mask_level = 0.5f;
	opt->mask_len = INT_MAX;
	
	opt->pri_ratio = 0.8f;
	opt->best_n = 5;

	opt->alt_drop = 0.15f;

	opt->a = 2, opt->b = 4;
	// opt->q = 4, opt->e = 2, opt->q2 = 24, opt->e2 = 1;
	// opt->sc_ambi = 1;
	// opt->zdrop = 400, opt->zdrop_inv = 200;
	// opt->end_bonus = -1;
	// opt->min_dp_max = opt->min_chaining_score * opt->a;
	// opt->min_ksw_len = 200;
	// opt->anchor_ext_len = 20, opt->anchor_ext_shift = 6;
	// opt->max_clip_ratio = 1.0f;
	opt->mini_batch_size = 500000000;
	// opt->max_sw_mat = 100000000;
	// opt->cap_kalloc = 1000000000;

	// opt->rank_min_len = 500;
	// opt->rank_frac = 0.9f;

	// opt->pe_ori = 0; // FF
	// opt->pe_bonus = 33;

	opt->step_size = 1; //read_seeding_step_size
	opt->min_events = 50;
	opt->max_num_chunk = 20;//max_num_chunks
	opt->min_chain_anchor = 10; //stop_mapping_min_num_anchors

	opt->min_bestmapq = 3;
	opt->min_mapq = 2;
	opt->min_bestmapq_ratio = 2, opt->min_meanmapq_ratio = 5;

	opt->min_bestchain_ratio = 5.0f; //stop_mapping_ratio
	opt->min_meanchain_ratio = 5.0f; //stop_mapping_mean_ratio

	opt->mini_batch_size = 500000000;

	//Default options for event detection.
	opt->window_length1 = 3;
    opt->window_length2 = 6;
    opt->threshold1 = 4.30265f;
    opt->threshold2 = 2.57058f;
    opt->peak_height = 1.0f;

	//TODO: RNA values:
	// opt->window_length1 = 7,
	// opt->window_length2 = 14,
	// opt->threshold1 = 2.5f,
	// opt->threshold2 = 9.0f,
	// opt->peak_height = 1.0f;

	//Sequence until parameters
	opt->t_threshold = 1.5f;
	opt->tn_samples = 5;
	opt->ttest_freq = 500;
	opt->tmin_reads = 500;
}