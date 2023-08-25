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

	opt->bp_per_sec = 450; //--bp-per-sec
	opt->sample_rate = 4000; //--sample-rate
	opt->chunk_size = 4000; //--chunk-size
	opt->sample_per_base = (float)opt->sample_rate / opt->bp_per_sec;

	//seeding
	opt->mid_occ_frac = 1e-2f; //--mid-occ-frac
	opt->min_mid_occ = 50; //--q-mid-occ [I1, I2]
	opt->max_mid_occ = 500000; //--q-mid-occ [I1, I2]
	// opt->q_mid_occ = 10;
	opt->q_occ_frac = 0.01f; //--q-occ-frac

	opt->max_max_occ = 32767;
	opt->occ_dist = 500;

	//chaining
	opt->bw = 500; //--bw
	opt->bw_long = 5000; //--bw-long
	opt->max_target_gap_length = 2500; //--max-target-gap
	opt->max_query_gap_length = 2500; //--max-query-gap
	opt->max_chain_iter = 200; //--max-iterations
	opt->max_num_skips = 5; //--max-skips
	opt->min_num_anchors = 2; //--min-anchors
	// opt->num_best_chains = 3; //--best-chains
	opt->min_chaining_score = 15; //--min-score
	opt->rmq_inner_dist = 1000; //--rmq-inner-dist
	opt->rmq_size_cap = 100000; //--rmq-size-cap
	opt->chain_gap_scale = 0.8f; //--chain-gap-scale
	opt->chain_skip_scale = 0.0f; //--chain-skip-scale

	opt->mask_level = 0.5f; //--primary-ratio
	opt->mask_len = INT_MAX; //--primary-length
	
	opt->pri_ratio = 0.8f;
	opt->best_n = 5;

	opt->alt_drop = 0.15f; //--alt-drop

	opt->w_bestq=0.1f; //--w-bestq
	opt->w_best2q=0.01f; //--w-best2q
	opt->w_best2c=0.4f; //--w-best2c
	opt->w_bestmq=0.01f; //--w-bestmq
	opt->w_bestmc=0.48f; //--w-bestmc
	opt->w_threshold = 0.6f; //--w-threshold

	opt->a = 2, opt->b = 4; //--chain-match-score, 
	// opt->q = 4, opt->e = 2, opt->q2 = 24, opt->e2 = 1;
	// opt->sc_ambi = 1;
	// opt->zdrop = 400, opt->zdrop_inv = 200;
	// opt->end_bonus = -1;
	// opt->min_dp_max = opt->min_chaining_score * opt->a;
	// opt->min_ksw_len = 200;
	// opt->anchor_ext_len = 20, opt->anchor_ext_shift = 6;
	// opt->max_clip_ratio = 1.0f;
	opt->mini_batch_size = 500000000; //-K
	// opt->max_sw_mat = 100000000;
	// opt->cap_kalloc = 1000000000;

	// opt->rank_min_len = 500;
	// opt->rank_frac = 0.9f;

	// opt->pe_ori = 0; // FF
	// opt->pe_bonus = 33;

	opt->step_size = 1; //read_seeding_step_size
	opt->min_events = 50; //--min-events
	opt->max_num_chunk = 15;//--max-chunks
	// opt->min_chain_anchor = 10;

	opt->min_bestmapq = 5; //--min-bestmapq
	opt->min_mapq = 2; //--min-mapq
	opt->min_bestmapq_ratio = 2.0f, opt->min_meanmapq_ratio = 5.0f; //--min-bestmapq-ratio, --min-meanmapq-ratio

	opt->min_bestchain_ratio = 2.5f; //--min-bestchain-ratio
	opt->min_meanchain_ratio = 2.5f; //--min-meanchain-ratio

	//Default options for event detection.
	opt->window_length1 = 3; //--seg-window-length1
    opt->window_length2 = 6; //--seg-window-length2
    opt->threshold1 = 4.30265f; //--seg-threshold1
    opt->threshold2 = 2.57058f; //--seg-threshold2
    opt->peak_height = 1.0f; //--seg-peak_height

	//TODO: RNA values:
	// opt->window_length1 = 7,
	// opt->window_length2 = 14,
	// opt->threshold1 = 2.5f,
	// opt->threshold2 = 9.0f,
	// opt->peak_height = 1.0f;

	//Sequence until parameters
	opt->t_threshold = 1.5f; //--threshold
	opt->tn_samples = 5; //--n-samples
	opt->ttest_freq = 500; //--test-frequency
	opt->tmin_reads = 500; //--min-reads
}