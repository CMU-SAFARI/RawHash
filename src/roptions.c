#include "roptions.h"
#include <limits.h>

void ri_idxopt_init(ri_idxopt_t *opt)
{
	memset(opt, 0, sizeof(ri_idxopt_t));
	opt->e = 8; opt->w = 0; opt->q = 4; opt->n = 0; opt->k = 6, opt->lev_col = 1;
	opt->b = 14;
	opt->diff = 0.35f;
	opt->mini_batch_size = 50000000;
	opt->batch_size = 4000000000ULL;

	opt->fine_min = -2.0f;
	opt->fine_max = 2.0f;
	opt->fine_range = 0.4;

	opt->window_length1 = 3; //--seg-window-length1
    opt->window_length2 = 9; //--seg-window-length2
    opt->threshold1 = 4.0f; //--seg-threshold1
    opt->threshold2 = 3.5f; //--seg-threshold2
    opt->peak_height = 0.4f; //--seg-peak_height

	// opt->window_length1 = 3; //--seg-window-length1
    // opt->window_length2 = 6; //--seg-window-length2
    // opt->threshold1 = 1.4f; //--seg-threshold1
    // opt->threshold2 = 9.0f; //--seg-threshold2
    // opt->peak_height = 0.2f; //--seg-peak_height

	opt->bp_per_sec = 450; //--bp-per-sec
	opt->sample_rate = 4000; //--sample-rate
	opt->sample_per_base = (float)opt->sample_rate / opt->bp_per_sec;
}

void ri_mapopt_init(ri_mapopt_t *opt)
{
	memset(opt, 0, sizeof(ri_mapopt_t));

	opt->bp_per_sec = 450; //--bp-per-sec
	opt->sample_rate = 4000; //--sample-rate
	opt->chunk_size = 4000; //--chunk-size
	opt->sample_per_base = (float)opt->sample_rate / opt->bp_per_sec;

	//seeding
	// opt->mid_occ_frac = 75e-4f; //--mid-occ-frac
	opt->mid_occ_frac = 1e-2f; //--mid-occ-frac
	// opt->mid_occ_frac = 5e-3f; //--mid-occ-frac
	opt->q_occ_frac = 1e-2f; //--q-occ-frac
	opt->min_mid_occ = 50; //--q-mid-occ [I1, I2]
	opt->max_mid_occ = 500000; //--q-mid-occ [I1, I2]

	opt->max_max_occ = 32767;
	opt->occ_dist = 500;

	//chaining
	opt->bw = 500; //--bw
	opt->bw_long = 0; //--bw-long
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
	
	opt->pri_ratio = 0.3f;
	opt->best_n = 0;

	opt->top_n_mean = 0; //--top-n-mean

	opt->alt_drop = 0.15f; //--alt-drop

	opt->w_bestmq=0.05f; //--w-bestmq
	opt->w_bestmc=0.6f; //--w-bestmc
	opt->w_bestq=0.35f; //--w-bestq
	opt->w_besta=0.2f; //--w-besta
	opt->w_bestma=0.2f; //--w-bestma
	// opt->w_best2c=0.1f; //--w-best2c
	opt->w_threshold = 0.45f; //--w-threshold

	opt->mini_batch_size = 500000000; //-K

	opt->step_size = 1;
	opt->min_events = 50; //--min-events
	opt->max_num_chunk = 10;//--max-chunks

	opt->min_mapq = 2; //--min-mapq

	//dtw
	opt->dtw_border_constraint = RI_M_DTW_BORDER_CONSTRAINT_SPARSE;
	opt->dtw_fill_method = RI_M_DTW_FILL_METHOD_BANDED;
	opt->dtw_band_radius_frac = 0.10f;
	opt->dtw_match_bonus = 0.4f;
	opt->dtw_min_score = 20.0f;

	//Default options for event detection.
	opt->window_length1 = 3; //--seg-window-length1
    opt->window_length2 = 9; //--seg-window-length2
    opt->threshold1 = 4.0f; //--seg-threshold1
    opt->threshold2 = 3.5f; //--seg-threshold2
    opt->peak_height = 0.4f; //--seg-peak_height

	// opt->window_length1 = 3; //--seg-window-length1
    // opt->window_length2 = 7; //--seg-window-length2
    // opt->threshold1 = 4.0f; //--seg-threshold1
    // opt->threshold2 = 3.0f; //--seg-threshold2
    // opt->peak_height = 0.4f; //--seg-peak_height

	// opt->window_length1 = 3; //--seg-window-length1
    // opt->window_length2 = 6; //--seg-window-length2
    // opt->threshold1 = 1.4f; //--seg-threshold1
    // opt->threshold2 = 9.0f; //--seg-threshold2
    // opt->peak_height = 0.2f; //--seg-peak_height

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