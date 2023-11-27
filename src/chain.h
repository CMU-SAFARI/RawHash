#ifndef RCHAIN_H
#define RCHAIN_H

#include <string.h> //for memset
#include "rutils.h"
#include "rindex.h"
#include "roptions.h"
#include "rseed.h"

#define RH_PARENT_UNSET   (-1)
#define RH_PARENT_TMP_PRI (-2)

#define MM_MAX_SEG       255

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint32_t capacity;                  // the capacity of cigar[]
	int32_t dp_score, dp_max, dp_max2;  // DP score; score of the max-scoring segment; score of the best alternate mappings
	uint32_t n_ambi:30, trans_strand:2; // number of ambiguous bases; transcript strand: 0 for unknown, 1 for +, 2 for -
	uint32_t n_cigar;                   // number of cigar operations in cigar[]
	uint32_t cigar[];
} mm_extra_t;

typedef struct {
	int32_t id;             // ID for internal uses (see also parent below)
	int32_t cnt;            // number of minimizers; if on the reverse strand
	int32_t rid;            // reference index; if this is an alignment from inversion rescue
	int32_t score;          // DP alignment score
	int32_t qs, qe, rs, re; // query start and end; reference start and end
	int32_t parent, subsc;  // parent==id if primary; best alternate mapping score
	int32_t as;             // offset in the a[] array (for internal uses only)
	int32_t mlen, blen;     // seeded exact match length; seeded alignment block length
	int32_t n_sub;          // number of suboptimal mappings
	int32_t score0;         // initial chaining score (before chain merging/spliting)
	uint32_t mapq:8, split:2, rev:1, inv:1, sam_pri:1, proper_frag:1, pe_thru:1, seg_split:1, seg_id:8, split_inv:1, is_alt:1, strand_retained:1, dummy:5;
	uint32_t hash;
	float div;
	mm_extra_t *p;
} mm_reg1_t;

typedef struct {
	int n_u, n_a;
	uint64_t *u;
	mm128_t *a;
} ri_seg_t;

// memory buffer for thread-local storage during mapping
typedef struct ri_tbuf_s {
	void *km;
	int rep_len, frag_gap;
}ri_tbuf_t;

/**
 * Create the chain from a list of anchors using the DP approach (more accurate but slower than RMQ)
 *
 * @param max_dist_t		max distance in the target (reference)
 * @param max_dist_q		max distance in the query (read)
 * @param bw				The absolute difference in gaps between two anchors: abs(target_gap - query_gap) t_gap = anchort0 - anchort1
 *							i.e., difference in diagonals
 * @param max_skip			The number of anchors we can skip when going back from an anchor to find its predecessor
 * @param max_iter			Max distance between an anchor its potential predecessor in the chain in terms their indices
 * @param min_cnt			Minimum number of anchors in a chain
 * @param min_sc    		Minimum score for a chain to be considered for backtracking
 * @param chn_pen_gap     	Gap penalty for the chain
 * @param chn_pen_skip    	Anchor skip penalty for the chain
 * @param n					Number of anchors
 * @param a					Anchors
 * 							a[].x: rev<<63 | tid<<32 | t_pos
 * 							a[].y: flags<<40 | q_span<<32 | q_pos
 * @param n_u_				Number of chains
 * @param u					Chains (output)
 * 							u[].x: score<<32 | #anchors (sum of lower 32 bits of u[] is the returned length of a[])
 * @param km				Memory pool
 * 
 * @return					Chains
 */
mm128_t *mg_lchain_dp(int max_dist_t,
                      int max_dist_q,
                      int bw,
                      int max_skip,
                      int max_iter,
                      int min_cnt,
                      int min_sc,
                      float chn_pen_gap,
                      float chn_pen_skip,
					//   int is_cdna, 
                    //   int n_seg, 
                      int64_t *n, 
                      mm128_t *a, 
					  mm128_t **_a,
                      int *n_u_,  
                      uint64_t **_u, 
                      void *km);
/**
 * Create the chain from a list of anchors using the Range Minimum Query (RMQ) technique (i.e., less accurate but faster than DP)
 *
 * @param max_dist			max distance in the target and query
 * @param max_dist_inner	
 * @param bw				The absolute difference in gaps between two anchors: abs(target_gap - query_gap) t_gap = anchort0 - anchort1
 *							i.e., difference in diagonals
 * @param max_skip			The number of anchors we can skip when going back from an anchor to find its predecessor
 * @param cap_rmq_size		
 * @param min_cnt			Minimum number of anchors in a chain
 * @param min_sc    		Minimum score for a chain to be considered for backtracking
 * @param chn_pen_gap     	Gap penalty for the chain
 * @param chn_pen_skip    	Anchor skip penalty for the chain
 * @param n					Number of anchors
 * @param a					Anchors
 * 							a[].x: rev<<63 | tid<<32 | t_pos
 * 							a[].y: flags<<40 | q_span<<32 | q_pos
 * @param n_u_				Number of chains
 * @param u					Chains (output)
 * 							u[]: score<<32 | #anchors
 * @param km				Memory pool
 * 
 * @return					Chains
 */
mm128_t *mg_lchain_rmq(int max_dist,
					   int max_dist_inner,
					   int bw,
					   int max_skip,
					   int cap_rmq_size,
					   int min_cnt,
					   int min_sc,
					   float chn_pen_gap,
					   float chn_pen_skip,
					   int64_t *n,
					   mm128_t *a,
					   mm128_t **_a,
					   int *n_u_,
					   uint64_t **_u,
					   void *km);

// mm_reg1_t *mm_gen_regs(void *km, uint32_t hash, int qlen, int n_u, uint64_t *u, mm128_t *a, int is_qstrand);
/*
* @brief 		Identify the regions from chains
*
* @param km 	pointer to the memory pool
* @param hash 	hash value
* @param qlen 	length of the query sequence
* @param n_u 	number of chains
* @param u 		Chains (output)
* 					u[]: score<<32 | #anchors
* @param a 		Anchors
* 					a[].x: rev<<63 	 | tid<<32 	  | t_pos
* 					a[].y: flags<<40 | q_span<<32 | q_pos
*
* @return 		pointer to the array of regions
*/
mm_reg1_t *mm_gen_regs(void *km, uint32_t hash, int qlen, int n_u, uint64_t *u, mm128_t *a);
void mm_mark_alt(const ri_idx_t *ri, int n, mm_reg1_t *r);
// void mm_split_reg(mm_reg1_t *r, mm_reg1_t *r2, int n, int qlen, mm128_t *a, int is_qstrand);
void mm_split_reg(mm_reg1_t *r, mm_reg1_t *r2, int n, int qlen, mm128_t *a);
// void mm_set_parent(void *km, float mask_level, int mask_len, int n, mm_reg1_t *r, int sub_diff, int hard_mask_level, 
void mm_set_parent(void *km, float mask_level, int mask_len, int n, mm_reg1_t *r, int hard_mask_level, float alt_diff_frac);
void mm_hit_sort(void *km, int *n_regs, mm_reg1_t *r, float alt_diff_frac);
int mm_set_sam_pri(int n, mm_reg1_t *r);
void mm_sync_regs(void *km, int n_regs, mm_reg1_t *regs);
// void mm_select_sub(void *km, float pri_ratio, int min_diff, int best_n, int check_strand, int min_strand_sc, int *n_, mm_reg1_t *r);
void mm_select_sub(void *km, float pri_ratio, int best_n, int check_strand, int min_strand_sc, int *n_, mm_reg1_t *r);
void mm_select_sub_multi(void *km, float pri_ratio, float pri1, float pri2, int max_gap_ref, int min_diff, int best_n, int n_segs, const int *qlens, int *n_, mm_reg1_t *r);
int mm_filter_strand_retained(int n_regs, mm_reg1_t *r);
void mm_filter_regs(const ri_mapopt_t *opt, int qlen, int *n_regs, mm_reg1_t *regs);
int mm_squeeze_a(void *km, int n_regs, mm_reg1_t *regs, mm128_t *a);
ri_seg_t *mm_seg_gen(void *km, uint32_t hash, int n_segs, const int *qlens, int n_regs0, const mm_reg1_t *regs0, int *n_regs, mm_reg1_t **regs, const mm128_t *a);
void mm_seg_free(void *km, int n_segs, ri_seg_t *segs);
// void mm_set_mapq(void *km, int n_regs, mm_reg1_t *regs, int min_chain_sc, int match_sc, int rep_len, int is_sr);
// void mm_set_mapq(void *km, int n_regs, mm_reg1_t *regs, int min_chain_sc, int match_sc, int rep_len);
void mm_set_mapq(void *km, int n_regs, mm_reg1_t *regs, int min_chain_sc, int rep_len);

#ifdef __cplusplus
}
#endif
#endif //RCHAIN_H