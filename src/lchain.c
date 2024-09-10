#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "kalloc.h"
#include "rutils.h"
#include "rseed.h"
#include "rsketch.h"
#include "krmq.h"
#include "chain.h"

#ifdef PROFILERH
double ri_sorttime = 0.0;
#endif

/*
* @brief 			Compute the log2 of a float number.
*
* @param x 			Float number
*
* @return 			Log2 of x
*/
static inline float mg_log2(float x) // NB: this doesn't work when x<2
{
	union { float f; uint32_t i; } z = { x };
	float log_2 = ((z.i >> 23) & 255) - 128;
	z.i &= ~(255 << 23);
	z.i += 127 << 23;
	log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
	return log_2;
}

/*
* @brief 			Backtrack from the end of a chain to the beginning.
*
* @param max_drop 	Maximum drop in score allowed between the last anchor (k) and its furthest predecessor in the chain
* @param z 			z[] keeps the anchors with acceptable scores.
* 						z[].x is the score
*						z[].y is the index of the anchor.
* @param f 			f[i] is the best score for anchor i
* @param p 			p[i] is the index of the best predecessor for anchor i
* @param t 			t[i] shows if a current anchor has been 'touched' (i.e., visited) so far in another chain backtracking.
* @param k 			Index of the anchor to backtrack from
*
* @return 			Index of the anchor that shows the beginning of the chain
*/
static int64_t mg_chain_bk_end(int32_t max_drop,
							   const mm128_t *z,
							   const int32_t *f,
							   const int64_t *p,
							   int32_t *t,
							   int64_t k)
{
	int64_t i = z[k].y, end_i = -1, max_i = i;
	int32_t max_s = 0;
	
	if(i < 0 || t[i] != 0) return i;

	do {
		// s is the score difference between the current anchor and the initial anchor (k)
		int32_t s;
		t[i] = 2;
		end_i = i = p[i];
		s = i < 0? z[k].x : (int32_t)z[k].x - f[i];
		// max_s is the maximum score difference between an anchor we processed so far and the initial anchor (k)
		if(s > max_s) max_s = s, max_i = i;
		// If the score difference between the current anchor and the initial anchor (k) is larger than the max_drop, we stop backtracking
		else if(max_s - s > max_drop) break;
	} while (i >= 0 && t[i] == 0); //do this until we reach a 'touched' anchor

	// reset modified t[]
	for (i = z[k].y; i >= 0 && i != end_i; i = p[i]) t[i] = 0;
	// return the index of the anchor that shows the beginning of the chain
	return max_i;
}

/**
 * @brief 			Backtrack from all possible chains to find the best chains.
 * 
 * @param km 		Pointer to memory pool.
 * @param n 		Number of anchors.
 * @param f 		f[i] is the best score for anchor i
 * @param p 		p[i] is the index of the best predecessor for anchor i
 * @param v 		v[] includes the backtracked anchors in the chains
 * @param t 		t[i] determines if anchor i has been 'touched' (i.e., visited/included) so far in any chain.
 * @param min_cnt 	Minimum number of anchors in a chain
 * @param min_sc  	Minimum score for a chain to be considered for backtracking
 * @param max_drop 	Maximum drop of a chain.
 * @param n_u_ 		Number of chains
 * @param n_v_ 		Number of anchors included in any chain
 * 
 * @return 			Pointer to the chains u[].
 * 					u[]: score<<32 | #anchors
 */
uint64_t *mg_chain_backtrack(void *km,
							 int64_t n,
							 const int32_t *f,
							 const int64_t *p,
							 int32_t *v,
							 int32_t *t,
							 int32_t min_cnt,
							 int32_t min_sc,
							 int32_t max_drop,
							 int32_t *n_u_,
							 int32_t *n_v_)
{
	mm128_t *z;
	uint64_t *u;
	int64_t i, k, n_z, n_v;
	int32_t n_u;

	*n_u_ = *n_v_ = 0;

	// n_z = # of anchors with acceptable scores.
	for (i = 0, n_z = 0; i < n; ++i) // precompute n_z
		if(f[i] >= min_sc) ++n_z;

	if(n_z == 0) return 0;

	KMALLOC(km, z, n_z);
	// z[] keeps the anchors with acceptable scores.
	// z[].x is the score, z[].y is the index of the anchor.
	for (i = 0, k = 0; i < n; ++i) // populate z[]
		if(f[i] >= min_sc) z[k].x = f[i], z[k++].y = i;

	// Sort z[] by score.
	#ifdef PROFILERH
	double csort_t = ri_realtime();
	#endif
	radix_sort_128x(z, z + n_z); //this sorting existed in the earlier version of RawHash as well.
	#ifdef PROFILERH
	ri_sorttime += ri_realtime() - csort_t;
	#endif

	// Filling t with 0s
	memset(t, 0, n * 4);

	// TODO: Do we need to precompute this? What is its memory overhead if we do not precompute it?
	// I will disable it for now.
	/*
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k){ // precompute n_u
		if(t[z[k].y] == 0){
			// n_v0 is the number of anchors we included in any chain so far
			// We keep it to calculate the number of anchors we are about to include in the current chain
			int64_t n_v0 = n_v, end_i;
			int32_t sc;
			// Find the beginning of the chain for the current anchor (k). end_i is the index to that anchor in the beginning of the chain
			end_i = mg_chain_bk_end(max_drop, z, f, p, t, k);
			// Backtrack from the end (k) to the beginning (end_i) and make all anchors in between 'touched' (i.e., visited)
			// n_v is the number of anchors we included in any chain so far
			for (i = z[k].y; i != end_i; i = p[i]) ++n_v, t[i] = 1;
			sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
			// If the score of the chain is acceptable, and we have enough anchors in the chain, increment the number of chains (n_u)
			if(sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt) ++n_u;
			else n_v = n_v0;
		}
	}
	KMALLOC(km, u, n_u);
	*/
	KMALLOC(km, u, n_z);
	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k){ // populate u[]
		if(t[z[k].y] == 0){
			// n_v0 is the number of anchors we included in any chain so far
			// We keep it to calculate the number of anchors we are about to include in the current chain
			int64_t n_v0 = n_v, end_i;
			int32_t sc;
			// Find the beginning of the chain for the current anchor (k). end_i is the index to that anchor in the beginning of the chain
			end_i = mg_chain_bk_end(max_drop, z, f, p, t, k);
			// Backtrack from the end (k) to the beginning (end_i) and make all anchors in between 'touched' (i.e., visited)
			// n_v is the number of anchors we included in any chain so far
			for (i = z[k].y; i != end_i; i = p[i]) v[n_v++] = i, t[i] = 1;
			
			sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
			// If the score of the chain is acceptable, and we have enough anchors in the chain, increment the number of chains (n_u)
			if(sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
				//TODO: Limited information stored for a chain to enable its backtracking. Store k and end_i as well? (or just k?)
				//Maybe v[] is enough to backtrack a chain?
				u[n_u++] = (uint64_t)sc << 32 | (n_v - n_v0);
			else n_v = n_v0;
		}
	}

	//Shrink u and copy the content of u to shrink_u.
	//TODO: this is enabled unless we precompute the size of u, which is commented out for now. Check in the profile
	// KMALLOC(km, shrink_u, n_u);
	// memcpy(shrink_u, u, n_u * 8);
	// ri_kfree(km, u); u = shrink_u;

	ri_kfree(km, z);
	assert(n_v < INT32_MAX);
	*n_u_ = n_u, *n_v_ = n_v;
	return u;
}

/*
* @brief 			Compact the anchors in a chain.
*
* @param km 		Pointer to memory pool.
* @param n_u 		Number of chains
* @param u 			Input: Chains
* 						u[]: score<<32 | #anchors in a chain
*					Output: Sorted chains by their target positions
* @param n_v 		Number of anchors included in any chain
* @param v 			Input: v[] includes the backtracked anchors in the chains
* 					Output: Deallocated at return.
* @param a 			Input: Anchors
* 						a[].x: rev<<63 	 | tid<<32 	  | t_pos
* 						a[].y: flags<<40 | q_span<<32 | q_pos
*					Output: Deallocated at return.
*
* @return 			The list of anchors in the chains (chains are sorted by their target positions). The structure is the same as a[]
*/
static mm128_t *compact_a(void *km,
						  int32_t n_u,
						  uint64_t *u,
						  int64_t *n_v,
						  int32_t *v,
						  mm128_t *a,
						  mm128_t **_a)
{
	// b is the list of anchors in the chains (chains are sorted by their target positions)
	mm128_t *b, *w;
	uint64_t *u2;
	int64_t i, j, k;

	// write the result to b[]
	if(*_a) ri_kfree(km, *_a);
	KMALLOC(km, b, *n_v);
	KMALLOC(km, *_a, *n_v);
	// KMALLOC(km, c, *n_v);
	for (i = 0, k = 0; i < n_u; ++i){
		//k0 is the index to the first anchor in the current chain
		//ni is the number of anchors in the current chain
		int32_t k0 = k, ni = (int32_t)u[i];
		//Copy the anchors in the current chain to b[]
		//Note that the anchors are copied in reverse order than they are in v[] (b is anchors, v are indices to anchors)
		//So the anchor with lowest index is at the beginning of the chain.
		for (j = 0; j < ni; ++j){
			b[k++] = a[v[k0 + (ni - j - 1)]];
			(*_a)[k-1] = b[k-1];
		}
	}
	ri_kfree(km, v);

	// sort u[] and a[] by the target position, such that adjacent chains may be joined
	KMALLOC(km, w, n_u);
	for (i = k = 0; i < n_u; ++i){
		// w[i].x: a[].x of the first anchor in the current chain: a[].x: rev<<63 | tid<<32 | t_pos
		// w[i].y: index of the first anchor in the current chain << 32 | chain_id
		w[i].x = b[k].x, w[i].y = (uint64_t)k<<32|i;
		k += (int32_t)u[i];
	}

	#ifdef PROFILERH
	double csort_t = ri_realtime();
	#endif
	radix_sort_128x(w, w + n_u);
	#ifdef PROFILERH
	ri_sorttime += ri_realtime() - csort_t;
	#endif
	
	KMALLOC(km, u2, n_u);
	for (i = k = 0; i < n_u; ++i){
		//j is the chain id, n is the number of anchors in chain j
		int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];
		// u2 is the sorted u[] by the target position
		u2[i] = u[j];
		// store sorted b[] (based on sorted u[] by the target position) in a[]
		memcpy(&a[k], &b[w[i].y>>32], n * sizeof(mm128_t));
		k += n;
	}
	// u is now sorted by the target position
	memcpy(u, u2, n_u * 8);
	// write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot.
	// b now includes only the anchors in the chains (i.e., sorted b or sorted v in reverse order)
	memcpy(b, a, k * sizeof(mm128_t));
	ri_kfree(km, a); ri_kfree(km, w); ri_kfree(km, u2);
	*n_v = k;
	return b;
}

/*
* @brief 			Compute the chaining score between two anchors i and j.
*
* @param ai 			Anchor i
* @param aj 			Anchor j
* @param max_dist_t 	Maximum distance in the target (reference)
* @param max_dist_q 	Maximum distance in the query (read)
* @param bw 			The absolute difference in gaps between two anchors: abs(target_gap - query_gap) t_gap = anchort0 - anchort1
*						i.e., difference in diagonals
* @param chn_pen_gap 	Gap penalty for the chain
* @param chn_pen_skip 	Anchor skip penalty for the chain
*
* @return 			Chaining score
*/
static inline int32_t compute_score(const mm128_t *ai,
                                	const mm128_t *aj,
                                	int32_t max_dist_t,
                                	int32_t max_dist_q,
                                	int32_t bw,
                                	float chn_pen_gap,
                                	float chn_pen_skip)
                                	// int is_cdna,
                                	// int n_seg)
{
    uint32_t span_mask = (1U<<RI_HASH_SHIFT)-1;

	//dq = distance query, dr = distance reference, dd = distance diagonal
	int32_t dq = (int32_t)ai->y - (int32_t)aj->y, dr, dd, dg, q_span, sc;
	// int32_t sidi = (ai->y & RI_SEED_SEG_MASK) >> RI_SEED_SEG_SHIFT;
	// int32_t sidj = (aj->y & RI_SEED_SEG_MASK) >> RI_SEED_SEG_SHIFT;

	if(dq <= 0 || dq > max_dist_q) return INT32_MIN;

	//Calculate the distance between two anchors in the reference
	dr = (int32_t)(ai->x - aj->x);
	// if(sidi == sidj && (dr == 0 || dq > max_dist_q)) return INT32_MIN;
	if(dr == 0 || dr > max_dist_t) return INT32_MIN;

	//Calculate the distance between two anchors in the diagonal
	dd = dr > dq? dr - dq : dq - dr;

	// if(sidi == sidj && dd > bw) return INT32_MIN;
	if(dd > bw || dr > max_dist_q) return INT32_MIN;
	// if(n_seg > 1 && !is_cdna && sidi == sidj && dr > max_dist_q) return INT32_MIN;
    // if(sidi == sidj && dr > max_dist_q) return INT32_MIN;
	// if(dr > max_dist_q) return INT32_MIN;

	// dg is the gap between two adjacent anchors (in query or target, whichever is shorter)
	dg = dr < dq? dr : dq;

	// TODO: currently the span is only determined by "e" (number of events concatanated in a seed)
	q_span = (aj->y>>RI_ID_SHIFT)&span_mask;

	// Matching bases. Consider the matching bases (dg) if it is larger than the span (q_span).
	sc = q_span < dg? q_span : dg;
	
	// Integrating penalties to the score
	if(dd || dg > q_span){
		float lin_pen, log_pen;
		lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
		log_pen = dd >= 1? mg_log2(dd + 1) : 0.0f; // mg_log2() only works for dd>=2
		// if(is_cdna || sidi != sidj){
        // if(sidi != sidj)
		// {
		// 	if(sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
		// 	else if(dr > dq || sidi != sidj) sc -= (int)(lin_pen < log_pen? lin_pen : log_pen); // deletion or jump between paired ends
		// 	else sc -= (int)(lin_pen + .5f * log_pen);
		// } else
		sc -= (int)(lin_pen + .5f * log_pen);
		// sc -= (int)(lin_pen);
	}

	return sc;
}

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
 * @param a					Anchors (input)
 * 							a[].x: rev<<63 | tid<<32 | t_pos
 * 							a[].y: flags<<40 | q_span<<32 | q_pos
 * @param _a				Anchors (output)
 * 							a[].x: rev<<63 | tid<<32 | t_pos
 * 							a[].y: flags<<40 | q_span<<32 | q_pos
 * @param n_u_				Number of chains
 * @param u					Chains (output)
 * 							u[]: score<<32 | #anchors
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
					//   float chn_rev_bump, 
                    //   int n_seg, 
                      int64_t *n, 
                      mm128_t *a, 
                      mm128_t **_a,
					  int *n_u_,  
                      uint64_t **_u, 
                      void *km,
					  double* profile_sort_time)
{	// TODO: make sure this works when *n has more than 32 bits

	#ifdef PROFILERH
	ri_sorttime = 0.0;
	#endif

	//t[i] provides the anchor index j whose predecessor (k) is i. So i -> k -> j (i.e., there is k between i and j)
	//f[i] is the best score for anchor i
	//p[i] is the index of the best predecessor for anchor i
	// n_u is the number of chains
	// mmax_f is the maximum score for an anchor within a band
	int32_t *f, *t, *v, n_u, n_v, mmax_f = 0, max_drop = bw;
	int64_t *p, i, j, max_ii, st = 0;
	// int64_t n_iter = 0;
	uint64_t *u;

	if(_u) *_u = 0, *n_u_ = 0;
	if(*n == 0 || a == 0){
		ri_kfree(km, a);
		if(*_a) ri_kfree(km, *_a);
		*_a = NULL;
		return 0;
	}
	if(max_dist_t < bw) max_dist_t = bw;
	// if(max_dist_q < bw && !is_cdna) max_dist_q = bw;
    if(max_dist_q < bw) max_dist_q = bw;
	// if(is_cdna) max_drop = INT32_MAX;
	KMALLOC(km, p, *n);
	KMALLOC(km, f, *n);
	KMALLOC(km, v, *n);
	RI_KCALLOC(km, t, *n);

    uint32_t span_mask = (1U<<RI_HASH_SHIFT)-1;

	// fill the score and backtrack arrays
	// *n is the number of anchors
	for(i = 0, max_ii = -1; i < *n; ++i){

		// max_j is the index of the best predecessor (with max score) for anchor i
		int64_t max_j = -1, end_j;

		//Chaining score of the best predecessor (with max score) for anchor i
		int32_t max_f = (a[i].y>>RI_ID_SHIFT)&span_mask, n_skip = 0;

		// Find a good starting point index (st) within the same target index, and the target pos difference is less than max_dist_t
		while(st < i && (a[i].x>>RI_ID_SHIFT != a[st].x>>RI_ID_SHIFT || a[i].x > a[st].x + max_dist_t)) ++st;

		// We can only traverse max_iter anchors to find the best predecessor for i
		if(i - st > max_iter) st = i - max_iter;

		// Iterate over all possible predecessors of i goind back from i to st
		for(j = i - 1; j >= st; --j){
			int32_t sc;
			sc = compute_score(&a[i], &a[j], max_dist_t, max_dist_q, bw, chn_pen_gap, chn_pen_skip);
			// ++n_iter;
			if(sc == INT32_MIN) continue;
			//Integrate the potential predecessor's (j) score with the score between anchors i and j
			sc += f[j];

			// Update the max chaining score we have seen so far for the current anchor if it is higher than the current max score
			if(sc > max_f){
				max_f = sc, max_j = j;
				if(n_skip > 0) --n_skip;
			}else if(t[j] == (int32_t)i){
				if(++n_skip > max_skip) break;
			}
			if(p[j] >= 0) t[p[j]] = i;
		}
		end_j = j;

		if(max_ii < 0 || a[i].x - a[max_ii].x > (int64_t)max_dist_t){
			int32_t max = INT32_MIN;
			max_ii = -1;

			// Find the predecessor within the band aka max_dist_t with the highest score
			// (not necessarily integrated with the score between anchors i and j)
			// max_ii is the index to that predecessor
			for (j = i - 1; j >= st; --j)
				if(max < f[j]) max = f[j], max_ii = j;
		}

		// We compute the score between the current anchor and the predecessor with highest score
		// If the score is higher than the current max score, we update the max score and the max_j index
		// This is probably already done in the above while loop. Not sure why we need to do it again. TODO: check later
		if(max_ii >= 0 && max_ii < end_j){
			int32_t tmp;
			tmp = compute_score(&a[i], &a[max_ii], max_dist_t, max_dist_q, bw, chn_pen_gap, chn_pen_skip);
			if(tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}

		// Update the score and the predecessor of anchor i based on the found best predecessor max_j
		f[i] = max_f, p[i] = max_j;
		// if(a[i].x >> 63) f[i] *= chn_rev_bump; //this is mainly for reverse complememting using signals, otherwise it has no effect (1.0)

		// v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		// v[i] is almost the same as f[i] but it can assign v[i] = v[max_j] where peak score of max_j can be better than f[i]
		v[i] = (max_j >= 0 && v[max_j] > max_f)? v[max_j] : max_f;

		// Update the max_ii index if the current anchor has a higher score than the previous max_ii (within the band aka max_dist_t)
		if(max_ii < 0 || (a[i].x - a[max_ii].x <= (int64_t)max_dist_t && f[max_ii] < f[i])) max_ii = i;
		if(mmax_f < max_f) mmax_f = max_f;
	}

	// Find all chains with their scores and the number of anchors in them
	u = mg_chain_backtrack(km, *n, f, p, v, t, min_cnt, min_sc, max_drop, &n_u, &n_v);
	
	// NB: note that u[] may not be sorted by score here
	*n_u_ = n_u, *_u = u;
	ri_kfree(km, p); ri_kfree(km, f); ri_kfree(km, t);
	if(n_u == 0){
		ri_kfree(km, a); ri_kfree(km, v);
		*n = 0; if(*_a) ri_kfree(km, *_a); *_a = NULL;
		return 0;
	}

	// Sort the chains by their target positions and return the anchors in the sorted chains
	// a is deallocated at return
	*n=n_v;
	
	#ifdef PROFILERH
	mm128_t* ret_val = compact_a(km, n_u, u, n, v, a, _a);
	if(profile_sort_time)*profile_sort_time = ri_sorttime;
	return ret_val;
	#else
		return compact_a(km, n_u, u, n, v, a, _a);
	#endif
}

typedef struct lc_elem_s {
	int32_t y;
	int64_t i;
	double pri;
	KRMQ_HEAD(struct lc_elem_s) head;
} lc_elem_t;

#define lc_elem_cmp(a, b) ((a)->y < (b)->y? -1 : (a)->y > (b)->y? 1 : ((a)->i > (b)->i) - ((a)->i < (b)->i))
#define lc_elem_lt2(a, b) ((a)->pri < (b)->pri)
KRMQ_INIT(lc_elem, lc_elem_t, head, lc_elem_cmp, lc_elem_lt2)

RI_KALLOC_POOL_INIT(rmq, lc_elem_t)

/*
* @brief 			Compute the chaining score between two anchors i and j.
*
* @param ai 			Anchor i
* @param aj 			Anchor j
* @param chn_pen_gap 	Gap penalty for the chain
* @param chn_pen_skip 	Anchor skip penalty for the chain
* @param exact 			Whether the score is exact
* @param width 			Width of the band
*
* @return 			Chaining score
*/
static inline int32_t comput_sc_simple(const mm128_t *ai,
									   const mm128_t *aj,
									   float chn_pen_gap,
									   float chn_pen_skip,
									   int32_t *exact,
									   int32_t *width)
{
	uint32_t span_mask = (1U<<RI_HASH_SHIFT)-1;

	int32_t dq = (int32_t)ai->y - (int32_t)aj->y, dr, dd, dg, q_span, sc;
	dr = (int32_t)(ai->x - aj->x);
	*width = dd = dr > dq? dr - dq : dq - dr;
	dg = dr < dq? dr : dq;
	q_span = (aj->y>>RI_ID_SHIFT)&span_mask;
	sc = q_span < dg? q_span : dg;
	if (exact) *exact = (dd == 0 && dg <= q_span);
	if (dd || dq > q_span) {
		float lin_pen, log_pen;
		lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
		log_pen = dd >= 1? mg_log2(dd + 1) : 0.0f; // mg_log2() only works for dd>=2
		sc -= (int)(lin_pen + .5f * log_pen);
	}
	return sc;
}

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
 * 							u[].x: score<<32 | #anchors (sum of lower 32 bits of u[] is the returned length of a[])
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
					   void *km)
{
	uint32_t span_mask = (1U<<RI_HASH_SHIFT)-1;

	int32_t *f,*t, *v, n_u, n_v, mmax_f = 0, max_rmq_size = 0, max_drop = bw;
	int64_t *p, i, i0, st = 0, st_inner = 0;
	uint64_t *u;
	lc_elem_t *root = 0, *root_inner = 0;
	void *mem_mp = 0;
	kmp_rmq_t *mp;

	if (_u) *_u = 0, *n_u_ = 0;
	if (*n == 0 || a == 0) {
		ri_kfree(km, a);
		if(*_a) ri_kfree(km, *_a);
		*_a = NULL;
		return 0;
	}
	if (max_dist < bw) max_dist = bw;
	if (max_dist_inner <= 0 || max_dist_inner >= max_dist) max_dist_inner = 0;
	KMALLOC(km, p, *n);
	KMALLOC(km, f, *n);
	RI_KCALLOC(km, t, *n);
	KMALLOC(km, v, *n);
	mem_mp = ri_km_init2(km, 0x10000);
	mp = kmp_init_rmq(mem_mp);

	// fill the score and backtrack arrays
	for (i = i0 = 0; i < *n; ++i) {
		int64_t max_j = -1;
		int32_t q_span = (a[i].y>>RI_ID_SHIFT)&span_mask, max_f = q_span;
		lc_elem_t s, *q, *r, lo, hi;
		// add in-range anchors
		if (i0 < i && a[i0].x != a[i].x) {
			int64_t j;
			for (j = i0; j < i; ++j) {
				q = kmp_alloc_rmq(mp);
				q->y = (int32_t)a[j].y, q->i = j, q->pri = -(f[j] + 0.5 * chn_pen_gap * ((int32_t)a[j].x + (int32_t)a[j].y));
				krmq_insert(lc_elem, &root, q, 0);
				if (max_dist_inner > 0) {
					r = kmp_alloc_rmq(mp);
					*r = *q;
					krmq_insert(lc_elem, &root_inner, r, 0);
				}
			}
			i0 = i;
		}
		// get rid of active chains out of range
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist || krmq_size(head, root) > cap_rmq_size)) {
			s.y = (int32_t)a[st].y, s.i = st;
			if ((q = krmq_find(lc_elem, root, &s, 0)) != 0) {
				q = krmq_erase(lc_elem, &root, q, 0);
				kmp_free_rmq(mp, q);
			}
			++st;
		}
		if (max_dist_inner > 0)  { // similar to the block above, but applied to the inner tree
			while (st_inner < i && (a[i].x>>32 != a[st_inner].x>>32 || 
				   a[i].x > a[st_inner].x + max_dist_inner || 
				   krmq_size(head, root_inner) > cap_rmq_size)) {
				s.y = (int32_t)a[st_inner].y, s.i = st_inner;
				if ((q = krmq_find(lc_elem, root_inner, &s, 0)) != 0) {
					q = krmq_erase(lc_elem, &root_inner, q, 0);
					kmp_free_rmq(mp, q);
				}
				++st_inner;
			}
		}
		// RMQ
		lo.i = INT32_MAX, lo.y = (int32_t)a[i].y - max_dist;
		hi.i = 0, hi.y = (int32_t)a[i].y;
		if ((q = krmq_rmq(lc_elem, root, &lo, &hi)) != 0) {
			int32_t sc, exact, width, n_skip = 0;
			int64_t j = q->i;
			assert(q->y >= lo.y && q->y <= hi.y);
			sc = f[j] + comput_sc_simple(&a[i], &a[j], chn_pen_gap, chn_pen_skip, &exact, &width);
			if (width <= bw && sc > max_f) max_f = sc, max_j = j;
			if (!exact && root_inner && (int32_t)a[i].y > 0) {
				lc_elem_t *lo, *hi;
				s.y = (int32_t)a[i].y - 1, s.i = *n;
				krmq_interval(lc_elem, root_inner, &s, &lo, &hi);
				if (lo) {
					const lc_elem_t *q;
					int32_t width, n_rmq_iter = 0;
					krmq_itr_t(lc_elem) itr;
					krmq_itr_find(lc_elem, root_inner, lo, &itr);
					while ((q = krmq_at(&itr)) != 0) {
						if (q->y < (int32_t)a[i].y - max_dist_inner) break;
						++n_rmq_iter;
						j = q->i;
						sc = f[j] + comput_sc_simple(&a[i], &a[j], chn_pen_gap, chn_pen_skip, 0, &width);
						if (width <= bw) {
							if (sc > max_f) {
								max_f = sc, max_j = j;
								if (n_skip > 0) --n_skip;
							} else if (t[j] == (int32_t)i) {
								if (++n_skip > max_skip)
									break;
							}
							if (p[j] >= 0) t[p[j]] = i;
						}
						if (!krmq_itr_prev(lc_elem, &itr)) break;
					}
				}
			}
		}
		// set max
		// assert(max_j < 0 || (a[max_j].x <= a[i].x && (int32_t)a[max_j].y <= (int32_t)a[i].y));

		// Update the score and the predecessor of anchor i based on the found best predecessor max_j
		f[i] = max_f, p[i] = max_j;

		// v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f;

		// Update mmax_f if the current anchor has a higher score than the previous mmax_f (within the band aka max_dist_t)
		if (mmax_f < max_f) mmax_f = max_f;
		if (max_rmq_size < krmq_size(head, root)) max_rmq_size = krmq_size(head, root);
	}
	ri_km_destroy(mem_mp);

	// Find all chains with their scores and the number of anchors in them
	u = mg_chain_backtrack(km, *n, f, p, v, t, min_cnt, min_sc, max_drop, &n_u, &n_v);
	// NB: note that u[] may not be sorted by score here
	*n_u_ = n_u, *_u = u;
	ri_kfree(km, p); ri_kfree(km, f); ri_kfree(km, t);
	if (n_u == 0) {
		ri_kfree(km, a); ri_kfree(km, v);
		*n = 0; if(*_a) ri_kfree(km, *_a); *_a = NULL;
		return 0;
	}

	// Sort the chains by their target positions and return the anchors in the sorted chains
	// a is deallocated at return
	*n = n_v;
	return compact_a(km, n_u, u, n, v, a, _a);
}
