#include "rseed.h"
#include "kalloc.h"
#include "rsketch.h"
#include "rutils.h"

#define MAX_MAX_HIGH_OCC 128

void ri_seed_select(int32_t n, ri_seed_t *a, int len, int max_occ, int max_max_occ, int dist)
{ // for high-occ minimizers, choose up to max_high_occ in each high-occ streak
	extern void ks_heapdown_uint64_t(size_t i, size_t n, uint64_t*);
	extern void ks_heapmake_uint64_t(size_t n, uint64_t*);
	int32_t i, last0, m;
	uint64_t b[MAX_MAX_HIGH_OCC]; // this is to avoid a heap allocation

	if (n == 0 || n == 1) return;
	for (i = m = 0; i < n; ++i)
		if (a[i].n > max_occ) ++m;
	if (m == 0) return; // no high-frequency k-mers; do nothing
	for (i = 0, last0 = -1; i <= n; ++i) {
		if (i == n || a[i].n <= max_occ) {
			if (i - last0 > 1) {
				int32_t ps = last0 < 0? 0 : (uint32_t)a[last0].q_pos>>1;
				int32_t pe = i == n? len : (uint32_t)a[i].q_pos>>1;
				int32_t j, k, st = last0 + 1, en = i;
				int32_t max_high_occ = (int32_t)((double)(pe - ps) / dist + .499);
				if (max_high_occ > 0) {
					if (max_high_occ > MAX_MAX_HIGH_OCC)
						max_high_occ = MAX_MAX_HIGH_OCC;
					for (j = st, k = 0; j < en && k < max_high_occ; ++j, ++k)
						b[k] = (uint64_t)a[j].n<<32 | j;
					ks_heapmake_uint64_t(k, b); // initialize the binomial heap
					for (; j < en; ++j) { // if there are more, choose top max_high_occ
						if (a[j].n < (int32_t)(b[0]>>32)) { // then update the heap
							b[0] = (uint64_t)a[j].n<<32 | j;
							ks_heapdown_uint64_t(0, k, b);
						}
					}
					for (j = 0; j < k; ++j) a[(uint32_t)b[j]].flt = 1;
				}
				for (j = st; j < en; ++j) a[j].flt ^= 1;
				for (j = st; j < en; ++j)
					if (a[j].n > max_max_occ)
						a[j].flt = 1;
			}
			last0 = i;
		}
	}
}

/**
 * @brief Collect all seeds from the given index and query vector.
 *
 * @param km Pointer to memory pool.
 * @param ri Pointer to the index.
 * @param riv Pointer to the query vector (vector of sketches).
 * @param n_seed_hits_ Pointer to the number of seeds collected.
 *
 * @return Pointer to the collected seeds.
 */
ri_seed_t *ri_seed_collect_all(void *km, const ri_idx_t *ri, const mm128_v *riv, int32_t *n_seed_hits_)
{
	ri_seed_t *seed_hits;
	size_t i;
	int32_t n_seed_hits;

    // uint32_t pos_mask = (1U<<31)-1;
    uint32_t span_mask = (1U<<RI_HASH_SHIFT)-1;
	seed_hits = (ri_seed_t*)ri_kmalloc(km, riv->n * sizeof(ri_seed_t));
	for (i = n_seed_hits = 0; i < riv->n; ++i) {
		const uint64_t *cr;
		ri_seed_t *q;
		mm128_t *p = &riv->a[i];
		uint32_t q_pos = (uint32_t)p->y, q_span = p->x & span_mask;
		int t;
		cr = ri_idx_get(ri, p->x>>RI_HASH_SHIFT, &t);
		if (t == 0) continue;
		q = &seed_hits[n_seed_hits++];
		q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> RI_ID_SHIFT;
		q->is_tandem = q->flt = 0;
		if (i > 0 && p->x>>RI_HASH_SHIFT == riv->a[i - 1].x>>RI_HASH_SHIFT) q->is_tandem = 1;
		if (i < riv->n - 1 && p->x>>RI_HASH_SHIFT == riv->a[i + 1].x>>RI_HASH_SHIFT) q->is_tandem = 1;
	}
	*n_seed_hits_ = n_seed_hits;
	return seed_hits;
}

/**
 * @brief Collects all seeds that match the given query vector and index.
 *
 * @param km Pointer to memory pool.
 * @param _n_sm Number of overall seed matches (not the positions).
//  * @param qlen Length of the query vector.
* @param max_occ Maximum number of occurrences.
//  * @param max_max_occ Maximum number of maximum occurrences.
//  * @param dist Distance between seeds.
 * @param ri Pointer to the index.
 * @param riv Pointer to the query vector.
 * @param n_sm_pos Overall number of positions retrieved from seed matches.
 * @param rep_len Pointer to the length of the repetition.
//  * @param n_seed_mini Pointer to the number of seed positions.
 * @param seed_pos Pointer to the seed positions. Needs to be freed.
 *
 * @return Pointer to the collected seeds. Needs to be freed.
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
                              int *rep_len)
                            //   int *n_seed_mini,
                            //   uint64_t **seed_mini)
{

    int rep_st = 0, rep_en = 0, n_seed_m, n_seed_m0;
	size_t i;
	ri_seed_t *seed_hits;
	// uint32_t mask = (1ULL<<31)-1;
	// *n_seed_mini = 0;
	// *seed_mini = (uint64_t*)ri_kmalloc(km, riv->n * sizeof(uint64_t));
	seed_hits = ri_seed_collect_all(km, ri, riv, &n_seed_m0);

	uint32_t n_flt = 0, n_seed_hits = 0;

    // if (dist > 0 && max_max_occ > max_occ) {
	// 	ri_seed_select(n_seed_m0, seed_hits, qlen, max_occ, max_max_occ, dist);
	// } else{
		for (i = 0; i < n_seed_m0; ++i){
			n_seed_hits += seed_hits[i].n;
        	if (seed_hits[i].n > max_occ) {seed_hits[i].flt = 1; n_flt += seed_hits[i].n;}
		}
	// }

	fprintf(stderr, "n_seed_hits: %u n_flt: %u ratio: %f\n", n_seed_hits, n_flt, (float)n_flt/n_seed_hits);
    for (i = 0, n_seed_m = 0, *rep_len = 0, *n_seed_pos = 0; i < n_seed_m0; ++i) {
	// for (i = 0, n_seed_m = 0, *n_seed_pos = 0; i < n_seed_m0; ++i) {
		ri_seed_t *q = &seed_hits[i];
		//fprintf(stderr, "X\t%d\t%d\t%d\n", q->q_pos>>RI_POS_SHIFT, q->n, q->flt);
		if (q->flt) {
			// int en = (q->q_pos >> RI_POS_SHIFT) + 1, st = en - q->q_span;
			int st = (q->q_pos >> RI_POS_SHIFT) + 1, en = st + q->q_span + 1;
			if (st > rep_en) {
				*rep_len += rep_en - rep_st;
				rep_st = st, rep_en = en;
			} else rep_en = en;
		} else {
			*n_seed_pos += q->n;
			// (*seed_mini)[(*n_seed_mini)++] = (uint64_t)q->q_span<<32 | q->q_pos>>RI_POS_SHIFT;
			seed_hits[n_seed_m++] = *q;
		}
	}
	*rep_len += rep_en - rep_st;
	*_n_seed_m = n_seed_m;
	return seed_hits;
}

void ri_seed_mz_flt(void *km, mm128_v *riv, int32_t q_occ_max, float q_occ_frac)
{
	mm128_t *a;
	size_t i, j, st;
	if (riv->n <= q_occ_max || q_occ_frac <= 0.0f || q_occ_max <= 0) return;
	KMALLOC(km, a, riv->n);
	for (i = 0; i < riv->n; ++i)
		a[i].x = riv->a[i].x, a[i].y = i;
	radix_sort_128x(a, a + riv->n);
	for (st = 0, i = 1; i <= riv->n; ++i) {
		if (i == riv->n || a[i].x != a[st].x) {
			int32_t cnt = i - st;
			if (cnt > riv->n * q_occ_frac){
				for (j = st; j < i; ++j)
					riv->a[a[j].y].x = 0;
			}
			st = i;
		}
	}
	ri_kfree(km, a);
	for (i = j = 0; i < riv->n; ++i)
		if (riv->a[i].x != 0)
			riv->a[j++] = riv->a[i];
	riv->n = j;
}