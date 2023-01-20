#include "rawindex.h"

#include <stdlib.h>
#include <assert.h>
#if defined(WIN32) || defined(_WIN32)
#include <io.h> // for open(2)
#else
#include <unistd.h>
#endif
#include <fcntl.h>
#include <stdio.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include <string.h>
#include <stdint.h>

#include <x86intrin.h> //TODO: insert ifdef here to include simd headers

#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#define mm_seq4_set(s, i, c) ((s)[(i)>>3] |= (uint32_t)(c) << (((i)&7)<<2))

double ri_realtime0;
int ri_verbose = 1;

//Some of the functions here are adopted from the Sigmap implementation (https://github.com/haowenz/sigmap/tree/c9a40483264c9514587a36555b5af48d3f054f6f). We have optimized the Sigmap implementation to work with the hash tables efficiently.

typedef struct { size_t n, m; char **a; } ri_char_v;

double ri_realtime(void)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

double ri_cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long ri_peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

unsigned char seq_nt4_table[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

typedef struct {
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	const float* pore_vals;
	mm_bseq_file_t* fp;
	ri_idx_t* ri;
} pipeline_t;

typedef struct {
    int n_seq;
	mm_bseq1_t* seq;
	mm128_v a;
} step_t;

static inline uint64_t hash64(uint64_t key, uint64_t mask){
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

static inline uint64_t MurmurHash3(uint64_t key, uint64_t mask) {
  key = (key^(key >> 33));
  key = (key*0xff51afd7ed558ccd);
  key = (key^(key >> 33));
  key = (key*0xc4ceb9fe1a85ec53);
  key = (key^(key >> 33));

  return key&mask;
}

//Quantize float in range [min,max] to 10-bit unsigned integer
static inline uint32_t quantize_float_uint32(const float x, const float min, const float max) {
	uint32_t a = (uint32_t)(((x-min)/(max-min))*127.0f);
	if(a>127) return 127;
	if(a<0) return 0;

	return a;
}

//https://stackoverflow.com/a/21673221
static inline __m256i movemask_inverse(const uint32_t hash_value) {
    __m256i vmask = _mm256_set1_epi32(hash_value);
    const __m256i shuffle = _mm256_setr_epi64x(0x0000000000000000,
      0x0101010101010101, 0x0202020202020202, 0x0303030303030303);
    vmask = _mm256_shuffle_epi8(vmask, shuffle);
    const __m256i bit_mask = _mm256_set1_epi64x(0x7fbfdfeff7fbfdfe);
    vmask = _mm256_or_si256(vmask, bit_mask);
    return _mm256_cmpeq_epi8(vmask, _mm256_set1_epi64x(-1));
}

static inline void calc_blend_simd(__m256i* blndcnt_lsb, __m256i* blndcnt_msb, __m256i* ma, __m256i* mb, uint64_t val, const uint64_t mask, const int bits) {

    (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*ma), (*mb), movemask_inverse(val&mask)));
    uint64_t blendVal = (uint64_t)_mm256_movemask_epi8((*blndcnt_lsb))&mask;
    if(bits > 32){
        (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb), _mm256_blendv_epi8((*ma), (*mb),movemask_inverse((val>>32)&mask)));
        blendVal |= ((uint64_t)_mm256_movemask_epi8((*blndcnt_msb)))<<32;
    }
}

static inline uint64_t calc_blend_rm_simd(__m256i* blndcnt_lsb, __m256i* blndcnt_msb, __m256i* ma, __m256i* mb, uint64_t val, uint64_t remval, const uint64_t mask, const int bits) {

    (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*ma), (*mb), movemask_inverse(val&mask)));
    uint64_t blendVal = (uint64_t)_mm256_movemask_epi8((*blndcnt_lsb))&mask;
    //removal of the oldest item
    (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*mb), (*ma), movemask_inverse(remval&mask)));

    if(bits > 32){
        (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb), _mm256_blendv_epi8((*ma), (*mb),movemask_inverse((val>>32)&mask)));
        blendVal |= ((uint64_t)_mm256_movemask_epi8((*blndcnt_msb)))<<32;

        (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb),_mm256_blendv_epi8((*mb), (*ma), movemask_inverse((remval>>32)&mask)));
    }
    
    return blendVal;
}

void ri_idx_stat(const ri_idx_t *ri)
{
	fprintf(stderr, "[M::%s] pore kmer size: %d; concatanated events: %d; Quantization model (most sig. q/least sig. l): %d/%d; w: %d; n: %d; #seq: %d\n", __func__, ri->k, ri->e, ri->q, ri->lq, ri->w, ri->n, ri->n_seq);
}

ri_idx_t* ri_idx_init(int b, int w, int e, int n, int q, int lq, int k){
	ri_idx_t* ri;
	ri = (ri_idx_t*)calloc(1, sizeof(ri_idx_t));
  	ri->b = b, ri->w = w; ri->e = e; ri->n = n; ri->q = q; ri->lq = lq, ri->k = k;
  	ri->B = (ri_idx_bucket_t*)calloc(1<<ri->b, sizeof(ri_idx_bucket_t));
  	ri->km = ri_km_init();

	return ri;
}

void ri_idx_destroy(ri_idx_t *ri){
	uint32_t i;
	if (ri == 0) return;
	if (ri->h) kh_destroy(str, (khash_t(str)*)ri->h);
	if (ri->B) {
		for (i = 0; i < 1U<<ri->b; ++i) {
			free(ri->B[i].p);
			free(ri->B[i].a.a);
			kh_destroy(idx, (idxhash_t*)ri->B[i].h);
		}
	}
	if(ri->km) ri_km_destroy(ri->km);
	else if(ri->n_seq && ri->seq){
		for (i = 0; i < ri->n_seq; ++i)
			free(ri->seq[i].name);
		free(ri->seq);
	}
	free(ri->B); free(ri);
}

void ri_idx_add(ri_idx_t *ri, int n, const mm128_t *a){
	int i, mask = (1<<ri->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &ri->B[a[i].x>>RI_HASH_SHIFT&mask].a;
		kv_push(mm128_t, 0, *p, a[i]);
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq); // read a mini-batch
		if (s->seq) {
			uint32_t old_m, m;
			assert((uint64_t)p->ri->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->ri->seq
			old_m = p->ri->n_seq, m = p->ri->n_seq + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m)
				p->ri->seq = (ri_idx_seq_t*)ri_krealloc(p->ri->km, p->ri->seq, m*sizeof(ri_idx_seq_t));
			for (i = 0; i < s->n_seq; ++i) {
				ri_idx_seq_t *seq = &p->ri->seq[p->ri->n_seq];
				seq->name = (char*)ri_kmalloc(p->ri->km, strlen(s->seq[i].name) + 1);
				strcpy(seq->name, s->seq[i].name);
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				p->sum_len += seq->len;
				s->seq[i].rid = p->ri->n_seq++;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			mm_bseq1_t* t = &s->seq[i];
			if (t->l_seq > 0){
				uint32_t s_len;
				float* s_values = (float*)calloc(t->l_seq, sizeof(float));

				ri_seq_to_sig(t->seq, t->l_seq, p->pore_vals, p->ri->k, 0, &s_len, s_values);
				ri_sketch(0, s_values, t->rid, 0, s_len, p->ri->w, p->ri->e, p->ri->n, p->ri->q, p->ri->lq, p->ri->k, &s->a);

				ri_seq_to_sig(t->seq, t->l_seq, p->pore_vals, p->ri->k, 1, &s_len, s_values);
				ri_sketch(0, s_values, t->rid, 1, s_len, p->ri->w, p->ri->e, p->ri->n, p->ri->q, p->ri->lq, p->ri->k, &s->a);

				free(s_values);
			}
			free(t->seq); free(t->name); 
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
        step_t *s = (step_t*)in;
		ri_idx_add(p->ri, s->a.n, s->a.a);
		ri_kfree(0, s->a.a); free(s);
	}
    return 0;
}

#include "ksort.h"

#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8) 

#define sort_key_64(x) (x)
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8)

KSORT_INIT_GENERIC(uint32_t)
KSORT_INIT_GENERIC(uint64_t)

static void worker_post(void *g, long i, int tid){
	int n, n_keys;
	size_t j, start_a, start_p;
	idxhash_t *h = 0;
	ri_idx_t *ri = (ri_idx_t*)g;
	ri_idx_bucket_t *b = &ri->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		//canf: make sure to shift b->a.a[j].x accordingly if we store more data there
		if (j == b->a.n || b->a.a[j].x>>RI_HASH_SHIFT != b->a.a[j-1].x>>RI_HASH_SHIFT) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>RI_HASH_SHIFT != b->a.a[j-1].x>>RI_HASH_SHIFT) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>RI_HASH_SHIFT>>ri->b<<1, &absent);
			assert(absent && j == start_a + n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == (int32_t)start_p);

	// deallocate and clear b->a
	ri_kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}
 
static void ri_idx_post(ri_idx_t *ri, int n_threads){
	kt_for(n_threads, worker_post, ri, 1<<ri->b);
}

void ri_idx_sort(ri_idx_t *ri, int n_threads){
	ri_idx_post(ri, n_threads);
}

const uint64_t *ri_idx_get(const ri_idx_t *ri, uint64_t hashval, int *n){

	int mask = (1<<ri->b) - 1;
	khint_t k;
	ri_idx_bucket_t *b = &ri->B[hashval&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, hashval>>ri->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

void ri_idx_dump(FILE* idx_file, const ri_idx_t* ri){

	uint32_t pars[8], i;
	// uint64_t sum_len = 0;

	pars[0] = ri->w, pars[1] = ri->e, pars[2] = ri->n, pars[3] = ri->q, pars[4] = ri->lq, pars[5] = ri->k, pars[6] = ri->n_seq, pars[7] = ri->flag;
	
	fwrite(RI_IDX_MAGIC, 1, RI_IDX_MAGIC_BYTE, idx_file);
	fwrite(pars, sizeof(uint32_t), 8, idx_file);

	for (i = 0; i < ri->n_seq; ++i) {
		if (ri->seq[i].name) {
			uint8_t l = strlen(ri->seq[i].name);
			fwrite(&l, 1, 1, idx_file);
			fwrite(ri->seq[i].name, 1, l, idx_file);
		} else {
			uint8_t l = 0;
			fwrite(&l, 1, 1, idx_file);
		}
		fwrite(&ri->seq[i].len, 4, 1, idx_file);
		// sum_len += ri->seq[i].len;
	}
	for (i = 0; i < 1U<<ri->b; ++i) {
		ri_idx_bucket_t *b = &ri->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, idx_file);
		fwrite(b->p, 8, b->n, idx_file);
		fwrite(&size, 4, 1, idx_file);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, idx_file);
		}
	}

	fflush(idx_file);
}

ri_idx_t* ri_idx_load(FILE* idx_file){

	ri_idx_t* ri;
	uint64_t sum_len = 0;

	char magic[RI_IDX_MAGIC_BYTE];
  	uint32_t i;

	if (fread(magic, 1, RI_IDX_MAGIC_BYTE, idx_file) != RI_IDX_MAGIC_BYTE) return 0;
	if (strncmp(magic, RI_IDX_MAGIC, RI_IDX_MAGIC_BYTE) != 0) return 0;
	int pars[8];
	fread(&pars[0], sizeof(int), 8, idx_file);

	ri = ri_idx_init(14, pars[0], pars[1], pars[2], pars[3], pars[4], pars[5]);
	ri->n_seq = pars[6];
	ri->flag = pars[7];
	ri->seq = (ri_idx_seq_t*)ri_kcalloc(ri->km, ri->n_seq, sizeof(ri_idx_seq_t));

	for (i = 0; i < ri->n_seq; ++i) {
		uint8_t l;
		ri_idx_seq_t *s = &ri->seq[i];
		fread(&l, 1, 1, idx_file);
		if (l) {
			s->name = (char*)ri_kmalloc(ri->km, l + 1);
			fread(s->name, 1, l, idx_file);
			s->name[l] = 0;
		}
		fread(&s->len, 4, 1, idx_file);
		s->offset = sum_len;
		sum_len += s->len;
	}
	for (i = 0; i < 1U<<ri->b; ++i) {
		ri_idx_bucket_t *b = &ri->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		fread(&b->n, 4, 1, idx_file);
		b->p = (uint64_t*)malloc(b->n * 8);
		fread(b->p, 8, b->n, idx_file);
		fread(&size, 4, 1, idx_file);
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
			fread(x, 8, 2, idx_file);
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}

	return ri;
}

ri_idx_reader_t* ri_idx_reader_open(const char *fn, const ri_idxopt_t *ipt, const char *fn_out)
{
	int64_t is_idx;
	ri_idx_reader_t *r;
	is_idx = ri_idx_is_idx(fn);
	if (is_idx < 0) return 0; // failed to open the index
	r = (ri_idx_reader_t*)calloc(1, sizeof(ri_idx_reader_t));
	r->is_idx = is_idx;
	if (ipt) r->opt = *ipt;
	else ri_idxopt_init(&r->opt);
	if (r->is_idx) {
		r->fp.idx = fopen(fn, "rb");
		r->idx_size = is_idx;
	} else r->fp.seq = mm_bseq_open(fn);
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
	return r;
}

void ri_idx_reader_close(ri_idx_reader_t* r){

	if (r->is_idx) fclose(r->fp.idx);
	else mm_bseq_close(r->fp.seq);
	if (r->fp_out) fclose(r->fp_out);
	free(r);
}

ri_idx_t* ri_idx_gen(mm_bseq_file_t* fp, float* pore_vals, int b, int w, int e, int n, int q, int lq, int k, int flag, int mini_batch_size, int n_threads, uint64_t batch_size)
{
	pipeline_t pl;
	if (fp == 0 || mm_bseq_eof(fp) || !pore_vals) return 0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = (uint64_t)mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.pore_vals = pore_vals;
	pl.ri = ri_idx_init(b, w, e, n, q, lq, k);

	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	// if (ri_verbose >= 3)
	// 	fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, ri_realtime() - ri_realtime0, cputime() / (ri_realtime() - ri_realtime0));

	ri_idx_post(pl.ri, n_threads);
	// if (ri_verbose >= 3)
	// 	fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, ri_realtime() - ri_realtime0, cputime() / (ri_realtime() - ri_realtime0));

	return pl.ri;
}

ri_idx_t* ri_idx_reader_read(ri_idx_reader_t* r, float* pore_vals, int n_threads){

	ri_idx_t *ri;
	if (r->is_idx) {
		ri = ri_idx_load(r->fp.idx);
	} else{
		ri = ri_idx_gen(r->fp.seq, pore_vals, r->opt.b, r->opt.w, r->opt.e, r->opt.n, r->opt.q, r->opt.lq, r->opt.k, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size);
	}
	if (ri) {
		if (r->fp_out) ri_idx_dump(r->fp_out, ri);
		ri->index = r->n_parts++;
	}
	return ri;
}

int ri_idx_reader_eof(const ri_idx_reader_t* r){
	return r->is_idx? (feof(r->fp.idx) || ftell(r->fp.idx) == r->idx_size) : mm_bseq_eof(r->fp.seq);
}

int64_t ri_idx_is_idx(const char* fn){

	int fd, is_idx = 0;
	int64_t ret, off_end;
	char magic[RI_IDX_MAGIC_BYTE];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
#ifdef WIN32
	if ((off_end = _lseeki64(fd, 0, SEEK_END)) >= RI_IDX_MAGIC_BYTE) {
		_lseeki64(fd, 0, SEEK_SET);
#else
	if ((off_end = lseek(fd, 0, SEEK_END)) >= RI_IDX_MAGIC_BYTE) {
		lseek(fd, 0, SEEK_SET);
#endif // WIN32
		ret = read(fd, magic, RI_IDX_MAGIC_BYTE);
		if (ret == RI_IDX_MAGIC_BYTE && strncmp(magic, RI_IDX_MAGIC, RI_IDX_MAGIC_BYTE) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

void ri_idxopt_init(ri_idxopt_t *opt)
{
	memset(opt, 0, sizeof(ri_idxopt_t));
	opt->e = 6; opt->w = 0; opt->q = 9; opt->lq = 3; opt->n = 0; opt->k = 6;
	opt->b = 14;
	opt->mini_batch_size = 50000000;
	opt->batch_size = 4000000000ULL;
}

#include <math.h>

//Find the outlier n-dimensional vector of float numbers in range [0,1] in a set of n-dimensional vectors
static inline float find_outlier(const float** x, const uint32_t n, const uint32_t m) {
	uint32_t outlier = 0;
	float max_dist = 0.0f;
	for(uint32_t i=0; i<m; i++) {
		float dist = 0.0f;
		for(uint32_t j=0; j<n; j++) {
			dist += (x[i][j]-x[outlier][j])*(x[i][j]-x[outlier][j]);
		}
		if(dist>max_dist) {
			max_dist = dist;
			outlier = i;
		}
	}
	// fprintf(stdout, "dist: %.9g %d\n", 10000*(max_dist/n*n), outlier);
	return max_dist;
}

void ri_seq_to_sig(const char *str, int len, const float* pore_vals, const int k, const int strand, uint32_t* s_len, float* s_values){

	int i, j, l, pos, n = 0;
	// uint64_t shift1 = 2 * (k - 1);
	uint64_t mask = (1ULL<<2*k) - 1, kmer = 0;
	double mean = 0, std_dev = 0, sum = 0, sum2 = 0, curval = 0;

	for (i = l = j = n = 0; i < len; ++i) {
		if(strand) pos = len - i -1;
		else pos = i;
		int c = seq_nt4_table[(uint8_t)str[pos]];
		if (c < 4) { // not an ambiguous base
			if(!strand) kmer = (kmer << 2 | c) & mask;    // forward k-mer
			// else kmer = (kmer >> 2) | (3ULL^c) << shift1; // reverse k-mer
			//TODO: this is currently based on the ordering in the original ordering in the sigmap implementation. Change later to above
			else kmer = ((kmer << 2) | (3ULL^c)) & mask; // reverse k-mer
		}else
      		kmer = (kmer << 2) & mask; //TODO: This is not the best approach. We are basically inserting 00 (A?) to kmer whenever c >= 4. Mask it instead

		if(i+1 < k) continue;

		curval = pore_vals[kmer];
		s_values[j++] = curval;
		sum += curval;
		sum2 += curval*curval;
	}

	mean = sum/j;
	std_dev = sqrt(sum2/j - (mean)*(mean));

	for(i = 0; i < j; ++i)
		s_values[i] = (s_values[i]-mean)/std_dev;

	*s_len = j;
}

void ri_sketch_blend_min(void *km, const float* s_values, uint32_t id, int strand, int len, int w, int e, int n, int q, int lq, int k, mm128_v *p){
	
}

void ri_sketch_blend(void *km, const float* s_values, uint32_t id, int strand, int len, int e, int n, int q, int lq, int k, mm128_v *p){

	uint64_t blendVal = 0;
    int blend_pos = 0;
	bool buffull = false;
	const uint64_t blendmask =  (1ULL<<28)-1;

	mm128_t blndBuf[n];
	memset(blndBuf, 0, n*sizeof(mm128_t));
    
    //SIMD-related variables
    __m256i ma = _mm256_set1_epi8(1);
    __m256i mb = _mm256_set1_epi8(-1);
    __m256i blndcnt_lsb = _mm256_set1_epi8(0);
    __m256i blndcnt_msb = _mm256_set1_epi8(0);

	int step = 1;//TODO: make this an argument
	uint32_t span = (k+e-1)*step; //for now single event is considered to span 6 bases.

	const uint32_t quant_bit = lq+2; 
	const uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1; //for least sig. 5 bits
	
	bool sigBufFull = false;
	uint32_t i, sigBufPos = 0, l_sigpos = 0; //last signal position
	uint32_t signal = 0, tmpQuantSignal = 0;
	uint64_t quantVal = 0, hashVal = 0;

	kv_resize(mm128_t, km, *p, p->n + len/step);

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

    for (i = 0; i < len; i += step) {
        if(i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) continue;

		l_sigpos = i;
		signal = *((uint32_t*)&s_values[i]);
		tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);

		mm128_t info = { UINT64_MAX, UINT64_MAX };

		quantVal = (quantVal<<quant_bit|tmpQuantSignal)&mask_events;
		hashVal = hash64(quantVal, mask);

		sigBuf[sigBufPos].y = id_shift | (uint32_t)i<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = true; sigBufPos = 0;}
		// sigBuf[sigBufPos].x = hashVal<<RI_HASH_SHIFT | span;

		if(!sigBufFull) continue;

		blndBuf[blend_pos].x = hash64(quantVal, mask);
		blndBuf[blend_pos].y = sigBuf[sigBufPos].y;

		if(++blend_pos == n) {buffull = true; blend_pos = 0;}

		if(buffull){
			blendVal = calc_blend_rm_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, hashVal, blndBuf[blend_pos].x, blendmask, 32);
			info.x = blendVal<<RI_HASH_SHIFT | span;
			info.y = blndBuf[blend_pos].y;
			kv_push(mm128_t, km, *p, info);
		}else{
			calc_blend_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, hashVal, blendmask, 32);
		}
    }
}

void ri_sketch_min(void *km, const float* s_values, uint32_t id, int strand, int len, int w, int e, int q, int lq, int k, mm128_v *p){
	
	int j, buf_pos, min_pos;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };

	int step = 1;//TODO: make this an argument
	uint32_t span = 6+e-1; //for now single event is considered to span 6 bases.

	assert(len > 0 && (w > 0 && w < 256) && (e > 0 && e <= 10));
	memset(buf, 0xff, w * 16);
	kv_resize(mm128_t, km, *p, p->n + (len/step)/w);
	
	const uint32_t quant_bit = lq+2; 
	const uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1; //for least sig. 5 bits

	bool sigBufFull = false;
	uint32_t i, l, sigBufPos = 0;
	uint32_t l_sigpos = 0; //last signal position
	uint32_t signal = 0, tmpQuantSignal = 0;
	uint64_t quantVal = 0, hashVal = 0;

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

    for (i = l = buf_pos = min_pos = 0; i < len; i += step) {
        if(i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) continue;

		l++;
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		l_sigpos = i;
		signal = *((uint32_t*)&s_values[i]);
		tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);
		// signal = quantize_float_uint32(s_values[i], -2.75f, 2.75f);
		// tmpQuantSignal = signal&mask_quant;

		quantVal = (quantVal<<quant_bit|tmpQuantSignal)&mask_events;

		sigBuf[sigBufPos].y = id_shift | (uint32_t)i<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = true; sigBufPos = 0;}
		sigBuf[sigBufPos].x = hash64(quantVal, mask)<<RI_HASH_SHIFT | span;

		if(!sigBufFull) continue;

		info.x = sigBuf[sigBufPos].x;
		info.y = sigBuf[sigBufPos].y;

		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + e - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + e && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + e - 1 && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest e-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + e - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
    }
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, km, *p, min);

}

void ri_sketch_reg(void *km, const float* s_values, uint32_t id, int strand, int len, int e, int q, int lq, int k, mm128_v *p){

	int step = 1;//TODO: make this an argument
	uint32_t span = (k+e-1)*step; //for now single event is considered to span 6 bases.
	
	uint32_t quant_bit = lq+2; 
	uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1; //for least sig. 5 bits

	bool sigBufFull = false;
	uint32_t i, sigBufPos = 0, l_sigpos = 0; //last signal position
	uint32_t signal = 0, tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	kv_resize(mm128_t, km, *p, p->n + len/step);

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

    for (i = 0; i < len; i += step) {
        if((i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) || s_values[i] == RI_MASK_SIGNAL) continue;

		l_sigpos = i;
		signal = *((uint32_t*)&s_values[i]);
		tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);
		// signal = quantize_float_uint32(s_values[i], -2.75f, 2.75f);
		// tmpQuantSignal = signal&mask_quant;

		sigBuf[sigBufPos].y = id_shift | (uint32_t)i<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = true; sigBufPos = 0;}

		quantVal = (quantVal<<quant_bit|tmpQuantSignal)&mask_events;
		sigBuf[sigBufPos].x = (hash64(quantVal, mask)<<RI_HASH_SHIFT) | span;
		
		if(!sigBufFull) continue;

		kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);

		// quantVal = (quantVal>>quant_bit<<quant_bit|(tmpQuantSignal-1))&mask_events;
		// sigBuf[sigBufPos].x = hash64(quantVal, mask)<<RI_HASH_SHIFT | span;
		// kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);

		// quantVal = (quantVal>>quant_bit<<quant_bit|(tmpQuantSignal+1))&mask_events;
		// sigBuf[sigBufPos].x = hash64(quantVal, mask)<<RI_HASH_SHIFT | span;
		// kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);
    }
}

void ri_sketch(void *km, const float* s_values, uint32_t id, int strand, int len, int w, int e, int n, int q, int lq, int k, mm128_v *p){

	assert(e > 1 && e < 10);

	if(w & n) ri_sketch_blend_min(km, s_values, id, strand, len, w, e, n, q, lq, k, p);
	else if(w) ri_sketch_min(km, s_values, id, strand, len, w, e, q, lq, k, p);
	else if(n) ri_sketch_blend(km, s_values, id, strand, len, e, n, q, lq, k, p);
	else ri_sketch_reg(km, s_values, id, strand, len, e, q, lq, k, p);
}

void ri_seed_mz_flt(void *km, mm128_v *mv, int32_t q_occ_max, float q_occ_frac)
{
	mm128_t *a;
	size_t i, j, st;
	if (mv->n <= q_occ_max || q_occ_frac <= 0.0f || q_occ_max <= 0) return;
	KMALLOC(km, a, mv->n);
	for (i = 0; i < mv->n; ++i)
		a[i].x = mv->a[i].x, a[i].y = i;
	radix_sort_128x(a, a + mv->n);
	for (st = 0, i = 1; i <= mv->n; ++i) {
		if (i == mv->n || a[i].x != a[st].x) {
			int32_t cnt = i - st;
			if (cnt > q_occ_max && cnt > mv->n * q_occ_frac)
				for (j = st; j < i; ++j)
					mv->a[a[j].y].x = 0;
			st = i;
		}
	}
	ri_kfree(km, a);
	for (i = j = 0; i < mv->n; ++i)
		if (mv->a[i].x != 0)
			mv->a[j++] = mv->a[i];
	mv->n = j;
}

void ri_mapopt_init(ri_mapopt_t *opt)
{
	memset(opt, 0, sizeof(ri_mapopt_t));

	opt->bp_per_sec = 450;
	opt->sample_rate = 4000;
	opt->chunk_size = 4000;

	//chaining
	opt->max_gap_length = 2000;
	opt->max_target_gap_length = 5000;
	opt->chaining_band_length = 5000;
	opt->max_num_skips = 25;
	opt->min_num_anchors = 2;
	opt->num_best_chains = 3;
	opt->min_chaining_score = 10.0f;

	opt->step_size = 1; //read_seeding_step_size
	opt->min_events = 50;
	opt->max_num_chunk = 30;//max_num_chunks
	opt->min_chain_anchor = 10; //stop_mapping_min_num_anchors
	opt->min_chain_anchor_out = 10;//output_mapping_min_num_anchors

	opt->min_bestmap_ratio = 1.4f; //stop_mapping_ratio
	opt->min_bestmap_ratio_out = 1.2f; //output_mapping_ratio

	opt->min_meanmap_ratio = 5; //stop_mapping_mean_ratio
	opt->min_meanmap_ratio_out = 5; //output_mapping_mean_ratio

	opt->mini_batch_size = 500000000;

	//Default options for event detection. TODO: Make it flexible so that we can change them to RNA values as well
	opt->window_length1 = 3;
    opt->window_length2 = 6;
    opt->threshold1 = 4.30265f;  // 4.60409f,//1.4f,
    opt->threshold2 = 2.57058f;  // 3.16927f,//9.0f,
    opt->peak_height = 1.0f;      // 0.2f

	opt->t_threshold = 1.5f;
	opt->tn_samples = 5;
	opt->ttest_freq = 500;
	opt->tmin_reads = 500;

	//TODO: RNA values:
	// opt->window_length1 = 7,
	// opt->window_length2 = 14,
	// opt->threshold1 = 2.5f,
	// opt->threshold2 = 9.0f,
	// opt->peak_height = 1.0f;
}

#include "kseq.h"
#include <zlib.h>
#include "hdf5_tools.hpp"
// #include <slow5/slow5.h>

typedef struct {
	uint32_t rid, l_sig;
	char *name;

	float dig, ran, offset;
	float* sig;
} ri_sig_t;

struct ri_sig_file_s {
	// gzFile fp;
	// kseq_t *ks;
	// ri_sig_t s;
	char** raw_path;
	char** ch_path;
	int num_read;
	int cur_read;
	hdf5_tools::File* fp;
};

struct ri_sig_file_s;
typedef struct ri_sig_file_s ri_sig_file_t;

struct ri_tbuf_s {
	void *km;
	int rep_len, frag_gap;
};

// memory buffer for thread-local storage during mapping
typedef struct ri_tbuf_s ri_tbuf_t;

typedef struct {
	int n_processed, n_threads, n_fp, cur_fp, n_f, cur_f;
	int64_t mini_batch_size;
	const ri_mapopt_t *opt;
	char **f;
	ri_sig_file_t *fp;
	const ri_idx_t *ri;
	const char **fn;
	uint32_t su_nreads, su_nestimations, ab_count, su_cur;
	float** su_estimations;
	uint32_t* su_c_estimations;
	int su_stop;
} pipeline_mt;

typedef struct {
	const pipeline_mt *p;
    int n_sig;
	ri_sig_t** sig;
	// int* n_reg, *seg_off, *n_seg, *rep_len, *frag_gap;
	ri_reg1_t** reg;
	ri_tbuf_t** buf;
} step_mt;

ri_tbuf_t *ri_tbuf_init(void)
{
	ri_tbuf_t *b;
	b = (ri_tbuf_t*)calloc(1, sizeof(ri_tbuf_t));
	b->km = ri_km_init();
	return b;
}

void ri_tbuf_destroy(ri_tbuf_t *b)
{
	if (b == 0) return;
	ri_km_destroy(b->km);
	free(b);
}

ri_sig_file_t *ri_sig_open(const char *fn)
{
	ri_sig_file_t *fp;

	hdf5_tools::File* fast5_file = new hdf5_tools::File();
	fast5_file->open(std::string(fn));
	// gzFile f;
	// f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (!fast5_file->is_open()) return 0;

	fp = (ri_sig_file_t*)calloc(1, sizeof(ri_sig_file_t));
	fp->fp = fast5_file;

	bool is_single = false;
	std::vector<std::string> fast5_file_groups = fast5_file->list_group("/");
	fp->num_read = fast5_file_groups.size();
	fp->ch_path = (char**)calloc(fp->num_read, sizeof(char*));
	fp->raw_path = (char**)calloc(fp->num_read, sizeof(char*));

	for (std::string &group : fast5_file_groups) {
		if (group == "Raw") {
			is_single = true;
			break;
		}
	}

	std::string raw_path;
	std::string ch_path;
	int i = 0;

	if (is_single) {
		ch_path = "/UniqueGlobalKey/channel_id";
		for (std::string &read : fast5_file->list_group("/Raw/Reads")) {
			raw_path = "/Raw/Reads/" + read;
			if(i == fp->num_read){
				fprintf(stdout, "ERROR: More reads than previously predicted (%d). Stopped reading the reads here.\n", fp->num_read);
				break;
			}
			fp->ch_path[i] = strdup(ch_path.c_str());
			fp->raw_path[i++] = strdup(raw_path.c_str());
		}
	} else {
		for (std::string &read : fast5_file_groups) {
			raw_path = "/" + read + "/Raw";
			ch_path = "/" + read + "/channel_id";
			fp->ch_path[i] = strdup(ch_path.c_str());
			fp->raw_path[i++] = strdup(raw_path.c_str());
		}
	}

	fp->num_read = i;
	fp->cur_read = 0;
	return fp;
}

void ri_sig_close(ri_sig_file_t *fp)
{
	if(!fp) return;
	// gzclose(fp->fp);
	fp->fp->close();
	for(int i = 0; i < fp->num_read; ++i){
		if(fp->ch_path[i])free(fp->ch_path[i]);
		if(fp->raw_path[i])free(fp->raw_path[i]);
	}
	free(fp);
}

static ri_sig_file_t *open_sig(const char *fn) //TODO: make this a part of the pipeline. Sequntially reading from many FAST5 files creates an overhead
{
	ri_sig_file_t *fp;
	int i, j;
	fp = (ri_sig_file_t*)calloc(1,sizeof(ri_sig_file_t));
	if ((fp = ri_sig_open(fn)) == 0) {
		fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i], strerror(errno));
		ri_sig_close(fp);
		return 0;
	}
	return fp;
}

static ri_sig_file_t **open_sigs(int n, const char **fn) //TODO: make this a part of the pipeline. Sequntially reading from many FAST5 files creates an overhead
{
	ri_sig_file_t **fp;
	int i, j;
	fp = (ri_sig_file_t**)calloc(n, sizeof(ri_sig_file_t*));
	for (i = 0; i < n; ++i) {
		if ((fp[i] = ri_sig_open(fn[i])) == 0) {
			fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i], strerror(errno));
	// 		if (mm_verbose >= 1)
	// 			fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i], strerror(errno));
			for (j = 0; j < i; ++j)
				ri_sig_close(fp[j]);
			free(fp);
			return 0;
		}
	}
	return fp;
}

typedef struct ri_detect_s {
	int DEF_PEAK_POS;
	float DEF_PEAK_VAL;
	float *sig;
	uint32_t s_len;
	float threshold;
	uint32_t window_length;
	uint32_t masked_to;
	int peak_pos;
	float peak_value;
	bool valid_peak;
};

typedef struct ri_detect_s ri_detect_t;

static inline void comp_prefix_prefixsq(const float *sig, uint32_t s_len, float* prefix_sum, float* prefix_sum_square) {
	
	assert(s_len > 0);

	prefix_sum[0] = 0.0f;
	prefix_sum_square[0] = 0.0f;
	for (uint32_t i = 0; i < s_len; ++i) {
		prefix_sum[i+1] = prefix_sum[i] + sig[i];
		prefix_sum_square[i+1] = prefix_sum_square[i] + sig[i]*sig[i];
	}
}

#include <float.h>
static inline float* comp_tstat(void *km, const float *prefix_sum, const float *prefix_sum_square, uint32_t s_len, uint32_t w_len) {
  
  	const float eta = FLT_MIN;

	// kvec_t(float) tstat = {0,0,0};
	// kv_resize(float, 0, tstat, s_len+1);
	// kv_pushp(float, 0, tstat, &s);

	float* tstat = (float*)ri_kcalloc(km, s_len+1, sizeof(float));
	// Quick return:
	// t-test not defined for number of points less than 2
	// need at least as many points as twice the window length
	if (s_len < 2*w_len || w_len < 2) return tstat;
	// fudge boundaries
	memset(tstat, 0, w_len*sizeof(float));
	// get to work on the rest
	for (uint32_t i = w_len; i <= s_len - w_len; ++i) {
		float sum1 = prefix_sum[i];
		float sumsq1 = prefix_sum_square[i];
		if (i > w_len) {
			sum1 -= prefix_sum[i - w_len];
			sumsq1 -= prefix_sum_square[i - w_len];
		}
		float sum2 = prefix_sum[i + w_len] - prefix_sum[i];
		float sumsq2 = prefix_sum_square[i + w_len] - prefix_sum_square[i];
		float mean1 = sum1 / w_len;
		float mean2 = sum2 / w_len;
		float combined_var = sumsq1 / w_len - mean1 * mean1 + sumsq2 / w_len - mean2 * mean2;
		// Prevent problem due to very small variances
		combined_var = fmaxf(combined_var, eta);
		// t-stat
		//  Formula is a simplified version of Student's t-statistic for the
		//  special case where there are two samples of equal size with
		//  differing variance
		const float delta_mean = mean2 - mean1;
		tstat[i] = fabs(delta_mean) / sqrt(combined_var / w_len);
	}
	// fudge boundaries
	memset(tstat+s_len-w_len+1, 0, (w_len)*sizeof(float));

	return tstat;
}

static inline uint32_t gen_peaks(ri_detect_t *short_detector, ri_detect_t *long_detector, const float peak_height, uint32_t* peaks) {
	
	assert(short_detector->s_len == long_detector->s_len);

	uint32_t curInd = 0;

	uint32_t ndetector = 2;
	ri_detect_t *detectors[ndetector];  // = {short_detector, long_detector};
	detectors[0] = short_detector;
	detectors[1] = long_detector;
	for (uint32_t i = 0; i < short_detector->s_len; i++) {
		for (uint32_t k = 0; k < ndetector; k++) {
			ri_detect_t *detector = detectors[k];
			// Carry on if we've been masked out
			if (detector->masked_to >= i) continue;

			float current_value = detector->sig[i];
			if (detector->peak_pos == detector->DEF_PEAK_POS) {
				// CASE 1: We've not yet recorded a maximum
				if (current_value < detector->peak_value) {
					// Either record a deeper minimum...
					detector->peak_value = current_value;
				} else if (current_value - detector->peak_value > peak_height) {  // TODO(Haowen): this might cause overflow, need to fix this
					// ...or we've seen a qualifying maximum
					detector->peak_value = current_value;
					detector->peak_pos = i;
					// otherwise, wait to rise high enough to be considered a peak
				}
			} else {
				// CASE 2: In an existing peak, waiting to see if it is good
				if (current_value > detector->peak_value) {
					// Update the peak
					detector->peak_value = current_value;
					detector->peak_pos = i;
				}
				// Dominate other tstat signals if we're going to fire at some point
				if (detector == short_detector) {
					if (detector->peak_value > detector->threshold) {
						long_detector->masked_to = detector->peak_pos + detector->window_length;
						long_detector->peak_pos = long_detector->DEF_PEAK_POS;
						long_detector->peak_value = long_detector->DEF_PEAK_VAL;
						long_detector->valid_peak = false;
					}
				}
				// Have we convinced ourselves we've seen a peak
				if (detector->peak_value - current_value > peak_height && detector->peak_value > detector->threshold) {
					detector->valid_peak = true;
				}
				// Finally, check the distance if this is a good peak
				if (detector->valid_peak && (i - detector->peak_pos) > detector->window_length / 2) {
					// Emit the boundary and reset
					peaks[curInd++] = detector->peak_pos;
					detector->peak_pos = detector->DEF_PEAK_POS;
					detector->peak_value = current_value;
					detector->valid_peak = false;
				}
			}
		}
	}

	return curInd;
}

static inline float* gen_events(void *km, const uint32_t *peaks, uint32_t peak_size, const float *prefix_sum, const float *prefix_sum_square, uint32_t s_len, uint32_t* n_events) {
	// Count number of events found
	uint32_t n_ev = 1;
	double mean = 0, std_dev = 0, sum = 0, sum2 = 0;

	for (uint32_t i = 1; i < peak_size; ++i)
		if (peaks[i] > 0 && peaks[i] < s_len)
			n_ev++;

	if(!n_ev) return 0;

	float* events = (float*)ri_kmalloc(km, n_ev*sizeof(float));

	// events[0] = (prefix_sum[peaks[0]])/peaks[0];

	float l_prefixsum = 0, l_peak = 0;

	// First event -- starts at zero
	// gen_event(0, peaks[0], prefix_sum, prefix_sum_square, s_len, events);
	// Other events -- peak[i-1] -> peak[i]
	for (uint32_t pi = 0; pi < n_ev - 1; pi++){
		// events[pi] = (prefix_sum[peaks[pi]] - prefix_sum[peaks[pi-1]]) / peaks[pi] - peaks[pi-1];
		events[pi] = (prefix_sum[peaks[pi]] - l_prefixsum)/(peaks[pi]-l_peak);
		sum += events[pi];
		sum2 += events[pi]*events[pi];
		l_prefixsum = prefix_sum[peaks[pi]];
		l_peak = peaks[pi];
	}
		// events[pi] = (prefix_sum[peaks[pi]] - prefix_sum[peaks[pi-1]]) / peaks[pi] - peaks[pi-1];
		// gen_event(peaks[pi - 1], peaks[pi], prefix_sum, prefix_sum_square, s_len, events+pi);

	// Last event -- ends at s_len
	events[n_ev-1] = (prefix_sum[s_len] - l_prefixsum)/(s_len-l_peak);
	sum += events[n_ev-1];
	sum2 += events[n_ev-1]*events[n_ev-1];

	// gen_event(peaks[n_ev - 2], s_len, prefix_sum, prefix_sum_square, s_len, events+n_ev-1);

	//normalization
	mean = sum/n_ev;
	std_dev = sqrt(sum2/n_ev - (mean)*(mean));

	for(uint32_t i = 0; i < n_ev; ++i){
		events[i] = (events[i]-mean)/std_dev;
	}

	(*n_events) = n_ev;
	return events;
}

static float* detect_events(void *km, uint32_t s_len, const float* sig, const ri_mapopt_t *opt, uint32_t *n) // kt_for() callback
{
	float* prefix_sum = (float*)ri_kcalloc(km, s_len+1, sizeof(float));
	float* prefix_sum_square = (float*)ri_kcalloc(km, s_len+1, sizeof(float));

	comp_prefix_prefixsq(sig, s_len, prefix_sum, prefix_sum_square);
	float* tstat1 = comp_tstat(km, prefix_sum, prefix_sum_square, s_len, opt->window_length1);
	float* tstat2 = comp_tstat(km, prefix_sum, prefix_sum_square, s_len, opt->window_length2);
	ri_detect_t short_detector = {.DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX, .sig = tstat1, .s_len = s_len, .threshold = opt->threshold1,
								  .window_length = opt->window_length1, .masked_to = 0, .peak_pos = -1, .peak_value = FLT_MAX, .valid_peak = false};

	ri_detect_t long_detector = {.DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX, .sig = tstat2, .s_len = s_len, .threshold = opt->threshold2,
								 .window_length = opt->window_length2, .masked_to = 0, .peak_pos = -1, .peak_value = FLT_MAX, .valid_peak = false};
	uint32_t* peaks = (uint32_t*)ri_kmalloc(km, s_len * sizeof(uint32_t));
	uint32_t n_peaks = gen_peaks(&short_detector, &long_detector, opt->peak_height, peaks);
	float* events = 0;
	if(n_peaks > 0) events = gen_events(km, peaks, n_peaks, prefix_sum, prefix_sum_square, s_len, n);
	ri_kfree(km, tstat1); ri_kfree(km, tstat2); ri_kfree(km, prefix_sum); ri_kfree(km, prefix_sum_square); ri_kfree(km, peaks);

	return events;
}

void comp_mapq(std::vector<ri_chain_t> &chains) {
  if (chains.size() == 1) {
    chains[0].mapq = 60;
    return;
  } else {
    int mapq = 40 * (1 - chains[1].score / chains[0].score);  // * std::min((size_t)1, chains[0].n_anchors / 20) * log(chains[0].score);
    if (mapq > 60) {
      mapq = 60;
    }
    if (mapq < 0) {
      mapq = 0;
    }
    chains[0].mapq = (uint8_t)mapq;
  }
}

void gen_primary_chains(void *km, std::vector<ri_chain_t> &chains) {

	std::sort(chains.begin(), chains.end(), std::greater<ri_chain_t>());
	std::vector<ri_chain_t> primary_chains;
	primary_chains.reserve(chains.size());
	primary_chains.emplace_back(chains[0]);
	for (uint32_t ci = 1; ci < chains.size(); ++ci) {
		bool is_primary = true;
		if (chains[ci].score < primary_chains.back().score / 3) {
			break;
		} else {
			for (uint32_t pi = 0; pi < primary_chains.size(); ++pi) {
				if (chains[ci].reference_sequence_index == primary_chains[pi].reference_sequence_index) {
					if (std::max(chains[ci].start_position, primary_chains[pi].start_position) <= std::min(chains[ci].end_position, primary_chains[pi].end_position)) {
						is_primary = false;
						break;
					}
				}
			}
		}
		if (is_primary) {
			primary_chains.emplace_back(chains[ci]);
		}else{
			if(chains[ci].anchors)ri_kfree(km, chains[ci].anchors);
		}
	}
	chains.swap(primary_chains);
}

void traceback_chains(void *km,
    int min_num_anchors, int strand, size_t chain_end_anchor_index,
    uint32_t chain_target_signal_index,
    const std::vector<float> &chaining_scores,
    const std::vector<size_t> &chaining_predecessors,
    const std::vector<std::vector<ri_anchor_t> > &anchors_fr,
    std::vector<bool> &anchor_is_used, std::vector<ri_chain_t> &chains) {

  if (!anchor_is_used[chain_end_anchor_index]) {
    std::vector<ri_anchor_t> anchors;
    anchors.reserve(100);
    bool stop_at_an_used_anchor = false;
    size_t chain_start_anchor_index = chain_end_anchor_index;
    // Add the end anchor
    anchors.push_back(anchors_fr[chain_target_signal_index][chain_start_anchor_index]);
    // The end anchor has an used predecessor
    if (chaining_predecessors[chain_start_anchor_index] != chain_start_anchor_index && anchor_is_used[chaining_predecessors[chain_start_anchor_index]]) {
      stop_at_an_used_anchor = true;
    }
    anchor_is_used[chain_start_anchor_index] = true;
    uint32_t chain_num_anchors = 1;
    while (chaining_predecessors[chain_start_anchor_index] != chain_start_anchor_index && !anchor_is_used[chaining_predecessors[chain_start_anchor_index]]) {
      chain_start_anchor_index = chaining_predecessors[chain_start_anchor_index];
      anchors.push_back(anchors_fr[chain_target_signal_index][chain_start_anchor_index]);
      if (chaining_predecessors[chain_start_anchor_index] != chain_start_anchor_index && anchor_is_used[chaining_predecessors[chain_start_anchor_index]]) {
        stop_at_an_used_anchor = true;
      }
      anchor_is_used[chain_start_anchor_index] = true;
      ++chain_num_anchors;
    }
    if (chain_num_anchors >= (uint32_t)min_num_anchors) {
      float adjusted_chaining_score = chaining_scores[chain_end_anchor_index];
      if (stop_at_an_used_anchor) {
        adjusted_chaining_score -= chaining_scores[chaining_predecessors[chain_start_anchor_index]];
      }
	  ri_anchor_t* a_anchors = (ri_anchor_t*)ri_kmalloc(km,anchors.size()*sizeof(ri_anchor_t));
	  std::copy(anchors.begin(), anchors.end(), a_anchors);
      chains.emplace_back(ri_chain_t{adjusted_chaining_score, chain_target_signal_index, 
		  				  anchors_fr[chain_target_signal_index][chain_start_anchor_index].target_position,
                          anchors_fr[chain_target_signal_index][chain_end_anchor_index].target_position,
                          chain_num_anchors, 0, strand, a_anchors});
    }
  }
}

bool compare(const std::pair<float, size_t> &left, const std::pair<float, size_t> &right) {
	if (left.first > right.first) return true;
	else if (left.first == right.first) return (left.second > right.second);
	else return false;
}

void gen_chains(void *km, const ri_idx_t *ri, const float* sig, const uint32_t l_sig, const uint32_t q_offset, const size_t n_seq, ri_reg1_t* reg, const ri_mapopt_t *opt){

	// Chaining parameters
	int max_gap_length = opt->max_gap_length;
	int max_target_gap_length = opt->max_target_gap_length;
	int chaining_band_length = opt->chaining_band_length;
	int max_num_skips = opt->max_num_skips;
	int min_num_anchors = opt->min_num_anchors;
	int num_best_chains = opt->num_best_chains;
	float min_chaining_score = opt->min_chaining_score;

	uint32_t mask = (1ULL<<31)-1;

	// kvec_t(ri_anchor_t)* anchors_fr[2] = {0};
	// kvec_t(ri_anchor_t)* anchors_f = anchors_fr;
	// kvec_t(ri_anchor_t)* anchors_r = anchors_fr+1;

	// anchors_f = (kvec_t(ri_anchor_t)*)calloc(n_seq, sizeof(kvec_t(ri_anchor_t)));
	// anchors_r = (kvec_t(ri_anchor_t)*)calloc(n_seq, sizeof(kvec_t(ri_anchor_t)));
	// memset(anchors_f, 0, n_seq*sizeof(kvec_t(ri_anchor_t)));
	// memset(anchors_r, 0, n_seq*sizeof(kvec_t(ri_anchor_t)));

	std::vector<std::vector<std::vector<ri_anchor_t> > > anchors_fr(2);
	anchors_fr[0] = std::vector<std::vector<ri_anchor_t> >(n_seq);
	anchors_fr[1] = std::vector<std::vector<ri_anchor_t> >(n_seq);
	std::vector<std::vector<ri_anchor_t> > &anchors_f = anchors_fr[0];
	std::vector<std::vector<ri_anchor_t> > &anchors_r = anchors_fr[1];

	// Get anchors in previous chains
	ri_chain_t* previous_chains = reg->chains;
	if (reg->n_chains > 0) {
		for (uint32_t c_ind = 0; c_ind < reg->n_chains; ++c_ind) {
			int strand = previous_chains[c_ind].strand;
			uint32_t reference_sequence_index = previous_chains[c_ind].reference_sequence_index;
			// anchors_fr[strand][reference_sequence_index].reserve(previous_chains[c_ind].n_anchors);
			anchors_fr[strand][reference_sequence_index].reserve(previous_chains[c_ind].n_anchors);
			// kv_resize(ri_anchor_t, 0, anchors_fr[strand][reference_sequence_index], anchors_fr[strand][reference_sequence_index].n + previous_chains[c_ind].n_anchors)
			for (uint32_t a_ind = 0; a_ind < previous_chains[c_ind].n_anchors; ++a_ind) {
				anchors_fr[strand][reference_sequence_index].emplace_back(previous_chains[c_ind].anchors[a_ind]);
				// kv_push(ri_anchor_t, 0, anchors_fr[strand][reference_sequence_index], previous_chains[c_ind].anchors[a_ind]);
			}
		}
	}

	if(reg->chains){
		for(int i = 0; i < reg->n_chains; ++i) if(reg->chains[i].anchors) ri_kfree(km, reg->chains[i].anchors);
		ri_kfree(km, reg->chains); reg->n_chains = 0;
	}

	uint32_t i, pi;
	mm128_v riv = {0,0,0};
	uint64_t hashVal, keyval;
	ri_sketch(km, sig, 0, 0, l_sig, ri->w, ri->e, ri->n, ri->q, ri->lq, ri->k, &riv);
	//   ri_seed_mz_flt(0, &riv, 1000, 0.01f);

	for (i = 0; i <= riv.n; ++i) {
		hashVal = riv.a[i].x>>RI_HASH_SHIFT;
		pi = (uint32_t)riv.a[i].y>>RI_POS_SHIFT;

		const uint64_t *cr;
		int t;
		cr = ri_idx_get(ri, hashVal, &t);

		if(t == 0) continue;

		for(int s = 0; s < t; ++s){
			keyval = cr[s];
			uint32_t t_ind = (uint32_t)(keyval>>RI_ID_SHIFT), target_signal_position = (uint32_t)(keyval>>RI_POS_SHIFT)&mask;
			anchors_fr[(keyval&1)][t_ind].emplace_back(ri_anchor_t{target_signal_position,  pi + q_offset});
			// kv_push(ri_anchor_t, 0, anchors_fr[(keyval&1)][t_ind], ri_anchor_t{target_signal_position, pi + q_offset})
		}
	}

	ri_kfree(km, riv.a);

	// Sort the anchors based on their occurrence on target signal
	for (size_t t_ind = 0; t_ind < n_seq; ++t_ind) {
		// std::sort(anchors_f[t_ind].a, anchors_f[t_ind].a + anchors_f[t_ind].n);
		// std::sort(anchors_r[t_ind].a, anchors_r[t_ind].a + anchors_r[t_ind].n);
		std::sort(anchors_f[t_ind].begin(), anchors_f[t_ind].end());
		std::sort(anchors_r[t_ind].begin(), anchors_r[t_ind].end());
	} 

	// Chaining DP done on each individual target signal
	float max_chaining_score = 0;  // std::numeric_limits<float>::min();
	std::vector<ri_chain_t> chains;
	for (size_t target_signal_index = 0; target_signal_index < n_seq; ++target_signal_index) {
		for (int strand_i = 0; strand_i < 2; ++strand_i) {
			std::vector<float> chaining_scores;
			chaining_scores.reserve(anchors_fr[strand_i][target_signal_index].size());
			std::vector<size_t> chaining_predecessors;
			chaining_predecessors.reserve(anchors_fr[strand_i][target_signal_index].size());
			std::vector<bool> anchor_is_used;
			anchor_is_used.reserve(anchors_fr[strand_i][target_signal_index].size());
			std::vector<std::pair<float, size_t> > end_anchor_index_chaining_scores;
			end_anchor_index_chaining_scores.reserve(10);
			for (size_t anchor_index = 0; anchor_index < anchors_fr[strand_i][target_signal_index].size(); ++anchor_index) {
				float distance_coefficient = 1; /*- 0.2 * anchors_fr[strand_i][target_signal_index][anchor_index].distance / search_radius;*/
				chaining_scores.emplace_back(distance_coefficient * ri->e);
				chaining_predecessors.emplace_back(anchor_index);
				anchor_is_used.push_back(false);
				int32_t current_anchor_target_position = anchors_fr[strand_i][target_signal_index][anchor_index].target_position;
				int32_t current_anchor_query_position = anchors_fr[strand_i][target_signal_index][anchor_index].query_position;
				int32_t start_anchor_index = 0;
				if (anchor_index > (size_t)chaining_band_length) start_anchor_index = anchor_index - chaining_band_length;

				int32_t previous_anchor_index = anchor_index - 1;
				int32_t num_skips = 0;
				for (; previous_anchor_index >= start_anchor_index; --previous_anchor_index) {
					int32_t previous_anchor_target_position = anchors_fr[strand_i][target_signal_index][previous_anchor_index].target_position;
					int32_t previous_anchor_query_position = anchors_fr[strand_i][target_signal_index][previous_anchor_index].query_position;

					if (previous_anchor_query_position == current_anchor_query_position) continue;
					if (previous_anchor_target_position == current_anchor_target_position) continue;
					if (previous_anchor_target_position + max_target_gap_length < current_anchor_target_position) break;

					int32_t target_position_diff = current_anchor_target_position - previous_anchor_target_position;
					assert(target_position_diff > 0);
					int32_t query_position_diff = current_anchor_query_position - previous_anchor_query_position;
					float current_chaining_score = 0;

					if (query_position_diff < 0) continue;
					
					float matching_dimensions = std::min(std::min(target_position_diff, query_position_diff), ri->e) * distance_coefficient;
					int gap_length = std::abs(target_position_diff - query_position_diff);
					float gap_scale = target_position_diff > 0? (float)query_position_diff / target_position_diff: 1;
					// float gap_cost = 0;
					// if (gap_length != 0) {
					if (gap_length < max_gap_length && gap_scale < 5 && gap_scale > 0.75) {
						current_chaining_score = chaining_scores[previous_anchor_index] + matching_dimensions;  // - gap_cost;
					}

					if (current_chaining_score > chaining_scores[anchor_index]) {
						chaining_scores[anchor_index] = current_chaining_score;
						chaining_predecessors[anchor_index] = previous_anchor_index;
						--num_skips;
					} else {
						++num_skips;
						if (num_skips > max_num_skips) break;
					}
				}
				// Update chain with max score
				if (chaining_scores[anchor_index] > max_chaining_score) {
					max_chaining_score = chaining_scores[anchor_index];
				}
				if (chaining_scores.back() >= min_chaining_score &&
					chaining_scores.back() > max_chaining_score / 2) {
					end_anchor_index_chaining_scores.emplace_back(chaining_scores.back(), anchor_index);
				}
			}
			// Sort a vector of <end anchor index, chaining score>
			std::sort(end_anchor_index_chaining_scores.begin(), end_anchor_index_chaining_scores.end(), compare);
			// Traceback all reg->chains from higest score to lowest
			for (size_t anchor_index = 0; anchor_index < end_anchor_index_chaining_scores.size() && anchor_index < (size_t)num_best_chains; ++anchor_index) {
				traceback_chains(km, min_num_anchors, strand_i, end_anchor_index_chaining_scores[anchor_index].second, target_signal_index, chaining_scores, 
								 chaining_predecessors, anchors_fr[strand_i], anchor_is_used, chains);

				if (chaining_scores[end_anchor_index_chaining_scores[anchor_index].second] < max_chaining_score / 2) break;
			}
		}
	}

	if (chains.size() > 0) {
		// Generate primary reg->chains
		gen_primary_chains(km, chains);
		// Compute MAPQ
		comp_mapq(chains);

		reg->n_chains = chains.size();
		reg->chains = (ri_chain_t*)ri_kmalloc(km, chains.size()*sizeof(ri_chain_t));
		std::copy(chains.begin(), chains.end(), reg->chains);
	}
}

//returns n_regs // currently we report one mapping
void ri_map_frag(const ri_idx_t *ri, const uint32_t s_len, const float *sig, ri_reg1_t* reg, ri_tbuf_t *b, const ri_mapopt_t *opt, const char *qname){
	
	uint32_t n_events = 0;
	// ri_reg1_t* reg = 0;
	uint32_t event_off = 0;

	float* events = detect_events(b->km, s_len, sig, opt, &n_events);

	if(n_events < opt->min_events) {
		if(events)ri_kfree(b->km, events);
		return;
	}

	gen_chains(b->km, ri, events, n_events, reg->offset, ri->n_seq, reg, opt);
	reg->offset += n_events;

	if(events)ri_kfree(b->km, events);
}

static void map_worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_mt *s = (step_mt*)_data; //s->sig and s->n_sig (signals read in this step and num of them)
	const ri_mapopt_t *opt = s->p->opt;
	double t = 0.0;
	ri_tbuf_t* b = s->buf[tid];
	ri_reg1_t* reg0 = s->reg[i];
	ri_sig_t* sig = s->sig[i];

	uint32_t l_chunk = opt->chunk_size;
	uint32_t max_chunk =  opt->max_num_chunk;
	uint32_t qlen = sig->l_sig;
	uint32_t s_qs, s_qe;

	uint32_t c_count = 0;

	t = ri_realtime();
	for (s_qs = c_count = 0; s_qs < qlen && c_count < max_chunk; s_qs += l_chunk, ++c_count) {
		s_qe = s_qs + l_chunk;
		if(s_qe > qlen) s_qe = qlen;

		ri_map_frag(s->p->ri, (const uint32_t)s_qe-s_qs, (const float*)&(sig->sig[s_qs]), reg0, b, opt, sig->name);
		
		if (reg0->n_chains >= 2) {
			if (reg0->chains[0].score / reg0->chains[1].score >= opt->min_bestmap_ratio) break;

			float mean_chain_score = 0;
			for (uint32_t c_ind = 0; c_ind < reg0->n_chains; ++c_ind)
				mean_chain_score += reg0->chains[c_ind].score;

			mean_chain_score /= reg0->n_chains;
			if (reg0->chains[0].score >= opt->min_meanmap_ratio * mean_chain_score) break;
		} else if (reg0->n_chains == 1 && reg0->chains[0].n_anchors >= opt->min_chain_anchor){
			break;
		}
	} double mapping_time = ri_realtime() - t;

	if (c_count > 0 && (s_qs >= qlen || c_count == max_chunk)) --c_count;

	float read_position_scale = ((float)(c_count+1)*l_chunk/reg0->offset) / ((float)opt->sample_rate/opt->bp_per_sec);
	// Save results in vector and output PAF
	
	ri_chain_s* chains = reg0->chains;

	if(!chains) {reg0->n_chains = 0;}

	uint32_t n_anchors0 = (reg0->n_chains)?chains[0].n_anchors:0;
	float mean_chain_score = 0;
	for (uint32_t c_ind = 0; c_ind < reg0->n_chains; ++c_ind)
		mean_chain_score += chains[c_ind].score;

	mean_chain_score /= reg0->n_chains;
	if (n_anchors0 > 0 && ((reg0->n_chains >= 2 && (chains[0].score / chains[1].score >= opt->min_bestmap_ratio_out || 
											   chains[0].score >= opt->min_meanmap_ratio_out * mean_chain_score)) ||
						(reg0->n_chains == 1 && n_anchors0 >= opt->min_chain_anchor_out))) {
		float anchor_ref_gap_avg_length = 0;
		float anchor_read_gap_avg_length = 0;
		for (size_t ai = 0; ai < n_anchors0; ++ai) {
			if (ai < n_anchors0 - 1) {
				anchor_ref_gap_avg_length += chains[0].anchors[ai].target_position - chains[0].anchors[ai + 1].target_position;
				anchor_read_gap_avg_length += chains[0].anchors[ai].query_position - chains[0].anchors[ai + 1].query_position;
			}
		}
		anchor_ref_gap_avg_length /= n_anchors0;
		anchor_read_gap_avg_length /= n_anchors0;
		std::string tags;
		tags.append("mt:f:" + std::to_string(mapping_time * 1000));
		tags.append("\tci:i:" + std::to_string(c_count + 1));
		tags.append("\tsl:i:" + std::to_string(qlen));
		tags.append("\tcm:i:" + std::to_string(n_anchors0));
		tags.append("\tnc:i:" + std::to_string(reg0->n_chains));
		tags.append("\ts1:f:" + std::to_string(chains[0].score));
		tags.append("\ts2:f:" + std::to_string(reg0->n_chains > 1 ? chains[1].score : 0));
		tags.append("\tsm:f:" + std::to_string(mean_chain_score));
		tags.append("\tat:f:" + std::to_string((float)anchor_ref_gap_avg_length));
		tags.append("\taq:f:" + std::to_string((float)anchor_read_gap_avg_length));

		reg0->read_id = sig->rid;
		reg0->ref_id = chains[0].reference_sequence_index;
		reg0->read_name = sig->name;
		reg0->read_length = qlen;
		reg0->read_start_position = (uint32_t)(read_position_scale*chains[0].anchors[n_anchors0-1].query_position);
		reg0->read_end_position = (uint32_t)(read_position_scale * chains[0].anchors[0].query_position);
		reg0->fragment_start_position = chains[0].strand?(uint32_t)(s->p->ri->seq[chains[0].reference_sequence_index].len+1-chains[0].end_position):chains[0].start_position;
		reg0->fragment_length = (uint32_t)(chains[0].end_position - chains[0].start_position + 1);
		reg0->mapq = chains[0].mapq;
		reg0->rev = chains[0].strand;
		reg0->mapped = 1; 
		reg0->tags = strdup(tags.c_str());
	} else {
		std::string tags;
		tags.append("mt:f:" + std::to_string(mapping_time * 1000));
		tags.append("\tci:i:" + std::to_string(c_count + 1));
		tags.append("\tsl:i:" + std::to_string(qlen));
		if (reg0->n_chains >= 1) {
			float anchor_ref_gap_avg_length = 0;
			float anchor_read_gap_avg_length = 0;
			for (size_t ai = 0; ai < n_anchors0; ++ai) {
				if (ai < n_anchors0 - 1) {
					anchor_ref_gap_avg_length += chains[0].anchors[ai].target_position - chains[0].anchors[ai + 1].target_position;
					anchor_read_gap_avg_length += chains[0].anchors[ai].query_position - chains[0].anchors[ai + 1].query_position;
				}
			}
			if(n_anchors0)anchor_ref_gap_avg_length /= n_anchors0;
			if(n_anchors0)anchor_read_gap_avg_length /= n_anchors0;
			tags.append("\tcm:i:" + std::to_string(n_anchors0));
			tags.append("\tnc:i:" + std::to_string(reg0->n_chains));
			tags.append("\ts1:f:" + std::to_string(chains[0].score));
			tags.append("\ts2:f:" + std::to_string(reg0->n_chains > 1 ? chains[1].score: 0));
			tags.append("\tsm:f:" + std::to_string(mean_chain_score));
			tags.append("\tat:f:" + std::to_string((float)anchor_ref_gap_avg_length));
			tags.append("\taq:f:" + std::to_string((float)anchor_read_gap_avg_length));
		}else {
			tags.append("\tcm:i:0");
			tags.append("\tnc:i:0");
			tags.append("\ts1:f:0");
			tags.append("\ts2:f:0");
			tags.append("\tsm:f:0");
			tags.append("\tat:f:0");
			tags.append("\taq:f:0");
		}

		reg0->read_id = sig->rid;
		reg0->ref_id = 0;
		reg0->read_name = sig->name;
		reg0->read_length = qlen;
		reg0->read_start_position = 0;
		reg0->read_end_position = 0;
		reg0->fragment_start_position = 0;
		reg0->fragment_length = 0;
		reg0->mapq = 0;
		reg0->rev = 0;
		reg0->mapped = 0;
		reg0->tags = strdup(tags.c_str());
	}

	for (uint32_t c_ind = 0; c_ind < reg0->n_chains; ++c_ind){
			if(reg0->chains[c_ind].anchors) ri_kfree(b->km, reg0->chains[c_ind].anchors);
	}
	ri_kfree(b->km, reg0->chains);

	if (b->km) {
		ri_km_stat_t kmst;
		ri_km_stat(b->km, &kmst);
		// assert(kmst.n_blocks == kmst.n_cores);
		ri_km_destroy(b->km);
		b->km = ri_km_init();
	}
}

//Check if the input const char* A is a directory
//Generated by GitHub Copilot
int is_dir(const char *A)
{
	struct stat st;
	if (stat(A, &st) == -1) return 0;
	return S_ISDIR(st.st_mode);
}

//Recursively find all files that ends with "fast5" under input directory const char *A
//Generated by GitHub Copilot
void find_fast5(const char *A, ri_char_v *fnames)
{
	if (!is_dir(A)) {
		//TODO: Add slow5 here later
		if (strstr(A, ".fast5")) {
			char** cur_fname;
			kv_pushp(char*, 0, *fnames, &cur_fname);
			(*cur_fname) = strdup(A);
		}
		return;
	}

	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir(A)) != NULL) {
		while ((ent = readdir(dir)) != NULL) {
			char *tmp = (char*)malloc(strlen(A) + strlen(ent->d_name) + 2);
			sprintf(tmp, "%s/%s", A, ent->d_name);
			if (is_dir(tmp)) {
				if (strcmp(ent->d_name, ".") && strcmp(ent->d_name, ".."))
					find_fast5(tmp, fnames);
			} else {
				if (strstr(ent->d_name, ".fast5")) {
					char** cur_fname;
					kv_pushp(char*, 0, *fnames, &cur_fname);
					(*cur_fname) = strdup(tmp);
				}
			}
			free(tmp);
		}
		closedir(dir);
	}
}

void ri_sig_read_sig(ri_sig_file_t* fp, ri_sig_t* s){

	assert(fp->cur_read < fp->num_read);

	s->name = 0;
	for (auto a : fp->fp->get_attr_map(fp->raw_path[fp->cur_read])) {
		if (a.first == "read_id") {
			s->name = strdup(a.second.c_str());
		}
	}

	assert(s->name);

	// float digitisation = 0, range = 0, offset = 0;
	for (auto a : fp->fp->get_attr_map(fp->ch_path[fp->cur_read])) {
		if (a.first == "channel_number") {
		// channel_idx = atoi(a.second.c_str()) - 1;
		} else if (a.first == "digitisation") {
			s->dig = atof(a.second.c_str());
		} else if (a.first == "range") {
			s->ran = atof(a.second.c_str());
		} else if (a.first == "offset") {
			s->offset = atof(a.second.c_str());
		}
	}

	std::string sig_path = std::string(fp->raw_path[fp->cur_read]) + "/Signal";
	std::vector<float> sig;
	fp->fp->read(sig_path, sig);
	// convert to pA
	uint32_t l_sig = 0;
	float scale = s->ran/s->dig;
	for (size_t i = 0; i < sig.size(); i++) {
		if ((sig[i] + s->offset) * scale > 30 && (sig[i] + s->offset) * scale < 200) {
			sig[l_sig] = (sig[i] + s->offset) * scale;
			++l_sig;
		}
	}

	s->sig = (float*)calloc(l_sig, sizeof(float));
	s->l_sig = l_sig;
	std::copy(sig.begin(), sig.begin() + l_sig, s->sig);
	fp->cur_read++;
}

ri_sig_t** ri_sig_read_frag(pipeline_mt *pl, int64_t chunk_size, int *n_){

	int i;
	int64_t size = 0;
	// int n_eof = 0; //at least one file with not end of file
	// kvec_t(ri_sig_t) a = {0,0,0};
	std::vector<ri_sig_t*> sigvec;
	*n_ = 0;
	if (pl->n_fp < 1) return 0;
	// if (a.m == 0) kv_resize(ri_sig_t, 0, a, 256);
	while (pl->fp) {
		// int n_read = 0;

		//Reading data in bulk if the buffer is emptied
		while(pl->fp && pl->fp->cur_read == pl->fp->num_read){
			ri_sig_close(pl->fp);
			if(pl->cur_f < pl->n_f){
				if((pl->fp = open_sig(pl->f[pl->cur_f++])) == 0) break;
				// ++n_read;
			}else if(pl->cur_fp < pl->n_fp){
				if(pl->f){
					for(int i = 0; i < pl->n_f; ++i) if(pl->f[i])free(pl->f[i]);
					free(pl->f);
				}
				pl->n_f = 0; pl->cur_f = 0;

				ri_char_v fnames = {0,0,0};
				kv_resize(char*, 0, fnames, 256);
				find_fast5(pl->fn[pl->cur_fp++], &fnames);
				pl->f =  fnames.a;
				if(!fnames.n || ((pl->fp = open_sig(pl->f[pl->cur_f++])) == 0)) break;
				pl->n_f = fnames.n;
				// ++n_read;
			}else {pl->fp = 0; break;}
		}

		if(!pl->fp || pl->fp->cur_read == pl->fp->num_read) break;

		// if (!n_read) break;
		
		ri_sig_t *s = (ri_sig_t*)calloc(1, sizeof(ri_sig_t));
		sigvec.push_back(s);
		ri_sig_read_sig(pl->fp, s);
		size += s->l_sig;

		if(size >= chunk_size) break;
	}

	ri_sig_t** a = 0;

	if(sigvec.size()){
		a = (ri_sig_t**)calloc(sigvec.size(), sizeof(ri_sig_t**));
		std::copy(sigvec.begin(), sigvec.end(), a);
	}
	
	*n_ = sigvec.size();
	return a;
}

static void *map_worker_pipeline(void *shared, int step, void *in)
{
	int i, k;
    pipeline_mt *p = (pipeline_mt*)shared;
    if (step == 0) { // step 0: read sequences
        step_mt *s;
        s = (step_mt*)calloc(1, sizeof(step_mt));
		s->sig = ri_sig_read_frag(p, p->mini_batch_size, &s->n_sig);
		if (s->n_sig && !p->su_stop) {
			s->p = p;
			for (i = 0; i < s->n_sig; ++i)
				(*(s->sig[i])).rid = p->n_processed++;
			s->buf = (ri_tbuf_t**)calloc(p->n_threads, sizeof(ri_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = ri_tbuf_init();
			s->reg = (ri_reg1_t**)calloc(s->n_sig, sizeof(ri_reg1_t*));
			for(i = 0; i < s->n_sig; ++i)
				s->reg[i] = (ri_reg1_t*)calloc(1, sizeof(ri_reg1_t));
			return s;
		} else if(s){
			if(s->sig) free(s->sig);
			free(s);
		}
    } else if (step == 1) { // step 1: detect events
		step_mt *s = (step_mt*)in;
		if(!p->su_stop) kt_for(p->n_threads, map_worker_for, in, s->n_sig);

		if(p->opt->flag & RI_M_SEQUENCEUNTIL && !p->su_stop){
			const ri_idx_t *ri = p->ri;
			for (k = 0; k < s->n_sig; ++k) {
				if(s->reg[k] && s->reg[k]->ref_id < ri->n_seq && s->reg[k]->read_name && s->reg[k]->mapped){
					ri_reg1_t* reg0 = s->reg[k];
					p->su_c_estimations[reg0->ref_id] += reg0->fragment_length;
					p->ab_count += reg0->fragment_length;
					p->su_nreads++;
					if(p->su_nreads > p->opt->tmin_reads && !(p->su_nreads%p->opt->ttest_freq)){
						//calculate abundance
						for(int ce = 0; ce < ri->n_seq; ++ce){
							p->su_estimations[p->su_cur][ce] =  (float)p->su_c_estimations[ce]/p->ab_count;
						}
						
						if(++p->su_cur >= p->opt->tn_samples) p->su_cur = 0;
						if(p->su_nestimations++ >= p->opt->tn_samples){
							if(find_outlier((const float**)p->su_estimations, ri->n_seq, p->opt->tn_samples) <= p->opt->t_threshold){
								//sending the stop signal.
								p->su_stop = k+1;
								fprintf(stderr, "[M::%s] Sequence Until is activated, stopping sequencing after processing %d mapped reads\n", __func__, p->su_nreads);
								break;
							}
						}
					}
				}
			}
		}
		return in;
    } else if (step == 2) { // step 2: output
		// void *km = 0;
        step_mt *s = (step_mt*)in;
		const ri_idx_t *ri = p->ri;
		for (k = 0; k < s->n_sig; ++k) {
			if(s->reg[k]){
				ri_reg1_t* reg0 = s->reg[k];
				char strand = reg0->rev?'-':'+';
				
				if(reg0->ref_id < ri->n_seq && reg0->read_name){
					if(reg0->mapped && (!p->su_stop || k < p->su_stop) ){
						fprintf(stdout, "%s\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%u\t%u\t%d\t%s\n", 
										reg0->read_name, reg0->read_length, reg0->read_start_position, reg0->read_end_position, 
										strand, ri->seq[reg0->ref_id].name, ri->seq[reg0->ref_id].len, reg0->fragment_start_position, reg0->fragment_start_position + reg0->fragment_length, reg0->read_end_position-reg0->read_start_position-1, reg0->fragment_length, reg0->mapq, reg0->tags);
					}else{
						fprintf(stdout, "%s\t%u\t*\t*\t*\t*\t*\t*\t*\t*\t*\t%d\t%s\n", reg0->read_name, reg0->read_length, reg0->mapq, reg0->tags);
					}
				}

				if(reg0->tags) free(reg0->tags);
				free(reg0);
			}
		}

		for (i = 0; i < p->n_threads; ++i) ri_tbuf_destroy(s->buf[i]);
		free(s->buf);

		if(s->reg)free(s->reg);

		for(int i = 0; i < s->n_sig; ++i){
			ri_sig_t *curS = s->sig[i];
			free(curS->sig);
			free(curS->name);
			free(curS);
		} if(s->sig)free(s->sig); free(s);
	}
    return 0;
}

int ri_map_file(const ri_idx_t *idx, const char *fn, const ri_mapopt_t *opt, int n_threads){
	return ri_map_file_frag(idx, 1, &fn, opt, n_threads);
}

int ri_map_file_frag(const ri_idx_t *idx, int n_segs, const char **fn, const ri_mapopt_t *opt, int n_threads){

	int i, pl_threads;
	pipeline_mt pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_mt));
	pl.n_fp = n_segs;
	pl.n_f = 0; pl.cur_f = 0;
	ri_char_v fnames = {0,0,0};
	kv_resize(char*, 0, fnames, 256);
	find_fast5(fn[0], &fnames);
	pl.f =  fnames.a;
	if(!fnames.n || ((pl.fp = open_sig(pl.f[0])) == 0)) return -1;
	if (pl.fp == 0) return -1;
	pl.fn = fn;
	pl.n_f = fnames.n;
	pl.cur_fp = 1;
	pl.cur_f = 1;
	pl.opt = opt, pl.ri = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	pl_threads = pl.n_threads == 1?1:2;
	pl.su_stop = 0;

	if(opt->flag & RI_M_SEQUENCEUNTIL){
		pl.su_nreads = 0;
		pl.su_nestimations = 0;
		pl.ab_count = 0;
		pl.su_cur = 0;
		pl.su_estimations = (float**)calloc(opt->tn_samples, sizeof(float*));
		for(int i = 0; i < opt->tn_samples; ++i) pl.su_estimations[i] = (float*)calloc(idx->n_seq, sizeof(float));

		pl.su_c_estimations = (uint32_t*)calloc(idx->n_seq, sizeof(uint32_t));
	}
	
	kt_pipeline(pl_threads, map_worker_pipeline, &pl, 3);

	if(opt->flag & RI_M_SEQUENCEUNTIL){
		// pl.su_nreads = 0;
		// pl.su_nestimations = 0;
		// pl.su_stop = 0;
		if(pl.su_estimations){
			for(int i = 0; i < opt->tn_samples; ++i) if(pl.su_estimations[i]) free(pl.su_estimations[i]);
			free(pl.su_estimations);
		}

		if(pl.su_c_estimations)free(pl.su_c_estimations);
	}

	return 0;
}