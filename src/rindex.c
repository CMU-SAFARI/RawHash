#include "rindex.h"
#include <assert.h>
#include <fcntl.h>
#include "rsketch.h"
#include "bseq.h"
#include "khash.h"
#include "rh_kvec.h"
#include "kthread.h"
#include "revent.h"

#if defined(WIN32) || defined(_WIN32)
#include <io.h> // for open(2)
#else
#include <unistd.h>
#endif

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

KHASH_MAP_INIT_STR(str, uint32_t)

#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#define mm_seq4_set(s, i, c) ((s)[(i)>>3] |= (uint32_t)(c) << (((i)&7)<<2))

double ri_realtime0;
int ri_verbose = 1;

typedef struct {
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	mm_bseq_file_t* fp;
	ri_idx_t* ri;

	int n_f, cur_f; //number of files, index of the current file (in sf)
	ri_sig_file_t *sfp; //pointer to the current fast5 file
	char **sf; //fast5 file names (multi if directory is given, single if single file is given)
} pipeline_t;

typedef struct {
    int n_seq;
	mm_bseq1_t* seq;
	ri_sig_t** sig;
	mm128_v a;
} step_t;

void ri_idx_stat(const ri_idx_t *ri)
{
	fprintf(stderr, "[M::%s] pore kmer size: %d; concatanated events: %d; quantization bits: %d; w: %d; n: %d; #seq: %d\n", __func__, ri->k, ri->e, ri->q, ri->w, ri->n, ri->n_seq);
}

ri_idx_t* ri_idx_init(float diff, int b, int w, int e, int n, int q, int k, float fine_min, float fine_max, float fine_range, int flag){
	ri_idx_t* ri;
	ri = (ri_idx_t*)calloc(1, sizeof(ri_idx_t));
  	ri->b = b, ri->w = w; ri->e = e; ri->n = n; ri->q = q; ri->k = k, ri->flag = flag;
	ri->diff = diff;
	ri->fine_min = fine_min, ri->fine_max = fine_max, ri->fine_range = fine_range;
	ri->seq = NULL; ri->sig = NULL; ri->F = NULL; ri->R = NULL; ri->f_l_sig = NULL; ri->r_l_sig = NULL; ri->pore = NULL; ri->h = NULL;
  	ri->B = (ri_idx_bucket_t*)calloc(1<<ri->b, sizeof(ri_idx_bucket_t));
  	ri->km = ri_km_init();

	return ri;
}

void ri_idx_destroy(ri_idx_t *ri){
	uint32_t i;
	if (!ri) return;
	if (ri->h) kh_destroy(str, (khash_t(str)*)ri->h);
	if (ri->B) {
		for (i = 0; i < 1U<<ri->b; ++i) {
			free(ri->B[i].p);
			free(ri->B[i].a.a);
			kh_destroy(idx, (idxhash_t*)ri->B[i].h);
		}
	}

	if(ri->km) ri_km_destroy(ri->km);
	free(ri->B); free(ri);
}

void ri_idx_add(ri_idx_t *ri, int n, const mm128_t *a){
	int i, mask = (1<<ri->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &ri->B[a[i].x>>RI_HASH_SHIFT&mask].a;
		rh_kv_push(mm128_t, 0, *p, a[i]);
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

		int seq_c = s->seq[s->n_seq-1].rid + 1;

		if(p->ri->flag&RI_I_STORE_SIG){
			p->ri->F = (float**)ri_krealloc(p->ri->km, p->ri->F,seq_c * sizeof(float*));
			p->ri->f_l_sig = (uint32_t*)ri_krealloc(p->ri->km, p->ri->f_l_sig,seq_c * sizeof(uint32_t));

			if(!(p->ri->flag&RI_I_REV_QUERY)){
				p->ri->R = (float**)ri_krealloc(p->ri->km, p->ri->R,seq_c * sizeof(float*));
				p->ri->r_l_sig = (uint32_t*)ri_krealloc(p->ri->km, p->ri->r_l_sig,seq_c * sizeof(uint32_t));
			}

			for (i = 0; i < s->n_seq; ++i) {
				mm_bseq1_t* t = &s->seq[i];
				if (t->l_seq > 0){
					uint32_t s_len;
					uint32_t r_id = s->seq[i].rid;
					p->ri->F[r_id] = (float*)ri_kcalloc(p->ri->km, t->l_seq, sizeof(float));
					float* s_values = p->ri->F[r_id];

					ri_seq_to_sig(t->seq, t->l_seq, p->ri->pore, p->ri->k, 0, &s_len, s_values);
					ri_sketch(0, s_values, t->rid, 0, s_len, p->ri->diff, p->ri->w, p->ri->e, p->ri->n, p->ri->q, p->ri->k, p->ri->fine_min, p->ri->fine_max, p->ri->fine_range, &s->a);
					p->ri->f_l_sig[r_id] = s_len;

					if(!(p->ri->flag&RI_I_REV_QUERY)){
						p->ri->R[r_id] = (float*)ri_kcalloc(p->ri->km, t->l_seq, sizeof(float));
						s_values = p->ri->R[r_id];
						ri_seq_to_sig(t->seq, t->l_seq, p->ri->pore, p->ri->k, 1, &s_len, s_values);
						ri_sketch(0, s_values, t->rid, 1, s_len, p->ri->diff, p->ri->w, p->ri->e, p->ri->n, p->ri->q, p->ri->k, p->ri->fine_min, p->ri->fine_max, p->ri->fine_range, &s->a);
						p->ri->r_l_sig[r_id] = s_len;
					}
				}
				free(t->seq); free(t->name);
			}
		}else{
			for (i = 0; i < s->n_seq; ++i) {
				mm_bseq1_t* t = &s->seq[i];
				if (t->l_seq > 0){
					uint32_t s_len;
					float* s_values = (float*)calloc(t->l_seq, sizeof(float));

					ri_seq_to_sig(t->seq, t->l_seq, p->ri->pore, p->ri->k, 0, &s_len, s_values);
					ri_sketch(0, s_values, t->rid, 0, s_len, p->ri->diff, p->ri->w, p->ri->e, p->ri->n, p->ri->q, p->ri->k, p->ri->fine_min, p->ri->fine_max, p->ri->fine_range, &s->a);

					if(!(p->ri->flag&RI_I_REV_QUERY)){
						ri_seq_to_sig(t->seq, t->l_seq, p->ri->pore, p->ri->k, 1, &s_len, s_values);
						ri_sketch(0, s_values, t->rid, 1, s_len, p->ri->diff, p->ri->w, p->ri->e, p->ri->n, p->ri->q, p->ri->k, p->ri->fine_min, p->ri->fine_max, p->ri->fine_range, &s->a);
					}

					free(s_values);
				}
				free(t->seq); free(t->name); 
			}
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

ri_sig_t** ri_sig_read(pipeline_t *pl,
					   int64_t chunk_size,
					   int *n_)
{
	int64_t size = 0;
	rhsig_v rsigv = {0,0,0};
	
	*n_ = 0;
	if (pl->n_f < 1) return 0;

	//Debugging for sweeping purposes
	// if(pl->su_nreads >= 2000) return 0;

	rh_kv_resize(ri_sig_t*, 0, rsigv, 4000);

	while (pl->sfp) {
		//Reading data in bulk if the buffer is emptied
		while(pl->sfp && pl->sfp->cur_read == pl->sfp->num_read){
			ri_sig_close(pl->sfp);
			if(pl->cur_f < pl->n_f){
				if((pl->sfp = open_sig(pl->sf[pl->cur_f++])) == 0) break;
			}else {pl->sfp = 0; break;}
		}

		if(!pl->sfp || pl->sfp->cur_read == pl->sfp->num_read) break;
		
		ri_sig_t *s = (ri_sig_t*)calloc(1, sizeof(ri_sig_t));
		rh_kv_push(ri_sig_t*, 0, rsigv, s);
		ri_read_sig(pl->sfp, s);
		size += s->l_sig;

		if(size >= chunk_size) break;
		
		//Debugging for sweeping purposes
		// pl->su_nreads++;
		// if(pl->su_nreads >= 2000) break;
	}

	ri_sig_t** a = 0;
	if(rsigv.n) a = rsigv.a;
	*n_ = rsigv.n;

	return a;
}

static void *worker_sig_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));

		s->sig = ri_sig_read(p, p->mini_batch_size, &s->n_seq); // read a mini-batch

		if (s->sig) {
			uint32_t old_m, m;
			assert((uint64_t)p->ri->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->ri->seq
			old_m = p->ri->n_seq, m = p->ri->n_seq + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m)
				p->ri->sig = (ri_sig_t*)ri_krealloc(p->ri->km, p->ri->sig, m*sizeof(ri_sig_t));
			for (i = 0; i < s->n_seq; ++i) {
				ri_sig_t *sig = &p->ri->sig[p->ri->n_seq];
				sig->name = (char*)ri_kmalloc(p->ri->km, strlen((*(s->sig[i])).name) + 1);
				strcpy(sig->name, (*(s->sig[i])).name);
				sig->l_sig = (*(s->sig[i])).l_sig;
				sig->offset = p->sum_len;
				p->sum_len += sig->l_sig;
				(*(s->sig[i])).rid = p->ri->n_seq++;
			}
		}else {
			free(s); s = 0;
		}
		return s;
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			ri_sig_t* t = s->sig[i];
			if (t->l_sig > 0){
				uint32_t s_len = 0;
				double s_sum = 0, s_std = 0;
				uint32_t n_events_sum = 0;
				float* s_values = detect_events(0, t->l_sig, t->sig, p->ri->window_length1, p->ri->window_length2, p->ri->threshold1, p->ri->threshold2, p->ri->peak_height, &s_sum, &s_std, &n_events_sum, &s_len);

				ri_sketch(0, s_values, t->rid, 0, s_len, p->ri->diff, p->ri->w, p->ri->e, p->ri->n, p->ri->q, p->ri->k, p->ri->fine_min, p->ri->fine_max, p->ri->fine_range, &s->a);

				if(s_values)free(s_values);
			}
			free(t->sig); free(t->name); 
		}
		free(s->sig); s->sig = 0;
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
        step_t *s = (step_t*)in;
		ri_idx_add(p->ri, s->a.n, s->a.a);
		ri_kfree(0, s->a.a); free(s);
	}
    return 0;
}

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

	uint32_t pars[7], i;

	pars[0] = ri->w, pars[1] = ri->e, pars[2] = ri->n, pars[3] = ri->q, pars[4] = ri->k, pars[5] = ri->n_seq, pars[6] = ri->flag;
	
	fwrite(RI_IDX_MAGIC, 1, RI_IDX_MAGIC_BYTE, idx_file);
	fwrite(pars, sizeof(uint32_t), 7, idx_file);
	fwrite(&ri->diff, sizeof(float), 1, idx_file);
	fwrite(&ri->fine_min, sizeof(float), 1, idx_file);
	fwrite(&ri->fine_max, sizeof(float), 1, idx_file);
	fwrite(&ri->fine_range, sizeof(float), 1, idx_file);
	fwrite(ri->pore, sizeof(ri_pore_t), 1, idx_file);
	fwrite(ri->pore->pore_vals, sizeof(float), ri->pore->n_pore_vals, idx_file);
	ri_kfree(ri->km, ri->pore->pore_vals); 
	ri->pore->pore_vals = 0;
	fwrite(ri->pore->pore_inds, sizeof(ri_porei_t), ri->pore->n_pore_vals, idx_file);
	ri_kfree(ri->km, ri->pore->pore_inds); 
	ri->pore->pore_inds = 0;
	ri_kfree(ri->km, ri->pore);
	// ri->pore = 0;

	for (i = 0; i < ri->n_seq; ++i) {
		if(ri->flag & RI_I_SIG_TARGET){
			if (ri->sig[i].name) {
				uint8_t l = strlen(ri->sig[i].name);
				fwrite(&l, 1, 1, idx_file);
				fwrite(ri->sig[i].name, 1, l, idx_file);
			} else {
				uint8_t l = 0;
				fwrite(&l, 1, 1, idx_file);
			}
			fwrite(&ri->sig[i].l_sig, 4, 1, idx_file);
		}else {
			if (ri->seq[i].name) {
				uint8_t l = strlen(ri->seq[i].name);
				fwrite(&l, 1, 1, idx_file);
				fwrite(ri->seq[i].name, 1, l, idx_file);
			} else {
				uint8_t l = 0;
				fwrite(&l, 1, 1, idx_file);
			}
			fwrite(&ri->seq[i].len, 4, 1, idx_file);
		}

		if(ri->flag & RI_I_STORE_SIG){
			fwrite(&(ri->f_l_sig[i]), 4, 1, idx_file);
			fwrite(ri->F[i], 4, ri->f_l_sig[i], idx_file);
			ri_kfree(ri->km, ri->F[i]);
			if(!(ri->flag&RI_I_REV_QUERY)){
				fwrite(&(ri->r_l_sig[i]), 4, 1, idx_file);
				fwrite(ri->R[i], 4, ri->r_l_sig[i], idx_file);
				ri_kfree(ri->km, ri->R[i]);
			}
		}
	}

	if(ri->flag & RI_I_STORE_SIG){
		ri_kfree(ri->km, ri->F);
		ri_kfree(ri->km, ri->f_l_sig);
		if(!(ri->flag&RI_I_REV_QUERY)){
			ri_kfree(ri->km, ri->R);
			ri_kfree(ri->km, ri->r_l_sig);
		}
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
	int pars[7];
	fread(&pars[0], sizeof(int), 7, idx_file);

	float diff, fine_min, fine_max, fine_range;
	fread(&diff, sizeof(float), 1, idx_file);
	fread(&fine_min, sizeof(float), 1, idx_file);
	fread(&fine_max, sizeof(float), 1, idx_file);
	fread(&fine_range, sizeof(float), 1, idx_file);

	ri = ri_idx_init(diff, 14, pars[0], pars[1], pars[2], pars[3], pars[4], fine_min, fine_max, fine_range, pars[6]);
	ri->n_seq = pars[5];
	if(ri->flag&RI_I_SIG_TARGET) ri->sig = (ri_sig_t*)ri_kcalloc(ri->km, ri->n_seq, sizeof(ri_sig_t));
	else ri->seq = (ri_idx_seq_t*)ri_kcalloc(ri->km, ri->n_seq, sizeof(ri_idx_seq_t));

	ri->pore = (ri_pore_t*)ri_kmalloc(ri->km, sizeof(ri_pore_t));
	fread(ri->pore, sizeof(ri_pore_t), 1, idx_file);
	ri->pore->pore_vals = (float*)ri_kmalloc(ri->km, ri->pore->n_pore_vals * sizeof(float));
	fread(ri->pore->pore_vals, sizeof(float), ri->pore->n_pore_vals, idx_file);
	ri->pore->pore_inds = (ri_porei_t*)ri_kmalloc(ri->km, ri->pore->n_pore_vals * sizeof(ri_porei_t));
	fread(ri->pore->pore_inds, sizeof(ri_porei_t), ri->pore->n_pore_vals, idx_file);

	if(ri->flag & RI_I_STORE_SIG){
		ri->F = (float**)ri_kcalloc(ri->km, ri->n_seq, sizeof(float*));
		ri->f_l_sig = (uint32_t*)ri_kcalloc(ri->km, ri->n_seq, sizeof(uint32_t));
		if(!(ri->flag&RI_I_REV_QUERY)){
			ri->R = (float**)ri_kcalloc(ri->km, ri->n_seq, sizeof(float*));
			ri->r_l_sig = (uint32_t*)ri_kcalloc(ri->km, ri->n_seq, sizeof(uint32_t));
		}
	}

	for (i = 0; i < ri->n_seq; ++i) {
		uint8_t l;
		if(ri->flag&RI_I_SIG_TARGET){
			ri_sig_t *s = &ri->sig[i];
			fread(&l, 1, 1, idx_file);
			if (l) {
				s->name = (char*)ri_kmalloc(ri->km, l + 1);
				fread(s->name, 1, l, idx_file);
				s->name[l] = 0;
			}
			fread(&s->l_sig, 4, 1, idx_file);
			s->offset = sum_len;
			sum_len += s->l_sig;
		}else {
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

		if(ri->flag & RI_I_STORE_SIG){
			fread(&(ri->f_l_sig[i]), 4, 1, idx_file);
			ri->F[i] = (float*)ri_kmalloc(ri->km, ri->f_l_sig[i] * sizeof(float));
			fread(ri->F[i], 4, ri->f_l_sig[i], idx_file);

			if(!(ri->flag&RI_I_REV_QUERY)){
				fread(&(ri->r_l_sig[i]), 4, 1, idx_file);
				ri->R[i] = (float*)ri_kmalloc(ri->km, ri->r_l_sig[i] * sizeof(float));
				fread(ri->R[i], 4, ri->r_l_sig[i], idx_file);
			}
		}
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
	} else if(r->opt.flag & RI_I_SIG_TARGET) {
		r->n_f = 0; r->cur_f = 0;
		ri_char_v fnames = {0,0,0};
		rh_kv_resize(char*, 0, fnames, 256);
		find_sfiles(fn, &fnames);
		r->sf =  fnames.a;
		if(!fnames.n || ((r->sfp = open_sig(r->sf[0])) == 0)) return 0;
		if (r->sfp == 0) return 0;
		r->n_f = fnames.n;
		r->cur_f = 1;
	}
	else r->fp.seq = mm_bseq_open(fn);
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
	return r;
}

void ri_idx_reader_close(ri_idx_reader_t* r){

	if (r->is_idx) fclose(r->fp.idx);
	else if(r->opt.flag & RI_I_SIG_TARGET && r->sfp) {
		if (r->sfp) ri_sig_close(r->sfp);
		if (r->sf) free(r->sf);
	} 
	else if(r->fp.seq) mm_bseq_close(r->fp.seq);
	else if(r->fp.seq) mm_bseq_close(r->fp.seq);
	if (r->fp_out) fclose(r->fp_out);
	free(r);
}

ri_idx_t* ri_idx_gen(mm_bseq_file_t* fp, ri_pore_t* pore, float diff, int b, int w, int e, int n, int q, int k, float fine_min, float fine_max, float fine_range, int flag, int mini_batch_size, int n_threads, uint64_t batch_size)
{

	if(flag&RI_I_SIG_TARGET) return 0;

	pipeline_t pl;
	if (fp == 0 || mm_bseq_eof(fp) || !pore || !pore->pore_vals) return 0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = (uint64_t)mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.ri = ri_idx_init(diff, b, w, e, n, q, k, fine_min, fine_max, fine_range, flag);
	
	pl.ri->pore = (ri_pore_t*)ri_kmalloc(pl.ri->km, sizeof(ri_pore_t));
	memcpy(pl.ri->pore, pore, sizeof(ri_pore_t));
	pl.ri->pore->pore_vals = (float*)ri_kmalloc(pl.ri->km, pore->n_pore_vals * sizeof(float));
	memcpy(pl.ri->pore->pore_vals, pore->pore_vals, pore->n_pore_vals * sizeof(float));
	pl.ri->pore->pore_inds = (ri_porei_t*)ri_kmalloc(pl.ri->km, pore->n_pore_vals * sizeof(ri_porei_t));
	memcpy(pl.ri->pore->pore_inds, pore->pore_inds, pore->n_pore_vals * sizeof(ri_porei_t));

	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	ri_idx_post(pl.ri, n_threads);

	return pl.ri;
}

ri_idx_t* ri_idx_siggen(ri_sig_file_t** fp, char **f, int &cur_f, int n_f, ri_pore_t* pore, float diff, int b, int w, int e, int n, int q, int k, float fine_min, float fine_max, float fine_range, uint32_t window_length1, uint32_t window_length2, float threshold1, float threshold2, float peak_height, int flag, int mini_batch_size, int n_threads, uint64_t batch_size)
{

	if(!(flag&RI_I_SIG_TARGET)) return 0;

	pipeline_t pl;
	if (fp == 0 || *fp == 0 || n_f <= 0 || cur_f > n_f || (cur_f == n_f && (*fp)->cur_read >= (*fp)->num_read) || !pore || !pore->pore_vals) return 0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = (uint64_t)mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.sfp = *fp;
	pl.sf = f;
	pl.n_f = n_f;
	pl.cur_f = cur_f;
	pl.ri = ri_idx_init(diff, b, w, e, n, q, k, fine_min, fine_max, fine_range, flag);

	pl.ri->pore = (ri_pore_t*)ri_kmalloc(pl.ri->km, sizeof(ri_pore_t));
	memcpy(pl.ri->pore, pore, sizeof(ri_pore_t));
	pl.ri->pore->pore_vals = (float*)ri_kmalloc(pl.ri->km, pore->n_pore_vals * sizeof(float));
	memcpy(pl.ri->pore->pore_vals, pore->pore_vals, pore->n_pore_vals * sizeof(float));
	pl.ri->pore->pore_inds = (ri_porei_t*)ri_kmalloc(pl.ri->km, pore->n_pore_vals * sizeof(ri_porei_t));
	memcpy(pl.ri->pore->pore_inds, pore->pore_inds, pore->n_pore_vals * sizeof(ri_porei_t));

	pl.ri->window_length1 = window_length1;
	pl.ri->window_length2 = window_length2;
	pl.ri->threshold1 = threshold1;
	pl.ri->threshold2 = threshold2;
	pl.ri->peak_height = peak_height;

	kt_pipeline(n_threads < 3? n_threads : 3, worker_sig_pipeline, &pl, 3);

	*fp = pl.sfp;
	cur_f = pl.cur_f;

	ri_idx_post(pl.ri, n_threads);

	return pl.ri;
}

ri_idx_t* ri_idx_reader_read(ri_idx_reader_t* r, ri_pore_t* pore, int n_threads){

	ri_idx_t *ri;
	if (r->is_idx) {
		ri = ri_idx_load(r->fp.idx);
	} else if(r->opt.flag&RI_I_SIG_TARGET) {
		ri = ri_idx_siggen(&(r->sfp), r->sf, r->cur_f, r->n_f, pore, r->opt.diff, r->opt.b, r->opt.w, r->opt.e, r->opt.n, r->opt.q, r->opt.k, r->opt.fine_min, r->opt.fine_max, r->opt.fine_range, r->opt.window_length1, r->opt.window_length2, r->opt.threshold1, r->opt.threshold2, r->opt.peak_height, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size);
	} else{
		ri = ri_idx_gen(r->fp.seq, pore, r->opt.diff, r->opt.b, r->opt.w, r->opt.e, r->opt.n, r->opt.q, r->opt.k, r->opt.fine_min, r->opt.fine_max, r->opt.fine_range, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size);
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

int32_t ri_idx_cal_max_occ(const ri_idx_t *ri, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return INT32_MAX;
	for (i = 0; i < 1<<ri->b; ++i)
		if (ri->B[i].h) n += kh_size((idxhash_t*)ri->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<ri->b; ++i) {
		idxhash_t *h = (idxhash_t*)ri->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

void ri_mapopt_update(ri_mapopt_t *opt, const ri_idx_t *ri)
{
	if (opt->mid_occ <= 0) {
		opt->mid_occ = ri_idx_cal_max_occ(ri, opt->mid_occ_frac);
		if (opt->mid_occ < opt->min_mid_occ)
			opt->mid_occ = opt->min_mid_occ;
		if (opt->max_mid_occ > opt->min_mid_occ && opt->mid_occ > opt->max_mid_occ)
			opt->mid_occ = opt->max_mid_occ;

		fprintf(stderr, "[M::%s::%.6f*%.6f] mid_occ = %d, min_mid_occ = %d, max_mid_occ = %d\n", __func__,
				opt->mid_occ_frac, opt->q_occ_frac, opt->mid_occ, opt->min_mid_occ, opt->max_mid_occ);
	}
	if (opt->bw_long < opt->bw) opt->bw_long = opt->bw;
}