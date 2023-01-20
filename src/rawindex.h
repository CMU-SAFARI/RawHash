#ifndef RAWINDEX_H
#define RAWINDEX_H

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include "kthread.h"
#include "kvec.h"
#include "khash.h"
#include "kalloc.h"
#include "bseq.h"
#include <tuple>
#include <algorithm>

#include <sys/time.h>

#include <sys/resource.h>
#include <sys/time.h>

extern double ri_realtime0;
extern int ri_verbose;

double ri_realtime(void);

double ri_cputime(void);

long ri_peakrss(void);

#define RI_I_NAIVE		0x1
#define RI_I_MIN		0x2
#define RI_I_BLEND		0x4
#define RI_I_SYNCMER	0x8

#define RI_M_SEQUENCEUNTIL	0x1

#define RI_IDX_MAGIC   "RI"
#define RI_IDX_MAGIC_BYTE 2

//A large float value can be used for masking if needed
#define RI_MASK_SIGNAL 3.402823466e+32F

#define LAST_SIG_DIFF 0.3F

//To store sketches in vectors
#define RI_HASH_SHIFT 6
#define RI_ID_SHIFT 32
#define RI_POS_SHIFT 1

// extern double ri_realtime0;
// extern int ri_verbose;

#ifdef __cplusplus
extern "C" {
#endif

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;

typedef struct {
	uint32_t n;
	uint32_t q_pos;
	uint32_t q_span:31, flt:1;
	uint32_t seg_id:31, is_tandem:1;
	const uint64_t *cr;
} ri_seed_t;

typedef struct ri_idx_bucket_s {
	mm128_v a;   // (seed, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for seeds appearing >1 times
	void *h;     // hash table indexing _p_ and seeds appearing once
} ri_idx_bucket_t;

typedef struct {
	char *name;      // name of the db sequence
	uint64_t offset; // offset in ri_idx_t::S
	uint32_t len;    // length
	// uint32_t is_alt;
} ri_idx_seq_t;

typedef struct {
	int32_t b, w, e, n, q, lq, k, flag;
	int32_t index;
	struct ri_idx_bucket_s *B; // index (hidden)

	void *km, *h;

	uint32_t n_seq;
	ri_idx_seq_t *seq;
	// uint32_t *S;

} ri_idx_t;

// indexing and mapping options
typedef struct {
	short b, w, e, n, q, lq, k, flag;
	int64_t mini_batch_size;
	uint64_t batch_size;
} ri_idxopt_t;

typedef struct {

	//ONT Device specific parameters
	uint32_t bp_per_sec;
	uint32_t sample_rate;
	uint32_t chunk_size;

	//Chaining parameters
	uint32_t min_events;
	uint32_t max_gap_length;
	uint32_t max_target_gap_length;
	uint32_t chaining_band_length;
	uint32_t max_num_skips;
	uint32_t min_num_anchors;
	uint32_t num_best_chains;
	float min_chaining_score;

	//Mapping parameters
	uint32_t step_size;
	uint32_t max_num_chunk;
	uint32_t min_chain_anchor;
	uint32_t min_chain_anchor_out;

	float min_bestmap_ratio;
	float min_bestmap_ratio_out;

	float min_meanmap_ratio;
	float min_meanmap_ratio_out;

	float t_threshold;
	uint32_t tn_samples;
	uint32_t ttest_freq;
	uint32_t tmin_reads;
	
	int64_t flag;    // see ri_F_* macros
	// int seed;
	// int sdust_thres; // score threshold for SDUST; 0 to disable

	// int max_qlen;    // max query length

	// int bw, bw_long; // bandwidth
	// int max_gap, max_gap_ref; // break a chain if there are no minimizers in a max_gap window
	// int max_frag_len;
	// int max_chain_skip, max_chain_iter;
	// int min_cnt;         // min number of minimizers on each chain
	// int min_chain_score; // min chaining score
	// float chain_gap_scale;
	// float chain_skip_scale;
	// int rmq_size_cap, rmq_inner_dist;
	// int rmq_rescue_size;
	// float rmq_rescue_ratio;

	// int best_n;      // top best_n chains are used

	// int rank_min_len;
	// float rank_frac;

	// int pe_ori, pe_bonus;

	// float mid_occ_frac;  // only used by ri_mapopt_update(); see below
	// float q_occ_frac;
	// int32_t min_mid_occ, max_mid_occ;
	// int32_t mid_occ;     // ignore seeds with occurrences above this threshold
	// int32_t max_occ, max_max_occ, occ_dist;
	int64_t mini_batch_size; // size of a batch of query bases to process in parallel
	// int64_t max_sw_mat;
	// int64_t cap_kalloc;

	// const char *split_prefix;


	//Event detector options
	uint32_t window_length1;
	uint32_t window_length2;
	float threshold1;
	float threshold2;
	float peak_height;
} ri_mapopt_t;

// index reader
typedef struct {
	int is_idx, n_parts;
	int64_t idx_size;
	ri_idxopt_t opt;
	FILE *fp_out;
	union {
		struct mm_bseq_file_s *seq;
		FILE *idx;
	} fp;
} ri_idx_reader_t;


// memory buffer for thread-local storage during mapping
typedef struct ri_tbuf_s ri_tbuf_t;

// double ri_realtime(void);

// double ri_cputime(void);

// long ri_peakrss(void);

KHASH_MAP_INIT_STR(str, uint32_t)
#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);

/*************
 * options    *
 *************/

void ri_idxopt_init(ri_idxopt_t *opt);

void ri_mapopt_init(ri_mapopt_t *opt);

/*************
 * index     *
 *************/

void ri_idx_stat(const ri_idx_t *ri);

/**
 * Initialize an index reader
 *
 * @param fn         index or fasta/fastq file name (this function tests the file type)
 * @param opt        indexing parameters
 * @param fn_out     if not NULL, write built index to this file
 *
 * @return an index reader on success; NULL if fail to open _fn_
 */
ri_idx_reader_t* ri_idx_reader_open(const char *fn, const ri_idxopt_t *ipt, const char *fn_out);

/**
 * Destroy/deallocate an index reader
 *
 * @param r          index reader
 */
void ri_idx_reader_close(ri_idx_reader_t* r);

int ri_idx_reader_eof(const ri_idx_reader_t* r);

/**
 * Check whether the file contains a rawindex index
 *
 * @param fn         file name
 *
 * @return the file size if fn is an index file; 0 if fn is not.
 */
int64_t ri_idx_is_idx(const char* fn);

/// @brief 
/// @param hashbits 
/// @param neighbors 
/// @return 
ri_idx_t* ri_idx_init(int b, int w, int e, int n, int q, int lq, int k);

ri_idx_t *ri_idx_reader_read(ri_idx_reader_t *r, float* pore_vals, int n_threads);

void ri_idx_add(ri_idx_t* ri, int n, const mm128_t* a);

void ri_idx_dump(FILE* idx_file, const ri_idx_t* ri);

ri_idx_t* ri_idx_load(FILE* fp);

/**
 * Destroy/deallocate an index
 *
 * @param r          minimap2 index
 */
void ri_idx_destroy(ri_idx_t* ri);

void ri_idx_sort(ri_idx_t* ri, int n_threads);

const uint64_t *ri_idx_get(const ri_idx_t *ri, uint64_t hashval, int *n);

/**************************
 * Sequence to Signal     *
 **************************/

void ri_seq_to_sig(const char *str, int len, const float* pore_vals, const int pore_kmer, const int strand, uint32_t* s_len, float* s_values);

/*************
 * sketch     *
 *************/

/**
 * Generate sketches (seeds) from signal values
 *
 * @param km       thread-local memory pool; using NULL falls back to malloc()
 * @param s_values signal values (float)
 * @param id       ID of the signals (e.g., ref ID); will be copied to the output $p array
 * @param strand   strand of signal values. Either 1 (forw) or 0 (reverse); will be copied to the output $p array
 * @param len      length of $s_values
 * @param w        If w>0, find a minimizer for every $w consecutive k-mers
 * @param e        Number of packed events in a single 32/64-bit hash value (i.e., creates the e-dimension of events as in Sigmap)
 * @param n        If n>0 Enables BLEND. Number of neighbors hash values to use with SimHash.
 * 				   Generates a single hash value from n hash values; will be stored as a hash value if enabled.
 * 				   Please see the BLEND paper for details of this hash value calculation.
 * @param p        seeds/sketches
 *                 p->a[i].x = hash value stored in a hash table
 *                 p->a[i].y = id<<32 | lastPos<<1 | strand
 *                 where lastPos is the position of the i-th signal,
 *                 and strand indicates whether the seed comes from the top or the bottom strand.
 *                 Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void ri_sketch(void *km, const float* s_values, uint32_t id, int strand, int len, int w, int e, int n, int q, int lq, int k, mm128_v *p);

void ri_seed_mz_flt(void *km, mm128_v *mv, int32_t q_occ_max, float q_occ_frac);

// ri_seed_t *ri_collect_matches(void *km, int *_n_m, int qlen, int max_occ, int max_max_occ, int dist, const ri_idx_t *mi, const mm128_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos, uint64_t **mini_pos);


/*************
 * mapping   *
 *************/

// struct ri_paf_s {
//   uint32_t read_id;
//   const char* read_name;
//   uint32_t read_length;
//   uint32_t read_start_position;
//   uint32_t read_end_position;
//   uint32_t fragment_start_position;
//   uint32_t fragment_length;
//   uint8_t mapq : 6, direction : 1, is_unique : 1;
//   const char* tags;
//   bool operator<(const ri_paf_s &m) const {
//     return std::tie(fragment_start_position, fragment_length, mapq, direction,
//                     is_unique, read_id, read_start_position,
//                     read_end_position) <
//            std::tie(m.fragment_start_position, m.fragment_length, m.mapq,
//                     m.direction, m.is_unique, m.read_id, m.read_start_position,
//                     m.read_end_position);
//   }
//   bool operator==(const ri_paf_s &m) const {
//     return std::tie(fragment_start_position) ==
//            std::tie(m.fragment_start_position);
//   }
// };

// typedef struct ri_paf_s ri_paf_t;

typedef struct ri_anchor_s{
  uint32_t target_position;
  uint32_t query_position;
  bool operator<(const ri_anchor_s &a) const {
    return std::tie(target_position, query_position) < std::tie(a.target_position, a.query_position);
  }
};

typedef struct ri_anchor_s ri_anchor_t;

typedef struct ri_chain_s{
  float score;
  uint32_t reference_sequence_index;
  uint32_t start_position;
  uint32_t end_position;
  uint32_t n_anchors;
  uint8_t mapq;
  int strand;
//   std::vector<ri_anchor_t> anchors;
  ri_anchor_t* anchors;
  bool operator>(const ri_chain_s &b) const {
    return std::tie(score, n_anchors, strand, reference_sequence_index, start_position, end_position) >
           std::tie(b.score, b.n_anchors, b.strand, b.reference_sequence_index, b.start_position, b.end_position);
  }
};

typedef struct ri_chain_s ri_chain_t;

typedef struct {
	// int32_t id;             // ID for internal uses (see also parent below)
	// int32_t cnt;            // number of minimizers; if on the reverse strand
	// int32_t rid;            // reference index; if this is an alignment from inversion rescue
	// int32_t score;          // DP alignment score
	// int32_t qs, qe, rs, re; // query start and end; reference start and end
	// int32_t parent, subsc;  // parent==id if primary; best alternate mapping score
	// int32_t as;             // offset in the a[] array (for internal uses only)
	// int32_t mlen, blen;     // seeded exact match length; seeded alignment block length
	// int32_t n_sub;          // number of suboptimal mappings
	// int32_t score0;         // initial chaining score (before chain merging/spliting)
	// uint32_t mapq:8, split:2, rev:1, inv:1, sam_pri:1, proper_frag:1, pe_thru:1, seg_split:1, seg_id:8, split_inv:1, is_alt:1, strand_retained:1, dummy:5;
	// uint32_t hash;
	// float div;

	uint32_t read_id;
	uint32_t ref_id;
	const char* read_name;
	uint32_t read_length;
	uint32_t read_start_position;
	uint32_t read_end_position;
	uint32_t fragment_start_position;
	uint32_t fragment_length;
	uint8_t mapq : 6, rev : 1, mapped : 1;
	char* tags;

	uint32_t offset;

	ri_chain_s* chains;
	uint32_t n_chains;
} ri_reg1_t;

// #include "kseq.h"
// #include <zlib.h>
// #include <slow5/slow5.h>
// #include "hdf5_tools.hpp"

// typedef struct {
// 	int l_seq, rid;
// 	char *name, *seq, *qual, *comment;
// } ri_sig_t;

// struct ri_sig_file_s {
// 	// gzFile fp;
// 	// kseq_t *ks;
// 	ri_sig_t s;
// 	hdf5_tools::File* fp;
// };

// struct ri_sig_file_s;
// typedef struct ri_sig_file_s ri_sig_file_t;

void ri_map_frag(const ri_idx_t *ri, const uint32_t s_len, const float *sig, ri_reg1_t* reg, ri_tbuf_t *b, const ri_mapopt_t *opt, const char *qname);

int ri_map_file(const ri_idx_t *idx, const char *fn, const ri_mapopt_t *opt, int n_threads);

int ri_map_file_frag(const ri_idx_t *idx, int n_segs, const char **fn, const ri_mapopt_t *opt, int n_threads);

#ifdef __cplusplus
}
#endif

#endif // RAWINDEX_H
