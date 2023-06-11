#ifndef RMAP_H
#define RMAP_H

#include "rawindex.h"
#include "rsig.h"
#include <tuple>
#include <vector>
#include <string>
#include <cassert>

#ifdef __cplusplus
extern "C" {
#endif

// typedef struct ri_seed_s{
// 	uint32_t n;
// 	uint32_t q_pos;
// 	uint32_t q_span:31, flt:1;
// 	uint32_t seg_id:31, is_tandem:1;
// 	const uint64_t *cr;
// } ri_seed_t;

typedef struct ri_anchor_s{
  uint32_t target_position;
  uint32_t query_position;
  bool operator<(const ri_anchor_s &a) const {
    return std::tie(target_position, query_position) < std::tie(a.target_position, a.query_position);
  }
} ri_anchor_t;

typedef struct ri_chain_s{
  float score;
  uint32_t reference_sequence_index;
  uint32_t start_position;
  uint32_t end_position;
  uint32_t n_anchors;
  uint8_t mapq;
  int strand;
  ri_anchor_t* anchors;
  bool operator>(const ri_chain_s &b) const {
    return std::tie(score, n_anchors, strand, reference_sequence_index, start_position, end_position) >
           std::tie(b.score, b.n_anchors, b.strand, b.reference_sequence_index, b.start_position, b.end_position);
  }
}ri_chain_t;

typedef struct ri_reg1_s{
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

typedef struct pipeline_ms{
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

// memory buffer for thread-local storage during mapping
typedef struct ri_tbuf_s {
	void *km;
	int rep_len, frag_gap;
}ri_tbuf_t;

typedef struct step_ms{
	const pipeline_mt *p;
    int n_sig;
	ri_sig_t** sig;
	// int* n_reg, *seg_off, *n_seg, *rep_len, *frag_gap;
	ri_reg1_t** reg;
	ri_tbuf_t** buf;
} step_mt;

/**
 * Map raw nanopore signals of a single read to a reference genome
 *
 * @param idx		rawindex (see rawindex.h)
 * @param fn		path to the signal files
 * @param opt		mapping options
 * @param n_threads	number of threads to use in mapping
 * 
 * @return		returns 0 if mapping is completed with no issues. -1, otherwise.
 */
int ri_map_file(const ri_idx_t *idx, const char *fn, const ri_mapopt_t *opt, int n_threads);

/**
 * Map raw nanopore signals of many reads to a reference genome
 *
 * @param idx		rawindex (see rawindex.h)
 * @param n_segs	number of signal files
 * @param fn		paths to the signal files
 * @param opt		mapping options
 * @param n_threads	number of threads to use in mapping
 * 
 * @return			returns 0 if mapping is completed with no issues. -1, otherwise.
 */
int ri_map_file_frag(const ri_idx_t *idx, int n_segs, const char **fn, const ri_mapopt_t *opt, int n_threads);

#ifdef __cplusplus
}
#endif
#endif //RMAP_H