#ifndef RMAP_H
#define RMAP_H

#include "rindex.h"
#include "rsig.h"
#include "chain.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ri_map_s{
	uint32_t c_id; //chain index
	uint32_t read_length;
	uint32_t ref_id;
	uint32_t read_start_position;
	uint32_t read_end_position;
	uint32_t fragment_start_position;
	uint32_t fragment_length;
	uint8_t mapq : 6, rev : 1, mapped : 1;
	char* tags;
} ri_map_t;

typedef struct ri_reg1_s{
	uint32_t read_id;
	// uint32_t ref_id;
	const char* read_name;
	// uint32_t read_start_position;
	// uint32_t read_end_position;
	// uint32_t fragment_start_position;
	// uint32_t fragment_length;
	// uint8_t mapq : 6, rev : 1, mapped : 1;
	// char* tags;

	ri_map_t* maps;
	uint32_t n_maps;

	uint32_t offset;
	float* events;
	mm128_t *prev_anchors;
	uint32_t n_prev_anchors;
	mm_reg1_t* creg; // This is for transition purposes.
	int n_cregs;
} ri_reg1_t;

typedef struct pipeline_ms{
	int n_processed, n_threads, io_n_threads, n_fp, cur_fp, n_f, cur_f;
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
 * @param idx		rindex (see rindex.h)
 * @param fn		path to the signal files
 * @param opt		mapping options
 * @param n_threads	number of threads to use in mapping
 * @param io_n_threads	number of threads to use in I/O operations (reading from a file)
 * 
 * @return		returns 0 if mapping is completed with no issues. -1, otherwise.
 */
int ri_map_file(const ri_idx_t *idx, const char *fn, const ri_mapopt_t *opt, int n_threads, int io_n_threads = 1);

/**
 * Map raw nanopore signals of many reads to a reference genome
 *
 * @param idx		rindex (see rindex.h)
 * @param n_segs	number of signal files
 * @param fn		paths to the signal files
 * @param opt		mapping options
 * @param n_threads	number of threads to use in mapping
 * @param io_n_threads	number of threads to use in I/O operations (reading from a file)
 * 
 * @return			returns 0 if mapping is completed with no issues. -1, otherwise.
 */
int ri_map_file_frag(const ri_idx_t *idx, int n_segs, const char **fn, const ri_mapopt_t *opt, int n_threads, int io_n_threads = 1);


// functions used in rmap.cpp and rawhash_wrapper.cpp

// init thread-local buffer
ri_tbuf_t *ri_tbuf_init(void);
// destroy thread-local buffer
void ri_tbuf_destroy(ri_tbuf_t *b);

// write out relevant mapping results to stdout, then free reg0
// was_mapped: if false, write out "*" for mapping info even when read was mapped, e.g. when sequence until terminates the pipeline)
void write_out_mappings_to_stdout(ri_reg1_t *reg0, const ri_idx_t *ri, bool was_mapped=true);

// frees mappings (works even if no mappings were set), which is not freed by "free_most_of_ri_reg1_t" and reg0 itself
void free_mappings_ri_reg1_t(ri_reg1_t* reg0);
// finally need to call free(reg0)

// map a single read to a reference, called in parallel
void map_worker_for(void *_data,
					long i,
					int tid);

#ifdef __cplusplus
}
#endif
#endif //RMAP_H