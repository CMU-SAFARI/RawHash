#ifndef RMAP_H
#define RMAP_H

#include "rindex.h"
#include "rsig.h"
#include "chain.h"
#include "tensorflow/lite/core/c/c_api.h"

#ifdef __cplusplus
extern "C" {
#endif

// info about mapping of a (partial) read to a reference
typedef struct ri_map_s{
	uint32_t c_id; //chain index
	// num_signals_per_event = processed_signal_length / num_events_in_processed_signal, per read averaged until reaching the mapping decision
	// num_signals_per_bp = 4000 / 450
	// num_bp_per_event = num_signals_per_event / num_signals_per_bp
	uint32_t read_length; // estimated length used for mapping (not full read length) in bps based on number of detected events (rather than signal length), num_events * num_bp_per_event
	uint32_t ref_id;
	uint32_t read_start_position; // mapping start position in the read, in terms of bps, num_events_before_mapping * num_bp_per_event
	uint32_t read_end_position; // mapping end position inclusive in the read, in terms of bps, num_events_until_mapping_end * num_bp_per_event
	uint32_t fragment_start_position; // start position of mapping on reference, in terms of samples (not basepairs)
	uint32_t fragment_length;
	uint8_t mapq : 6, rev : 1, mapped : 1;
	char* tags;
} ri_map_t;

// mapping of (partial) read chunks of one read to a reference, can be used to iteratively add new chunks to check if it maps (reusing mapping attempts from previous chunks of read)
typedef struct ri_reg1_s{
	uint32_t read_id;
	// uint32_t ref_id;
	const char* read_name; // assigned for outputting the mapping only, not owned, so pointer should remain valid during the lifetime of this object
	// uint32_t read_start_position;
	// uint32_t read_end_position;
	// uint32_t fragment_start_position;
	// uint32_t fragment_length;
	// uint8_t mapq : 6, rev : 1, mapped : 1;
	// char* tags;

	ri_map_t* maps;
	uint32_t n_maps; // number of mapping chains

	uint32_t offset; // number of events processed from this read so far
	float* events;
	mm128_t *prev_anchors;
	uint32_t n_prev_anchors;
	mm_reg1_t* creg; // This is for transition purposes.
	int n_cregs;
} ri_reg1_t;

// memory buffer for thread-local storage during mapping
typedef struct ri_tbuf_s {
	void *km;
	int rep_len, frag_gap; // todo4: these seem to be unused in map_worker_for
	// for ML predictions
	TfLiteInterpreter* interpreter;
	TfLiteTensor* input_tensor;
}ri_tbuf_t;

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

	TfLiteModel* map_model;
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
 * 
 * @return		returns 0 if mapping is completed with no issues. -1, otherwise.
 */
int ri_map_file(const ri_idx_t *idx, const char *fn, const ri_mapopt_t *opt, int n_threads);

/**
 * Map raw nanopore signals of many reads to a reference genome
 *
 * @param idx		rindex (see rindex.h)
 * @param n_segs	number of signal files
 * @param fn		paths to the signal files
 * @param opt		mapping options
 * @param n_threads	number of threads to use in mapping
 *
 * @return			returns 0 if mapping is completed with no issues. -1, otherwise.
 */
int ri_map_file_frag(const ri_idx_t *idx, int n_segs, const char **fn, const ri_mapopt_t *opt, int n_threads);


// functions used in rmap.cpp and rawhash_ruclient.cpp

// init thread-local buffer
ri_tbuf_t *ri_tbuf_init(void);
// destroy thread-local buffer
void ri_tbuf_destroy(ri_tbuf_t *b);
// init mapping info
void init_reg1_t(ri_reg1_t* reg0);

// write out relevant mapping results to stdout, then free reg0
// was_mapped: if false, write out "*" for mapping info even when read was mapped, e.g. when sequence until terminates the pipeline)
// frees things not freed by "free_most_of_ri_reg1_t" and reg0 itself
void write_out_mappings_and_free(ri_reg1_t *reg0, const ri_idx_t *ri, bool was_mapped=true);

// compute tag field and mapping info, for later writing out mapping results to stdout
// this allows to call free_most_of_ri_reg1_t() after this function while guaranteeing it can still be written out
// num_chunks: number of processed chunks
// processed_len: length of processed signal in terms of samples
// qlen: total length of read (all received raw signals)
void compute_tag_and_mapping_info(uint32_t c_count, uint32_t processed_len, ri_reg1_t *reg0, const ri_mapopt_t *opt, double mapping_time, uint32_t qlen, const ri_idx_t* ri);

// if no mapping was found (because several chains mapped, but none with sufficient quality), but the first chain still has sufficient quality (as if only one chain was found), then accept it as a mapping
void try_mapping_if_none_found(ri_reg1_t *reg0, const ri_mapopt_t *opt);

// does not free "maps" field which is freed after outputting
void free_most_of_ri_reg1_t(void* km, ri_reg1_t* reg0);

// map new chunk to a reference, reusing mapping computations for previous chunks
// return: true if mapping was successful
bool continue_mapping_with_new_chunk(const ri_idx_t *ri, // reference index
	const uint32_t s_len, // chunk length
	const float *sig, // raw signal for chunk
	ri_reg1_t* reg0, // to store mapping results, reuse previous computations
	ri_tbuf_t *b,
	const ri_mapopt_t *opt,
	const char *qname, // read id
	double* mean_sum, // for raw signal normalization, reused from old chunks and updated with new chunk
	double* std_dev_sum, // see mean_sum arg
	uint32_t* n_events_sum // number of events in all seen chunks (not including current chunk)
);
#ifdef __cplusplus
}
#endif
#endif //RMAP_H