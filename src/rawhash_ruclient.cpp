#ifdef RUCLIENT_ENABLED
#include "rawhash_ruclient.hpp"

#include "kalloc.h"


// convert to float and remove outliers
// we receive int16 (checked in main.cpp) which we convert to float similar to ri_read_sig_slow5
// todo1: into function
// todo1: check this function works
float* convert_to_float_signal_with_outliers_removed(int16_t* int_raw_signal, uint32_t orig_length, uint32_t* new_l_sig, double range, uint32_t digitisation, double offset) {
	float *sigF = (float*)malloc(orig_length * sizeof(float));
	uint32_t l_sig = 0;
	float pa = 0.0f;
	float scale = range/float(digitisation);
	
	for(int i = 0; i < orig_length; ++i){
		pa = (int_raw_signal[i]+offset)*scale;
		if (pa > 30.0f && pa < 200.0f) {
			sigF[l_sig++] = pa;
		}
	}
    if (l_sig == 0) {
        // realloc may not work when size is 0 bytes, depending on the implementation
        free(sigF);
        sigF = NULL;
    } else {
        sigF = (float*)realloc(sigF, l_sig * sizeof(float));
	    assert(sigF != NULL);
    }
	*new_l_sig = l_sig;
	return sigF;
}

RawHashDecisionMaker::RawHashDecisionMaker(const ri_idx_t *ri, const ri_mapopt_t *opt, SingleChannelCalibration calibration): ri(ri), opt(opt), calibration(calibration) {
	b = ri_tbuf_init();
}
RawHashDecisionMaker::~RawHashDecisionMaker() {
	// free buffer and reg
	ri_tbuf_destroy(b);
}
Decision RawHashDecisionMaker::decide(ReadIdentifier const& read_ident, ChunkType const& chunk, uint32_t chunk_idx) {
	// check if it maps
	if (chunk_idx == 0) {
		mean_sum = 0;
		std_dev_sum = 0;
		n_events_sum = 0;
		c_count = 0;
		mapping_time = 0;
		qlen = 0; // received signal length
		reg0 = (ri_reg1_t*)malloc(sizeof(ri_reg1_t)); // todo: calloc and push to queue
		// start new read
		init_reg1_t(reg0);
	}

	bool unmapped;
	if (c_count < opt->max_num_chunk) {
		double t = ri_realtime();

		uint32_t l_sig = 0; // actual chunk length after removing outliers
		float* signal_with_outliers_removed = convert_to_float_signal_with_outliers_removed(
			(int16_t*)chunk.raw_data.get(), // from bytes to int16
			chunk.chunk_length(), &l_sig,
			calibration.range, calibration.digitisation, calibration.offset
		);

		unmapped = !continue_mapping_with_new_chunk(ri, l_sig, signal_with_outliers_removed, reg0, b, opt, read_ident.read_id.c_str(), &mean_sum, &std_dev_sum, &n_events_sum);
		mapping_time += ri_realtime() - t;
	
		c_count += 1;
		qlen += l_sig;
	} else {
		try_mapping_if_none_found(reg0, opt);
		unmapped = (reg0->n_maps == 0);
	}
	

	if (unmapped) {
		return Decision::STOP_RECEIVING;
	} else {
		return Decision::REJECT;
	}
}

void RawHashDecisionMaker::mark_final_decision(ReadIdentifier const& read_ident, Decision const& decision) {
	// write to file

	#ifdef PROFILERH
	ri_maptime += mapping_time;
	#endif

	// todo4: c_count or c_count + 1?
	// passing 0 for the full read length as well because we only know for sure at the end of sequencing (e.g. for stop_receiving)
	reg0->read_name = read_ident.read_id.c_str();
	compute_tag_and_mapping_info(c_count, qlen, reg0, opt, mapping_time, 0, ri); // need to use strdup when not using cast and free read_id at the end
	free_most_of_ri_reg1_t(b->km, reg0);
	km_destroy_and_recreate(&(b->km));

	// todo: write to queue instead which writes out reads
	write_out_mappings_to_stdout(reg0, ri);
	free_mappings_ri_reg1_t(reg0);
	free(reg0);
}

// // map reads coming from readuntil
// void map_reads_from_ru() {
// 	// see map_worker_pipeline for what is freed and how (memory allocator)
// 	// reuse buffer and reg for each read
// 	ri_tbuf_t* b = ri_tbuf_init();
// 	ri_reg1_t* reg = (ri_reg1_t*)malloc(sizeof(ri_reg1_t));

// todo: cannot reuse reg0 due to write out

// todo: write_out_mappings_and_free(reg0, ri);
// }

// map reads from a single channel to a single reference
// class ChannelDataMapper {
// 	public:
// 		ChannelDataMapper(const ri_idx_t *ri, const ri_mapopt_t *opt): ri(ri), opt(opt) {
// 			b = ri_tbuf_init();
// 		}

// 		// called by thread, tries to map current read in channel
// 		void map_reads() {

// 			{
// 				ri_reg1_t* reg0 = (ri_reg1_t*)malloc(sizeof(ri_reg1_t)); // todo: calloc and push to queue
// 				// start new read
// 				init_reg1_t(reg0);

// 				uint32_t c_count = 0;
// 				double mean_sum = 0, std_dev_sum = 0;
// 				uint32_t n_events_sum = 0;
// 				double t = ri_realtime();

// 				// todo
// 				uint32_t l_chunk = 10;
// 				ri_sig_t* sig = NULL;
// 				uint32_t qlen = 100;
				
// 				// process chunks

// 				double mapping_time = ri_realtime() - t;
// 				#ifdef PROFILERH
// 				ri_maptime += mapping_time;
// 				#endif

// 				try_mapping_if_none_found(reg0, opt);
// 				compute_tag_and_mapping_info(c_count, c_count * l_chunk, reg0, opt, mapping_time, qlen, sig->name, ri);
// 				free_most_of_ri_reg1_t(b->km, reg0);
// 				km_destroy_and_recreate(&(b->km));

// 				// todo: write to queue instead which writes out reads
// 				write_out_mappings_and_free(reg0, ri);

// 				// wait until read has ended
// 			}

// 		}

// 		~ChannelDataMapper() {
// 			// free buffer and reg
// 			ri_tbuf_destroy(b);
// 		}
// 	private:
// 		ri_tbuf_t* b = NULL; // thread-local buffer
// 		const ri_idx_t *ri; // reference index
// 		const ri_mapopt_t *opt; // mapping options
// };

#endif