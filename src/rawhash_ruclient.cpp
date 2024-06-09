#ifdef RUCLIENT_ENABLED
#include "rawhash_ruclient.hpp"

#include "kalloc.h"

#ifdef PROFILERH
// total counts across all reads
double ri_maptime = 0.0;
#endif


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

Decision ChannelReadMapper::try_mapping_new_chunk(ReadIdentifier const& read_ident, ChunkType const& chunk, uint32_t chunk_idx) {
	assert(c_count == chunk_idx);
	assert(c_count < opt->max_num_chunk);
	
	double t = ri_realtime();
	uint32_t l_sig = 0; // actual chunk length after removing outliers
	float* signal_with_outliers_removed = convert_to_float_signal_with_outliers_removed(
		(int16_t*)chunk.raw_data.get(), // from bytes to int16, todo1: use .release()
		chunk.chunk_length(), &l_sig,
		calibration.range, calibration.digitisation, calibration.offset
	);

	bool mapped = continue_mapping_with_new_chunk(ri, l_sig, signal_with_outliers_removed, reg0, b, opt, read_ident.read_id.c_str(), &mean_sum, &std_dev_sum, &n_events_sum);
	mapping_time += ri_realtime() - t;

	c_count += 1;
	qlen += l_sig;

	// try mapping with less stringent criteria
	if (c_count == opt->max_num_chunk) {
		try_mapping_if_none_found(reg0, opt);
		mapped = (reg0->n_maps != 0);
	}

	// mapping strategy
	if ((c_count == opt->max_num_chunk) && !mapped) {
		return Decision::STOP_RECEIVING;
	} else if (!mapped) {
		return Decision::UNDECIDED;
	} else {
		return Decision::REJECT;
	}
}

ChannelReadMapper::ChannelReadMapper(const ri_idx_t *ri, const ri_mapopt_t *opt, 
		SingleChannelCalibration calibration, ri_tbuf_t* b, std::string read_id) :
		ri(ri), opt(opt), calibration(calibration), b(b), read_id(read_id) {
	
	mean_sum = 0;
	std_dev_sum = 0;
	n_events_sum = 0;
	c_count = 0;
	mapping_time = 0;
	qlen = 0; // received signal length
	reg0 = (ri_reg1_t*)malloc(sizeof(ri_reg1_t)); // todo: calloc and push to queue
	init_reg1_t(reg0);
}

void ChannelReadMapper::populate_mappings() {
	#ifdef PROFILERH
	ri_maptime += mapping_time;
	#endif
	// todo4: c_count or c_count + 1?
	// passing 0 for the full read length as well because we only know for sure at the end of sequencing (e.g. for stop_receiving)
	reg0->read_name = read_id.c_str(); // todo4: does not work if put into constructor
	compute_tag_and_mapping_info(c_count, qlen, reg0, opt, mapping_time, 0, ri); // need to use strdup when not using cast and free read_id at the end

	free_most_of_ri_reg1_t(b->km, reg0);
	km_destroy_and_recreate(&(b->km));
}

void ChannelReadMapper::write_out_mappings() {

	// todo: write to queue instead which writes out reads
	write_out_mappings_to_stdout(reg0, ri);
}

void ChannelReadMapper::free_mappings() {
	if (reg0 != NULL) {
		free_mappings_ri_reg1_t(reg0);
		free(reg0);
		reg0 = NULL;
	}
}

ChannelReadMapper::~ChannelReadMapper() {
	// free_mappings(); // todo4: why can't this be called again?? if(reg0 != NULL) should prevent it appropriately
}

RawHashDecisionMaker::RawHashDecisionMaker(const ri_idx_t *ri, const ri_mapopt_t *opt, SingleChannelCalibration calibration): 
	b(ri_tbuf_init()), ri(ri), opt(opt), calibration(calibration), channel_read_mapper(ri, opt, calibration, b, "") {
}
RawHashDecisionMaker::~RawHashDecisionMaker() {
	// free buffer and reg
	ri_tbuf_destroy(b);
}

Decision RawHashDecisionMaker::decide(ReadIdentifier const& read_ident, ChunkType const& chunk, uint32_t chunk_idx) {
	// check if it maps
	if (chunk_idx == 0) {
		channel_read_mapper = ChannelReadMapper(ri, opt, calibration, b, read_ident.read_id);
	}

	return channel_read_mapper.try_mapping_new_chunk(read_ident, chunk, chunk_idx);
}

void RawHashDecisionMaker::mark_final_decision(ReadIdentifier const& read_ident, Decision const& decision) {
	// channel_read_mapper.mark_final_decision(read_ident, decision);
	channel_read_mapper.populate_mappings();
	channel_read_mapper.write_out_mappings();
	channel_read_mapper.free_mappings();
}

#endif