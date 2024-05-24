#pragma once
#ifdef RUCLIENT_ENABLED

#include "rmap.h"

#include "ru_client/ru_method.hpp"
using namespace ru_client;
// #include "ru_client/ont_device_client.hpp"
// typedef DeviceServiceClient::Calibration Calibration;

// calibration of a single channel
struct SingleChannelCalibration {
	double range;
	uint32_t digitisation;
	double offset;
};

struct Alignment {
	char* contig;
	uint32_t ref_start;
	uint32_t ref_end;
	bool pos_strand;
};

typedef DecisionMaker::ChunkType ChunkType;

/**
 * @brief Maps read in current channel, progressively receiving new chunks
 * 
 */
class ChannelReadMapper {
	// todo: rename to RawHashSingleReadMapper

public:
// todo: calibration by reference
	ChannelReadMapper(const ri_idx_t *ri, const ri_mapopt_t *opt, SingleChannelCalibration calibration, ri_tbuf_t* b, std::string read_id);
	Decision try_mapping_new_chunk(ReadIdentifier const& read_ident, ChunkType const& chunk, uint32_t chunk_idx);

	// populate mappings from reg0, and free irrelevant stuff
	void populate_mappings();
	// write out to stdout
	void write_out_mappings();
	// destroy immediately to avoid doing it in the destructor (when going to new read where we want to minimize latency), only destroys if reg0 was set, i.e. at least one chunk received
	void free_mappings();

	~ChannelReadMapper();

	double mapping_time;
private:
	const ri_idx_t *ri; // reference index
	const ri_mapopt_t *opt; // mapping options
	SingleChannelCalibration calibration;
	ri_tbuf_t* b = NULL; // thread-local buffer, not owned by this class!

	// read-specific
	std::string read_id; // since reg0->read_name is const char*, we manage it explicitly
	ri_reg1_t* reg0;
	double mean_sum, std_dev_sum;
	uint32_t n_events_sum;
	uint32_t c_count;
	uint32_t qlen; // number of raw signals seen for this read
};

/**
 * @brief Decision maker taking decisions for reads on a single channel
 * 
 */
class RawHashDecisionMaker : public DecisionMaker {
public:
	RawHashDecisionMaker(const ri_idx_t *ri, const ri_mapopt_t *opt, SingleChannelCalibration calibration);
	virtual ~RawHashDecisionMaker();
	virtual Decision decide(ReadIdentifier const& read_ident, ChunkType const& chunk, uint32_t chunk_idx) override;

	virtual void mark_final_decision(ReadIdentifier const& read_ident, Decision const& decision);
private:
	ri_tbuf_t* b = NULL; // thread-local buffer
	const ri_idx_t *ri; // reference index
	const ri_mapopt_t *opt; // mapping options
	SingleChannelCalibration calibration;

	// per read
	ChannelReadMapper channel_read_mapper;
};

#endif