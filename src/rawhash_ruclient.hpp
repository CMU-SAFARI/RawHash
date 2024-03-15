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
	ri_reg1_t* reg0;
	double mean_sum, std_dev_sum;
	uint32_t n_events_sum;
	uint32_t c_count;
	double mapping_time;
	uint32_t qlen; // number of raw signals seen for this read
};

#endif