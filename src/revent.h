#ifndef REVENT_H
#define REVENT_H

#include <stdint.h>
#include "roptions.h"

#ifdef __cplusplus
extern "C" {
#endif

float* normalize_signal(void *km,
						const float* sig,
						const uint32_t s_len,
						double* mean_sum,
						double* std_dev_sum,
						uint32_t* n_events_sum,
						uint32_t* n_sig);
/**
 * Detects events from signals
 *
 * @param km	thread-local memory pool; using NULL falls back to malloc()
 * @param s_len	length of $sig
 * @param sig	signal values
 * @param opt	mapping options @TODO: Should be decoupled from the mapping options
 * @param n		number of events
 * 
 * @return		list of event values of length $n
 */
float* detect_events(void *km,
					 const uint32_t s_len,
					 const float* sig,
					 const uint32_t window_length1,
					 const uint32_t window_length2,
					 const float threshold1,
					 const float threshold2,
					 const float peak_height,
					 double* mean_sum,
					 double* std_dev_sum,
					 uint32_t* n_events_sum,
					 uint32_t *n_events);

#ifdef __cplusplus
}
#endif
#endif //REVENT_H