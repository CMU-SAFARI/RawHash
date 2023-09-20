#ifndef REVENT_H
#define REVENT_H

#include <stdint.h>
#include "roptions.h"

#ifdef __cplusplus
extern "C" {
#endif

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
float* detect_events(void *km, uint32_t s_len, const float* sig, uint32_t window_length1, uint32_t window_length2, float threshold1, float threshold2, float peak_height, uint32_t *n);

#ifdef __cplusplus
}
#endif
#endif //REVENT_H