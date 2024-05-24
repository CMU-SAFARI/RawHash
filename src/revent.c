#include "revent.h"
#include "kalloc.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include "rutils.h"

//Some of the functions here are adopted from the Sigmap implementation (https://github.com/haowenz/sigmap/tree/c9a40483264c9514587a36555b5af48d3f054f6f). We have optimized the Sigmap implementation to work with the hash tables efficiently.

typedef struct ri_detect_s {
	int DEF_PEAK_POS;
	float DEF_PEAK_VAL;
	float *sig;
	uint32_t s_len;
	float threshold;
	uint32_t window_length;
	uint32_t masked_to;
	int peak_pos;
	float peak_value;
	int valid_peak;
}ri_detect_t;

static inline void comp_prefix_prefixsq(const float *sig,
										const uint32_t s_len,
										float* prefix_sum,
										float* prefix_sum_square)
{
	assert(s_len > 0);

	prefix_sum[0] = 0.0f;
	prefix_sum_square[0] = 0.0f;
	for (uint32_t i = 0; i < s_len; ++i) {
		prefix_sum[i+1] = prefix_sum[i] + sig[i];
		prefix_sum_square[i+1] = prefix_sum_square[i] + sig[i]*sig[i];
	}
}

static inline float* comp_tstat(void *km,
							 	const float *prefix_sum,
								const float *prefix_sum_square,
								const uint32_t s_len,
								const uint32_t w_len)
{
  	const float eta = FLT_MIN;
	float* tstat = (float*)ri_kcalloc(km, s_len+1, sizeof(float));
	if (s_len < 2*w_len || w_len < 2) return tstat;
	memset(tstat, 0, w_len*sizeof(float));
	
	for (uint32_t i = w_len; i <= s_len - w_len; ++i) {
		float sum1 = prefix_sum[i];
		float sumsq1 = prefix_sum_square[i];
		if (i > w_len) {
			sum1 -= prefix_sum[i - w_len];
			sumsq1 -= prefix_sum_square[i - w_len];
		}
		float sum2 = prefix_sum[i + w_len] - prefix_sum[i];
		float sumsq2 = prefix_sum_square[i + w_len] - prefix_sum_square[i];
		float mean1 = sum1 / w_len;
		float mean2 = sum2 / w_len;
		float combined_var = (sumsq1/w_len - mean1*mean1 + sumsq2/w_len - mean2*mean2)/w_len;
		// Prevent problem due to very small variances
		combined_var = fmaxf(combined_var, eta);
		// t-stat
		//  Formula is a simplified version of Student's t-statistic for the
		//  special case where there are two samples of equal size with
		//  differing variance
		const float delta_mean = mean2 - mean1;
		tstat[i] = fabs(delta_mean) / sqrt(combined_var);
	}
	// fudge boundaries
	memset(tstat+s_len-w_len+1, 0, (w_len)*sizeof(float));

	return tstat;
}

// static inline float calculate_adaptive_peak_height(const float *prefix_sum, const float *prefix_sum_square, uint32_t current_index, uint32_t window_length, float base_peak_height) {
//     // Ensure we don't go beyond signal bounds
//     uint32_t start_index = current_index > window_length ? current_index - window_length : 0;
//     uint32_t end_index = current_index + window_length;

//     float sum = prefix_sum[end_index] - prefix_sum[start_index];
//     float sumsq = prefix_sum_square[end_index] - prefix_sum_square[start_index];
//     float mean = sum / (end_index - start_index);
//     float variance = (sumsq / (end_index - start_index)) - (mean * mean);
//     float stddev = sqrtf(variance);

//     // Example adaptive strategy: Increase peak height in high-variance regions
//     return base_peak_height * (1 + stddev);
// }

static inline uint32_t gen_peaks(ri_detect_t **detectors,
								 const uint32_t n_detectors,
								 const float peak_height,
								 const float *prefix_sum,
								 const float *prefix_sum_square,
								 uint32_t* peaks) {

	uint32_t curInd = 0;
	for (uint32_t i = 0; i < detectors[0]->s_len; i++) {
		for (uint32_t k = 0; k < n_detectors; k++) {
			ri_detect_t *detector = detectors[k];
			if (detector->masked_to >= i) continue;

			float current_value = detector->sig[i];
			// float adaptive_peak_height = calculate_adaptive_peak_height(prefix_sum, prefix_sum_square, i, detector->window_length, peak_height);

			if (detector->peak_pos == detector->DEF_PEAK_POS) {
				// CASE 1: We've not yet recorded any maximum
				if (current_value < detector->peak_value) { // A deeper minimum:
					detector->peak_value = current_value;
				} else if (current_value - detector->peak_value > peak_height) {
					// ...or a qualifying maximum:
					detector->peak_value = current_value;
					detector->peak_pos = i;
					// otherwise, wait to rise high enough to be considered a peak
				}
			} else {
				// CASE 2: In an existing peak, waiting to see if it is good
				if (current_value > detector->peak_value) {
					// Update the peak
					detector->peak_value = current_value;
					detector->peak_pos = i;
				}
				// Tell other detectors no need to check for a peak until a certain point
				if (detector->peak_value > detector->threshold) {
					for(int n_d = k+1; n_d < n_detectors; n_d++){
						detectors[n_d]->masked_to = detector->peak_pos + detectors[0]->window_length;
						detectors[n_d]->peak_pos = detectors[n_d]->DEF_PEAK_POS;
						detectors[n_d]->peak_value = detectors[n_d]->DEF_PEAK_VAL;
						detectors[n_d]->valid_peak = 0;
					}
				}
				// There is a good peak
				if (detector->peak_value - current_value > peak_height && 
					detector->peak_value > detector->threshold) {
					detector->valid_peak = 1;
				}
				// Check if we are now further away from the current peak
				if (detector->valid_peak && (i - detector->peak_pos) > detector->window_length / 2) {
					peaks[curInd++] = detector->peak_pos;
					detector->peak_pos = detector->DEF_PEAK_POS;
					detector->peak_value = current_value;
					detector->valid_peak = 0;
				}
			}
		}
	}

	return curInd;
}

int compare_floats(const void* a, const void* b) {
    const float* da = (const float*) a;
    const float* db = (const float*) b;
    return (*da > *db) - (*da < *db);
}

float calculate_mean_of_filtered_segment(float* segment,
										 const uint32_t segment_length)
{
    // Calculate median and IQR
    qsort(segment, segment_length, sizeof(float), compare_floats); // Assuming compare_floats is already defined
    float q1 = segment[segment_length / 4];
    float q3 = segment[3 * segment_length / 4];
    float iqr = q3 - q1;
    float lower_bound = q1 - iqr;
    float upper_bound = q3 + iqr;

    float sum = 0.0;
    uint32_t count = 0;
    for (uint32_t i = 0; i < segment_length; i++) {
        if (segment[i] >= lower_bound && segment[i] <= upper_bound) {
            sum += segment[i];
            ++count;
        }
    }

    // Return the mean of the filtered segment
    return count > 0 ? sum / count : 0; // Ensure we don't divide by zero
}

/**
 * @brief Generates events from peaks, prefix sums and s_len.
 * 
 * @param km Pointer to memory manager.
 * @param peaks Array of peak positions.
 * @param peak_size Size of peaks array.
 * @param prefix_sum Array of prefix sums.
 * @param s_len Length of the signal.
 * @param n_events Pointer to the number of events generated.
 * @return float* Pointer to the array of generated events.
 */
static inline float* gen_events(void *km,
								float* sig,
								const uint32_t *peaks,
								const uint32_t peak_size,
								const uint32_t s_len,
								uint32_t* n_events)
{
	uint32_t n_ev = 0;

	for (uint32_t pi = 0; pi < peak_size; ++pi)
		if (peaks[pi] > 0 && peaks[pi] < s_len) n_ev++;

	float* events = (float*)ri_kmalloc(km, n_ev*sizeof(float));

	uint32_t start_idx = 0, segment_length = 0;

	for (uint32_t pi = 0, i = 0; pi < peak_size && i < n_ev; pi++){
		if (!(peaks[pi] > 0 && peaks[pi] < s_len)) continue;

    	segment_length = peaks[pi] - start_idx;
		events[i++] = calculate_mean_of_filtered_segment(sig + start_idx, segment_length);
		start_idx = peaks[pi];
	}

	(*n_events) = n_ev;
	return events;
}

// sig: events
// return values:
// n_sig: will be assigned the number of events detected in the signal
float* normalize_signal(void *km,
						const float* sig,
						const uint32_t s_len,
						double* mean_sum,
						double* std_dev_sum,
						uint32_t* n_events_sum,
						uint32_t* n_sig)
{
	double sum = (*mean_sum), sum2 = (*std_dev_sum);
	double mean = 0, std_dev = 0;
	float* events = (float*)ri_kcalloc(km, s_len, sizeof(float));

	for (uint32_t i = 0; i < s_len; ++i) {
		sum += sig[i];
		sum2 += sig[i]*sig[i];
	}

	(*n_events_sum) += s_len;
	(*mean_sum) = sum;
	(*std_dev_sum) = sum2;

	mean = sum/(*n_events_sum);
	std_dev = sqrt(sum2/(*n_events_sum) - (mean)*(mean));

	float norm_val = 0;
	int k = 0;
	for(uint32_t i = 0; i < s_len; ++i){
		norm_val = (sig[i]-mean)/std_dev;
		if(norm_val < 3 && norm_val > -3) events[k++] = norm_val;
	}

	(*n_sig) = k;

	return events;
}

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
					 uint32_t* n_events)
{
	float* prefix_sum = (float*)ri_kcalloc(km, s_len+1, sizeof(float));
	float* prefix_sum_square = (float*)ri_kcalloc(km, s_len+1, sizeof(float));

	//Normalize the signal
	uint32_t n_signals = 0;
	(*n_events) = 0;
	float* norm_signals = normalize_signal(km, sig, s_len, mean_sum, std_dev_sum, n_events_sum, &n_signals);
	if(n_signals == 0) return 0;
	comp_prefix_prefixsq(norm_signals, n_signals, prefix_sum, prefix_sum_square);
	
	float* tstat1 = comp_tstat(km, prefix_sum, prefix_sum_square, n_signals, window_length1);
	float* tstat2 = comp_tstat(km, prefix_sum, prefix_sum_square, n_signals, window_length2);
	ri_detect_t short_detector = {.DEF_PEAK_POS = -1,
								  .DEF_PEAK_VAL = FLT_MAX,
								  .sig = tstat1,
								  .s_len = n_signals,
								  .threshold = threshold1,
								  .window_length = window_length1,
								  .masked_to = 0,
								  .peak_pos = -1,
								  .peak_value = FLT_MAX,
								  .valid_peak = 0};

	ri_detect_t long_detector = {.DEF_PEAK_POS = -1,
								 .DEF_PEAK_VAL = FLT_MAX,
								 .sig = tstat2,
								 .s_len = n_signals,
								 .threshold = threshold2,
								 .window_length = window_length2,
								 .masked_to = 0,
								 .peak_pos = -1,
								 .peak_value = FLT_MAX,
								 .valid_peak = 0};

	uint32_t* peaks = (uint32_t*)ri_kmalloc(km, n_signals * sizeof(uint32_t));
	ri_detect_t *detectors[2] = {&short_detector, &long_detector};
	uint32_t n_peaks = gen_peaks(detectors, 2, peak_height, prefix_sum, prefix_sum_square, peaks);
	ri_kfree(km, tstat1); ri_kfree(km, tstat2); ri_kfree(km, prefix_sum); ri_kfree(km, prefix_sum_square);

	float* events = 0;
	if(n_peaks > 0) events = gen_events(km, norm_signals, peaks, n_peaks, n_signals, n_events);
	ri_kfree(km, norm_signals); ri_kfree(km, peaks);

	return events;
}
