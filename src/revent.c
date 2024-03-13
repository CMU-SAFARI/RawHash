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

static inline void comp_prefix_prefixsq(const float *sig, uint32_t s_len, float* prefix_sum, float* prefix_sum_square) {
	
	assert(s_len > 0);

	prefix_sum[0] = 0.0f;
	prefix_sum_square[0] = 0.0f;
	for (uint32_t i = 0; i < s_len; ++i) {
		prefix_sum[i+1] = prefix_sum[i] + sig[i];
		prefix_sum_square[i+1] = prefix_sum_square[i] + sig[i]*sig[i];
	}
}

static inline float* comp_tstat(void *km, const float *prefix_sum, const float *prefix_sum_square, uint32_t s_len, uint32_t w_len) {
  
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
		float combined_var = sumsq1 / w_len - mean1 * mean1 + sumsq2 / w_len - mean2 * mean2;
		// Prevent problem due to very small variances
		combined_var = fmaxf(combined_var, eta);
		// t-stat
		//  Formula is a simplified version of Student's t-statistic for the
		//  special case where there are two samples of equal size with
		//  differing variance
		const float delta_mean = mean2 - mean1;
		tstat[i] = fabs(delta_mean) / sqrt(combined_var / w_len);
	}
	// fudge boundaries
	memset(tstat+s_len-w_len+1, 0, (w_len)*sizeof(float));

	return tstat;
}

static inline uint32_t gen_peaks(ri_detect_t *short_detector,
								 ri_detect_t *long_detector,
								 const float peak_height,
								 uint32_t* peaks) {
	
	assert(short_detector->s_len == long_detector->s_len);

	uint32_t curInd = 0;

	// uint32_t ndetector = 2;
	ri_detect_t *detectors[2];  // = {short_detector, long_detector};
	detectors[0] = short_detector;
	detectors[1] = long_detector;
	for (uint32_t i = 0; i < short_detector->s_len; i++) {
		for (uint32_t k = 0; k < 2; k++) {
			ri_detect_t *detector = detectors[k];
			// Carry on if we've been masked out
			if (detector->masked_to >= i) continue;

			float current_value = detector->sig[i];
			if (detector->peak_pos == detector->DEF_PEAK_POS) {
				// CASE 1: We've not yet recorded a maximum
				if (current_value < detector->peak_value) {
					// Either record a deeper minimum...
					detector->peak_value = current_value;
				} else if (current_value - detector->peak_value > peak_height) {
					// ...or we've seen a qualifying maximum
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
				// Dominate other tstat signals if we're going to fire at some point
				if (detector == short_detector) {
					if (detector->peak_value > detector->threshold) {
						long_detector->masked_to = detector->peak_pos + detector->window_length;
						long_detector->peak_pos = long_detector->DEF_PEAK_POS;
						long_detector->peak_value = long_detector->DEF_PEAK_VAL;
						long_detector->valid_peak = 0;
					}
				}
				// Have we convinced ourselves we've seen a peak
				if (detector->peak_value - current_value > peak_height && detector->peak_value > detector->threshold) {
					detector->valid_peak = 1;
				}
				// Finally, check the distance if this is a good peak
				if (detector->valid_peak && (i - detector->peak_pos) > detector->window_length / 2) {
					// Emit the boundary and reset
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

/**
 * @brief Generates events from peaks, prefix sums and s_len.
 * 
 * @param km Pointer to memory manager.
 * @param peaks Array of peak positions.
 * @param peak_size Size of peaks array.
 * @param prefix_sum Array of prefix sums.
 * @param prefix_sum_square Array of prefix sums squared.
 * @param s_len Length of the signal.
 * @param n_events Pointer to the number of events generated.
 * @return float* Pointer to the array of generated events.
 */
static inline float* gen_events(void *km,
								const uint32_t *peaks,
								const uint32_t peak_size,
								const float *prefix_sum,
								const float *prefix_sum_square,
								const uint32_t s_len,
								uint32_t* n_events)
{
	uint32_t n_ev = 1;

	for (uint32_t i = 1; i < peak_size; ++i)
		if (peaks[i] > 0 && peaks[i] < s_len) n_ev++;

	float* events = (float*)ri_kmalloc(km, n_ev*sizeof(float));
	float l_prefixsum = 0, l_peak = 0;

	for (uint32_t pi = 0; pi < n_ev - 1; pi++){
		events[pi] = (prefix_sum[peaks[pi]] - l_prefixsum)/(peaks[pi]-l_peak);
		l_prefixsum = prefix_sum[peaks[pi]];
		l_peak = peaks[pi];
	}

	events[n_ev-1] = (prefix_sum[s_len] - l_prefixsum)/(s_len-l_peak);

	(*n_events) = n_ev;
	return events;
}

static inline float* normalize_signal(void *km,
									  const float* sig,
									  uint32_t s_len,
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
					 uint32_t s_len,
					 const float* sig,
					 uint32_t window_length1,
					 uint32_t window_length2,
					 float threshold1,
					 float threshold2,
					 float peak_height,
					 double* mean_sum,
					 double* std_dev_sum,
					 uint32_t* n_events_sum,
					 uint32_t *n)
{
	float* prefix_sum = (float*)ri_kcalloc(km, s_len+1, sizeof(float));
	float* prefix_sum_square = (float*)ri_kcalloc(km, s_len+1, sizeof(float));
	float* events = 0;

	//Normalize the signal
	events = normalize_signal(km, sig, s_len, mean_sum, std_dev_sum, n_events_sum, n);
	comp_prefix_prefixsq(events, (*n), prefix_sum, prefix_sum_square);
	
	float* tstat1 = comp_tstat(km, prefix_sum, prefix_sum_square, (*n), window_length1);
	float* tstat2 = comp_tstat(km, prefix_sum, prefix_sum_square, (*n), window_length2);
	ri_detect_t short_detector = {.DEF_PEAK_POS = -1,
								  .DEF_PEAK_VAL = FLT_MAX,
								  .sig = tstat1,
								  .s_len = (*n),
								  .threshold = threshold1,
								  .window_length = window_length1,
								  .masked_to = 0,
								  .peak_pos = -1,
								  .peak_value = FLT_MAX,
								  .valid_peak = 0};

	ri_detect_t long_detector = {.DEF_PEAK_POS = -1,
								 .DEF_PEAK_VAL = FLT_MAX,
								 .sig = tstat2,
								 .s_len = (*n),
								 .threshold = threshold2,
								 .window_length = window_length2,
								 .masked_to = 0,
								 .peak_pos = -1,
								 .peak_value = FLT_MAX,
								 .valid_peak = 0};

	uint32_t* peaks = (uint32_t*)ri_kmalloc(km, (*n) * sizeof(uint32_t));
	uint32_t n_peaks = gen_peaks(&short_detector, &long_detector, peak_height, peaks);
	if(n_peaks > 0) events = gen_events(km, peaks, n_peaks, prefix_sum, prefix_sum_square, (*n), n);
	else {ri_kfree(km, events); events = 0;}
	ri_kfree(km, tstat1); ri_kfree(km, tstat2); ri_kfree(km, prefix_sum); ri_kfree(km, prefix_sum_square); ri_kfree(km, peaks);

	return events;
}
