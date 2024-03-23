#include "rsketch.h"
#include <assert.h>
#include "rh_kvec.h"
#include <math.h>
#include <float.h>

static inline uint64_t hash64(uint64_t key, uint64_t mask){
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

uint32_t dynamic_quantize(float signal,
						  float fine_min,
						  float fine_max,
						  float fine_range,
						  uint32_t n_buckets)
{
    // Total range for normalization
    float minVal = -3.0, maxVal = 3.0;
    float range = maxVal - minVal;
	float coarse_coef1 = (1-fine_range)/2;
	float coarse_coef2 = fine_range + coarse_coef1;

    // Normalize the signal to [0, 1]
    float normalized = (signal - minVal) / range;

	float a = (fine_min - minVal) / range;
	float b = (fine_max - minVal) / range;

    // Conditional quantization based on the segment
    float quantized = fine_max;
    if (signal >= fine_min && signal <= fine_max) {
        // Within [fine_min, fine_max], map to a sub-range [a, b] in [0, 1],
		//then scale to [0, fine_range] for finer granularity
        quantized = fine_range * ((normalized - a) / (b - a));
    }
	else {
        // Outside [fine_min, fine_max], split the rest of [0, 1] into two and map accordingly
        if (normalized < 0.5) quantized = fine_range + coarse_coef1 * normalized; // Coarser granularity
        else quantized = coarse_coef2 + coarse_coef1 * normalized; // Coarser granularity
    }

    // Map the quantized value back to the range [0, 2^n_buckets - 1]
    uint32_t quantizedValue = (uint32_t)(quantized * (n_buckets-1));

    return quantizedValue;
}

void ri_sketch_min(void *km,
				   const float* s_values,
				   uint32_t id,
				   int strand,
				   int len,
				   float diff,
				   int w,
				   int e,
				   uint32_t quant_bit,
				   int k,
				   float fine_min,
				   float fine_max,
				   float fine_range,
				   mm128_v *p)
{
	assert(len > 0 && (w > 0 && w < 256) && e*quant_bit <= (64-RI_HASH_SHIFT));
	
	int j, buf_pos, min_pos;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };

	uint32_t span = k+e-1;
	
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_quant_bit = (1ULL<<quant_bit)-1;

	uint32_t n_buckets = 1UL<<quant_bit;

	memset(buf, 0xff, w * 16);
	rh_kv_resize(mm128_t, km, *p, p->n + len/w);

	int sigBufFull = 0;
	uint32_t f_pos, l, sigBufPos = 0;
	uint32_t l_sigpos = 0; //last signal position
	uint32_t f_tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

    for (f_pos = l = buf_pos = min_pos = 0; f_pos < len; ++f_pos) {
        if(f_pos > 0 && fabs(s_values[f_pos] - s_values[l_sigpos]) < diff) continue;

		l++;
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		l_sigpos = f_pos;

		f_tmpQuantSignal = dynamic_quantize(s_values[f_pos], fine_min, fine_max, fine_range, n_buckets)&mask_quant_bit;

		quantVal = (quantVal<<quant_bit|f_tmpQuantSignal)&mask_events;

		sigBuf[sigBufPos].y = id_shift | (uint32_t)f_pos<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}
		sigBuf[sigBufPos].x = hash64(quantVal, mask)<<RI_HASH_SHIFT | span;

		if(!sigBufFull) continue;

		info.x = sigBuf[sigBufPos].x;
		info.y = sigBuf[sigBufPos].y;

		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + e - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) rh_kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) rh_kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + e && min.x != UINT64_MAX) rh_kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + e - 1 && min.x != UINT64_MAX) rh_kv_push(mm128_t, km, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest e-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + e - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) rh_kv_push(mm128_t, km, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) rh_kv_push(mm128_t, km, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
    }
	if (min.x != UINT64_MAX)
		rh_kv_push(mm128_t, km, *p, min);
}

void ri_sketch_reg(void *km,
				   const float* s_values,
				   uint32_t id,
				   int strand,
				   int len,
				   float diff,
				   int e,
				   uint32_t quant_bit,
				   int k,
				   float fine_min,
				   float fine_max,
				   float fine_range,
				   mm128_v *p){

	assert(len > 0 && (uint32_t)e*quant_bit <= 64-RI_HASH_SHIFT);

	uint32_t span = k+e-1;
	
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_quant_bit = (1ULL<<quant_bit)-1;

	uint32_t n_buckets = 1UL<<quant_bit;

	int sigBufFull = 0;
	uint32_t f_pos, sigBufPos = 0, l_sigpos = 0; //last signal position
	uint32_t f_tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	rh_kv_resize(mm128_t, km, *p, p->n + len);

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

    for (f_pos = 0; f_pos < len; ++f_pos) {
        if((f_pos > 0 && fabs(s_values[f_pos] - s_values[l_sigpos]) < diff)) continue;

		l_sigpos = f_pos;

		f_tmpQuantSignal = dynamic_quantize(s_values[f_pos], fine_min, fine_max, fine_range, n_buckets)&mask_quant_bit;

		sigBuf[sigBufPos].y = id_shift | (uint32_t)f_pos<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}

		quantVal = (quantVal<<quant_bit|f_tmpQuantSignal)&mask_events;
		sigBuf[sigBufPos].x = (hash64(quantVal, mask)<<RI_HASH_SHIFT) | span;
		
		if(!sigBufFull) continue;

		rh_kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);
    }
}

void ri_sketch_reg_rev(void *km,
					   const float* s_values,
					   uint32_t id,
					   int strand,
					   int len,
					   float diff,
					   int e,
					   uint32_t quant_bit,
					   int k,
					   float fine_min,
				   	   float fine_max,
				   	   float fine_range,
					   mm128_v *p,
					   TfLiteInterpreter* interpreter,
					   TfLiteTensor* input_tensor)
{
	assert(interpreter != NULL);
	assert(input_tensor != NULL);

	assert(len > 0 && e*quant_bit <= 64-RI_HASH_SHIFT);

	uint32_t span = k+e-1;
	
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_quant_bit = (1ULL<<quant_bit)-1;

	uint32_t n_buckets = 1UL<<quant_bit;

	//For reverse complementing
	const int rev_span = 11; //TODO: make this an argument (see rindex.c)
	const int mid_rev = rev_span/2; 

	int sigBufFull = 0;
	uint32_t f_pos, r_pos = len-1-mid_rev, sigBufPos = 0, l_sigpos = 0; //last signal position
	uint32_t f_tmpQuantSignal = 0, r_tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	rh_kv_resize(mm128_t, km, *p, p->n + len);

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));
	
	float* revBuf = (float*)malloc(rev_span*sizeof(float));
	memset(revBuf, 0, rev_span*sizeof(float));
	float* out_results = (float*)malloc(32*sizeof(float));
	memset(out_results, 0, 32*sizeof(float));

    for (f_pos = 0; f_pos < len; ++f_pos) {
        if((f_pos > 0 && fabs(s_values[f_pos] - s_values[l_sigpos]) < diff)) continue;

		l_sigpos = f_pos;

		f_tmpQuantSignal = dynamic_quantize(s_values[f_pos], fine_min, fine_max, fine_range, n_buckets)&mask_quant_bit;
		
		if(f_pos+1 < rev_span) {revBuf[f_pos] = f_tmpQuantSignal; continue;}
		
		revBuf[rev_span-1] = f_tmpQuantSignal;
		
		//Do the inference here
		TfLiteTensorCopyFromBuffer(input_tensor, revBuf, sizeof(float)*rev_span);
		TfLiteInterpreterInvoke(interpreter);
		const TfLiteTensor* output_tensor = TfLiteInterpreterGetOutputTensor(interpreter, 0);
		if(!output_tensor) fprintf(stderr, "Error: output_tensor is NULL\n");
		else {
			// fprintf(stderr, "size: %u %u ", TfLiteTensorByteSize(output_tensor), sizeof(float));
			TfLiteTensorCopyToBuffer(output_tensor, out_results, 32*sizeof(float));
			r_tmpQuantSignal = 0;
			float maxVal = out_results[0];
			for(int revout = 1; revout < 32; revout++) {
				if(out_results[revout] > maxVal) {
					maxVal = out_results[revout];
					r_tmpQuantSignal = revout;
				}
			}
		}

		for(int i = 0; i < rev_span-1; i++) revBuf[i] = revBuf[i+1]; //TODO: Do it more efficiently

		sigBuf[sigBufPos].y = id_shift | (uint32_t)r_pos--<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}

		quantVal = (quantVal<<quant_bit|r_tmpQuantSignal)&mask_events;
		sigBuf[sigBufPos].x = (hash64(quantVal, mask)<<RI_HASH_SHIFT) | span;

		if(!sigBufFull) continue;

		rh_kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);
    }

	free(out_results);
	free(revBuf);
}

void ri_sketch(void *km,
               const float* s_values,
               uint32_t id,
               int strand,
               int len,
               float diff,
               int w,
               int e,
               int n,
               uint32_t quant_bit,
               int k,
               float fine_min,
               float fine_max,
               float fine_range,
               mm128_v *p)
{
	if(w) ri_sketch_min(km, s_values, id, strand, len, diff, w, e, quant_bit, k, fine_min, fine_max, fine_range, p);
	else ri_sketch_reg(km, s_values, id, strand, len, diff, e, quant_bit, k, fine_min, fine_max, fine_range, p);
}

void ri_sketch_rev(void *km,
				   const float* s_values,
				   uint32_t id,
				   int strand,
				   int len,
				   float diff,
				   int w,
				   int e,
				   int n,
				   uint32_t quant_bit,
				   int k,
				   float fine_min,
               	   float fine_max,
               	   float fine_range,
				   mm128_v *p,
				   TfLiteInterpreter* interpreter,
				   TfLiteTensor* input_tensor)
{
	// if(w) ri_sketch_min(km, s_values, id, strand, len, diff, w, e, q, lq, k, p);
	ri_sketch_reg_rev(km, s_values, id, strand, len, diff, e, quant_bit, k, fine_min, fine_max, fine_range, p, interpreter, input_tensor);
}

#ifdef TRAIN_REVERSE
void train_reverse(const float* f_values,
                   const float* r_values,
                   int len,
                   float diff,
                   int rev_span,
                   uint32_t quant_bit,
                   float fine_min,
                   float fine_max,
                   float fine_range) {
							
	const int mid_rev = rev_span/2;

	int sigBufFull = 0;
	uint32_t f_signal = 0, f_tmpQuantSignal = 0, r_signal = 0, r_tmpQuantSignal = 0, mask_quant_bit = (1ULL<<quant_bit)-1;

	uint32_t revBuf[rev_span];
	memset(revBuf, 0, rev_span*sizeof(uint32_t));
	int posBuf[rev_span];
	memset(posBuf, 0, rev_span*sizeof(int));

	uint32_t n_buckets = 1UL<<quant_bit, l_sigpos = 0;

	int f_pos, cnt, r_pos = len-1-mid_rev;
	for (cnt = 0, f_pos = 0; f_pos < len; ++f_pos) {
		if((f_pos > 0 && fabs(f_values[f_pos] - f_values[l_sigpos]) < diff)) continue;

		l_sigpos = f_pos;

		f_tmpQuantSignal = dynamic_quantize(f_values[f_pos], fine_min, fine_max, fine_range, n_buckets)&mask_quant_bit;

		if(posBuf[rev_span-1] == 0) {
			revBuf[cnt] = f_tmpQuantSignal;
			posBuf[cnt++] = f_pos;
			continue;
		}

		revBuf[rev_span-1] = f_tmpQuantSignal;
		posBuf[rev_span-1] = f_pos;

		r_pos = len-1-posBuf[mid_rev];
		r_tmpQuantSignal = dynamic_quantize(r_values[r_pos], fine_min, fine_max, fine_range, n_buckets)&mask_quant_bit;

		for(int i = 0; i < rev_span; ++i){
			if(posBuf[i] == posBuf[mid_rev]-mid_rev+i)
				fprintf(stderr, "\t%d", revBuf[i]+1);
			else
				fprintf(stderr, "\t0");
		}
		for(uint32_t bit = 0; bit < n_buckets; ++bit){
			if(bit == r_tmpQuantSignal) fprintf(stderr, "\t1");
			else fprintf(stderr, "\t0");
		} fprintf(stderr, "\n");

		for(int i = 0; i < rev_span; ++i){
			revBuf[i] = revBuf[i+1];
			posBuf[i] = posBuf[i+1];
		}
	}
}
#endif