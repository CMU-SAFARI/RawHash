#include "rsketch.h"
#include <assert.h>
#include "rh_kvec.h"
#include <math.h>
#include <float.h>

#define TRAIN_MAP

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

void ri_sketch_min(void *km, const float* s_values, uint32_t id, int strand, int len, float diff, int w, int e, int q, int lq, int k, mm128_v *p){

	const uint32_t quant_bit = lq+2;
	assert(len > 0 && (w > 0 && w < 256) && e*quant_bit <= 64);
	
	int j, buf_pos, min_pos;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };

	// int step = 1;//TODO: make this an argument
	// uint32_t span = (k+e-1)*step;
	uint32_t span = k+e-1;
	
	const uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1; //for least sig. 5 bits

	memset(buf, 0xff, w * 16);
	// rh_kv_resize(mm128_t, km, *p, p->n + (len/step)/w);
	rh_kv_resize(mm128_t, km, *p, p->n + len/w);

	int sigBufFull = 0;
	uint32_t i, l, sigBufPos = 0;
	uint32_t l_sigpos = 0; //last signal position
	uint32_t signal = 0, tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

	// for (i = l = buf_pos = min_pos = 0; i < len; i += step) {
    for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
        if(i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < diff) continue;

		l++;
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		l_sigpos = i;
		signal = *((uint32_t*)&s_values[i]);
		tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);

		quantVal = (quantVal<<quant_bit|tmpQuantSignal)&mask_events;

		sigBuf[sigBufPos].y = id_shift | (uint32_t)i<<RI_POS_SHIFT | strand;
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

void ri_sketch_reg(void *km, const float* s_values, uint32_t id, int strand, int len, float diff, int e, int q, int lq, int k, mm128_v *p){

	uint32_t quant_bit = lq+2;
	assert(len > 0 && e*quant_bit <= 64);

	uint32_t span = k+e-1;
	
	uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1;

	int sigBufFull = 0;
	uint32_t i, sigBufPos = 0, l_sigpos = 0; //last signal position
	uint32_t signal = 0, f_tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	rh_kv_resize(mm128_t, km, *p, p->n + len);

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

	int f_pos;
    for (f_pos = 0; f_pos < len; ++f_pos) {
        if((f_pos > 0 && fabs(s_values[f_pos] - s_values[l_sigpos]) < diff)) continue;

		l_sigpos = f_pos;
		signal = *((uint32_t*)&s_values[f_pos]);
		f_tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);

		sigBuf[sigBufPos].y = id_shift | (uint32_t)f_pos<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}

		quantVal = (quantVal<<quant_bit|f_tmpQuantSignal)&mask_events;
		sigBuf[sigBufPos].x = (hash64(quantVal, mask)<<RI_HASH_SHIFT) | span;
		
		if(!sigBufFull) continue;

		rh_kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);
    }
}

void ri_sketch_reg_rev(void *km, const float* s_values, uint32_t id, int strand, int len, float diff, int e, int q, int lq, int k, mm128_v *p, TfLiteInterpreter* interpreter, TfLiteTensor* input_tensor){

	uint32_t quant_bit = lq+2;
	assert(len > 0 && e*quant_bit <= 64);

	uint32_t span = k+e-1;
	
	uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1;

	int sigBufFull = 0;
	uint32_t i, sigBufPos = 0, l_sigpos = 0; //last signal position
	uint32_t signal = 0, f_tmpQuantSignal = 0, r_tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	rh_kv_resize(mm128_t, km, *p, p->n + len);

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

	//For reverse complementing
	const int rev_span = 11; //TODO: make this an argument (see rindex.c)
	const int mid_rev = rev_span/2; 
	float* revBuf = (float*)malloc(rev_span*sizeof(float));
	memset(revBuf, 0, rev_span*sizeof(float));
	float* out_results = (float*)malloc(32*sizeof(float));
	memset(out_results, 0, 32*sizeof(float));

	int f_pos, r_pos = len-1-mid_rev;
    for (f_pos = 0; f_pos < len; ++f_pos) {
        if((f_pos > 0 && fabs(s_values[f_pos] - s_values[l_sigpos]) < diff)) continue;

		l_sigpos = f_pos;
		signal = *((uint32_t*)&s_values[f_pos]);
		f_tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);
		
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

void ri_sketch(void *km, const float* s_values, uint32_t id, int strand, int len, float diff, int w, int e, int n, int q, int lq, int k, mm128_v *p){

	if(w) ri_sketch_min(km, s_values, id, strand, len, diff, w, e, q, lq, k, p);
	else ri_sketch_reg(km, s_values, id, strand, len, diff, e, q, lq, k, p);
}

void ri_sketch_rev(void *km, const float* s_values, uint32_t id, int strand, int len, float diff, int w, int e, int n, int q, int lq, int k, mm128_v *p, TfLiteInterpreter* interpreter, TfLiteTensor* input_tensor){

	// if(w) ri_sketch_min(km, s_values, id, strand, len, diff, w, e, q, lq, k, p);
	ri_sketch_reg_rev(km, s_values, id, strand, len, diff, e, q, lq, k, p, interpreter, input_tensor);
}