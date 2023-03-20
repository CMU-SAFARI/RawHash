#include "rsketch.h"
#include <assert.h>
#include "kvec.h"
#include <math.h>

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

// static inline uint64_t MurmurHash3(uint64_t key, uint64_t mask) {
//   key = (key^(key >> 33));
//   key = (key*0xff51afd7ed558ccd);
//   key = (key^(key >> 33));
//   key = (key*0xc4ceb9fe1a85ec53);
//   key = (key^(key >> 33));

//   return key&mask;
// }

//Quantize float in range [min,max] to 10-bit unsigned integer
static inline uint32_t quantize_float_uint32(const float x, const float min, const float max) {
	uint32_t a = (uint32_t)(((x-min)/(max-min))*127.0f);
	if(a>127) return 127;
	if(a<0) return 0;

	return a;
}

// #include <x86intrin.h> //TODO: insert ifdef here to include simd headers

//https://stackoverflow.com/a/21673221
// static inline __m256i movemask_inverse(const uint32_t hash_value) {
//     __m256i vmask = _mm256_set1_epi32(hash_value);
//     const __m256i shuffle = _mm256_setr_epi64x(0x0000000000000000,
//       0x0101010101010101, 0x0202020202020202, 0x0303030303030303);
//     vmask = _mm256_shuffle_epi8(vmask, shuffle);
//     const __m256i bit_mask = _mm256_set1_epi64x(0x7fbfdfeff7fbfdfe);
//     vmask = _mm256_or_si256(vmask, bit_mask);
//     return _mm256_cmpeq_epi8(vmask, _mm256_set1_epi64x(-1));
// }

// static inline void calc_blend_simd(__m256i* blndcnt_lsb, __m256i* blndcnt_msb, __m256i* ma, __m256i* mb, uint64_t val, const uint64_t mask, const int bits) {

//     (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*ma), (*mb), movemask_inverse(val&mask)));
//     uint64_t blendVal = (uint64_t)_mm256_movemask_epi8((*blndcnt_lsb))&mask;
//     if(bits > 32){
//         (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb), _mm256_blendv_epi8((*ma), (*mb),movemask_inverse((val>>32)&mask)));
//         blendVal |= ((uint64_t)_mm256_movemask_epi8((*blndcnt_msb)))<<32;
//     }
// }

// static inline uint64_t calc_blend_rm_simd(__m256i* blndcnt_lsb, __m256i* blndcnt_msb, __m256i* ma, __m256i* mb, uint64_t val, uint64_t remval, const uint64_t mask, const int bits) {

//     (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*ma), (*mb), movemask_inverse(val&mask)));
//     uint64_t blendVal = (uint64_t)_mm256_movemask_epi8((*blndcnt_lsb))&mask;
//     //removal of the oldest item
//     (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*mb), (*ma), movemask_inverse(remval&mask)));

//     if(bits > 32){
//         (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb), _mm256_blendv_epi8((*ma), (*mb),movemask_inverse((val>>32)&mask)));
//         blendVal |= ((uint64_t)_mm256_movemask_epi8((*blndcnt_msb)))<<32;

//         (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb),_mm256_blendv_epi8((*mb), (*ma), movemask_inverse((remval>>32)&mask)));
//     }
    
//     return blendVal;
// }

// void ri_sketch_blend_min(void *km, const float* s_values, uint32_t id, int strand, int len, int w, int e, int n, int q, int lq, int k, mm128_v *p){	
// }

// void ri_sketch_blend(void *km, const float* s_values, uint32_t id, int strand, int len, int e, int n, int q, int lq, int k, mm128_v *p){

// 	uint64_t blendVal = 0;
//     int blend_pos = 0;
// 	bool buffull = false;
// 	const uint64_t blendmask =  (1ULL<<28)-1;

// 	mm128_t blndBuf[n];
// 	memset(blndBuf, 0, n*sizeof(mm128_t));
    
//     //SIMD-related variables
//     __m256i ma = _mm256_set1_epi8(1);
//     __m256i mb = _mm256_set1_epi8(-1);
//     __m256i blndcnt_lsb = _mm256_set1_epi8(0);
//     __m256i blndcnt_msb = _mm256_set1_epi8(0);

// 	int step = 1;//TODO: make this an argument
// 	uint32_t span = (k+e-1)*step; //for now single event is considered to span 6 bases.

// 	const uint32_t quant_bit = lq+2; 
// 	const uint32_t shift_r = 32-q;
// 	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1; //for least sig. 5 bits
	
// 	int sigBufFull = 0;
// 	uint32_t i, sigBufPos = 0, l_sigpos = 0; //last signal position
// 	uint32_t signal = 0, tmpQuantSignal = 0;
// 	uint64_t quantVal = 0, hashVal = 0;

// 	kv_resize(mm128_t, km, *p, p->n + len/step);

// 	mm128_t sigBuf[e];
// 	memset(sigBuf, 0, e*sizeof(mm128_t));

//     for (i = 0; i < len; i += step) {
//         if(i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) continue;

// 		l_sigpos = i;
// 		signal = *((uint32_t*)&s_values[i]);
// 		tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);

// 		mm128_t info = { UINT64_MAX, UINT64_MAX };

// 		quantVal = (quantVal<<quant_bit|tmpQuantSignal)&mask_events;
// 		hashVal = hash64(quantVal, mask);

// 		sigBuf[sigBufPos].y = id_shift | (uint32_t)i<<RI_POS_SHIFT | strand;
// 		if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}
// 		// sigBuf[sigBufPos].x = hashVal<<RI_HASH_SHIFT | span;

// 		if(!sigBufFull) continue;

// 		blndBuf[blend_pos].x = hash64(quantVal, mask);
// 		blndBuf[blend_pos].y = sigBuf[sigBufPos].y;

// 		if(++blend_pos == n) {buffull = true; blend_pos = 0;}

// 		if(buffull){
// 			blendVal = calc_blend_rm_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, hashVal, blndBuf[blend_pos].x, blendmask, 32);
// 			info.x = blendVal<<RI_HASH_SHIFT | span;
// 			info.y = blndBuf[blend_pos].y;
// 			kv_push(mm128_t, km, *p, info);
// 		}else{
// 			calc_blend_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, hashVal, blendmask, 32);
// 		}
//     }
// }

void ri_sketch_min(void *km, const float* s_values, uint32_t id, int strand, int len, int w, int e, int q, int lq, int k, mm128_v *p){
	
	int j, buf_pos, min_pos;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };

	int step = 1;//TODO: make this an argument
	uint32_t span = 6+e-1; //for now single event is considered to span 6 bases.

	assert(len > 0 && (w > 0 && w < 256) && (e > 0 && e <= 10));
	memset(buf, 0xff, w * 16);
	kv_resize(mm128_t, km, *p, p->n + (len/step)/w);
	
	const uint32_t quant_bit = lq+2; 
	const uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1; //for least sig. 5 bits

	int sigBufFull = 0;
	uint32_t i, l, sigBufPos = 0;
	uint32_t l_sigpos = 0; //last signal position
	uint32_t signal = 0, tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

    for (i = l = buf_pos = min_pos = 0; i < len; i += step) {
        if(i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) continue;

		l++;
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		l_sigpos = i;
		signal = *((uint32_t*)&s_values[i]);
		tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);
		// signal = quantize_float_uint32(s_values[i], -2.75f, 2.75f);
		// tmpQuantSignal = signal&mask_quant;

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
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + e && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + e - 1 && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest e-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + e - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
    }
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, km, *p, min);

}

void ri_sketch_reg(void *km, const float* s_values, uint32_t id, int strand, int len, int e, int q, int lq, int k, mm128_v *p){

	int step = 1;//TODO: make this an argument
	uint32_t span = (k+e-1)*step;
	
	uint32_t quant_bit = lq+2; 
	uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1; //for least sig. 5 bits

	int sigBufFull = 0;
	uint32_t i, sigBufPos = 0, l_sigpos = 0; //last signal position
	uint32_t signal = 0, tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	kv_resize(mm128_t, km, *p, p->n + len/step);

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

    for (i = 0; i < len; i += step) {
        if((i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) || s_values[i] == RI_MASK_SIGNAL) continue;

		l_sigpos = i;
		signal = *((uint32_t*)&s_values[i]);
		tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);
		// signal = quantize_float_uint32(s_values[i], -2.75f, 2.75f);
		// tmpQuantSignal = signal&mask_quant;

		sigBuf[sigBufPos].y = id_shift | (uint32_t)i<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}

		quantVal = (quantVal<<quant_bit|tmpQuantSignal)&mask_events;
		sigBuf[sigBufPos].x = (hash64(quantVal, mask)<<RI_HASH_SHIFT) | span;
		
		if(!sigBufFull) continue;

		kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);

		// quantVal = (quantVal>>quant_bit<<quant_bit|(tmpQuantSignal-1))&mask_events;
		// sigBuf[sigBufPos].x = hash64(quantVal, mask)<<RI_HASH_SHIFT | span;
		// kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);

		// quantVal = (quantVal>>quant_bit<<quant_bit|(tmpQuantSignal+1))&mask_events;
		// sigBuf[sigBufPos].x = hash64(quantVal, mask)<<RI_HASH_SHIFT | span;
		// kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);
    }
}

void ri_sketch(void *km, const float* s_values, uint32_t id, int strand, int len, int w, int e, int n, int q, int lq, int k, mm128_v *p){

	assert(e > 1 && e < 10);

	// if(w & n) ri_sketch_blend_min(km, s_values, id, strand, len, w, e, n, q, lq, k, p);
	// else if(n) ri_sketch_blend(km, s_values, id, strand, len, e, n, q, lq, k, p);
	if(w) ri_sketch_min(km, s_values, id, strand, len, w, e, q, lq, k, p);
	else ri_sketch_reg(km, s_values, id, strand, len, e, q, lq, k, p);
}