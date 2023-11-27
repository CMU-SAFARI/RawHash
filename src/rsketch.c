#include "rsketch.h"
#include <assert.h>
#include "rh_kvec.h"
#include <math.h>
#include <string.h>
#include <x86intrin.h>

static inline uint64_t hash64(uint64_t key, uint64_t mask){
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
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

typedef struct { // a simplified version of kdq
    int front, count;
    int a[32];
} tiny_queue_t;

typedef struct { uint64_t x, y, sk_x, sk_y; uint32_t i; uint16_t k; } b304_t;
typedef struct { uint64_t x, y; uint32_t i; uint16_t k; } b176_t;
typedef struct { uint64_t x; uint32_t i; uint16_t k; } b112_t;

// #include <x86intrin.h> //TODO: insert ifdef here to include simd headers

//https://stackoverflow.com/a/21673221
static inline __m256i movemask_inverse(const uint32_t hash_value) {
    __m256i vmask = _mm256_set1_epi32(hash_value);
    const __m256i shuffle = _mm256_setr_epi64x(0x0000000000000000,
      0x0101010101010101, 0x0202020202020202, 0x0303030303030303);
    vmask = _mm256_shuffle_epi8(vmask, shuffle);
    const __m256i bit_mask = _mm256_set1_epi64x(0x7fbfdfeff7fbfdfe);
    vmask = _mm256_or_si256(vmask, bit_mask);
    return _mm256_cmpeq_epi8(vmask, _mm256_set1_epi64x(-1));
}

static inline void calc_blend_simd(__m256i* blndcnt_lsb, __m256i* blndcnt_msb, __m256i* ma, __m256i* mb, uint64_t val, const uint64_t mask, const int bits) {

    (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*ma), (*mb), movemask_inverse(val&mask)));
    uint64_t blendVal = (uint64_t)_mm256_movemask_epi8((*blndcnt_lsb))&mask;
    if(bits > 32){
        (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb), _mm256_blendv_epi8((*ma), (*mb),movemask_inverse((val>>32)&mask)));
        blendVal |= ((uint64_t)_mm256_movemask_epi8((*blndcnt_msb)))<<32;
    }
}

static inline uint64_t calc_blend_rm_simd(__m256i* blndcnt_lsb, __m256i* blndcnt_msb, __m256i* ma, __m256i* mb, uint64_t val, uint64_t remval, const uint64_t mask, const int bits) {

    (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*ma), (*mb), movemask_inverse(val&mask)));
    uint64_t blendVal = (uint64_t)_mm256_movemask_epi8((*blndcnt_lsb))&mask;
    //removal of the oldest item
    (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*mb), (*ma), movemask_inverse(remval&mask)));

    if(bits > 32){
        (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb), _mm256_blendv_epi8((*ma), (*mb),movemask_inverse((val>>32)&mask)));
        blendVal |= ((uint64_t)_mm256_movemask_epi8((*blndcnt_msb)))<<32;

        (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb),_mm256_blendv_epi8((*mb), (*ma), movemask_inverse((remval>>32)&mask)));
    }
    
    return blendVal;
}

void ri_sketch_blend_min(void *km, const float *s_values, uint32_t id, int strand, int len, int w, int e, int n, int q, int lq, int k, mm128_v *p)
{
	assert(len > 0 && (w > 0 && w+k < 8192) && (k > 0 && k <= 28) && (e > 0 && e+k < 8192) && (e <= 56));
	int j, buf_pos, min_pos;
	mm128_t buf[256], min = {UINT64_MAX, UINT64_MAX};

	int step = 1;//TODO: make this an argument
	uint32_t span = (k+e-1)*step;

	// SIMD-related variables
	__m256i ma = _mm256_set1_epi8(1);
	__m256i mb = _mm256_set1_epi8(-1);
	__m256i blndcnt_lsb = _mm256_set1_epi8(0);
	__m256i blndcnt_msb = _mm256_set1_epi8(0);

	const uint32_t quant_bit = lq + 2;
	const uint32_t shift_r = 32 - q;
	const uint64_t id_shift = (uint64_t)id << RI_ID_SHIFT, mask = (1ULL << 32) - 1, mask_events = (1ULL << (q * e)) - 1, mask_l_quant = (1UL << lq) - 1, mask_q_quant = (1UL << quant_bit) - 1; // for least sig. 5 bits
	assert(e * quant_bit <= 64); 
	memset(buf, 0xff, w * 16);
	rh_kv_resize(mm128_t, km, *p, p->n + len / w);

	int sigBufFull = 0;
	int blendBufFull = 0; 

	uint32_t i, l, sigBufPos = 0;
	uint32_t l_sigpos = 0; // last signal position
	uint32_t signal = 0, tmpQuantSignal = 0;
	uint64_t quantVal = 0, hashVal = 0; 
	uint64_t blendVal = 0;
	int blend_pos = 0;
	const uint64_t blendmask =  (1ULL<<28)-1;
	
	mm128_t blendBuf[n];
	memset(blendBuf, 0, n*sizeof(mm128_t));

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e * sizeof(mm128_t));

	for (i = l = buf_pos = min_pos = 0; i < len; i+= step)
	{
		if (i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) continue;
		l++;
		mm128_t info = {UINT64_MAX, UINT64_MAX};
		l_sigpos = i;
		signal = *((uint32_t *)&s_values[i]);
		tmpQuantSignal = signal >> 30 << lq | ((signal >> shift_r) & mask_l_quant);

		quantVal = (quantVal << quant_bit | tmpQuantSignal) & mask_events;
		hashVal = hash64(quantVal, mask);

		sigBuf[sigBufPos].y = id_shift | (uint32_t)i<< RI_POS_SHIFT | strand;
		if (++sigBufPos == e) { sigBufFull = 1; sigBufPos = 0;}
		//sigBuf[sigBufPos].x = hash64(quantVal, mask) << RI_HASH_SHIFT | span; 

		if (!sigBufFull) continue;

		blendBuf[blend_pos].x = hashVal; 
		blendBuf[blend_pos].y = sigBuf[sigBufPos].y;
		if(++blend_pos == n){blendBufFull = 1; blend_pos = 0;}

		if(blendBufFull) {
			blendVal = calc_blend_rm_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, hashVal, blendBuf[blend_pos].x, mask, 32);
			info.x = blendVal<<RI_HASH_SHIFT | span;
			info.y = blendBuf[blend_pos].y;
		} else {
			calc_blend_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, hashVal, mask, 32);
			continue;
		}

		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + e - 1 && min.x != UINT64_MAX)
		{ // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y)
					rh_kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y)
					rh_kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info.x <= min.x)
		{ // a new minimum; then write the old min
			if (l >= w + e && min.x != UINT64_MAX) rh_kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
		}
		else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + e - 1 && min.x != UINT64_MAX) { 
				//fprintf(stderr, "min 1: %u  \n", min);
				rh_kv_push(mm128_t, km, *p, min);
				}
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x)
					min = buf[j], min_pos = j; // >= is important s.t. min is always the closest e-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x)
					min = buf[j], min_pos = j;
			if (l >= w + e - 1 && min.x != UINT64_MAX)
			{									  // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) {
						//fprintf(stderr, "buf[j] 1 : %u \n", buf[j]);
						rh_kv_push(mm128_t, km, *p, buf[j]);
					}
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) {
						//fprintf(stderr, "buf[j] 2 :%u \n", buf[j]);
						rh_kv_push(mm128_t, km, *p, buf[j]);
					}
			}
		}
		if (++buf_pos == w)
			buf_pos = 0;
	}
	if (min.x != UINT64_MAX) {
		//fprintf(stderr, "min 2: %u \n", min);	
		rh_kv_push(mm128_t, km, *p, min);
	}
}



void ri_sketch_blend_hybrid(void *km, const float* s_values, uint32_t id, int strand, int len, int e, int n, int q, int lq, int k, mm128_v *p) {
	//for the hybrid method, choose m > 0 to enable concatenating some blend values. 
	// set m = len, if only concatenated quantized values are wanted. 
	//fprintf(stderr, "we are inside"); 
	uint64_t blendVal = 0;
	const uint64_t blendmask = (1ULL<<28)-1;
	int m = 1; 
	//fprintf(stderr, "m is : %d", m);
    mm128_t blendBuf[n]; 
	memset(blendBuf, 0, n*sizeof(mm128_t)); 
	
    //SIMD-related variables
    __m256i ma = _mm256_set1_epi8(1);
    __m256i mb = _mm256_set1_epi8(-1);
    __m256i blndcnt_lsb = _mm256_set1_epi8(0);
    __m256i blndcnt_msb = _mm256_set1_epi8(0);

	int step = 1;//TODO: make this an argument
	uint32_t span = (k+e-1)*step; //for now single event is considered to span 6 bases.

	const uint32_t quant_bit = lq+2;
	const uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id << RI_ID_SHIFT, mask = (1ULL << 32) - 1, mask_events = (1ULL << (q * e)) - 1, mask_l_quant = (1UL << lq) - 1, mask_q_quant = (1UL << quant_bit) - 1; // for least sig. 5 bits

	int sigBufFull = 0, blendBufFull = 0, blendPos = 0, counter = 0, limit = 4;
	uint32_t i, sigBufPos = 0, l_sigpos = 0; //last signal position
	uint32_t signal = 0, shortQuantVal = 0;
	uint64_t longQuantVal = 0;

	rh_kv_resize(mm128_t, km, *p, p->n + len/step);

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));
	
    for (i = 0; i < len; i += step) {
        if(i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) continue;
		
		l_sigpos = i;
		signal = *((uint32_t*)&s_values[i]);
		
		shortQuantVal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);
		longQuantVal = (longQuantVal << quant_bit | shortQuantVal) & mask_events;// we concatenate these until we have n many 
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		sigBuf[sigBufPos].y = id_shift | (uint32_t)i << RI_POS_SHIFT | strand;
		//blendBuf[blendPos].y = id_shift | (uint32_t)i << RI_POS_SHIFT | strand;
		
		//blendBuf[blendPos].x = hash64(longQuantVal, mask); 
		
		//if(++blendPos == n){blendBufFull = 1; blendPos = 0;};
		if(++sigBufPos == e) {sigBufFull =1; sigBufPos = 0;}
		if(!sigBufFull) continue; 

		blendBuf[blendPos].x = hash64(longQuantVal, mask);
		blendBuf[blendPos].y = sigBuf[sigBufPos].y; 

		if(++blendPos == n){blendBufFull = 1; blendPos = 0;};

		if (blendBufFull)
		{
			blendVal = calc_blend_rm_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, hash64(longQuantVal, mask), blendBuf[blendPos].x, blendmask, 32);
			//sigBuf[sigBufPos].y = blendBuf[blendPos].y;
			
			//if (++sigBufPos == e) { sigBufFull = 1; sigBufPos = 0; }
			//if (!sigBufFull) continue;
			info.y = blendBuf[blendPos].y;
			if(i%m != 0) { // this is the regular rawhash2 implementation
				info.x = hash64(longQuantVal, mask)<<RI_HASH_SHIFT | span;
				rh_kv_push(mm128_t, km, *p, info);
			}
			else {
				info.x = blendVal<<RI_HASH_SHIFT | span;
				rh_kv_push(mm128_t, km, *p, info); // here, we push the blended quantVal 
			}
		}else {
			calc_blend_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, hash64(longQuantVal, mask), blendmask, 32);
		}
	}
}
// void ri_sketch_blend_quant(void *km, const float* s_values, uint32_t id, int strand, int len, int e, int q, int lq, int k, mm128_v *p) {

// 	uint64_t blendVal = 0;
// 	const uint64_t blendmask =  (1ULL<<28)-1;
    
//     //SIMD-related variables
//     __m256i ma = _mm256_set1_epi8(1);
//     __m256i mb = _mm256_set1_epi8(-1);
//     __m256i blndcnt_lsb = _mm256_set1_epi8(0);
//     __m256i blndcnt_msb = _mm256_set1_epi8(0);

// 	int step = 1;//TODO: make this an argument
// 	uint32_t span = (k+e-1)*step; //for now single event is considered to span 6 bases.

// 	const uint32_t quant_bit = lq+2; 
// 	const uint32_t shift_r = 32-q;
// 	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(q*e))-1, mask_l_quant = (1UL<<lq)-1; //for least sig. 5 bits
	
// 	int sigBufFull = 0;
// 	uint32_t i, sigBufPos = 0, l_sigpos = 0; //last signal position
// 	uint32_t signal = 0, tmpQuantSignal = 0;
// 	uint64_t quantVal = 0, hashVal = 0;

// 	rh_kv_resize(mm128_t, km, *p, p->n + len/step);

// 	mm128_t sigBuf[e];
// 	memset(sigBuf, 0, e*sizeof(mm128_t));

//     for (i = 0; i < len; i += step) {
//         if(i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) continue;

// 		l_sigpos = i;
// 		signal = *((uint32_t*)&s_values[i]);
// 		tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);
// 		//tmpQuantSignal = (signal >> shift_r);

// 		mm128_t info = { UINT64_MAX, UINT64_MAX };

// 		quantVal = (quantVal<<quant_bit|tmpQuantSignal)&mask_events;

// 		sigBuf[sigBufPos].x = quantVal; 
// 		sigBuf[sigBufPos].y = id_shift | (uint32_t)i<<RI_POS_SHIFT | strand;
// 		if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}

// 		if(!sigBufFull) continue;
		
// 		if(sigBufFull){
// 			blendVal = calc_blend_rm_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, quantVal, sigBuf[sigBufPos].x, blendmask, 32);  
// 			info.x = blendVal<<RI_HASH_SHIFT | span;
// 			info.y = sigBuf[sigBufPos].y;
// 			rh_kv_push(mm128_t, km, *p, info);
// 		}else{
// 			calc_blend_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, quantVal, blendmask, 32);
// 		}
//     }
// }

void ri_sketch_blend_quant(void *km, const float* s_values, uint32_t id, int strand, int len, int e, int n, int q, int lq, int k, mm128_v *p) {
 
	int blendConcat = 0, f=1; 
	uint64_t blendVal = 0;
	const uint64_t blendmask = (1ULL<<28)-1;
    mm128_t blendBuf[n]; 
	memset(blendBuf, 0, n*sizeof(mm128_t)); 
	
    //SIMD-related variables
    __m256i ma = _mm256_set1_epi8(1);
    __m256i mb = _mm256_set1_epi8(-1);
    __m256i blndcnt_lsb = _mm256_set1_epi8(0);
    __m256i blndcnt_msb = _mm256_set1_epi8(0);

	int step = 1;//TODO: make this an argument
	uint32_t span = (k+e-1)*step; //for now single event is considered to span 6 bases.

	const uint32_t quant_bit = lq+2;
	const uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id << RI_ID_SHIFT, mask = (1ULL << 32) - 1, mask_events = (1ULL << (q * e)) - 1, mask_l_quant = (1UL << lq) - 1, mask_q_quant = (1UL << quant_bit) - 1; // for least sig. 5 bits

	int sigBufFull = 0, blendBufFull = 0, blendPos = 0, counter = 0, limit = 4;
	uint32_t i, sigBufPos = 0, l_sigpos = 0; //last signal position
	uint32_t signal = 0, shortQuantVal = 0;
	uint64_t longQuantVal = 0, concatBlend = 0;

	rh_kv_resize(mm128_t, km, *p, p->n + len/step);

	mm128_t sigBuf[e];
	
	memset(sigBuf, 0, e*sizeof(mm128_t));
    for (i = 0; i < len; i += step) {
        if(i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) continue;
		
		l_sigpos = i;
		signal = *((uint32_t*)&s_values[i]);
		shortQuantVal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);

		mm128_t info = { UINT64_MAX, UINT64_MAX };
		// check the size of shortQuantVal here
		longQuantVal = (longQuantVal << quant_bit | hash64(shortQuantVal, mask)) & mask_events;// we concatenate these until we have n many 
		sigBuf[sigBufPos].y = id_shift | (uint32_t)i << RI_POS_SHIFT | strand;
		if(++sigBufPos == n) {sigBufFull = 1; sigBufPos = 0;}

		if(!sigBufFull) continue; 

		blendBuf[blendPos].x = longQuantVal; 
		blendBuf[blendPos].y = sigBuf[sigBufPos].y;

		if(++blendPos == n){blendBufFull = 1; blendPos = 0;};

		if (blendBufFull)
		{
			//if(f) {
			blendVal = calc_blend_rm_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, longQuantVal, blendBuf[blendPos].x, mask_q_quant, 32);
			info.x = blendVal<<RI_HASH_SHIFT |span; 
			info.y = blendBuf[blendPos].y; 
			//}
			/*else {
				blendConcat++; 
				blendVal = calc_blend_rm_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, shortQuantVal, blendBuf[blendPos].x, mask_q_quant, 32);
				concatBlend = (concatBlend << quant_bit | hash64(blendVal, mask)) & mask_events;
				if(blendConcat == 3) {
					info.x = concatBlend << RI_HASH_SHIFT | span;
					info.y = blendBuf[blendPos].y;
				}
			}*/
			rh_kv_push(mm128_t, km, *p, info);
		}else {
			calc_blend_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, longQuantVal, mask_q_quant, 32);
		}
	}
}

void ri_sketch_blend(void *km, const float* s_values, uint32_t id, int strand, int len, int e, int n, int q, int lq, int k, mm128_v *p){

	uint64_t blendVal = 0;
    int blend_pos = 0;
	bool buffull = false;
	const uint64_t blendmask =  (1ULL<<28)-1;

	mm128_t blndBuf[n];
	memset(blndBuf, 0, n*sizeof(mm128_t));
    
    //SIMD-related variables
    __m256i ma = _mm256_set1_epi8(1);
    __m256i mb = _mm256_set1_epi8(-1);
    __m256i blndcnt_lsb = _mm256_set1_epi8(0);
    __m256i blndcnt_msb = _mm256_set1_epi8(0);

	int step = 1;//TODO: make this an argument
	uint32_t span = (k+e-1)*step; //for now single event is considered to span 6 bases.

	const uint32_t quant_bit = lq+2; 
	const uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1; //for least sig. 5 bits
	
	int sigBufFull = 0;
	uint32_t i, sigBufPos = 0, l_sigpos = 0; //last signal position
	uint32_t signal = 0, tmpQuantSignal = 0;
	uint64_t quantVal = 0, hashVal = 0;

	rh_kv_resize(mm128_t, km, *p, p->n + len/step);

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

    for (i = 0; i < len; i += step) {
        if(i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) continue;

		l_sigpos = i;
		signal = *((uint32_t*)&s_values[i]);
		tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);

		mm128_t info = { UINT64_MAX, UINT64_MAX };

		quantVal = (quantVal<<quant_bit|tmpQuantSignal)&mask_events;
		hashVal = hash64(quantVal, mask); 

		sigBuf[sigBufPos].y = id_shift | (uint32_t)i<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}
		// sigBuf[sigBufPos].x = hashVal<<RI_HASH_SHIFT | span;

		if(!sigBufFull) continue;

		blndBuf[blend_pos].x = hash64(quantVal, mask);
		blndBuf[blend_pos].y = sigBuf[sigBufPos].y;

		if(++blend_pos == n) {buffull = true; blend_pos = 0;}

		if(buffull){
			blendVal = calc_blend_rm_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, hashVal, blndBuf[blend_pos].x, blendmask, 32);
			info.x = blendVal<<RI_HASH_SHIFT | span;
			info.y = blndBuf[blend_pos].y;
			rh_kv_push(mm128_t, km, *p, info);
		}else{
			calc_blend_simd(&blndcnt_lsb, &blndcnt_msb, &ma, &mb, hashVal, blendmask, 32);
		}
    }
}

void ri_sketch_min(void *km, const float* s_values, uint32_t id, int strand, int len, int w, int e, int q, int lq, int k, mm128_v *p){
	
	int j, buf_pos, min_pos;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };

	// int step = 1;//TODO: make this an argument
	// uint32_t span = (k+e-1)*step;
	uint32_t span = k+e-1;
	const uint32_t quant_bit = lq+2;
	const uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1; //for least sig. 5 bits

	assert(len > 0 && (w > 0 && w < 256) && e*quant_bit <= 64);

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
        if(i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF) continue;

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

void ri_sketch_reg(void *km, const float* s_values, uint32_t id, int strand, int len, int e, int q, int lq, int k, mm128_v *p){

	// int step = 1;//TODO: make this an argument
	// uint32_t span = (k+e-1)*step;
	uint32_t span = k+e-1;
	
	uint32_t quant_bit = lq+2; 
	uint32_t shift_r = 32-q;
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_l_quant = (1UL<<lq)-1; //for least sig. 5 bits

	int sigBufFull = 0;
	uint32_t i, sigBufPos = 0, l_sigpos = 0; //last signal position
	uint32_t signal = 0, tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	assert(len > 0 && e*quant_bit <= 64);

	// rh_kv_resize(mm128_t, km, *p, p->n + len/step);
	rh_kv_resize(mm128_t, km, *p, p->n + len);

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

	// for (i = 0; i < len; i += step) {
    for (i = 0; i < len; ++i) {
        if((i > 0 && fabs(s_values[i] - s_values[l_sigpos]) < LAST_SIG_DIFF)) continue;

		l_sigpos = i;
		signal = *((uint32_t*)&s_values[i]);
		tmpQuantSignal = signal>>30<<lq | ((signal>>shift_r)&mask_l_quant);

		sigBuf[sigBufPos].y = id_shift | (uint32_t)i<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}

		quantVal = (quantVal<<quant_bit|tmpQuantSignal)&mask_events;
		sigBuf[sigBufPos].x = (hash64(quantVal, mask)<<RI_HASH_SHIFT) | span;
		
		if(!sigBufFull) continue;

		rh_kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);

		// quantVal = (quantVal>>quant_bit<<quant_bit|(tmpQuantSignal-1))&mask_events;
		// sigBuf[sigBufPos].x = hash64(quantVal, mask)<<RI_HASH_SHIFT | span;
		// rh_kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);

		// quantVal = (quantVal>>quant_bit<<quant_bit|(tmpQuantSignal+1))&mask_events;
		// sigBuf[sigBufPos].x = hash64(quantVal, mask)<<RI_HASH_SHIFT | span;
		// rh_kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);
    }
}

void ri_sketch(void *km, const float* s_values, uint32_t id, int strand, int len, int w, int e, int n, int q, int lq, int k, mm128_v *p){

	assert(e >= 1 && e < 10);
	//int f = 1; // decides whether or not to use the hybrid method or Blend_quant
	//if(w & n) ri_sketch_blend_min(km, s_values, id, strand, len, w, e, n, q, lq, k, p);
	//else if(n) ri_sketch_blend_quant(km, s_values, id, strand, len, e, n, q, lq, k, p);
	//q = 11;
	//if(w) ri_sketch_min(km, s_values, id, strand, len, w, e, q, lq, k, p);
	//ri_sketch_blend(km, s_values, id, strand, len, e, n, q, lq, k, p);
	ri_sketch_blend_quant(km, s_values, id, strand, len, e, n, q, lq, k, p);
	//if(n) ri_sketch_blend_hybrid(km, s_values, id, strand, len, e, n, q, lq, k, p);
	//else if(n) ri_sketch_blend_hybrid(km, s_values, id, strand, len, e, n, q, lq, k, p);
	//ri_sketch_reg(km, s_values, id, strand, len, e, q, lq, k, p);
}