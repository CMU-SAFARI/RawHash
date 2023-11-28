#include <iostream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <immintrin.h>
#include <iomanip>
#include <limits>
#include "dtw.h"

using namespace std;

#define DISTANCE(A, B) abs((A)-(B))
#define DISTANCE16(A, B) _mm512_abs_ps(_mm512_sub_ps((A), (B)))

//check the maximum available float SIMD width
#if defined(__AVX512F__)
	#define F32SIMD_WIDTH 16
#elif defined(__AVX__)
	#define F32SIMD_WIDTH 8
#elif defined(__SSE__)
	#define F32SIMD_WIDTH 4
#else
	#define F32SIMD_WIDTH 1
#endif

//#define DEBUG

#define PRINTDP {							\
	for(int i = 0; i < dp.size(); i++){		\
			std::cout << std::fixed;		\
			std::cout << std::setprecision(2);	\
			std::cout << dp[i]<< " ";		\
		}									\
		std::cout << std::endl;				\
	}

float DTW_global(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element){
	std::vector<float> dp(a_length, 1e10);
    dp[0] = DISTANCE(a_values[0], b_values[0]);

    for(uint32_t j = 1; j < a_length; j++){
        dp[j] = dp[j-1] + DISTANCE(a_values[j], b_values[0]);
    }

	for(uint32_t i = 1; i < b_length; i++){
		float old_left = dp[0];
		dp[0] = dp[0]+DISTANCE(a_values[0], b_values[i]);
		for(uint32_t j = 1; j < a_length; j++){
			float top = dp[j-1];
			float left = dp[j];
			float topleft = old_left;
			float center = std::min(
							std::min(top, left),
							topleft
						) + DISTANCE(a_values[j], b_values[i]);
			dp[j] = center;
			old_left = left;
		}
	}
	if(exclude_last_element){
		return dp[a_length-1] - DISTANCE(a_values[a_length-1], b_values[b_length-1]);
	}
	else{
		return dp[a_length-1];
	}
}

//512 bit SIMD vector print macro
#define PRINT512(VEC) {					\
	std::cout << #VEC << ": ";			\
	for(int i = 0; i < 16; i++){		\
		std::cout << (VEC)[i] << " ";	\
	}									\
	std::cout << std::endl;				\
}


float DTW_global_slow(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element){
	assert(a_length > 0 && b_length > 0);
    
    vector<vector<float>> dp(a_length, vector<float>(b_length));
    
    dp[0][0] = DISTANCE(a_values[0], b_values[0]);
    
    for(uint32_t i = 1; i < a_length; i++){
        dp[i][0] = dp[i-1][0] + DISTANCE(a_values[i], b_values[0]);
    }
    for(uint32_t j = 1; j < b_length; j++){
        dp[0][j] = dp[0][j-1] + DISTANCE(a_values[0], b_values[j]);
    }

    for(uint32_t i = 1; i < a_length; i++){
        for(uint32_t j = 1; j < b_length; j++){
            float best_in = min(min(dp[i-1][j], dp[i][j-1]), dp[i-1][j-1]);
            dp[i][j] = best_in + DISTANCE(a_values[i], b_values[j]);
        }
    }

	if(exclude_last_element){
		return dp[a_length-1][b_length-1] - DISTANCE(a_values[a_length-1], b_values[b_length-1]);
	}
	else{
		return dp[a_length-1][b_length-1];
	}
}

float DTW_global_diagonalbanded(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, int band_radius, bool exclude_last_element){
	assert(a_length > 0 && b_length > 0);
	assert(a_length < std::numeric_limits<int>::max());
	assert(b_length < std::numeric_limits<int>::max());
	assert(band_radius >= 0);

	const int band_width = band_radius*2+1;
	vector<float> dp(band_width, 1e10);
	const int dp_center = band_radius;

	float prev = 0.0f;
	for(int row_offset = 0; row_offset <= min(band_radius, (int)b_length-1); row_offset++){
		int j = row_offset;
		float cur = prev + DISTANCE(a_values[0], b_values[j]);
		dp[dp_center+row_offset] = cur;
		prev = cur;
	}

	for(int i = 1; i < (int)a_length; i++){
		const int center_row = i;
		const int row_offset_start = max(-(int)band_radius, -center_row);
		const int row_offset_end = min(band_radius, (int)b_length - center_row - 1);

		float top = 1e10;
		for(int row_offset = row_offset_start; row_offset <= row_offset_end; row_offset++){
			int j = center_row + row_offset;
			float topleft = dp[dp_center+row_offset]; //for the first few iterations this might go oob, make sure to init to 1e10
			float left = row_offset==band_radius?1e10:dp[dp_center+row_offset+1];
			float center = min(
							min(top, left),
							topleft
						) + DISTANCE(a_values[i], b_values[j]);
			dp[dp_center+row_offset] = center;
			top = center;
		}
	}

	const int center_row = a_length-1;
	const int row_offset_start = max(-(int)band_radius, -center_row);
	const int row_offset_end = min(band_radius, (int)b_length - center_row - 1);
	const int desired_j = b_length-1;
	const int desired_row_offset = desired_j - center_row;

	if(row_offset_start > desired_row_offset ||
	   row_offset_end   < desired_row_offset     ){
		//diagonal band does not cover the desired element
		//i.e., the bottom right corner of the matrix
		return 1e10;
	}
	else{
		float res = dp[dp_center+desired_row_offset];
		if(exclude_last_element){
			return res - DISTANCE(a_values[a_length-1], b_values[b_length-1]);
		}
		else{
			return res;
		}
	}
}

float DTW_global_slantedbanded(const float* a_values, uint32_t a_length, const float* b_values, uint32_t b_length, int band_radius, bool exclude_last_element){
	assert(a_length > 0 && b_length > 0);
	assert(a_length < std::numeric_limits<int>::max());
	assert(b_length < std::numeric_limits<int>::max());
	assert(band_radius >= 0);

	#ifdef DEBUG
		vector<vector<float>> debug_matrix(a_length, vector<float>(b_length, 0.0f));
	#endif

	//make sure a is the longer sequence
	if(a_length < b_length){
		//swap
		const float* tmp_values = a_values;
		const uint32_t tmp_length = a_length;
		a_values = b_values;
		a_length = b_length;
		b_values = tmp_values;
		b_length = tmp_length;
	}

	const int band_width = band_radius*2+1;
	vector<float> dp(band_width, 1e10);
	const int dp_center = band_radius;

	float prev = 0.0f;
	for(int row_offset = 0; row_offset <= min(band_radius, (int)b_length-1); row_offset++){
		int j = row_offset;
		float cur = prev + DISTANCE(a_values[0], b_values[j]);
		dp[dp_center+row_offset] = cur;
		prev = cur;

		#ifdef DEBUG
			debug_matrix[0][j] = cur;
		#endif
	}

	int center_row = 0;
	for(int i = 1; i < (int)a_length; i++){
		//increment center_row to follow the slope from top left to bottom right
		int next_row = center_row+1;
		//float next_slope = next_row/(float)i;
		//float target_slope = b_length/(float)a_length;	
		//floating point logic reformulated and implemented as integers for performance
		int64_t next_slope = next_row*(int64_t)a_length;
		int64_t target_slope = b_length*(int64_t)i;
		bool increment_center_row = false;
		if(next_slope <= target_slope){
			center_row++;
			increment_center_row = true;
		}

		const int row_offset_start = max(-(int)band_radius, -center_row);
		const int row_offset_end = min(band_radius, (int)b_length - center_row - 1);

		float top = 1e10;
		float topleft = increment_center_row && (center_row + row_offset_start > 0)?dp[dp_center+row_offset_start]:1e10;
		for(int row_offset = row_offset_start; row_offset <= row_offset_end; row_offset++){
			int j = center_row + row_offset;
			//float topleft = dp[dp_center+row_offset]; //for the first few iterations this might go oob, make sure to init to 1e10
			float left;
			if(increment_center_row){
				left = row_offset==band_radius?1e10:dp[dp_center+row_offset+1];
			}
			else{
				left = dp[dp_center+row_offset];
			}
			float center = min(
							min(top, left),
							topleft
						) + DISTANCE(a_values[i], b_values[j]);
			dp[dp_center+row_offset] = center;
			top = center;
			topleft = left;

			#ifdef DEBUG
				if(i < (int)a_length && j < (int)b_length)
					debug_matrix[i][j] = center;
			#endif
		}
	}

	#ifdef DEBUG
		std::cout << "debug_matrix for DTW_global_slantedbanded:" << std::endl;
		for(int j = 0; j < b_length; j++){
			for(int i = 0; i < a_length; i++){
				std::cout << std::fixed;
				std::cout << std::setprecision(2);
				std::cout << debug_matrix[i][j] << "\t";
			}
			std::cout << std::endl;
		}
	#endif

	const int desired_j = b_length-1;
	const int desired_row_offset = desired_j - center_row;
	float res = dp[dp_center+desired_row_offset];
	
	if(exclude_last_element){
		return res - DISTANCE(a_values[a_length-1], b_values[b_length-1]);
	}
	else{
		return res;
	}
}

float DTW_global_slantedbanded_antidiagonalwise(const float* __restrict__ a_values, uint32_t a_length, const float* __restrict__ b_values, uint32_t b_length, int band_radius, bool exclude_last_element){
	assert(a_length > 0 && b_length > 0);
	assert(a_length < std::numeric_limits<int>::max());
	assert(b_length < std::numeric_limits<int>::max());
	assert(band_radius >= 0);

	#ifdef DEBUG
		vector<vector<float>> debug_matrix(a_length, vector<float>(b_length, 0.0f));
	#endif

	//make sure a is the longer sequence
	if(a_length < b_length){
		//swap
		const float* tmp_values = a_values;
		const uint32_t tmp_length = a_length;
		a_values = b_values;
		a_length = b_length;
		b_values = tmp_values;
		b_length = tmp_length;
	}

	//math.ceil((1-m/n)*band_radius)
	//math.ceil((n-m)/n*band_radius)
	//math.ceil((n-m)*band_radius/n)
	//((n-m)*band_radius + n - 1) / n
	int extra_band_radius_due_to_slanting =	((a_length-b_length)*band_radius + a_length - 1) / a_length;

	band_radius += extra_band_radius_due_to_slanting;
	const int primary_antidiagonal_length = band_radius + (band_radius % 2 == 0 ? 1 : 0);
    const int secondary_antidiagonal_length = band_radius + (band_radius % 2 == 1 ? 1 : 0);
	const bool primary_larger = primary_antidiagonal_length > secondary_antidiagonal_length;

	const int dpsize = max(primary_antidiagonal_length, secondary_antidiagonal_length);
	float dp_storage[dpsize*3]; //maybe this needs to be malloc's or new'd instead
	float* __restrict__ dp0 = dp_storage;
	float* __restrict__ dp1 = dp_storage + dpsize;
	float* __restrict__ dp2 = dp_storage + dpsize*2;
	for(int i = 0; i < dpsize; i++){ //initialize to 1e10 to simplify (literal) corner cases
		dp0[i] = 1e10;
		dp1[i] = 1e10;
		dp2[i] = 1e10;
	}

	int center_row = 0;
	{ //iteration 0
		int iteration = 0;
		int center_column = iteration;
		int antidiagonal_start_i = center_column + primary_antidiagonal_length/2;
		int antidiagonal_start_j = center_row - primary_antidiagonal_length/2;

		//this entire loop could be replaced by initializing the top left corner to DISTANCE(a_values[0], b_values[0])
		//for(int antidiagonal_offset = 0; antidiagonal_offset < primary_antidiagonal_length; antidiagonal_offset++){
		{
			int antidiagonal_offset = primary_antidiagonal_length/2;
			int i = antidiagonal_start_i - antidiagonal_offset;
			int j = antidiagonal_start_j + antidiagonal_offset;
			if(j >= 0 && j < (int)b_length && i >= 0 && i < (int)a_length){
				if(primary_larger)
					dp2[antidiagonal_offset] = DISTANCE(a_values[i], b_values[j]);
				else
					dp2[antidiagonal_offset+1] = DISTANCE(a_values[i], b_values[j]);
				
				#ifdef DEBUG
					if(primary_larger)
						debug_matrix[i][j] = dp2[antidiagonal_offset];
					else
						debug_matrix[i][j] = dp2[antidiagonal_offset+1];
				#endif
			}
		}
		float *tmp = dp0;
		dp0 = dp1;
		dp1 = dp2;
		dp2 = tmp;
	}

	bool previous_increment_center_row = false;
	for(int iteration = 1; (uint32_t)iteration < a_length; iteration++){
		int center_column = iteration;
		int next_row = center_row+1;
		int64_t next_slope = next_row*(int64_t)a_length;
		int64_t target_slope = b_length*(int64_t)center_column;
		bool increment_center_row = false;
		if(next_slope <= target_slope){
			center_row++;
			increment_center_row = true;
		}

		if(increment_center_row){
			int antidiagonal_start_i = center_column + secondary_antidiagonal_length/2 - 1;
			int antidiagonal_start_j = center_row - secondary_antidiagonal_length/2;

			int antidiagonal_offset_start = max(max(0, antidiagonal_start_i-(int)a_length+1), -antidiagonal_start_j);
			int antidiagonal_offset_end = min(min(secondary_antidiagonal_length, (int)antidiagonal_start_i+1), (int)b_length-antidiagonal_start_j);
			if(primary_larger){
				for(int antidiagonal_offset = antidiagonal_offset_start; antidiagonal_offset < antidiagonal_offset_end; antidiagonal_offset++){
					int i = antidiagonal_start_i - antidiagonal_offset;
					int j = antidiagonal_start_j + antidiagonal_offset;

					//calculate dp entry
					float top = dp1[antidiagonal_offset];
					float topleft = dp0[antidiagonal_offset];
					float left = dp1[antidiagonal_offset+1];
					float center = min(
									min(top, left),
									topleft
								) + DISTANCE(a_values[i], b_values[j]);
					dp2[antidiagonal_offset] = center;
					
					#ifdef DEBUG
						debug_matrix[i][j] = center;
					#endif
				}
			}
			else{
				for(int antidiagonal_offset = antidiagonal_offset_start; antidiagonal_offset < antidiagonal_offset_end; antidiagonal_offset++){
					int i = antidiagonal_start_i - antidiagonal_offset;
					int j = antidiagonal_start_j + antidiagonal_offset;

					bool is_first = antidiagonal_offset==0;
					bool is_last = antidiagonal_offset==secondary_antidiagonal_length-1;
					float top = is_first?1e10:dp1[antidiagonal_offset];
					float topleft = is_first && !previous_increment_center_row ?
						1e10 : dp0[antidiagonal_offset]; //when the secondary is larger, topleft is only available if dp0 was a secondary one
					float left = is_last?1e10:dp1[antidiagonal_offset+1];
					float center = min(
									min(top, left),
									topleft
								) + DISTANCE(a_values[i], b_values[j]);
					dp2[antidiagonal_offset] = center;
					
					#ifdef DEBUG
						debug_matrix[i][j] = center;
					#endif
				}
			}

			float *tmp = dp0;
			dp0 = dp1;
			dp1 = dp2;
			dp2 = tmp;
		}

		int antidiagonal_start_i = center_column + primary_antidiagonal_length/2;
		int antidiagonal_start_j = center_row - primary_antidiagonal_length/2;

		int antidiagonal_offset_start = max(max(0, antidiagonal_start_i-(int)a_length+1), -antidiagonal_start_j);
		int antidiagonal_offset_end = min(min(primary_antidiagonal_length, (int)antidiagonal_start_i+1), (int)b_length-antidiagonal_start_j);

		if(primary_larger){
			for(int antidiagonal_offset = antidiagonal_offset_start; antidiagonal_offset < antidiagonal_offset_end; antidiagonal_offset++){
				int i = antidiagonal_start_i - antidiagonal_offset;
				int j = antidiagonal_start_j + antidiagonal_offset;
			
				float top, topleft, left;
				if(increment_center_row){
					bool is_first = antidiagonal_offset==0;
					bool is_last = antidiagonal_offset==primary_antidiagonal_length-1;
					top = is_first?1e10:dp1[antidiagonal_offset-1]; //the first element of a primary antidiagonal never has a top
					topleft = dp0[antidiagonal_offset]; //all elements have a topleft when going down
					left = is_last?1e10:dp1[antidiagonal_offset]; //the last element of a primary antidiagonal does not have a left when going down
				}
				else{
					bool is_first = antidiagonal_offset==0;
					//bool is_last = antidiagonal_offset==primary_antidiagonal_length-1;
					top = is_first?1e10:dp1[antidiagonal_offset-1]; //the first element of a primary antidiagonal never has a top
					topleft = is_first?1e10:dp0[antidiagonal_offset-1]; //the first element of a primary antidiagonal does not have a topleft when not going down
					left = dp1[antidiagonal_offset]; //all elements have a left when not going down
				}
			
				float center = min(
								min(top, left),
								topleft
							) + DISTANCE(a_values[i], b_values[j]);
				dp2[antidiagonal_offset] = center;
									
				#ifdef DEBUG
					debug_matrix[i][j] = center;
				#endif
			}
		}
		else{
			for(int antidiagonal_offset = antidiagonal_offset_start; antidiagonal_offset < antidiagonal_offset_end; antidiagonal_offset++){
				int i = antidiagonal_start_i - antidiagonal_offset;
				int j = antidiagonal_start_j + antidiagonal_offset;
				
				//to simplify the code, accesses to primary anti diagonal will be starting at dp0[1] instead of dp0[0]
				float top, topleft, left;
				if(increment_center_row){
					top = dp1[antidiagonal_offset]; //the first element of a primary antidiagonal always has a top when going down
					topleft = dp0[antidiagonal_offset+1]; //all elements have a topleft when going down. +1 due to the simplification (see comment above)
					left = dp1[antidiagonal_offset+1]; //all elements have a left when going down
				}
				else{
					bool is_first = antidiagonal_offset==0;
					//bool is_last = antidiagonal_offset==primary_antidiagonal_length-1;
					top = is_first?1e10:dp1[antidiagonal_offset]; //the first element of a primary antidiagonal never has a top. No -1 due to the simplification (see comment above)
					topleft = is_first && !previous_increment_center_row ? 
						1e10:dp0[antidiagonal_offset]; //the first element of a primary antidiagonal does not have a topleft when not previously going down
					left = dp1[antidiagonal_offset+1]; //all elements have a left when not going down
				}
				
				float center = min(
								min(top, left),
								topleft
							) + DISTANCE(a_values[i], b_values[j]);
				dp2[antidiagonal_offset+1] = center; //+1 due to the simplification (see comment above)
					
				#ifdef DEBUG
					debug_matrix[i][j] = center;
				#endif
			}
		}

		float *tmp = dp0;
		dp0 = dp1;
		dp1 = dp2;
		dp2 = tmp;
		previous_increment_center_row = increment_center_row;
	}

	#ifdef DEBUG
		cout << "debug_matrix for DTW_global_slantedbanded_antidiagonalwise:" << endl;
		//set numberic precision to 2 decimal places
		cout << fixed << setprecision(2);
		for(int j = 0; j < (int)b_length; j++){
			for(int i = 0; i < (int)a_length; i++){
				cout << debug_matrix[i][j] << "\t";
			}
			cout << endl;
		}
	#endif

	float res;
	if(primary_larger){
		res = dp1[primary_antidiagonal_length/2];
	}
	else{
		res = dp1[primary_antidiagonal_length/2+1]; //+1 due to the simplification (see comment above)
	}

	if(exclude_last_element){
		return res - DISTANCE(a_values[a_length-1], b_values[b_length-1]);
	}
	else{
		return res;
	}
}

/*
 * a is aligned fully (globally) to the best matching substring of b (i.e., b is not aligned globally)
 * a is typically the shorter sequence
 */
float DTW_semiglobal(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element){
	assert(a_length >= 1);
	assert(b_length >= 1);
	
	std::vector<float> dp(a_length, 1e10);
	float best = 1e10;
	for(uint32_t i = 0; i < b_length; i++){
		float old_left = dp[0];
		dp[0] = DISTANCE(a_values[0], b_values[i]);

		for(uint32_t j = 1; j < a_length; j++){
			float top = dp[j-1];
			float left = dp[j];
			float topleft = old_left;
			float center = std::min(
							std::min(top, left),
							topleft
						) + DISTANCE(a_values[j], b_values[i]);
			dp[j] = center;
			old_left = left;
		}
		best = std::min(best, dp[a_length-1]);
	}
	return best;
}

/*
 * a is aligned fully (globally) to the best matching substring of b (i.e., b is not aligned globally)
 * a is typically the shorter sequence
 */
float DTW_semiglobal_slow(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element){
	assert(a_length > 0 && b_length > 0);
    
    vector<vector<float>> dp(a_length, vector<float>(b_length));
    
    dp[0][0] = DISTANCE(a_values[0], b_values[0]);
    
    for(uint32_t i = 1; i < a_length; i++){
        dp[i][0] = dp[i-1][0] + DISTANCE(a_values[i], b_values[0]);
    }
    for(uint32_t j = 1; j < b_length; j++){
        dp[0][j] = DISTANCE(a_values[0], b_values[j]);
    }

    for(uint32_t i = 1; i < a_length; i++){
        for(uint32_t j = 1; j < b_length; j++){
            float best_in = min(min(dp[i-1][j], dp[i][j-1]), dp[i-1][j-1]);
            dp[i][j] = best_in + DISTANCE(a_values[i], b_values[j]);
        }
    }

    float best = dp[a_length-1][0];
	uint32_t best_j = 0;
    for(uint32_t j = 1; j < b_length; j++){
        //best = min(best, dp[a_length-1][j]);
		if(dp[a_length-1][j] < best){
			best = dp[a_length-1][j];
			best_j = j;
		}
    }

	if(exclude_last_element){
		return best - DISTANCE(a_values[a_length-1], b_values[best_j]);
	}
	else{
	    return best;
	}
}

dtw_result DTW_global_tb(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element){
	assert(a_length > 0 && b_length > 0);
    
    vector<vector<float>> dp(a_length, vector<float>(b_length));
    
    dp[0][0] = DISTANCE(a_values[0], b_values[0]);
    
    for(uint32_t i = 1; i < a_length; i++){
        dp[i][0] = dp[i-1][0] + DISTANCE(a_values[i], b_values[0]);
    }
    for(uint32_t j = 1; j < b_length; j++){
        dp[0][j] = dp[0][j-1] + DISTANCE(a_values[0], b_values[j]);
    }

    for(uint32_t i = 1; i < a_length; i++){
        for(uint32_t j = 1; j < b_length; j++){
            float best_in = min(min(dp[i-1][j], dp[i][j-1]), dp[i-1][j-1]);
            dp[i][j] = best_in + DISTANCE(a_values[i], b_values[j]);
        }
    }

	uint32_t i = a_length-1;
	uint32_t j = b_length-1;
	vector<alignment_element> reverse_alignment;
	reverse_alignment.push_back(
		(alignment_element){
			(position_pair){i, j},
			DISTANCE(a_values[i], b_values[j])
		}
	);
	while(i > 0 || j > 0){
		if(i==0){
			j--;
		}
		else if(j==0){
			i--;
		}
		else{
			float left = dp[i-1][j];
			float top = dp[i][j-1];
			float topleft = dp[i-1][j-1];

			if(left < min(top, topleft)){
				i--;
			}
			else if(top < min(left, topleft)){
				j--;
			}
			else{
				i--;
				j--;
			}
		}

		alignment_element ae;
		ae.position.i = i;
		ae.position.j = j;
		ae.difference = DISTANCE(a_values[i], b_values[j]);
		reverse_alignment.push_back(ae);
	}

	vector<alignment_element> alignment;
    std::copy(reverse_alignment.rbegin(), reverse_alignment.rend(), std::back_inserter(alignment));

	if(exclude_last_element){
		alignment.pop_back();
		float score = dp[a_length-1][b_length-1] - DISTANCE(a_values[a_length-1], b_values[b_length-1]);
		return (dtw_result){score, alignment};
	}
	else{
		return (dtw_result){dp[a_length-1][b_length-1], alignment};
	}
}

dtw_result DTW_semiglobal_tb(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element){
	assert(a_length > 0 && b_length > 0);
    
    vector<vector<float>> dp(a_length, vector<float>(b_length));
    
    dp[0][0] = DISTANCE(a_values[0], b_values[0]);
    
    for(uint32_t i = 1; i < a_length; i++){
        dp[i][0] = dp[i-1][0] + DISTANCE(a_values[i], b_values[0]);
    }
    for(uint32_t j = 1; j < b_length; j++){
        dp[0][j] = DISTANCE(a_values[0], b_values[j]);
    }

    for(uint32_t i = 1; i < a_length; i++){
        for(uint32_t j = 1; j < b_length; j++){
            float best_in = min(min(dp[i-1][j], dp[i][j-1]), dp[i-1][j-1]);
            dp[i][j] = best_in + DISTANCE(a_values[i], b_values[j]);
        }
    }

    float best = dp[a_length-1][0];
	uint32_t best_j = 0;
    for(uint32_t j = 1; j < b_length; j++){
		if(dp[a_length-1][j] < best){
			best = dp[a_length-1][j];
			best_j = j;
		}
	}

	uint32_t i = a_length-1;
	uint32_t j = best_j;
	vector<alignment_element> reverse_alignment;
	reverse_alignment.push_back(
		(alignment_element){
			(position_pair){i, j},
			DISTANCE(a_values[i], b_values[j])
		}
	);
	
	while(i > 0){
		if(i==0){
			j--;
		}
		else if(j==0){
			i--;
		}
		else{
			float left = dp[i-1][j];
			float top = dp[i][j-1];
			float topleft = dp[i-1][j-1];

			if(left < min(top, topleft)){
				i--;
			}
			else if(top < min(left, topleft)){
				j--;
			}
			else{
				i--;
				j--;
			}
		}

		alignment_element ae;
		ae.position.i = i;
		ae.position.j = j;
		ae.difference = DISTANCE(a_values[i], b_values[j]);
		reverse_alignment.push_back(ae);
	}

	
	vector<alignment_element> alignment;
    std::copy(reverse_alignment.rbegin(), reverse_alignment.rend(), std::back_inserter(alignment));

	if(exclude_last_element){
		const alignment_element last = alignment.back();
		float score = dp[a_length-1][best_j] - last.difference;
		alignment.pop_back();
		return (dtw_result){score, alignment};
	}
	else{
	    return (dtw_result){dp[a_length-1][best_j], alignment};
	}
}
