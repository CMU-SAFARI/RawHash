#include "sequence_until.h"

//Find the outlier n-dimensional vector of float numbers in range [0,1] in a set of n-dimensional vectors
float find_outlier(const float** x, const uint32_t n, const uint32_t m) {
	uint32_t outlier = 0;
	float max_dist = 0.0f;
	for(uint32_t i=0; i<m; i++) {
		float dist = 0.0f;
		for(uint32_t j=0; j<n; j++) {
			dist += (x[i][j]-x[outlier][j])*(x[i][j]-x[outlier][j]);
		}
		if(dist>max_dist) {
			max_dist = dist;
			outlier = i;
		}
	}
	return max_dist;
}