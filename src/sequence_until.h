#ifndef SEQUENCE_UNTIL_H
#define SEQUENCE_UNTIL_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Finds the outlier point (i.e., the point that is further away from all other points) among n many points in a m-dimensional space
 * For relative abundance estimation, this translates to the distance of an estimation that is furthest away from all other most recent n many estimations
 *
 * @param x		x[point][dim]. List of $dim-dimensional points in space.
 * 				For relative abundance estimation this would be x[most recent $n estimations][$m many genomes to estimate the relative abundance].
 * @param n		number of points in the space.
 * 				For relative abundance estimations, this would be the number of most recent relative abundance estimations.
 * @param m		length of the dimension. Dimension of each point ($x[point]) should be the same for all 0 <= point < $n.
 * 				For relative abundance estimation, this would be the number of genomes whose abundance is estimated across all other genomes.
 * 
 * @return		the distance of $x[i] to all other points such that x[i] is the furthest point among all other x[y] where 0 <= y < n and y != i; 0 if all points are at the same position
 */
float find_outlier(const float** x, const uint32_t n, const uint32_t m);

#ifdef __cplusplus
}
#endif
#endif //SEQUENCE_UNTIL_H