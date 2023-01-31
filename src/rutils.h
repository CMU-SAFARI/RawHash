#ifndef RUTILS_H
#define RUTILS_H

#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

extern double ri_realtime0;
extern int ri_verbose;
extern unsigned char seq_nt4_table[256];

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; char **a; } ri_char_v;

void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);

double ri_realtime(void);
double ri_cputime(void);
long ri_peakrss(void);

#ifdef __cplusplus
}
#endif
#endif //RUTILS_H