#include "rutils.h"
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

unsigned char seq_nt4_table[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
	
double ri_realtime(void)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

double ri_cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long ri_peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

char* strsep(char** stringp, const char* delim) {
    char* start = *stringp;
    char* p;

    p = (start != NULL) ? strpbrk(start, delim) : NULL;

    if (p == NULL) {
        *stringp = NULL;
    } else {
        *p = '\0';
        *stringp = p + 1;
    }

    return start;
}

uint32_t rev_complement(uint32_t x, const short k) {
    uint32_t y = 0;
    for (short i = 0; i < k; i++) {
        y = (y<<2) | ((x&3)^3);
        x >>= 2;
    }
    return y;
}

int compare_value_index_pairs(const void* a, const void* b) {
    float diff = ((ri_porem_t*)a)->pore_val - ((ri_porem_t*)b)->pore_val;
    return (diff < 0) ? -1 : (diff > 0);
}

ri_porem_t* create_sorted_pairs(const float* pore_vals, uint32_t n_pore_vals, const short k) {
    double sum = 0, sum2 = 0, mean, std_dev;

    // First pass: Calculate mean and standard deviation
    for (uint32_t i = 0; i < n_pore_vals; i++) {
        sum += pore_vals[i];
        sum2 += pore_vals[i] * pore_vals[i];
    }
    mean = sum / n_pore_vals;
    std_dev = sqrt(sum2 / n_pore_vals - mean * mean);

    // Allocate memory for pairs
    ri_porem_t* pairs = (ri_porem_t*)malloc(n_pore_vals * sizeof(ri_porem_t));

    // Second pass: Normalize values and create pairs
    for (uint32_t i = 0; i < n_pore_vals; i++) {
        pairs[i].pore_val = (pore_vals[i] - mean) / std_dev; // Normalization
        pairs[i].ind = i;
        pairs[i].rev_ind = rev_complement(i, k);
    }

    // Sort pairs based on normalized pore_val
    qsort(pairs, n_pore_vals, sizeof(ri_porem_t), compare_value_index_pairs);
    return pairs;
}

void load_pore(const char* fpore, const short k, const short lev_col, float** pore_vals, ri_porem_t** pore_inds){
	FILE* fp = fopen(fpore, "r");
	if(fp == NULL){
		fprintf(stderr, "Error: cannot open file %s\n", fpore);
		return;
	}

	(*pore_vals) = (float*)malloc(sizeof(float) * (1U<<(2*k)));
	char line[1024];
	char* token;
	int i = 0;
	while(fgets(line, sizeof(line), fp) != NULL){
		if(!strncmp(line, "kmer", 4)) continue;
    	char* rest = line;
    	int j = 0;
		while((token = strsep(&rest, "\t")) != NULL){
			if(j++ == lev_col){
				float value;
				if (sscanf(token, "%f", &value) == 1) {
					(*pore_vals)[i] = value;
				} else {
					fprintf(stderr, "Error: cannot convert '%s' to float\n", token);
					free(pore_vals); pore_vals = NULL;
					fclose(fp);
					return;
				}
				break;
			}
		}
		i++;
	}
	fclose(fp);

    (*pore_inds) = create_sorted_pairs(*pore_vals, (1U<<(2*k)), k);
}

// #define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8) 

// #define sort_key_64(x) (x)
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8)

KSORT_INIT_GENERIC(uint32_t)
KSORT_INIT_GENERIC(uint64_t)