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

#ifdef __linux__
#include <sys/resource.h>
#include <sys/time.h>
void liftrlimit()
{
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
}
#else
void liftrlimit() {}
#endif

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
    float diff = ((ri_porei_t*)a)->pore_val - ((ri_porei_t*)b)->pore_val;
    return (diff < 0)?-1:(diff > 0);
}

ri_porei_t* create_sorted_pairs(const ri_pore_t* pore) {
    double sum = 0, sum2 = 0, mean, std_dev;

    // First pass: Calculate mean and standard deviation
    for (uint32_t i = 0; i < pore->n_pore_vals; i++) {
        sum += pore->pore_vals[i];
        sum2 += pore->pore_vals[i] * pore->pore_vals[i];
    }
    mean = sum / pore->n_pore_vals;
    std_dev = sqrt(sum2 / pore->n_pore_vals - mean * mean);

    // Allocate memory for pairs
    ri_porei_t* pairs = (ri_porei_t*)malloc(pore->n_pore_vals * sizeof(ri_porei_t));

    // Second pass: Normalize values and create pairs
    for (uint32_t i = 0; i < pore->n_pore_vals; i++) {
        pairs[i].pore_val = (pore->pore_vals[i] - mean) / std_dev; // Normalization
        pairs[i].ind = i;
        pairs[i].rev_ind = rev_complement(i, pore->k);
    }

    // Sort pairs based on normalized pore_val
    qsort(pairs, pore->n_pore_vals, sizeof(ri_porei_t), compare_value_index_pairs);
    return pairs;
}

void load_pore(const char* fpore, const short k, const short lev_col, ri_pore_t* pore){
	FILE* fp = fopen(fpore, "r");
	if(fp == NULL){
		fprintf(stderr, "Error: cannot open file %s\n", fpore);
		return;
	}

	pore->pore_vals = (float*)malloc(sizeof(float) * (1U<<(2*k)));
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
					pore->pore_vals[i] = value;
				} else {
					fprintf(stderr, "Error: cannot convert '%s' to float\n", token);
					free(pore->pore_vals); pore->pore_vals = NULL;
					fclose(fp);
					return;
				}
				break;
			}
		}
		i++;
	}
	fclose(fp);

    pore->n_pore_vals = 1U<<(2*k);
    pore->k = k;
    pore->pore_inds = create_sorted_pairs(pore);
}

// Function to swap two elements
void swap(float *a, float *b) {
    float t = *a; *a = *b; *b = t;
}

// Partition function for QuickSelect
int partition(float arr[], int low, int high) {
    float pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

// QuickSelect function to find Median
float quickSelect(float arr[], int low, int high, int k) {
    if (low == high) return arr[low];

    int pivotIndex = partition(arr, low, high);

    if (k == pivotIndex) return arr[k];
    else if (k < pivotIndex) return quickSelect(arr, low, pivotIndex - 1, k);
    else return quickSelect(arr, pivotIndex + 1, high, k);
}

// Function to find median using QuickSelect
float findMedian(float arr[], int n) {
    if (n % 2 != 0) return quickSelect(arr, 0, n - 1, n / 2);
    else return 0.5 * (quickSelect(arr, 0, n - 1, n / 2 - 1) + quickSelect(arr, 0, n - 1, n / 2));
}

float calculateMAD(float* arr, int n, float median) {
    float* deviations = (float*)malloc(n * sizeof(float));
    for (int i = 0; i < n; ++i) {
        deviations[i] = fabs(arr[i] - median);
    }
    float mad = findMedian(deviations, n);
    free(deviations);
    return mad;
}

// #define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8) 

// #define sort_key_64(x) (x)
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8)

KSORT_INIT_GENERIC(uint32_t)
KSORT_INIT_GENERIC(uint64_t)