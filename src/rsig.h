#ifndef RSIG_H
#define RSIG_H

#include "rutils.h"
#ifndef NHDF5RH
#include "hdf5_tools.hpp"
#endif
#ifndef NPOD5RH
#include "pod5_format/c_api.h"
#endif
#ifndef NSLOW5RH
#include <slow5/slow5.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ri_sig_s{
	uint32_t rid, l_sig; //read id and length of the signal values
	char *name; //name of the read

	// float dig, ran, offset; //digitalisation, range, offset
	float* sig; //signal values of a read
} ri_sig_t;

typedef struct ri_sig_file_s {
	// gzFile fp;
	// kseq_t *ks;
	// ri_sig_t s;
	int num_read; //Number of reads
	int cur_read; //Number of processed reads by RawHash (shows the id of the next read to process)
	
	//HDF5-related
	char** raw_path; //List of paths to raw values
	char** ch_path; //List of paths to channel
	#ifndef NHDF5RH
	hdf5_tools::File* fp; //FAST5 file pointer
	#endif
	
	//POD5-related
	unsigned long int pod5_batch_count;
	unsigned long int pod5_batch;
	unsigned long int pod5_row_count;
	unsigned long int pod5_row;
	#ifndef NPOD5RH
	Pod5ReadRecordBatch_t* batch;
	
	Pod5FileReader_t* pp; //POD5 file pointer
	#endif
} ri_sig_file_t;

/**
 * Closes the signal file.
 *
 * @param fp	a struct that includes the file pointer to a signal file
 * 				Allocated memories in the struct are destroyed here.
 */
void ri_sig_close(ri_sig_file_t *fp);

/**
 * Opens the signal file (e.g., a FAST5 file)
 *
 * @param fn	path to the signal file
 * 
 * @return		a struct that includes the file pointer to the opened signal file
 * 				Returned struct (and its variables) is allocated in this function.
 */
ri_sig_file_t *open_sig(const char *fn);

/**
 * Opens all the signal files (e.g., FAST5 files)
 *
 * @param n		number of files
 * @param fn	list of paths to the files
 * 
 * @return		List of structs that include the file pointers to each opened signal file
 * 				Returned structs (and their variables) are allocated in this function.
 */
ri_sig_file_t **open_sigs(int n, const char **fn);

/**
 * Converts the sequence into its expected event values
 *
 * @param str		sequence to convert to its expected event values
 * @param len		length of the $str
 * @param pore_vals	expected event values for each possible k-mer
 * @param pore_kmer	k-mer size of each event
 * @param strand	directin of strand. 1 is forward, 0 is reverse direction.
 * @param s_len		length of $s_values
 * @param s_values	expected event values of each k-mer in $str
 */
void ri_seq_to_sig(const char *str, int len, const float* pore_vals, const int pore_kmer, const int strand, uint32_t* s_len, float* s_values);

/**
 * Reads the entire signal values of the next read from a file
 *
 * @param fp	file pointer to the signal file (i.e., either FAST5 or SLOW5)
 * @param s		attribute of the read and the signal values.
 * 				$s->name = name of the read
 * 				$s->sig = signal values
 * 				$s->l_sig = number of signal values
 */
void ri_read_sig(ri_sig_file_t* fp, ri_sig_t* s);

/**
 * Recursively find all files that ends with "fast5" under input directory const char *A
 *
 * @param A			path to a directory where the fast5 files are searched
 * @param fnames	list of fast5 or pod5 files
 * 					fnames->a = List of file names
 * 					fnames->n = Number of fast5 files
 */
void find_sfiles(const char *A, ri_char_v *fnames);

#ifdef __cplusplus
}
#endif
#endif //RSIG_H