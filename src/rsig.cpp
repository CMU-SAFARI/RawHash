#include "rsig.h"
#include "kvec.h"
#include <math.h>
#include <sys/stat.h>
#include <dirent.h>

void ri_seq_to_sig(const char *str, int len, const float* pore_vals, const int k, const int strand, uint32_t* s_len, float* s_values){

	int i, j, l, pos, n = 0;
	// uint64_t shift1 = 2 * (k - 1);
	uint64_t mask = (1ULL<<2*k) - 1, kmer = 0;
	double mean = 0, std_dev = 0, sum = 0, sum2 = 0, curval = 0;

	for (i = l = j = n = 0; i < len; ++i) {
		if(strand) pos = len - i -1;
		else pos = i;
		int c = seq_nt4_table[(uint8_t)str[pos]];
		if (c < 4) { // not an ambiguous base
			if(!strand) kmer = (kmer << 2 | c) & mask;    // forward k-mer
			// else kmer = (kmer >> 2) | (3ULL^c) << shift1; // reverse k-mer
			//TODO: this is currently based on the ordering in the original ordering in the sigmap implementation. Change later to above
			else kmer = ((kmer << 2) | (3ULL^c)) & mask; // reverse k-mer
		}else
      		kmer = (kmer << 2) & mask; //TODO: This is not the best approach. We are basically inserting 00 (A?) to kmer whenever c >= 4. Mask it instead

		if(i+1 < k) continue;

		curval = pore_vals[kmer];
		s_values[j++] = curval;
		sum += curval;
		sum2 += curval*curval;
	}

	mean = sum/j;
	std_dev = sqrt(sum2/j - (mean)*(mean));

	for(i = 0; i < j; ++i)
		s_values[i] = (s_values[i]-mean)/std_dev;

	*s_len = j;
}

static inline ri_sig_file_t *ri_sig_open_fast5(const char *fn)
{
	ri_sig_file_t *fp;

	hdf5_tools::File* fast5_file = new hdf5_tools::File();
	fast5_file->open(std::string(fn));
	// gzFile f;
	// f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (!fast5_file->is_open()) return 0;

	fp = (ri_sig_file_t*)calloc(1, sizeof(ri_sig_file_t));
	fp->pp = NULL;
	fp->fp = fast5_file;

	bool is_single = false;
	std::vector<std::string> fast5_file_groups = fast5_file->list_group("/");
	fp->num_read = fast5_file_groups.size();
	fp->ch_path = (char**)calloc(fp->num_read, sizeof(char*));
	fp->raw_path = (char**)calloc(fp->num_read, sizeof(char*));

	for (std::string &group : fast5_file_groups) {
		if (group == "Raw") {
			is_single = true;
			break;
		}
	}

	std::string raw_path;
	std::string ch_path;
	int i = 0;

	if (is_single) {
		ch_path = "/UniqueGlobalKey/channel_id";
		for (std::string &read : fast5_file->list_group("/Raw/Reads")) {
			raw_path = "/Raw/Reads/" + read;
			if(i == fp->num_read){
				fprintf(stderr, "ERROR: More reads than previously predicted (%d). Stopped reading the reads here.\n", fp->num_read);
				break;
			}
			fp->ch_path[i] = strdup(ch_path.c_str());
			fp->raw_path[i++] = strdup(raw_path.c_str());
		}
	} else {
		for (std::string &read : fast5_file_groups) {
			raw_path = "/" + read + "/Raw";
			ch_path = "/" + read + "/channel_id";
			fp->ch_path[i] = strdup(ch_path.c_str());
			fp->raw_path[i++] = strdup(raw_path.c_str());
		}
	}

	fp->num_read = i;
	fp->cur_read = 0;
	return fp;
}

static inline ri_sig_file_t *ri_sig_open_pod5(const char *fn){
	pod5_init();

	ri_sig_file_t *fp;

    // Open the file ready for walking:
    Pod5FileReader_t* pod5_file = pod5_open_file(fn);

	if (!pod5_file) {
		fprintf(stderr, "ERROR: Failed to open file (%s) %s\n", fn, pod5_get_error_string());
		return 0;
    }

	long unsigned int batch_count = 0;
    if (pod5_get_read_batch_count(&batch_count, pod5_file) != POD5_OK) {
        fprintf(stderr, "Failed to query batch count: %s\n", pod5_get_error_string());
		pod5_close_and_free_reader(pod5_file);
		return 0;
    }

	fp = (ri_sig_file_t*)calloc(1, sizeof(ri_sig_file_t));
	fp->pp = pod5_file;
	fp->fp = NULL;

	fp->num_read = batch_count; //num_read is the number of batches for pod5
	fp->cur_read = 0; //cur_read is cur_batch for pod5

	Pod5ReadRecordBatch_t* batch = NULL;
	if(pod5_get_read_batch(&batch, pod5_file, fp->cur_read) != POD5_OK){
		fprintf(stderr, "Failed to get batch: %s\n", pod5_get_error_string());
		pod5_close_and_free_reader(pod5_file);
		free(fp);
		return 0;
	}

	long unsigned int batch_row_count = 0;
	if(pod5_get_read_batch_row_count(&batch_row_count, batch) != POD5_OK) {
		fprintf(stderr, "Failed to get batch row count\n");
		pod5_close_and_free_reader(pod5_file);
		pod5_free_read_batch(batch);
		free(fp);
		return 0;
	}

	fp->pod5_row_count = batch_row_count;
	fp->pod5_row = 0;
	
	fp->batch = batch;
	fp->num_read = batch_count;

	return fp;
}

static inline ri_sig_file_t *ri_sig_open_slow5(const char *fn){
	//TODO: complete here
	return 0;
}

ri_sig_file_t *ri_sig_open(const char *fn){
	if (strstr(fn, ".fast5")) {
		return ri_sig_open_fast5(fn);
	} else if (strstr(fn, ".pod5")) {
		return ri_sig_open_pod5(fn);
	} else if (strstr(fn, ".slow5") || strstr(fn, ".blow5")) {
		// TODO: complete here
	}

	return 0;
}

void ri_sig_close(ri_sig_file_t *fp)
{
	if(!fp) return;
	// gzclose(fp->fp);

	if(fp->fp){
		fp->fp->close();
		for(int i = 0; i < fp->num_read; ++i){
			if(fp->ch_path[i])free(fp->ch_path[i]);
			if(fp->raw_path[i])free(fp->raw_path[i]);
		}
	}else if(fp->pp){
		pod5_close_and_free_reader(fp->pp);
	}
	free(fp);
}

ri_sig_file_t *open_sig(const char *fn) //TODO: make this a part of the pipeline. Sequntially reading from many FAST5 files creates an overhead
{
	ri_sig_file_t *fp;
	fp = (ri_sig_file_t*)calloc(1,sizeof(ri_sig_file_t));
	if ((fp = ri_sig_open(fn)) == 0) {
		fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn, strerror(errno));
		ri_sig_close(fp);
		return 0;
	}
	return fp;
}

ri_sig_file_t **open_sigs(int n, const char **fn) //TODO: make this a part of the pipeline. Sequntially reading from many FAST5 files creates an overhead
{
	ri_sig_file_t **fp;
	int i, j;
	fp = (ri_sig_file_t**)calloc(n, sizeof(ri_sig_file_t*));
	for (i = 0; i < n; ++i) {
		if ((fp[i] = ri_sig_open(fn[i])) == 0) {
			fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i], strerror(errno));
			for (j = 0; j < i; ++j) ri_sig_close(fp[j]);
			free(fp);
			return 0;
		}
	}
	return fp;
}

//Check if the input const char* A is a directory
//Generated by GitHub Copilot
int is_dir(const char *A)
{
	struct stat st;
	if (stat(A, &st) == -1) return 0;
	return S_ISDIR(st.st_mode);
}

//Recursively find all files that ends with "fast5" or "pod5" under input directory const char *A
//Generated by GitHub Copilot
void find_sfiles(const char *A, ri_char_v *fnames)
{
	if (!is_dir(A)) {
		//TODO: Add slow5 here later
		if (strstr(A, ".fast5") || strstr(A, ".pod5")) {
			char** cur_fname;
			kv_pushp(char*, 0, *fnames, &cur_fname);
			(*cur_fname) = strdup(A);
		}
		return;
	}

	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir(A)) != NULL) {
		while ((ent = readdir(dir)) != NULL) {
			char *tmp = (char*)malloc(strlen(A) + strlen(ent->d_name) + 2);
			sprintf(tmp, "%s/%s", A, ent->d_name);
			if (is_dir(tmp)) {
				if (strcmp(ent->d_name, ".") && strcmp(ent->d_name, ".."))
					find_sfiles(tmp, fnames);
			} else {
				if (strstr(ent->d_name, ".fast5") || strstr(ent->d_name, ".pod5")) {
					char** cur_fname;
					kv_pushp(char*, 0, *fnames, &cur_fname);
					(*cur_fname) = strdup(tmp);
				}
			}
			free(tmp);
		}
		closedir(dir);
	}
}

static inline void ri_read_sig_fast5(ri_sig_file_t* fp, ri_sig_t* s){

	if(fp->cur_read >= fp->num_read) return;
	
	s->name = 0;
	for (auto a : fp->fp->get_attr_map(fp->raw_path[fp->cur_read])) {
		if (a.first == "read_id") {
			s->name = strdup(a.second.c_str());
		}
	}

	assert(s->name);

	float dig = 0, ran = 0, offset = 0;

	// float digitisation = 0, range = 0, offset = 0;
	for (auto a : fp->fp->get_attr_map(fp->ch_path[fp->cur_read])) {
		if (a.first == "channel_number") {
		// channel_idx = atoi(a.second.c_str()) - 1;
		} else if (a.first == "digitisation") {
			dig = atof(a.second.c_str());
		} else if (a.first == "range") {
			ran = atof(a.second.c_str());
		} else if (a.first == "offset") {
			offset = atof(a.second.c_str());
		}
	}

	std::string sig_path = std::string(fp->raw_path[fp->cur_read]) + "/Signal";
	std::vector<float> sig;
	fp->fp->read(sig_path, sig);
	// convert to pA
	uint32_t l_sig = 0;
	float scale = ran/dig;
	float pa = 0;
	for (size_t i = 0; i < sig.size(); i++) {
		pa = (sig[i]+offset)*scale;
		if (pa > 30.0f && pa < 200.0f) {
			sig[l_sig++] = pa;
		}
	}

	s->sig = (float*)calloc(l_sig, sizeof(float));
	s->l_sig = l_sig;
	std::copy(sig.begin(), sig.begin() + l_sig, s->sig);
	fp->cur_read++;
}

static inline void ri_read_sig_pod5(ri_sig_file_t* fp, ri_sig_t* s){

	if(fp->cur_read >= fp->num_read) return;
	if(fp->pod5_row >= fp->pod5_row_count) return;

	pod5_init();

    // Open the file ready for walking:
    Pod5FileReader_t* pod5_file = fp->pp;

    if (!pod5_file) {
        fprintf(stderr, "Failed to open file %s\n", pod5_get_error_string());
		return;
    }

	uint16_t read_table_version = 0;
	ReadBatchRowInfo_t read_data;
	if (pod5_get_read_batch_row_info_data(fp->batch, fp->pod5_row, READ_BATCH_ROW_INFO_VERSION, &read_data, &read_table_version) != POD5_OK) {
		fprintf(stderr, "Failed to get read %lu\n", fp->pod5_row);
		return;
	}

	//Retrieve global information for the run
	RunInfoDictData_t* run_info_data;
	if (pod5_get_run_info(fp->batch, read_data.run_info, &run_info_data) != POD5_OK) {
		fprintf(stderr, "Failed to get Run Info %lu %s\n", fp->pod5_row, pod5_get_error_string());
		pod5_release_run_info(run_info_data);
		return;
	}

	// Retrieves the digitisation and range for the channel
	// CalibrationExtraData_t calibration_extra_info;
	// pod5_get_calibration_extra_info(fp->batch, fp->pod5_row, &calibration_extra_info);
	// float scale = calibration_extra_info.range / (float)calibration_extra_info.digitisation;

	//This is same as read_data.num_samples
	// unsigned long int sample_count = 0;
	// pod5_get_read_complete_sample_count(fp->pp, fp->batch, fp->pod5_row, &sample_count);

	int16_t *sig = (int16_t*)malloc(read_data.num_samples * sizeof(int16_t));
	float *sigF = (float*)malloc(read_data.num_samples * sizeof(float));
	
	if (pod5_get_read_complete_signal(pod5_file, fp->batch, fp->pod5_row, read_data.num_samples, sig) != POD5_OK) {
		fprintf(stderr, "Failed to get read %lu signal: %s\n", fp->pod5_row, pod5_get_error_string());
		pod5_release_run_info(run_info_data);
		if(sig) free(sig);
		return;
	}

	char read_id_tmp[37];
	pod5_format_read_id(read_data.read_id, read_id_tmp);

	uint32_t l_sig = 0;
	float pa = 0.0f;
	for(uint64_t i = 0; i < read_data.num_samples; i++){
		pa = (sig[i]+read_data.calibration_offset)*read_data.calibration_scale;
		if (pa > 30 && pa < 200) {
			sigF[l_sig++] = pa;
		}
	}

	free(sig);
	s->sig = (float*)calloc(l_sig, sizeof(float));
	s->l_sig = l_sig;
	memcpy(s->sig, sigF, l_sig*sizeof(float));
	s->name = strdup(read_id_tmp);
	free(sigF);

	// fprintf(stderr, "%s %lu\n", read_id_tmp, s->l_sig);
	// for(int i = 0; i < s->l_sig; i++){
	// 	fprintf(stderr, "%f ", s->sig[i]);
	// }
	// fprintf(stderr, "\n");

	pod5_release_run_info(run_info_data);

	fp->pod5_row++;
	if(fp->pod5_row >= fp->pod5_row_count){

		if (pod5_free_read_batch(fp->batch) != POD5_OK) {
            fprintf(stderr, "Failed to release batch\n");
			// pod5_close_and_free_reader(fp->pp);
			return;
        }

		fp->cur_read++;
		if(fp->cur_read < fp->num_read){
			Pod5ReadRecordBatch_t* batch = NULL;
			if(pod5_get_read_batch(&batch, fp->pp, fp->cur_read) != POD5_OK){
				fprintf(stderr, "Failed to get batch: %s\n", pod5_get_error_string());
				return;
			}

			long unsigned int batch_row_count = 0;
			if(pod5_get_read_batch_row_count(&batch_row_count, batch) != POD5_OK) {
				fprintf(stderr, "Failed to get batch row count\n");
				return;
			}

			fp->pod5_row_count = batch_row_count;
			fp->pod5_row = 0;
			
			fp->batch = batch;
		}
	}
}

void ri_read_sig(ri_sig_file_t* fp, ri_sig_t* s){

	assert(fp->cur_read < fp->num_read);

	if(fp->fp) ri_read_sig_fast5(fp, s);
	else if(fp->pp) ri_read_sig_pod5(fp, s);
}