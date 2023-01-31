#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "rawhash.h"
#include "ketopt.h"
#include "pore_model.h" //TODO, remove the CPP dependency

#define RI_VERSION "0.91"

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

static ko_longopt_t long_options[] = {
	{ (char*)"min-events",			ko_required_argument, 300 },
	{ (char*)"max-gap",			ko_required_argument, 301 },
	{ (char*)"max-target-gap",		ko_required_argument, 302 },
	{ (char*)"max-chains",			ko_required_argument, 303 },
	{ (char*)"min-anchors",		ko_required_argument, 304 },
	{ (char*)"best-chains",		ko_required_argument, 305 },
	{ (char*)"min-score",			ko_required_argument, 306 },
	{ (char*)"max-chunks",			ko_required_argument, 307 },
	{ (char*)"stop-min-anchor",	ko_required_argument, 308 },
	{ (char*)"map-min-anchor",		ko_required_argument, 309 },
	{ (char*)"stop-best-ratio",	ko_required_argument, 310 },
	{ (char*)"map-best-ratio",		ko_required_argument, 311 },
	{ (char*)"stop-mean-ratio",	ko_required_argument, 312 },
	{ (char*)"map-mean-ratio",		ko_required_argument, 313 },
	{ (char*)"bp-per-sec",			ko_required_argument, 314 },
	{ (char*)"sample-rate",		ko_required_argument, 315 },
	{ (char*)"chunk-size",			ko_required_argument, 316 },
	{ (char*)"version",			ko_no_argument, 317 },
	{ (char*)"threshold",			ko_required_argument, 318 },
	{ (char*)"n-samples",			ko_required_argument, 319 },
	{ (char*)"test-frequency",		ko_required_argument, 320 },
	{ (char*)"min-reads",			ko_required_argument, 321 },
	{ (char*)"sequence-until",     ko_no_argument,       322 },
	{ 0, 0, 0 }
};

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	return (int64_t)(x + .499);
}

static inline void yes_or_no(ri_mapopt_t *opt, int64_t flag, int long_idx, const char *arg, int yes_to_set)
{
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
}

int ri_set_opt(const char *preset, ri_idxopt_t *io, ri_mapopt_t *mo)
{
	if (preset == 0) {
		ri_idxopt_init(io);
		ri_mapopt_init(mo);
	} else if (strcmp(preset, "sensitive") == 0) {
		io->e = 6; io->q = 9; io->lq = 3; io->w = 0; io->n = 0;
	} else if (strcmp(preset, "fast") == 0) {
		io->e = 7; io->q = 9; io->lq = 3; io->w = 0; io->n = 0; mo->mini_batch_size = 750000000;
		// mo->min_bestmap_ratio = 1.2f; mo->min_meanmap_ratio = 3; mo->min_meanmap_ratio_out = 3;
	} else if (strcmp(preset, "faster") == 0) {
		io->e = 7; io->q = 9; io->lq = 3; io->w = 5; io->n = 0; mo->mini_batch_size = 1000000000;
		// mo->min_bestmap_ratio = 1.2f; mo->min_meanmap_ratio = 3; mo->min_meanmap_ratio_out = 3;
	} else if (strcmp(preset, "viral") == 0) {
		io->e = 5; io->q = 9; io->lq = 3; io->w = 0; io->n = 0;
	} else if (strcmp(preset, "sequence-until") == 0) {
		io->e = 7; io->q = 9; io->lq = 3; io->w = 0; io->n = 0; mo->mini_batch_size = 750000000;
	} else return -1;
	return 0;
}

int main(int argc, char *argv[])
{
	const char *opt_str = "k:d:p:e:q:l:w:n:o:t:K:x:";
	ketopt_t o = KETOPT_INIT;
	ri_mapopt_t opt;
  	ri_idxopt_t ipt;
	int c, n_threads = 3;
	// int n_parts;
	char *fnw = 0, *fpore = 0;
	FILE *fp_help = stderr;
	ri_idx_reader_t *idx_rdr;
	ri_idx_t *ri;

	ri_verbose = 3;
	liftrlimit();
	ri_realtime0 = ri_realtime();
	ri_set_opt(0, &ipt, &opt);

	// test command line options and apply option -x/preset first
	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'x') {
			if (ri_set_opt(o.arg, &ipt, &opt) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", o.arg);
				return 1;
			}
		} else if (c == ':') {
			fprintf(stderr, "[ERROR] missing option argument\n");
			return 1;
		} else if (c == '?') {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[o.i - 1]);
			return 1;
		}
	}
	o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'd') fnw = o.arg;
		else if (c == 'p') fpore = o.arg;
		else if (c == 'k') ipt.k = atoi(o.arg);
		else if (c == 'e') ipt.e = atoi(o.arg);
		else if (c == 'q') ipt.q = atoi(o.arg);
		else if (c == 'l') ipt.lq = atoi(o.arg);
		else if (c == 'w') ipt.w = atoi(o.arg);
		else if (c == 'n') ipt.n = atoi(o.arg);
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'v') ri_verbose = atoi(o.arg);
		else if (c == 'K') {opt.mini_batch_size = mm_parse_num(o.arg);}
		else if (c == 'h') fp_help = stdout;
		else if (c == 'o') {
			if (strcmp(o.arg, "-") != 0) {
				if (freopen(o.arg, "wb", stdout) == NULL) {
					fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m: %s\n", o.arg, strerror(errno));
					exit(1);
				}
			}
		}
		else if (c == 300) opt.min_events = atoi(o.arg); // --min-events
		else if (c == 301) opt.max_gap_length = atoi(o.arg);// --max-gap
		else if (c == 302) opt.max_target_gap_length = atoi(o.arg);// --max-target-gap
		else if (c == 303) opt.chaining_band_length = atoi(o.arg);// --max-chains
		else if (c == 304) opt.min_num_anchors = atoi(o.arg);// --min-anchors
		else if (c == 305) opt.num_best_chains = atoi(o.arg);// --best-chains
		else if (c == 306) opt.min_chaining_score = atof(o.arg);// --min-score
		else if (c == 307) opt.max_num_chunk = atoi(o.arg);// --max-chunks
		else if (c == 308) opt.min_chain_anchor = atoi(o.arg);// --stop-min-anchor
		else if (c == 309) opt.min_chain_anchor_out = atoi(o.arg);// --map-min-anchor
		else if (c == 310) opt.min_bestmap_ratio = atof(o.arg);// --stop-best-ratio
		else if (c == 311) opt.min_bestmap_ratio_out = atoi(o.arg);// --map-best-ratio
		else if (c == 312) opt.min_meanmap_ratio = atoi(o.arg);// --stop-mean-ratio
		else if (c == 313) opt.min_meanmap_ratio_out = atoi(o.arg);// --map-mean-ratio
		else if (c == 314) opt.bp_per_sec = atoi(o.arg);// --bp-per-sec
		else if (c == 315) opt.sample_rate = atoi(o.arg);// --sample-rate
		else if (c == 316) opt.chunk_size = atoi(o.arg);// --chunk-size
		else if (c == 317) opt.min_events = atoi(o.arg);// --version
		else if (c == 318) opt.t_threshold = atof(o.arg);// --threshold
		else if (c == 319) opt.tn_samples = atoi(o.arg);// --n-samples
		else if (c == 320) opt.ttest_freq = atoi(o.arg);// --test-frequency
		else if (c == 321) opt.tmin_reads = atoi(o.arg);// --min-reads
		else if (c == 322) opt.flag |= RI_M_SEQUENCEUNTIL;// --sequence-until
		else if (c == 'V') {puts(RI_VERSION); return 0;}
	}

	if (argc == o.ind || fp_help == stdout) {
		fprintf(fp_help, "Usage: rawhash [options] <target.fa>|<target.idx> [query.fast5] [...]\n");
		fprintf(fp_help, "Options:\n");
		fprintf(fp_help, "  Indexing:\n");
		fprintf(fp_help, "    -d FILE      [Strongly recommended to create before mapping] dump index to FILE [].\n");
		fprintf(fp_help, "    -p FILE      pore model FILE [].\n");
		fprintf(fp_help, "    -k INT       size of the k-mers in the pore model [%d]. This is usually 6 for R9.4 and 9 for R10\n", ipt.k);
		fprintf(fp_help, "    -e INT       number of events concatanated in a single hash (usually no larger than 10). Also applies during mapping [%d]\n", ipt.e);
		fprintf(fp_help, "    -q INT       most significant bits of signal values to process [%d]. Signal values are assumed to be in the IEEE standard for floating-point arithmetic\n", ipt.q);
		fprintf(fp_help, "    -l INT       least significant bits of the q bits to quantize along with the most signficant 2 bits of the q bits [%d]\n", ipt.lq);
		fprintf(fp_help, "    -w INT       minimizer window size [%d]. Enables minimizer-based seeding in indexing and mapping (may reduce accuracy but improves the performance and memory space efficiency)\n", ipt.w);
		fprintf(fp_help, "    -n NUM       number of consecutive seeds to use for BLEND-based seeding [%d]. Enables the BLEND mechanism (may improve accuracy but reduces the performance at the moment)\n", ipt.n);
		fprintf(fp_help, "\n  Chaining:\n");
		fprintf(fp_help, "    --min-events INT     minimum number of INT events in a chunk to start chain enlongation [%d]\n", opt.min_events);
		fprintf(fp_help, "    --max-gap INT     maximum INT gap length in a chain [%d]\n", opt.max_gap_length);
		fprintf(fp_help, "    --max-target-gap INT     maximum INT gap length in a chain of the target genome [%d]\n", opt.max_target_gap_length);
		fprintf(fp_help, "    --max-chains INT     maximum INT number of chains kept [%d]\n", opt.max_num_skips);
		fprintf(fp_help, "    --min-anchors INT     minimum INT number of anchors to create a chain [%d]\n", opt.min_num_anchors);
		fprintf(fp_help, "    --best-chains INT     best INT chains to keep in the next iteration of the next chunk [%d]\n", opt.num_best_chains);
		fprintf(fp_help, "    --min-score FLOAT     minumum chaining score FLOAT [%g]\n", opt.min_chaining_score);
		fprintf(fp_help, "\n  Mapping:\n");
		fprintf(fp_help, "    --max-chunks INT     maximum INT number of chunks to try for mapping per read [%u]\n", opt.max_num_chunk);
		fprintf(fp_help, "    --stop-min-anchor INT     stop chain enlongation if there is only one chain and the chain has minimum INT number of anchors [%u]\n", opt.min_chain_anchor);
		fprintf(fp_help, "    --map-min-anchor INT     map the read if the best chain has minimum INT number of anchors [%u]\n", opt.min_chain_anchor_out);
		fprintf(fp_help, "    --stop-best-ratio FLOAT     stop chain enlongation if the ratio between the best and the second-best chain scores is => FLOAT [%g]\n", opt.min_bestmap_ratio);
		fprintf(fp_help, "    --map-best-ratio FLOAT     map the read if the ratio between the best and the second-best chain scores is >= FLOAT [%g]\n", opt.min_bestmap_ratio_out);
		fprintf(fp_help, "    --stop-mean-ratio FLOAT     stop chain enlongation if the ratio between the best chain score and the mean chain score is >= FLOAT [%g]\n", opt.min_bestmap_ratio);
		fprintf(fp_help, "    --map-mean-ratio FLOAT     map the read if the ratio between the best chain score and the mean chain score is >= FLOAT [%g]\n", opt.min_meanmap_ratio_out);
		fprintf(fp_help, "\n  ONT Device:\n");
		fprintf(fp_help, "    --bp-per-sec INT     DNA molecules transiting through the pore (bp per second) [%u]\n", opt.bp_per_sec);
		fprintf(fp_help, "    --sample-rate INT     current sample rate in Hz [%u]\n", opt.sample_rate);
		fprintf(fp_help, "    --chunk-size INT     current samples in a single chunk (by default set to the amount of signals sampled in 1 second) [%u]\n", opt.chunk_size);
		fprintf(fp_help, "\n  Sequence Until Parameters:\n");
		fprintf(fp_help, "    --sequence-until - Activates Sequence Until and performs real-time relative abundance calculations. The computation will stop as soon as an estimation with high confidence is reached without processing further reads from the set.\n");
		fprintf(fp_help, "    --threshold FLOAT     outliers are determined if cross-correlation distance > FLOAT [%g]. Sequencing will stop if there are no outliers in the sample of estimations.\n", opt.t_threshold);
		fprintf(fp_help, "    --n-samples INT     New estimation is tested against INT many previous estimations [%u]\n", opt.tn_samples);
		fprintf(fp_help, "    --test-frequency INT     Make a new estimation after every INT reads [%u]\n", opt.ttest_freq);
		fprintf(fp_help, "    --min-reads INT     Minimum number of reads to sequence before making the first estimation [%u]\n", opt.tmin_reads);
		fprintf(fp_help, "\n  Input/Output:\n");
		fprintf(fp_help, "    -o FILE      output mappings to FILE [stdout]\n");
		fprintf(fp_help, "    -t INT       number of threads [%d]\n", n_threads);
		fprintf(fp_help, "    -K NUM       minibatch size for mapping [500M]. Increasing this value may increase thread utilization. If there are many larger FAST5 files, it is recommended to keep this value between 500M - 5G to use less memory while utilizing threads nicely.\n");
//		fprintf(fp_help, "    -v INT       verbose level [%d]\n", ri_verbose);
		fprintf(fp_help, "    --version    show version number\n");
		fprintf(fp_help, "\n  Preset:\n");
		fprintf(fp_help, "    -x STR       preset (always applied before other options) []\n");
		fprintf(fp_help, "                 - sensitive - Enables sensitive mapping. Suitable when low recall is needed or when working with small genomes < 50M.\n");
		fprintf(fp_help, "                 - fast - Enables fast mapping with slightly reduced accuracy. Suitable when reads are mapped to genomes > 50M\n");
		fprintf(fp_help, "                 - faster - Enables faster mapping than -x fast and reduced memory space usage for indexing with slightly reduced accuracy. Suitable when many reads are mapped to very large genomes > 3Gb\n");
		fprintf(fp_help, "                 - viral - Enables accurate mapping to very small genomes such as viral genomes.\n");;
		// fprintf(fp_help, "\nSee `man ./rawhash.1' for detailed description of these and other advanced command-line options.\n");
		return fp_help == stdout? 0 : 1;
	}

	if(ipt.w && ipt.n){
		fprintf(stderr, "[ERROR] minimizer window 'w' ('%d') and BLEND 'neighbor' ('%d') values cannot be set together. At least one of them must be zero to enable one of the seeding options: %s\n", ipt.w, ipt.n, strerror(errno));
		return 1;
	}

	idx_rdr = ri_idx_reader_open(argv[o.ind], &ipt, fnw);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", argv[o.ind], strerror(errno));
		return 1;
	}

	if (!idx_rdr->is_idx && fnw == 0 && argc - o.ind < 2) {
		fprintf(stderr, "[ERROR] missing input: please specify a query FAST5/SLOW5 file(s) to map or option -d to store the index in a file before running the mapping\n");
		ri_idx_reader_close(idx_rdr);
		return 1;
	}

	float* pore_vals = 0;
	if(!idx_rdr->is_idx && fpore == 0){
		fprintf(stderr, "[ERROR] missing input: please specify a pore model file with -p when generating the index from a sequence file\n");
		ri_idx_reader_close(idx_rdr);
		return 1;
	}else if(!idx_rdr->is_idx && fpore){
		sigmap::PoreModel pore_model;
  		pore_model.Load(std::string(fpore));

		if(pore_model.GetKmerSize() <= 4){
			fprintf(stderr, "[ERROR] cannot parse the k-mer pore model file. Please see the example k-mer model files provided in the rawhash repository.\n");
			ri_idx_reader_close(idx_rdr);
			return 1;
		}
		pore_vals = (float*)calloc(1U<<(pore_model.GetKmerSize()*2), sizeof(float));
		for(uint32_t i = 0; i < 1U<<(pore_model.GetKmerSize()*2); ++i)
			pore_vals[i] = pore_model.pore_models_[i].level_mean;
	}

	while ((ri = ri_idx_reader_read(idx_rdr, pore_vals, n_threads)) != 0) {
		int ret;
		if (ri_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, ri_realtime() - ri_realtime0, ri_cputime() / (ri_realtime() - ri_realtime0), ri->n_seq);
		if (ri_verbose >= 3) ri_idx_stat(ri);
		if (argc - (o.ind + 1) == 0) {
			ri_idx_destroy(ri);
			continue; // no query files
		}
		ret = 0;
		// if (!(opt.flag & MM_F_FRAG_MODE)) { //TODO: enable frag mode directly from options
		// for (i = o.ind + 1; i < argc; ++i) {
		// 	ret = ri_map_file(ri, argv[i], &opt, n_threads);
		// 	if (ret < 0) break;
		// }
		// }
		// else { //TODO: enable frag mode directly from options
			ret = ri_map_file_frag(ri, argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &opt, n_threads);
		// }
		ri_idx_destroy(ri);
		if (ret < 0) {
			fprintf(stderr, "ERROR: failed to map the query file\n");
			exit(EXIT_FAILURE);
		}
	}
	// n_parts = idx_rdr->n_parts;
	ri_idx_reader_close(idx_rdr);
	if(pore_vals)free(pore_vals);

	if (fflush(stdout) == EOF) {
		perror("[ERROR] failed to write the results");
		exit(EXIT_FAILURE);
	}

	if (ri_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, RI_VERSION);
		// fprintf(stderr, "[M::%s] CMD:", __func__);
		// for (i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, ri_realtime() - ri_realtime0, ri_cputime(), ri_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
