#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "rawhash.h"
#include "ketopt.h"

#define RH_VERSION "2.0"

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
	{ (char*)"level_column",     	ko_required_argument, 300 },
	{ (char*)"q-mid-occ",     		ko_required_argument, 301 },
	{ (char*)"q-occ-frac",     		ko_required_argument, 302 },
	{ (char*)"min-events",			ko_required_argument, 303 },
	{ (char*)"bw",					ko_required_argument, 304 },
	{ (char*)"max-target-gap",		ko_required_argument, 305 },
	{ (char*)"max-query-gap",		ko_required_argument, 306 },
	{ (char*)"min-anchors",			ko_required_argument, 307 },
	{ (char*)"min-score",			ko_required_argument, 308 },
	{ (char*)"chain-gap-scale",		ko_required_argument, 309 },
	{ (char*)"chain-skip-scale",	ko_required_argument, 310 },
	{ (char*)"chain-match-score",	ko_required_argument, 311 },
	{ (char*)"primary-ratio",		ko_required_argument, 312 },
	{ (char*)"primary-length",		ko_required_argument, 313 },
	{ (char*)"max-skips",			ko_required_argument, 314 },
	{ (char*)"max-iterations",		ko_required_argument, 315 },
	{ (char*)"rmq",					ko_no_argument, 	  316 },
	{ (char*)"rmq-inner-dist",		ko_required_argument, 317 },
	{ (char*)"rmq-size-cap",		ko_required_argument, 318 },
	{ (char*)"bw-long",				ko_required_argument, 319 },
	{ (char*)"max-chunks",			ko_required_argument, 320 },
	{ (char*)"min-mapq",			ko_required_argument, 321 },
	{ (char*)"min-bestmapq",		ko_required_argument, 322 },
	{ (char*)"min-bestmapq-ratio",	ko_required_argument, 323 },
	{ (char*)"min-bestchain-ratio",	ko_required_argument, 324 },
	{ (char*)"min-meanmapq-ratio",	ko_required_argument, 325 },
	{ (char*)"min-meanchain-ratio",	ko_required_argument, 326 },
	{ (char*)"bp-per-sec",			ko_required_argument, 327 },
	{ (char*)"sample-rate",			ko_required_argument, 328 },
	{ (char*)"chunk-size",			ko_required_argument, 329 },
	{ (char*)"seg-window-length1",	ko_required_argument, 330 },
	{ (char*)"seg-window-length2",	ko_required_argument, 331 },
	{ (char*)"seg-threshold1",		ko_required_argument, 332 },
	{ (char*)"seg-threshold2",		ko_required_argument, 333 },
	{ (char*)"seg-peak_height",		ko_required_argument, 334 },
	{ (char*)"sequence-until",     	ko_no_argument,       335 },
	{ (char*)"threshold",			ko_required_argument, 336 },
	{ (char*)"n-samples",			ko_required_argument, 337 },
	{ (char*)"test-frequency",		ko_required_argument, 338 },
	{ (char*)"min-reads",			ko_required_argument, 339 },
	{ (char*)"mid-occ-frac",		ko_required_argument, 340 },
	{ (char*)"version",				ko_no_argument, 	  341 },
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
	} else if (strcmp(preset, "r10sensitive") == 0) {
		io->e = 8; io->q = 8; io->lq = 2; io->w = 0; io->n = 0;
		mo->window_length1 = 5; mo->window_length2 = 12;
		mo->threshold1 = 4.20265f; mo->threshold2 = 3.67058f;  
		mo->peak_height = 0.5f;
	} else if (strcmp(preset, "fast") == 0) {
		io->e = 8; io->q = 9; io->lq = 3; io->w = 0; io->n = 0; mo->mini_batch_size = 750000000;
		// mo->min_bestchain_ratio = 1.2f; mo->min_meanchain_ratio = 3; mo->min_meanchain_ratio_out = 3;
		mo->max_num_chunk = 20;
	} else if (strcmp(preset, "faster") == 0) {
		io->e = 7; io->q = 9; io->lq = 3; io->w = 5; io->n = 0; mo->mini_batch_size = 1000000000;
		// mo->min_bestchain_ratio = 1.2f; mo->min_meanchain_ratio = 3; mo->min_meanchain_ratio_out = 3;
	} else if (strcmp(preset, "viral") == 0) {
		io->e = 5; io->q = 9; io->lq = 3; io->w = 0; io->n = 0;
		mo->min_num_anchors = 2; mo->min_chaining_score = 10;
		// mo->mini_batch_size = 1000000000;
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
	char *fnw = 0, *fpore = 0, *s;
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
		else if (c == 300) ipt.lev_col = atoi(o.arg);// --level_column
		else if (c == 301) { //--q-mid-occ
			opt.min_mid_occ = strtol(o.arg, &s, 10); // min
			if (*s == ',') opt.max_mid_occ = strtol(s + 1, &s, 10); //max
			// opt.q_mid_occ = atoi(o.arg);// --q-mid-occ
		}
		else if (c == 302) opt.q_occ_frac = atof(o.arg);// --q-occ-frac
		else if (c == 303) opt.min_events = atoi(o.arg); // --min-events
		else if (c == 304) opt.bw = atoi(o.arg);// --bw
		else if (c == 305) opt.max_target_gap_length = atoi(o.arg);// --max-target-gap
		else if (c == 306) opt.max_query_gap_length = atoi(o.arg);// --max-query-gap
		else if (c == 307) opt.min_num_anchors = atoi(o.arg);// --min-anchors
		else if (c == 308) opt.min_chaining_score = atoi(o.arg);// --min-score
		else if (c == 308) opt.chain_gap_scale = atof(o.arg);// --chain-gap-scale
		else if (c == 310) opt.chain_skip_scale = atof(o.arg);// --chain-skip-scale
		else if (c == 311) opt.a = atoi(o.arg);// --chain-match-score
		else if (c == 312) opt.mask_level = atof(o.arg);// --primary-ratio
		else if (c == 313) opt.mask_len = atoi(o.arg);// --primary-length
		else if (c == 314) opt.max_num_skips = atoi(o.arg);// --max-skips
		else if (c == 315) opt.max_chain_iter = atoi(o.arg);// --max-iterations
		else if (c == 316) opt.flag |= RI_M_RMQ; // --rmq
		else if (c == 317) opt.rmq_inner_dist = atoi(o.arg); // --rmq-inner-dist
		else if (c == 318) opt.rmq_size_cap = atoi(o.arg); // --rmq-size-cap
		else if (c == 319) opt.bw_long = atoi(o.arg);// --bw-long
		else if (c == 320) opt.max_num_chunk = atoi(o.arg);// --max-chunks
		else if (c == 321) opt.min_mapq = atoi(o.arg);// --min-mapq
		else if (c == 322) opt.min_bestmapq = atoi(o.arg);// --min-bestmapq
		else if (c == 323) opt.min_bestmapq_ratio = atof(o.arg);// --min-bestmapq-ratio
		else if (c == 324) opt.min_bestchain_ratio = atof(o.arg);// --min-bestchain-ratio
		else if (c == 325) opt.min_meanmapq_ratio = atof(o.arg);// --min-meanmapq-ratio
		else if (c == 326) opt.min_meanchain_ratio = atof(o.arg);// --min-meanchain-ratio
		else if (c == 327) {opt.bp_per_sec = atoi(o.arg); opt.sample_per_base = (float)opt.sample_rate / opt.bp_per_sec;}// --bp-per-sec
		else if (c == 328) {opt.sample_rate = atoi(o.arg); opt.sample_per_base = (float)opt.sample_rate / opt.bp_per_sec;}// --sample-rate
		else if (c == 329) opt.chunk_size = atoi(o.arg);// --chunk-size
		else if (c == 330) opt.window_length1 = atoi(o.arg);// --seg-window-length1
		else if (c == 331) opt.window_length2 = atoi(o.arg);// --seg-window-length2
		else if (c == 332) opt.threshold1 = atof(o.arg);// --seg-threshold1
		else if (c == 333) opt.threshold2 = atof(o.arg);// --seg-threshold2
		else if (c == 334) opt.peak_height = atof(o.arg);// --seg-peak_height
		else if (c == 335) opt.flag |= RI_M_SEQUENCEUNTIL;// --sequence-until
		else if (c == 336) opt.t_threshold = atof(o.arg);// --threshold
		else if (c == 337) opt.tn_samples = atoi(o.arg);// --n-samples
		else if (c == 338) opt.ttest_freq = atoi(o.arg);// --test-frequency
		else if (c == 339) opt.tmin_reads = atoi(o.arg);// --min-reads
		else if (c == 340) opt.mid_occ_frac = atof(o.arg);// --mid-occ-frac
		else if (c == 341) {puts(RH_VERSION); return 0;}// --version
		else if (c == 'V') {puts(RH_VERSION); return 0;}
	}

	if (argc == o.ind || fp_help == stdout) {
		fprintf(fp_help, "Usage: rawhash [options] <target.fa>|<target.idx> [query.fast5] [...]\n");
		fprintf(fp_help, "Options:\n");
		
		fprintf(fp_help, "  K-mer (pore) Model:\n");
		fprintf(fp_help, "    -p FILE      pore model FILE [].\n");
		fprintf(fp_help, "    -k INT       size of the k-mers in the pore model [%d]. This is usually 6 for R9.4 and 9 for R10\n", ipt.k);
		fprintf(fp_help, "    --level_column INT       0-based column index where the mean values are stored in the pore file [%d]. This is usually 1 for both R9.4 and R10\n", ipt.lev_col);
		
		fprintf(fp_help, "\n  Indexing:\n");
		fprintf(fp_help, "    -d FILE     [Strongly recommended to create before mapping] dump index to FILE [].\n");
		
		fprintf(fp_help, "    -e INT     number of events concatanated in a single hash (usually no larger than 10). Also applies during mapping [%d]\n", ipt.e);
		fprintf(fp_help, "    -q INT     most significant bits of signal values to process [%d]. Signal values are assumed to be in the IEEE standard for floating-point arithmetic\n", ipt.q);
		fprintf(fp_help, "    -l INT     least significant bits of the q bits to quantize along with the most signficant 2 bits of the q bits [%d]\n", ipt.lq);
		fprintf(fp_help, "    -w INT     minimizer window size [%d]. Enables minimizer-based seeding in indexing and mapping (may reduce accuracy but improves the performance and memory space efficiency)\n", ipt.w);
		fprintf(fp_help, "    -n NUM     number of consecutive seeds to use for BLEND-based seeding [%d]. Enables the BLEND mechanism (may improve accuracy but reduces the performance at the moment)\n", ipt.n);
		
		fprintf(fp_help, "\n  Seeding:\n");
		fprintf(fp_help, "    --q-mid-occ INT1[,INT2]     Lower and upper bounds of k-mer occurrences [%d, %d]. The final k-mer occurrence threshold is max{INT1, min{INT2, --q-occ-frac}}. This option prevents excessively small or large -f estimated from the input reference.\n", opt.min_mid_occ, opt.max_mid_occ);
		fprintf(fp_help, "    --q-occ-frac FLOAT     Discard a query seed if its occurrence is higher than FLOAT fraction of all query seeds [%g]. Set 0 to disable. [Note: Both --q-mid-occ and --q-occ-frac should be met for a seed to be discarded]\n", opt.q_occ_frac);
		
		fprintf(fp_help, "\n  Chaining Parameters:\n");
		fprintf(fp_help, "    --min-events INT     minimum number of INT events in a chunk to start chain elongation [%d]\n", opt.min_events);
		fprintf(fp_help, "    --bw INT     maximum INT gap length in a chain [%d]\n", opt.bw);
		fprintf(fp_help, "    --max-target-gap INT     maximum INT target gap length in a chain [%d]\n", opt.max_target_gap_length);
		fprintf(fp_help, "    --max-query-gap INT     maximum INT query gap length in a chain [%d]\n", opt.max_query_gap_length);
		fprintf(fp_help, "    --min-anchors INT     chain is discarded if it contains less than INT number of anchors [%d]\n", opt.min_num_anchors);
		// fprintf(fp_help, "    --best-chains INT     best INT chains to keep in the next iteration of the next chunk [%d]\n", opt.num_best_chains);
		fprintf(fp_help, "    --min-score INT     chain is discarded if its score is < INT [%d]\n", opt.min_chaining_score);
		fprintf(fp_help, "    --chain-gap-scale FLOAT     [Advanced] Determines [chain gap penalty] = FLOAT * 0.01 * e  [%g]\n", opt.chain_gap_scale);
		fprintf(fp_help, "    --chain-skip-scale FLOAT     [Advanced] Determines [chain skip penalty] = FLOAT * 0.01 * e  [%g]\n", opt.chain_skip_scale);
		fprintf(fp_help, "    --chain-match-score INT     [Advanced] Match score (Used in MAPQ and Primary chain identification) [%d]\n", opt.a);
		fprintf(fp_help, "    --primary-ratio FLOAT     [Advanced] The chain is primary if its region ratio uncovered by other chains is larger than FLOAT [%g]\n", opt.mask_level);
		fprintf(fp_help, "    --primary-length INT     [Advanced] The chain is primary if its region length uncovered by other chains is larger than INT [%d]\n", opt.mask_len);
		fprintf(fp_help, "    --max-skips INT     [Advanced] stop looking for a predecessor for an anchor if the best predecessor is not updated after INT many iterations [%d]\n", opt.max_num_skips);
		fprintf(fp_help, "    --max-iterations INT     [Advanced] maximum INT number predecessor anchors to check to calculate the best score for an anchor [%d]\n", opt.max_chain_iter);
		fprintf(fp_help, "    --rmq     [Advanced] Uses RMQ-based chaining. Faster but less accurate than default (DP)\n");
		fprintf(fp_help, "    --rmq-inner-dist INT     [Advanced] RMQ inner distance [%d]\n", opt.rmq_inner_dist);
		fprintf(fp_help, "    --rmq-size-cap INT     [Advanced] RMQ cap size [%d]\n", opt.rmq_size_cap);
		fprintf(fp_help, "    --bw-long INT     [Advanced] maximum long INT gap length in a chain [Used in RMQ] [%d]\n", opt.bw_long);
		
		fprintf(fp_help, "\n  Mapping Decisions (Mapping and sequencing is stopped after taking any of these decisions):\n");
		fprintf(fp_help, "    --max-chunks INT     stop mapping (read not mapped) after sequencing INT number of chunks [%u]\n", opt.max_num_chunk);
		fprintf(fp_help, "    --min-mapq INT     map the read if there is only one chain and its MAPQ > INT [%d]\n", opt.min_mapq);
		fprintf(fp_help, "    --min-bestmapq INT     map the read if the best chain MAPQ > INT and the second-best chain MAPQ == 0 [%d]\n", opt.min_bestmapq);
		fprintf(fp_help, "    --min-bestmapq-ratio FLOAT     map the read if the MAPQ ratio between the best and the second-best chain is larger than or equal to FLOAT [%g]\n", opt.min_bestmapq_ratio);
		fprintf(fp_help, "    --min-bestchain-ratio FLOAT     map the read if the chain score ratio between the best and the second-best chain is larger than or equal to FLOAT [%g]\n", opt.min_bestchain_ratio);
		fprintf(fp_help, "    --min-meanmapq-ratio FLOAT     map the read if the MAPQ ratio between the best and all chains is larger than or equal to FLOAT [%g]\n", opt.min_meanmapq_ratio);
		fprintf(fp_help, "    --min-meanchain-ratio FLOAT     map the read if the chain score ratio between the best and all chains is larger than or equal to FLOAT [%g]\n", opt.min_meanchain_ratio);
		
		fprintf(fp_help, "\n  Nanopore Parameters:\n");
		fprintf(fp_help, "    --bp-per-sec INT     DNA molecules transiting through the pore (bp per second) [%u]\n", opt.bp_per_sec);
		fprintf(fp_help, "    --sample-rate INT     current sample rate in Hz [%u]\n", opt.sample_rate);
		fprintf(fp_help, "    --chunk-size INT     current samples in a single chunk (by default set to the amount of signals sampled in 1 second) [%u]\n", opt.chunk_size);
		
		fprintf(fp_help, "    --seg-window-length1 INT     [Advanced] First window length in segmentation [%u]\n", opt.window_length1);
		fprintf(fp_help, "    --seg-window-length2 INT     [Advanced] Second window length in segmentation [%u]\n", opt.window_length2);
		fprintf(fp_help, "    --seg-threshold1 FLOAT     [Advanced] Peak value threshold for the first window in segmentation [%g]\n", opt.threshold1);
		fprintf(fp_help, "    --seg-threshold2 FLOAT     [Advanced] Peak value threshold for the first window in segmentation [%g]\n", opt.threshold2);
		fprintf(fp_help, "    --seg-peak_height FLOAT     [Advanced] Peak height than the current signal to confirm the peak point in segmentation [%g]\n", opt.peak_height);

		fprintf(fp_help, "\n  Sequence Until Parameters:\n");
		fprintf(fp_help, "    --sequence-until     Activates Sequence Until and performs real-time relative abundance calculations. The computation will stop as soon as an estimation with high confidence is reached without processing further reads from the set.\n");
		fprintf(fp_help, "    --threshold FLOAT     outliers are determined if cross-correlation distance > FLOAT [%g]. Sequencing will stop if there are no outliers in the sample of estimations.\n", opt.t_threshold);
		fprintf(fp_help, "    --n-samples INT     New estimation is tested against INT many previous estimations [%u]\n", opt.tn_samples);
		fprintf(fp_help, "    --test-frequency INT     Make a new estimation after every INT reads [%u]\n", opt.ttest_freq);
		fprintf(fp_help, "    --min-reads INT     Minimum number of reads to sequence before making the first estimation [%u]\n", opt.tmin_reads);
		
		fprintf(fp_help, "\n  Input/Output:\n");
		fprintf(fp_help, "    -o FILE     output mappings to FILE [stdout]\n");
		fprintf(fp_help, "    -t INT      number of threads [%d]\n", n_threads);
		fprintf(fp_help, "    -K NUM      minibatch size for mapping [500M]. Increasing this value may increase thread utilization. If there are many larger FAST5 files, it is recommended to keep this value between 500M - 5G to use less memory while utilizing threads nicely.\n");
//		fprintf(fp_help, "    -v INT     verbose level [%d]\n", ri_verbose);
		fprintf(fp_help, "    --version     show version number\n");
		
		fprintf(fp_help, "\n  Preset:\n");
		fprintf(fp_help, "    -x STR     preset (always applied before other options) []\n");
		fprintf(fp_help, "                 - sensitive     Enables sensitive mapping (for R9.4 model). Suitable when low recall is needed or when working with small genomes < 500M.\n");
		fprintf(fp_help, "                 - r10sensitive     Enables sensitive mapping for R10. Currently under development and may produce slighlty less accurate results than R9.4.\n");
		fprintf(fp_help, "                 - fast     Enables fast mapping with slightly reduced accuracy. Suitable when reads are mapped to genomes > 500M and < 5Gb\n");
		fprintf(fp_help, "                 - faster     Enables faster mapping than '-x fast' and reduced memory space usage for indexing with slightly reduced accuracy. This mechanism uses the minimizer sketching technique and should be used when '-x fast' cannot meet the real-time requirements for a particular genome (e.g., for very large genomes > 5Gb)\n");
		fprintf(fp_help, "                 - viral     Enables accurate mapping to very small genomes such as viral genomes.\n");;
		
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
		load_pore(fpore, ipt.k, ipt.lev_col, &pore_vals);
		if(!pore_vals){
			fprintf(stderr, "[ERROR] cannot parse the k-mer pore model file. Please see the example k-mer model files provided in the rawhash repository.\n");
			ri_idx_reader_close(idx_rdr);
			return 1;
		}
	}

	while ((ri = ri_idx_reader_read(idx_rdr, pore_vals, n_threads)) != 0) {
		int ret;
		if (ri_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, ri_realtime() - ri_realtime0, ri_cputime() / (ri_realtime() - ri_realtime0), ri->n_seq);
		if (argc != o.ind + 1) ri_mapopt_update(&opt, ri);
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
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, RH_VERSION);
		// fprintf(stderr, "[M::%s] CMD:", __func__);
		// for (i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, ri_realtime() - ri_realtime0, ri_cputime(), ri_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
