
#include <string.h>
#include <errno.h>

#include "rawhash_wrapper.hpp"
#include "rawhash.h"
#include "cli_parsing.cpp"

RawHashMapper::RawHashMapper(int argc, char *argv[]) {
    
    // copied from main()
    liftrlimit();
    ri_realtime0 = ri_realtime();
    
    CLIParsedArgs parsed_args = parse_args(argc, argv);
    ri_mapopt_t& opt = parsed_args.opt;
    ri_idxopt_t& ipt = parsed_args.ipt;
    int& n_threads = parsed_args.n_threads;
    // int n_parts;
    char*& idx_out_filename = parsed_args.idx_out_filename;
    char*& fpore = parsed_args.fpore;
    // FILE*& fp_help = parsed_args.fp_help;
    int& ru_server_port = parsed_args.ru_server_port;
    ketopt_t& o = parsed_args.o;

    ri_idx_reader_t *idx_rdr;
    idx_rdr = ri_idx_reader_open(argv[o.ind], &ipt, idx_out_filename);
    if (idx_rdr == 0) {
        fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", argv[o.ind], strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (!idx_rdr->is_idx && idx_out_filename == 0 && argc - o.ind < 2 && !(ipt.flag&RI_I_OUT_QUANTIZE)) {
        fprintf(stderr, "[ERROR] missing input: please specify a query FAST5/SLOW5/POD5 file(s) to map or option -d to store the index in a file before running the mapping\n");
        ri_idx_reader_close(idx_rdr);
        exit(EXIT_FAILURE);
    }

    // load the pore model
    ri_pore_t pore;
    pore.pore_vals = NULL;
    pore.pore_inds = NULL;
    pore.max_val = -5000.0;
    pore.min_val = 5000.0;
    if(!(ipt.flag&RI_I_OUT_QUANTIZE)) {
        if(!idx_rdr->is_idx && fpore == 0) {
            fprintf(stderr, "[ERROR] missing input: please specify a pore model file with -p when generating the index from a sequence file\n");
            ri_idx_reader_close(idx_rdr);
            exit(EXIT_FAILURE);
        } else if(!idx_rdr->is_idx && fpore) {
            load_pore(fpore, ipt.k, ipt.lev_col, &pore);
            if(!pore.pore_vals){
                fprintf(stderr, "[ERROR] cannot parse the k-mer pore model file. Please see the example k-mer model files provided in the RawHash repository.\n");
                ri_idx_reader_close(idx_rdr);
                exit(EXIT_FAILURE);
            }
        }
    }


    ri_idx_t* ri = ri_idx_reader_read(idx_rdr, &pore, n_threads);
    if (ri == 0) {
        fprintf(stderr, "[ERROR] failed to read the index\n");
        ri_idx_reader_close(idx_rdr);
        exit(EXIT_FAILURE);
    }

    if (ri_verbose >= 3)
        fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
                __func__, ri_realtime() - ri_realtime0, ri_cputime() / (ri_realtime() - ri_realtime0), ri->n_seq);
    
    ri_mapopt_update(&opt, ri);
    if (ri_verbose >= 3) ri_idx_stat(ri);

    _ri = ri;
}

RawHashMapper::~RawHashMapper() {
    fprintf(stderr, "[M::%s] closing the index\n", __func__);
    ri_idx_destroy((ri_idx_t*)_ri);
}

void RawHashMapper::info() {
    fprintf(stderr, "[M::%s] index info\n", __func__); // todo
}
