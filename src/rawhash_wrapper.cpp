
#include <string.h>
#include <errno.h>

#include "rawhash_wrapper.hpp"
#include "rawhash.h"
#include "rmap.h"
#include "cli_parsing.cpp"

RawHashMapper::RawHashMapper(int argc, char *argv[]) {
    std::cout << "RawHashMapper constructor" << std::endl;
    
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
    _opt = new ri_mapopt_t;
    *((ri_mapopt_t*)_opt) = opt; // copy constructor
}

RawHashMapper::~RawHashMapper() {
    fprintf(stderr, "[M::%s] closing the index\n", __func__);
    ri_idx_destroy((ri_idx_t*)_ri);
    delete _opt;
}

void RawHashMapper::idx_info() const {
    fprintf(stderr, "[M::%s] index info\n", __func__); // todo
    ri_idx_stat((ri_idx_t*)_ri);
}

std::vector<Alignment> RawHashMapper::map(float* raw_signal, int signal_length) {
    // print
    // fprintf(stderr, "[M::%s] mapping\n", __func__);

    if (ri_verbose >= 3) {
        for (int i = 0; (i < signal_length) && (i <= 10); i++) {
            fprintf(stderr, "%f ", raw_signal[i]);
        }
        fprintf(stderr, "...");
    }

    ri_idx_t* ri = (ri_idx_t*)_ri;
    pipeline_mt pipeline = {
        // .n_processed = 0,
        // .n_threads = 1,
        // .n_fp = 0,
        // .cur_fp = 0,
        // .n_f = 0,
        // .cur_f = 0,
        // .mini_batch_size = 0,
        .opt = (ri_mapopt_t*)_opt,
        // .f = NULL,
        // .fp = NULL,
        .ri = ri,
        // .fn = NULL,
        // .su_nreads = 0,
        // .su_nestimations = 0,
        // .ab_count = 0,
        // .su_cur = 0,
        // .su_estimations = NULL,
        // .su_c_estimations = NULL,
        // .su_stop = 0,
        .map_model = NULL
    };
    ri_tbuf_t *buf = ri_tbuf_init();
    char* read_name = "read123";
    ri_sig_t signal = {
        .rid = 123,
        .l_sig = signal_length,
        .name = read_name,

        .offset = 0, // todo
        .sig = raw_signal,
    };
    ri_sig_t* signal_ptr = &signal; // since &&signal not working since &signal is temporary
    ri_reg1_t reg0_noptr; // will be initialized
    ri_reg1_t* reg0 = &reg0_noptr;
    step_mt data = {
        .p = &pipeline,
        .n_sig = 1,
        .sig = (ri_sig_t**)(&signal_ptr),
        .reg = (ri_reg1_t**)(&reg0),
        .buf = &buf
    };
    map_worker_for(&data, 0, 0);

    // write out
    if (ri_verbose >= 3) {
        write_out_mappings_to_stdout(reg0, ri);
    }
    std::vector<Alignment> alignments;
    for (int m = 0; m < reg0->n_maps; ++m) {
        if (reg0->maps[m].ref_id < ri->n_seq) {
            // mapped
            alignments.push_back(Alignment {
                .contig = (ri->flag & RI_I_SIG_TARGET) ? ri->sig[reg0->maps[m].ref_id].name : ri->seq[reg0->maps[m].ref_id].name,
                .ref_start = reg0->maps[m].fragment_start_position,
                .ref_end = reg0->maps[m].fragment_start_position + reg0->maps[m].fragment_length,
                .is_pos_strand = not reg0->maps[m].rev
            });
        }
    }
    // dummy alignment
    // alignments.push_back(Alignment {
    //     .contig = "fakealignment_chr123", .ref_start = 123, .ref_end = 456, .is_pos_strand = true
    // });
    free_mappings_ri_reg1_t(reg0);

    ri_tbuf_destroy(buf);

    return alignments;
}