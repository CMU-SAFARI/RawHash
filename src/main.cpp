#include <string.h>
#include <errno.h>

#include "cli_parsing.cpp"
#include "rawhash.h"

#ifdef RUCLIENT_ENABLED
	#include "rawhash_ruclient.hpp"

	#include "ru_client/ru_method.hpp"
	#include "ru_client/grpc_utils.hpp"
	#include "ru_client/data_client.hpp"
	#include "ru_client/ont_device_client.hpp"
	#include "ru_client/utils.hpp" // for is_big_endian()

	#include "ru_client/spdlog_dropin.hpp"
	namespace spdlog {
		namespace level {

			using namespace ru_client::internal::level;
		}
	}
	namespace {
		void log(const std::string& message, ru_client::internal::level::spdlog_level level = spdlog::level::info) {
			// std::cout << "Sev" << level << ": " << message << std::endl;
			fprintf(stderr, "%s\n", ("Sev" + std::to_string(level) + ": " + message).c_str());
		}
	}
	using namespace ru_client;
#endif

int main(int argc, char *argv[])
{
	liftrlimit();
	ri_realtime0 = ri_realtime();
	
	CLIParsedArgs parsed_args = parse_args(argc, argv);
	ri_mapopt_t& opt = parsed_args.opt;
  	ri_idxopt_t& ipt = parsed_args.ipt;
	int& n_threads = parsed_args.n_threads;
	// int n_parts;
	char*& idx_out_filename = parsed_args.idx_out_filename;
	char*& fpore = parsed_args.fpore;
	FILE*& fp_help = parsed_args.fp_help;
	int& ru_server_port = parsed_args.ru_server_port;
	ketopt_t& o = parsed_args.o;

	ri_idx_reader_t *idx_rdr;
	ri_idx_t *ri;
	idx_rdr = ri_idx_reader_open(argv[o.ind], &ipt, idx_out_filename);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", argv[o.ind], strerror(errno));
		exit(EXIT_FAILURE);
	}

	if (!idx_rdr->is_idx && idx_out_filename == 0 && argc - o.ind < 2) {
		fprintf(stderr, "[ERROR] missing input: please specify a query FAST5/SLOW5 file(s) to map or option -d to store the index in a file before running the mapping\n");
		ri_idx_reader_close(idx_rdr);
		exit(EXIT_FAILURE);
	}

	ri_pore_t pore;
	pore.pore_vals = NULL;
	pore.pore_inds = NULL;
	pore.max_val = -5000.0;
	pore.min_val = 5000.0;
	if(!idx_rdr->is_idx && fpore == 0){
		fprintf(stderr, "[ERROR] missing input: please specify a pore model file with -p when generating the index from a sequence file\n");
		ri_idx_reader_close(idx_rdr);
		exit(EXIT_FAILURE);
	}else if(!idx_rdr->is_idx && fpore){
		load_pore(fpore, ipt.k, ipt.lev_col, &pore);
		if(!pore.pore_vals){
			fprintf(stderr, "[ERROR] cannot parse the k-mer pore model file. Please see the example k-mer model files provided in the RawHash repository.\n");
			ri_idx_reader_close(idx_rdr);
			exit(EXIT_FAILURE);
		}
	}

	#ifdef RUCLIENT_ENABLED
	// todo1: remove
	int& port = ru_server_port;
	if (port == 0) {
        std::ifstream port_file("/home/mmordig/rawhash_project/ru_python/example_run/server_run/ont_device_server_port.txt");
        if (!port_file.is_open()) {
            log("Could not open port file", spdlog::level::err);
            return EXIT_FAILURE;
        }
        port_file >> port;
    }
    if (port <= 0) {
        log("Invalid port " + std::to_string(port), spdlog::level::err);
        return EXIT_FAILURE;
    }
	log("Using port " + std::to_string(port));
	#endif

	if (ru_server_port == -1) {
		// offline processing
		while ((ri = ri_idx_reader_read(idx_rdr, &pore, n_threads)) != 0) {
			int ret;
			if (ri_verbose >= 3)
				fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
						__func__, ri_realtime() - ri_realtime0, ri_cputime() / (ri_realtime() - ri_realtime0), ri->n_seq);
			if (argc != o.ind + 1) ri_mapopt_update(&opt, ri); // only update the index if it is used for querying later, todo4: can be moved down to before "ret = 0"
			if (ri_verbose >= 3) ri_idx_stat(ri);
			if (argc - (o.ind + 1) == 0) {
				fprintf(stderr, "[INFO] no files to query index on, just created the index\n");
				ri_idx_destroy(ri);
				continue; // no query files, just creating the index
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
	} else {
		// online processing and decision making by connecting to a real device
		#ifndef RUCLIENT_ENABLED
			fprintf(stderr, "[ERROR] The --ru-server-port option is not enabled in this build\n");
			exit(EXIT_FAILURE);
		#else

		fprintf(stderr, "Connecting to the RU device on port %d\n", ru_server_port);
		
		if (argc - (o.ind + 1) > 0) {
			fprintf(stderr, "[ERROR] The --ru-server-port option cannot be used with query files\n");
			exit(EXIT_FAILURE);
		}

		// todo4: using while loop to detect only one index is read because it is unclear what happens if there are several because ri_idx_reader_read does not increment file pointer when loading idx
		bool read_one = false;
		while ((ri = ri_idx_reader_read(idx_rdr, &pore, n_threads)) != 0) {
			if (read_one) {
				fprintf(stderr, "[ERROR] Only one index file can be read when using the --ru-server-port option\n");
				exit(EXIT_FAILURE);
			}
			read_one = true;

			if (ri_verbose >= 3)
				fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
						__func__, ri_realtime() - ri_realtime0, ri_cputime() / (ri_realtime() - ri_realtime0), ri->n_seq);
			ri_mapopt_update(&opt, ri);
			if (ri_verbose >= 3) ri_idx_stat(ri);

			auto grpc_channel = setup_grpc_channel({.port = ru_server_port});
    		AcquisitionServiceClient acquisition_client(grpc_channel);
			
			DeviceServiceClient device_client(grpc_channel);
        	uint32_t n_channels = device_client.get_nb_channels();

			fprintf(stderr, "Connected to the device with %d channels\n", n_channels);

			DataServiceClient data_client(grpc_channel);

			auto data_type = data_client.get_data_types()["calibrated"];
			assert((data_type.number_type == "int"));
			assert((data_type.n_bytes == sizeof(int16_t)));
			assert((data_type.big_endian == is_big_endian()));

			auto calibrations = device_client.get_calibration(1, n_channels);
        	log("Received calibrations: " + calibrations.to_string(), spdlog::level::info);

			auto start_acquisition = [&acquisition_client] {
				if (!acquisition_client.start_sequencing()) {
					fprintf(stderr, "Failed to start sequencing, may already be running");
				}
			};

			int first_channel = 1;
    		int last_channel = n_channels;
			ReportNumSamplesBehindCallback report_num_samples_behind_cb(acquisition_client, 10); // todo1: n_channels instead of 10

			auto create_decision_maker = [&report_num_samples_behind_cb, &ri, &opt, &calibrations, &first_channel](uint32_t channel) -> std::shared_ptr<DecisionMaker> {
            	std::unique_ptr<DecisionMaker> dec_maker;
				dec_maker = std::make_unique<RawHashDecisionMaker>(ri, &opt, SingleChannelCalibration {
					.range = calibrations.ranges[channel-first_channel],
					.digitisation = calibrations.digitisation,
					.offset = calibrations.offsets[channel-first_channel]
				});
				return std::make_unique<CombinedDecisionMaker>(std::vector<std::shared_ptr<DecisionMaker>> {
					std::make_shared<ReportChunkSampleStartDecisionMaker>(report_num_samples_behind_cb),
					std::move(dec_maker)
				});
			};

			DataServiceClient::ReadUntilSetup setup = {
				.first_channel = (uint32_t)first_channel, .last_channel = (uint32_t)last_channel, .calibrated = true, 
				.sample_minimum_chunk_size = 100, .max_unblock_read_length_samples = 1000, 
				.accepted_first_chunk_classifications = {} // all
			};
			if (!data_client.start_read_until(setup)) {
				log("Failed to start read-until", spdlog::level::err);
				exit(EXIT_FAILURE);
			}

			start_acquisition();
			uint32_t num_monitored_channels = (uint32_t)(last_channel - first_channel + 1);
			auto decision_maker_stats = start_decision_maker_threads(data_client, create_decision_maker, first_channel, last_channel);

			// todo1: remove
			int run_time = 2000; // todo
        	wait_until_elapsed_or_interrupted(run_time);
			log("Stopping readuntil");
			if (!data_client.stop_read_until()) {
				log("Failed to stop read-until", spdlog::level::err);
				exit(EXIT_FAILURE);
			} 
			// signals threads to stop

			log("Waiting for decision maker threads to finish");
			DecisionStats decision_stats;
			decision_stats.num_channels = 0;
			for (size_t i(0); i < decision_maker_stats.size(); ++i) {
				decision_stats += decision_maker_stats[i].get();
				log("Joined ru thread for channel " + std::to_string(i+1) + " / " + std::to_string(decision_maker_stats.size()));
			}
			log("Decision maker threads finished");
			log("Decision stats: " + decision_stats.to_string());

			if ((uint32_t)decision_stats.num_channels != num_monitored_channels) {
				log("Not all channels received data: " + std::to_string(decision_stats.num_channels) + " / " + std::to_string(num_monitored_channels), spdlog::level::warn);
			}

			if (!acquisition_client.stop_sequencing()) {
				log("Failed to stop sequencing, may already be stopped", spdlog::level::warn);
			}

			ri_idx_destroy(ri);
		}
		#endif
	}

	// n_parts = idx_rdr->n_parts;
	ri_idx_reader_close(idx_rdr);
	if(pore.pore_vals)free(pore.pore_vals);
	if(pore.pore_inds)free(pore.pore_inds);

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
	exit(EXIT_SUCCESS);
}
