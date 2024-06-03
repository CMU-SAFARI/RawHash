#pragma once
#include <cstdint>
#include <vector>

struct Alignment {
	const char* contig;
	uint32_t ref_start;
	uint32_t ref_end; // exclusive
	bool is_pos_strand;
};

class RawHashMapper {
    /**
     * @brief Wrapper that can be used from other code (as an alternative to the CLI)
     * 
     * This API is inefficient (no threading) and only for experimenting 
     * with RawHash from Python without having to load the index repeatedly.
     * 
     * Also mixes malloc/free and new/delete
     */
public:
    RawHashMapper(int argc, char *argv[]);
    ~RawHashMapper();

    /*
    * Map the reads
    * Same as ri_map_file_frag -> map_worker_pipeline -> map_worker_for, 
    * except it reads from memory rather than the file
    * and does not perform the pipeline step
    */
    std::vector<Alignment> map(float* raw_signal, int signal_length);

    /**
     * Get info about the index
    */
    void idx_info() const;

private:
    // ri_idx_t *ri;
    void *_ri; // ri_idx_t*
    void* _opt; // ri_mapopt_t*
};