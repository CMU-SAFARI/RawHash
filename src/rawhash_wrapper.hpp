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
    /**
     * @brief Use the same arguments as the CLI
     * 
     * Including the leading program name (which will be ignored)
     * Note that the option "-d index" (to dump the index) is ignored, so the index should be 
     * dumped with the rawhash CLI, then it can be used for mapping with this function.
     * 
     * @param argc 
     * @param argv 
     */
    RawHashMapper(int argc, char *argv[]);
    ~RawHashMapper();

    /*
    * Map the raw read, returning the alignments
    *
    * Same as ri_map_file_frag -> map_worker_pipeline -> map_worker_for, 
    * except it reads from memory rather than the file
    * and does not perform pipeling or threading
    */
    std::vector<Alignment> map(float* raw_signal, int signal_length);

    /**
     * Get info about the index
    */
    void print_idx_info() const;

private:
    void *_ri; // ri_idx_t*
    void* _opt; // ri_mapopt_t*
};