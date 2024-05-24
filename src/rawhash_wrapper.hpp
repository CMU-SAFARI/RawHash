#pragma once

class RawHashMapper {
    /**
     * @brief Wrapper that can be used from other code (as an alternative to the CLI)
     * 
     */
public:
    RawHashMapper(int argc, char *argv[]);
    ~RawHashMapper();
    void info();

// private:
    // ri_idx_t *ri;
    void *_ri;
};