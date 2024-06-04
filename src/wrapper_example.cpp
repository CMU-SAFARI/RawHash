#include "rawhash_wrapper.hpp"

int main(int argc, char *argv[]) {
    RawHashMapper mapper(argc, argv);
    mapper.print_idx_info();
    return 0;
}