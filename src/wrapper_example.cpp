#include "rawhash_wrapper.hpp"

int main(int argc, char *argv[]) {
    RawHashMapper mapper(argc, argv);
    mapper.idx_info();
    return 0;
}