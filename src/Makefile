CPPFLAGS=-std=c++11 -Wall -O3 -fopenmp -march=native -mavx2 -msse4.1 -I${HDF5_INCLUDE_DIR} -I${SLOW5_INCLUDE_DIR}
# For debugging
# CPPFLAGS=-std=c++11 -Wall -O0 -g -pg -fopenmp -march=native -mavx2 -msse4.1 -I${HDF5_INCLUDE_DIR} -I${SLOW5_INCLUDE_DIR}

OBJS= kthread.o kalloc.o bseq.o pore_model.o rawindex.o main.o

WORKDIR = $(shell pwd)
HDF5_DIR ?= ${WORKDIR}/../extern/hdf5/build
HDF5_INCLUDE_DIR ?= ${HDF5_DIR}/include
HDF5_LIB_DIR ?= ${HDF5_DIR}/lib
HDF5_LIB ?= hdf5
SLOW5_DIR ?= ${WORKDIR}/../extern/slow5lib/
SLOW5_INCLUDE_DIR ?= ${SLOW5_DIR}/include
SLOW5_LIB_DIR ?= ${SLOW5_DIR}/lib

LDFLAGS=${HDF5_LIB_DIR}/lib${HDF5_LIB}.a ${SLOW5_LIB_DIR}/libslow5.a -lm -lz -ldl -g

PROG=rawhash

.PHONY: all clean
# .SUFFIXES:.cpp .o

all: hdf5 slow5 check_hdf5 check_slow5 $(PROG)
subset: check_slow5 check_hdf5 $(PROG)

check_hdf5:
	@[ -f "${HDF5_INCLUDE_DIR}/H5pubconf.h" ] || { echo "HDF5 headers not found" >&2; exit 1; }
	@[ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.so" ] || [ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.a" ] || { echo "HDF5 library not found" >&2; exit 1; }

check_slow5:
	@[ -f "${SLOW5_INCLUDE_DIR}/slow5/slow5.h" ] || { echo "SLOW5 headers not found" >&2; exit 1; }
	@[ -f "${SLOW5_LIB_DIR}/libslow5.so" ] || [ -f "${SLOW5_LIB_DIR}/libslow5.a" ] || { echo "SLOW5 library not found" >&2; exit 1; }

hdf5:
	cd ../extern/hdf5;\
  	mkdir build;\
	./configure --enable-threadsafe --disable-hl --prefix="${HDF5_DIR}";\
	make -j;\
	make install

slow5:
	make -C ${SLOW5_DIR}

$(PROG): $(OBJS)
	${CXX} $(CPPFLAGS) $(OBJS) -o $(PROG) $(LDFLAGS)

%.o: %.cpp
	${CXX} $(CPPFLAGS) -c $< -o $@

clean:
	rm -fr *.o $(PROG) *~