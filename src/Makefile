# CPPFLAGS=-std=c++11 -Wall -O3 -fopenmp -march=native -mavx2 -msse4.1 -I${HDF5_INCLUDE_DIR} -I${SLOW5_INCLUDE_DIR}
#wo SIMD
CFLAGS=-Wall -O3 -Wc++-compat -DHAVE_KALLOC
CPPFLAGS=-std=c++11 -Wall -O3 -fopenmp -march=native -DHAVE_KALLOC
# For debugging
# CPPFLAGS=--std=c++11 -Wall -O2 -g -pg -fopenmp -march=native -DHAVE_KALLOC
# CFLAGS=-Wall -O2 -g -pg -Wc++-compat -DHAVE_KALLOC

OBJS= kthread.o kalloc.o bseq.o roptions.o sequence_until.o rutils.o pore_model.o rsig.o revent.o rsketch.o rawindex.o rmap.o main.o

CMAKE_CXX_COMPILER ?= g++
CMAKE_CXX_COMPILER_VERSION ?= $(shell $(CMAKE_CXX_COMPILER) -dumpversion)
CMAKE_SYSTEM_PROCESSOR ?= $(shell uname -m)
CMAKE_SYSTEM_NAME ?= $(shell uname -s)

# WORKDIR = $(CURDIR)/../extern
WORKDIR = $(shell pwd)/../extern
# WORKDIR +=

POD5_VERSION = 0.1.8
POD5_URLDIR = pod5-$(POD5_VERSION)-$(CMAKE_SYSTEM_NAME)
POD5_REPO = https://github.com/nanoporetech/pod5-file-format
POD5_INCLUDE_DIR = $(WORKDIR)/$(POD5_URLDIR)/include

ifeq ($(CMAKE_SYSTEM_NAME),Linux)
    LIB_DIR = lib64
    ifeq ($(shell echo $(CMAKE_SYSTEM_PROCESSOR) | grep -E '^(aarch64|arm)'),)
        POD5_URL = $(POD5_REPO)/releases/download/$(POD5_VERSION)/lib_pod5-$(POD5_VERSION)-linux-x64.tar.gz
    else ifeq ($(CMAKE_CXX_COMPILER_ID),GNU)
        ifeq ($(shell expr $(CMAKE_CXX_COMPILER_VERSION) \< 9.0),1)
            POD5_VERSION = 0.1.8
            POD5_URLDIR = pod5-$(POD5_VERSION)-$(CMAKE_SYSTEM_NAME)
            POD5_INCLUDE_DIR = $(WORKDIR)/$(POD5_URLDIR)/include
            POD5_URL = $(POD5_REPO)/releases/download/$(POD5_VERSION)/lib_pod5-$(POD5_VERSION)-linux-gcc7-arm64.tar.gz
            LIB_DIR = lib
        endif
    else ifeq ($(shell echo $(CMAKE_CXX_COMPILER) | grep -E 'clang\+\+$'),)
        ifeq ($(shell expr $(CMAKE_CXX_COMPILER_VERSION) \< 11.0),1)
            POD5_VERSION = 0.1.8
            POD5_URLDIR = pod5-$(POD5_VERSION)-$(CMAKE_SYSTEM_NAME)
            POD5_INCLUDE_DIR = $(WORKDIR)/$(POD5_URLDIR)/include
            POD5_URL = $(POD5_REPO)/releases/download/$(POD5_VERSION)/lib_pod5-$(POD5_VERSION)-linux-clang11-arm64.tar.gz
        endif
    endif
    POD5_LIBRARIES = $(WORKDIR)/$(POD5_URLDIR)/$(LIB_DIR)/libpod5_format.a \
                     $(WORKDIR)/$(POD5_URLDIR)/$(LIB_DIR)/libarrow.a \
                     $(WORKDIR)/$(POD5_URLDIR)/$(LIB_DIR)/libjemalloc_pic.a
endif

ifeq ($(CMAKE_SYSTEM_NAME),Darwin)
    POD5_URL = $(POD5_REPO)/releases/download/$(POD5_VERSION)/lib_pod5-$(POD5_VERSION)-osx-11.0-arm64.tar.gz
    POD5_LIBRARIES = $(WORKDIR)/$(POD5_URLDIR)/lib/libpod5_format.a \
                     $(WORKDIR)/$(POD5_URLDIR)/lib/libarrow.a
endif

ifeq ($(CMAKE_SYSTEM_NAME),Windows_NT)
    POD5_URL = $(POD5_REPO)/releases/download/$(POD5_VERSION)/lib_pod5-$(POD5_VERSION)-win
endif

# WORKDIR = $(shell pwd)
HDF5_DIR ?= ${WORKDIR}/hdf5/build
HDF5_INCLUDE_DIR ?= ${HDF5_DIR}/include
HDF5_LIB_DIR ?= ${HDF5_DIR}/lib
HDF5_LIB ?= hdf5

SLOW5_DIR ?= ${WORKDIR}/slow5lib/
SLOW5_INCLUDE_DIR ?= ${SLOW5_DIR}/include
SLOW5_LIB_DIR ?= ${SLOW5_DIR}/lib

# POD5_DIR ?= ${WORKDIR}/pod5/
# POD5_INCLUDE_DIR ?= ${POD5_DIR}/include
# POD5_LIB_DIR ?= ${POD5_DIR}/lib64

INCLUDES=-I${HDF5_INCLUDE_DIR} -I${SLOW5_INCLUDE_DIR} -I${POD5_INCLUDE_DIR} 
LIBS=-L${WORKDIR}/zstd/lib/ ${HDF5_LIB_DIR}/lib${HDF5_LIB}.a ${SLOW5_LIB_DIR}/libslow5.a

# LIBS+=${POD5_LIB_DIR}/libarrow.a ${POD5_LIB_DIR}/libjemalloc_pic.a ${POD5_LIB_DIR}/libpod5_format.a
# POD5_LIB_DIR = ${POD5_DIR}/lib
# LIBS+=${POD5_DIR}/lib/libpod5_format.so
LIBS+=${POD5_LIBRARIES}
LIBS+=-lm -lz -ldl -lzstd

PROG=rawhash

ifneq ($(aarch64),)
	arm_neon=1
endif

ifneq ($(arm_neon),) # if arm_neon is defined
ifeq ($(aarch64),)	#if aarch64 is not defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -mfpu=neon -fsigned-char
else				#if aarch64 is defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -fsigned-char
endif
endif

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

ifneq ($(tsan),)
	CFLAGS+=-fsanitize=thread
	LIBS+=-fsanitize=thread
endif

.PHONY: all clean

all: zstd hdf5 slow5 pod5 check_hdf5 check_slow5 check_pod5 $(PROG)
subset: check_zstd check_pod5 check_slow5 check_hdf5 $(PROG)

check_zstd:
	@[ -f "${WORKDIR}/zstd/lib/libzstd.so" ] || { echo "ZSTD library not found" >&2; exit 1; }

check_hdf5:
	@[ -f "${HDF5_INCLUDE_DIR}/H5pubconf.h" ] || { echo "HDF5 headers not found" >&2; exit 1; }
	@[ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.so" ] || [ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.a" ] || { echo "HDF5 library not found" >&2; exit 1; }

check_slow5:
	@[ -f "${SLOW5_INCLUDE_DIR}/slow5/slow5.h" ] || { echo "SLOW5 headers not found" >&2; exit 1; }
	@[ -f "${SLOW5_LIB_DIR}/libslow5.so" ] || [ -f "${SLOW5_LIB_DIR}/libslow5.a" ] || { echo "SLOW5 library not found" >&2; exit 1; }

check_pod5:
	@[ -f "${POD5_INCLUDE_DIR}/pod5_format/c_api.h" ] || { echo "POD5 headers not found" >&2; exit 1; }
	@[ -f "$(WORKDIR)/$(POD5_URLDIR)/$(LIB_DIR)/libpod5_format.a" ] || { echo "POD5 library not found" >&2; exit 1; }

zstd:
	cd ${WORKDIR}/zstd && make -j

hdf5:
	cd ${WORKDIR}/hdf5;\
  	mkdir build;\
	./configure --enable-threadsafe --disable-hl --prefix="${HDF5_DIR}";\
	make -j;\
	make install

slow5:
	make -C ${SLOW5_DIR}

pod5: check_zstd
	cd $(WORKDIR) && mkdir -p $(POD5_URLDIR); cd $(POD5_URLDIR); wget -qO- $(POD5_URL) | tar -xzv;

$(PROG): $(OBJS)
	${CXX} $(CPPFLAGS) $(OBJS) -o $(PROG) $(LIBS)

.SUFFIXES:
.SUFFIXES:.cpp .o

%.o: %.cpp
	${CXX} -c $(CPPFLAGS) $(INCLUDES) $< -o $@

.SUFFIXES:
.SUFFIXES:.c .o

%.o: %.c
	${CC} -c $(CFLAGS) $(INCLUDES) $< -o $@

clean:
	rm -fr *.o $(PROG) *~

rsketch.o: rutils.h kvec.h
rsig.o: hdf5_tools.hpp kvec.h
rmap.o: rawindex.h rsig.h kthread.h kvec.h rutils.h rsketch.h revent.h sequence_until.h
revent.o: roptions.h kalloc.h
rawindex.o: roptions.h rutils.h rsketch.h rsig.h bseq.h khash.h kvec.h kthread.h
main:o rawhash.h ketopt.h pore_model.h
