#pass DEBUG=1 to make to enable debugging
ifdef DEBUG
	CFLAGS=-Wall -O2 -fsanitize=address -g -Wc++-compat -fopenmp -march=native -pthread -DHAVE_KALLOC
	CPPFLAGS=--std=c++11 -Wall -O2 -fsanitize=address -g -fopenmp -march=native -pthread -DHAVE_KALLOC
else
	CFLAGS=-Wall -O3 -Wc++-compat -fopenmp -march=native -pthread -DHAVE_KALLOC
	CPPFLAGS=-std=c++11 -Wall -O3 -fopenmp -march=native -pthread -DHAVE_KALLOC
endif

#pass PROFILE=1 to make to enable profiling
ifdef PROFILE
	CFLAGS+=-g -fno-omit-frame-pointer -march=native -DPROFILERH=1
	CPPFLAGS+=-g -fno-omit-frame-pointer -march=native -DPROFILERH=1
endif

OBJS= kthread.o kalloc.o bseq.o roptions.o sequence_until.o rutils.o rsig.o revent.o rsketch.o rindex.o lchain.o rseed.o rmap.o hit.o main.o

CXX_COMPILER_VERSION ?= $(shell $(CXX) -dumpversion)
SYSTEM_PROCESSOR ?= $(shell uname -m)
SYSTEM_NAME ?= $(shell uname -s)

WORKDIR = $(shell pwd)/../extern

POD5_VERSION = 0.2.2
POD5_URLDIR = pod5-$(POD5_VERSION)-$(SYSTEM_NAME)
POD5_REPO = https://github.com/nanoporetech/pod5-file-format
ifdef POD5_INCLUDE_DIR
	POD5_DOWNLOAD = 0
endif
POD5_INCLUDE_DIR ?= $(WORKDIR)/$(POD5_URLDIR)/include

ifeq ($(SYSTEM_NAME),Linux)
    POD5_LIB ?= lib64
    ifeq ($(shell echo $(SYSTEM_PROCESSOR) | grep -E '^(aarch64|arm)'),)
        POD5_URL = $(POD5_REPO)/releases/download/$(POD5_VERSION)/lib_pod5-$(POD5_VERSION)-linux-x64.tar.gz
    else ifeq ($(shell expr $(CXX_COMPILER_VERSION) \< 9.0),1)
        POD5_URLDIR = pod5-$(POD5_VERSION)-$(SYSTEM_NAME)
        POD5_URL = $(POD5_REPO)/releases/download/$(POD5_VERSION)/lib_pod5-$(POD5_VERSION)-linux-gcc7-arm64.tar.gz
        POD5_LIB ?= lib
    else
        POD5_URLDIR = pod5-$(POD5_VERSION)-$(SYSTEM_NAME)
        POD5_URL = $(POD5_REPO)/releases/download/$(POD5_VERSION)/lib_pod5-$(POD5_VERSION)-linux-arm64.tar.gz
    endif
	POD5_LIB_DIR ?= $(WORKDIR)/$(POD5_URLDIR)/$(POD5_LIB)
    POD5_LIBRARIES ?= $(POD5_LIB_DIR)/libpod5_format.a \
                     $(POD5_LIB_DIR)/libarrow.a \
                     $(POD5_LIB_DIR)/libjemalloc_pic.a
endif

ifeq ($(SYSTEM_NAME),Darwin)
    POD5_URL = $(POD5_REPO)/releases/download/$(POD5_VERSION)/lib_pod5-$(POD5_VERSION)-osx-11.0-arm64.tar.gz
	POD5_LIB_DIR ?= $(WORKDIR)/$(POD5_URLDIR)/lib
    POD5_LIBRARIES ?= $(POD5_LIB_DIR)/libpod5_format.a \
                     $(POD5_LIB_DIR)/libarrow.a
endif

ifeq ($(SYSTEM_NAME),Windows_NT)
    POD5_URL = $(POD5_REPO)/releases/download/$(POD5_VERSION)/lib_pod5-$(POD5_VERSION)-win
endif

HDF5_DIR ?= ${WORKDIR}/hdf5/
HDF5_BUILD_DIR ?= ${HDF5_DIR}/build
ifdef HDF5_INCLUDE_DIR
	HDF5_COMPILE = 0
endif
HDF5_INCLUDE_DIR ?= ${HDF5_BUILD_DIR}/include
HDF5_LIB_DIR ?= ${HDF5_BUILD_DIR}/lib
HDF5_LIB ?= hdf5

SLOW5_DIR ?= ${WORKDIR}/slow5lib/
ifdef SLOW5_INCLUDE_DIR
	SLOW5_COMPILE = 0
endif
SLOW5_INCLUDE_DIR ?= ${SLOW5_DIR}/include
SLOW5_LIB_DIR ?= ${SLOW5_DIR}/lib

#pass NOPOD5=1 to make to disable compiling with POD5
ifeq ($(NOPOD5),1)
	CFLAGS+=-DNPOD5RH=1
	CPPFLAGS+=-DNPOD5RH=1
else
	INCLUDES+=-I${POD5_INCLUDE_DIR}
	LIBS+=-L${WORKDIR}/zstd/lib/ ${POD5_LIBRARIES} -lzstd
endif

#pass NOHDF5=1 to make to disable compiling with HDF5
ifeq ($(NOHDF5),1)
	CFLAGS+=-DNHDF5RH=1
	CPPFLAGS+=-DNHDF5RH=1
else
	INCLUDES+=-I${HDF5_INCLUDE_DIR}
	LIBS+=${HDF5_LIB_DIR}/lib${HDF5_LIB}.a 
endif

#pass NOSLOW5=1 to make to disable compiling with SLOW5
ifeq ($(NOSLOW5),1)
	CFLAGS+=-DNSLOW5RH=1
	CPPFLAGS+=-DNSLOW5RH=1
else
	INCLUDES+=-I${SLOW5_INCLUDE_DIR}
	LIBS+=${SLOW5_LIB_DIR}/libslow5.a
endif

LIBS+=-lm -lz -ldl 

PROG=rawhash2

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

.PHONY: all subset clean help

all: zstd hdf5 slow5 pod5 check_hdf5 check_slow5 check_pod5 $(PROG) ## Build RawHash with HDF5, SLOW5, and POD5
subset: check_zstd check_pod5 check_slow5 check_hdf5 $(PROG) ## Build RawHash without recompiling HDF5, SLOW5, and POD5

help: ##Show help
	@echo "Set PROFILE=1 to enable profiling"
	@echo "Set DEBUG=1 to enable debugging"
	@echo
	@echo "Set NOHDF5=1 to disable compiling with HDF5 (note that at least one of HDF5, SLOW5, and POD5 must be enabled)"
	@echo "Set NOSLOW5=1 to disable compiling with SLOW5 (note that at least one of HDF5, SLOW5, and POD5 must be enabled)"
	@echo "Set NOPOD5=1 to disable compiling with POD5 (note that at least one of HDF5, SLOW5, and POD5 must be enabled)"
	@echo
	@echo "use \"make subset\" to prevent compiling HDF5, SLOW5, and POD5 libraries from scratch if they are already built"
	@echo "use \"make clean\" to remove all generated files"
	
check_zstd:
	@if [ "$(NOPOD5)" != "1" ]; then \
		[ -f "${WORKDIR}/zstd/lib/libzstd.so" ] || { echo "ZSTD library not found" >&2; exit 1; }; \
	fi

check_hdf5:
	@if [ "$(NOHDF5)" != "1" ]; then \
		[ -f "${HDF5_INCLUDE_DIR}/H5pubconf.h" ] || { echo "HDF5 headers not found" >&2; exit 1; }; \
		[ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.so" ] || [ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.a" ] || { echo "HDF5 library not found" >&2; exit 1; }; \
	fi

check_slow5:
	@if [ "$(NOSLOW5)" != "1" ]; then \
		[ -f "${SLOW5_INCLUDE_DIR}/slow5/slow5.h" ] || { echo "SLOW5 headers not found" >&2; exit 1; }; \
		[ -f "${SLOW5_LIB_DIR}/libslow5.so" ] || [ -f "${SLOW5_LIB_DIR}/libslow5.a" ] || { echo "SLOW5 library not found" >&2; exit 1; }; \
	fi

check_pod5:
	@if [ "$(NOPOD5)" != "1" ]; then \
		[ -f "${POD5_INCLUDE_DIR}/pod5_format/c_api.h" ] || { echo "POD5 headers not found" >&2; exit 1; }; \
		[ -f "$(POD5_LIB_DIR)/libpod5_format.a" ] || { echo "POD5 library not found" >&2; exit 1; }; \
	fi

zstd:
	@if [ "$(NOPOD5)" != "1" ]; then \
		cd ${WORKDIR}/zstd && make -j; \
	fi

hdf5:
	@if [ "$(NOHDF5)" != "1" ] && [ "$(HDF5_COMPILE)" != "0" ]; then \
		cd ${HDF5_DIR}; \
		mkdir -p build; \
		./configure --enable-threadsafe --disable-hl --prefix="${HDF5_BUILD_DIR}"; \
		make -j; \
		make install; \
	fi

slow5:
	@if [ "$(NOSLOW5)" != "1" ] && [ "$(SLOW5_COMPILE)" != "0" ]; then \
		make -C ${SLOW5_DIR}; \
	fi

pod5: check_zstd
	@if [ "$(NOPOD5)" != "1" ] && [ "$(POD5_DOWNLOAD)" != "0" ]; then \
		cd $(WORKDIR) && mkdir -p $(POD5_URLDIR); cd $(POD5_URLDIR); wget -qO- $(POD5_URL) | tar -xzv; \
	fi

ifeq ($(NOHDF5),1)
FASTDEPCOM=${CC}
FASTDEPCOM_FLAGS=$(CFLAGS)
else
FASTDEPCOM=${CXX}
FASTDEPCOM_FLAGS=$(CPPFLAGS)
endif

# $(PROG): $(OBJS)
# 	${FASTDEPCOM} $(FASTDEPCOM_FLAGS) $(OBJS) -o $(PROG) $(LIBS)

$(PROG): $(OBJS)
	${CXX} $(CPPFLAGS) $(OBJS) -o $(PROG) $(LIBS) -lstdc++

rsig.o: rsig.c
	$(FASTDEPCOM) -c $(FASTDEPCOM_FLAGS) $(INCLUDES) $< -o $@

main.o: main.cpp
	$(FASTDEPCOM) -c $(FASTDEPCOM_FLAGS) $(INCLUDES) $< -o $@

rmap.o: rmap.c
	$(FASTDEPCOM) -c $(FASTDEPCOM_FLAGS) $(INCLUDES) $< -o $@

rindex.o: rindex.c
	$(FASTDEPCOM) -c $(FASTDEPCOM_FLAGS) $(INCLUDES) $< -o $@

.SUFFIXES:
.SUFFIXES:.cpp .o

%.o: %.cpp
	${CXX} -c $(CPPFLAGS) $(INCLUDES) $< -o $@

.SUFFIXES:
.SUFFIXES:.c .o

%.o: %.c
	${CC} -c $(CFLAGS) $(INCLUDES) $< -o $@

clean: ## Remove all generated files
	rm -fr *.o $(PROG) *~

rsketch.o: rutils.h rh_kvec.h
lchain.o: kalloc.h rutils.h rseed.h rsketch.h chain.h krmq.h
hit.o: chain.h kalloc.h khash.h
rsig.o: hdf5_tools.hpp rh_kvec.h
rseed.o: rsketch.h kalloc.h rutils.h rindex.h
hit.o: rmap.h kalloc.h khash.h
rmap.o: rindex.h rsig.h kthread.h rh_kvec.h rutils.h rsketch.h revent.h sequence_until.h
revent.o: roptions.h kalloc.h
rindex.o: roptions.h rutils.h rsketch.h rsig.h bseq.h khash.h rh_kvec.h kthread.h
main:o rawhash.h ketopt.h rutils.h
