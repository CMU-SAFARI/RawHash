.PHONY: all subset clean help

all:rawhash2

HOSTNAME := $(shell hostname)
ifeq ($(HOSTNAME),mikado1)
binary_name = rawhash_pybinding
else
binary_name = rawhash2
endif

help: ##Show help
	+$(MAKE) -C src help
	
rawhash2:
	@if [ ! -e bin ] ; then mkdir -p ./bin/ ; fi
	+$(MAKE) -C src
#	mv ./src/rawhash2 ./bin/
	mv ./src/$(binary_name) ./bin/

subset:
	@if [ ! -e bin ] ; then mkdir -p ./bin/ ; fi
	+$(MAKE) -C src subset
	mv ./src/rawhash2 ./bin/
	
clean:
	rm -rf bin/
	+$(MAKE) clean -C ./src/
	
