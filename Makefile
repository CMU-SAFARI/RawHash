.PHONY: all subset clean help

all:rawhash2

help: ##Show help
	+$(MAKE) -C src help
	
rawhash2:
	@if [ ! -e bin ] ; then mkdir -p ./bin/ ; fi
	+$(MAKE) -C src
	mv ./src/rawhash2 ./bin/

subset:
	@if [ ! -e bin ] ; then mkdir -p ./bin/ ; fi
	+$(MAKE) -C src subset
	mv ./src/rawhash2 ./bin/
	
clean:
	rm -rf bin/
	+$(MAKE) clean -C ./src/
	
