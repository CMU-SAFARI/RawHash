.PHONY: all subset clean help

all:rawhash

help: ##Show help
	+$(MAKE) -C src help
	
rawhash:
	@if [ ! -e bin ] ; then mkdir -p ./bin/ ; fi
	+$(MAKE) -C src
	mv ./src/rawhash ./bin/

subset:
	@if [ ! -e bin ] ; then mkdir -p ./bin/ ; fi
	+$(MAKE) -C src subset
	mv ./src/rawhash ./bin/
	
clean:
	rm -rf bin/
	+$(MAKE) clean -C ./src/
	