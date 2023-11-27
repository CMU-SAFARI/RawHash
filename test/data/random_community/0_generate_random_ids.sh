#!/bin/bash

#Note: The multi fast5 files should be converted to single fast5 files and stored in "single_fast5_files" folders inside each dataset. We do this to make sure the reads in a multi-fast5 file are not processed sequentially, but rather to simulate the environment where each read from any sample can be sequenced at any random order. For D4 and D3, the datasets are already in single fast5 files. For further details regarding each read name extracted from each fast5 file, please refer to the read_ids.txt file.

ln -s
echo '#!/bin/bash' > 1_symbolic_links.sh;
echo >> 1_symbolic_links.sh;
echo 'mkdir -p fast5_files' >> 1_symbolic_links.sh;
echo >> 1_symbolic_links.sh;
echo 'cd fast5_files' >> 1_symbolic_links.sh;
echo >> 1_symbolic_links.sh;
shuf read_ids.txt | awk 'BEGIN{id=0}{print "ln -s "$0" read"count++".fast5"}' >> 1_symbolic_links.sh
