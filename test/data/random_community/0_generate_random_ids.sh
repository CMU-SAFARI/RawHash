#!/bin/bash

echo '#!/bin/bash' > 1_symbolic_links.sh;
echo >> 1_symbolic_links.sh;
echo 'mkdir -p fast5_files' >> 1_symbolic_links.sh;
echo >> 1_symbolic_links.sh;
echo 'cd fast5_files' >> 1_symbolic_links.sh;
echo >> 1_symbolic_links.sh;
shuf -n50000 read_ids.txt | awk 'BEGIN{id=0}{print "ln -s "$0" read"count++".fast5"}' >> 1_symbolic_links.sh

