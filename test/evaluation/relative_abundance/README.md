# Relative Abundance Estimation

We assume your current directory points to this directory ([`relative_abundance`](./)).

## RawHash2

Run RawHash2 to map the raw nanopore signals to the corresponding reference genome. This will create the `rawhash2` directory including a `PAF` file that stores the mapping output from RawHash2.

**Important:** RawHash2 will require around 150GB of peak memory during this run.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_rawhash2.sh 32
```

## RawHash

Run RawHash to map the raw nanopore signals to the corresponding reference genome. This will create the `rawhash` directory including a `PAF` file that stores the mapping output from RawHash.

**Important:** RawHash will require around 150GB of peak memory during this run.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_rawhash.sh 32
```

## UNCALLED

Run UNCALLED to map the raw nanopore signals to the corresponding reference genome. This will create the `uncalled` directory including a `PAF` file that stores the mapping output from UNCALLED.

**Important:** UNCALLED will require around 50GB of peak memory during this run.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_uncalled.sh 32
```

## Sigmap

Run Sigmap to map the raw nanopore signals to the corresponding reference genome. This will create the `sigmap` directory including a `PAF` file that stores the mapping output from Sigmap.

**Important:** Sigmap will require around 490GB of peak memory during this run.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_sigmap.sh 32
```

## Ground truth mapping (minimap2)

Run minimap2 to map the basecalled sequences to the corresponding reference genome. This will create the `true_mappings.paf` file that includes the ground truth mapping output from minimap2.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_minimap2.sh 32
```

Note that this ground truth mapping is used to calculate the precision, recall, and F-1 score. For estimating the relative abundance, we use the actual number of bases of the basecalled sequences used in this particular use case. To calculate the number of bases, we suggest using `seqkit stat` on each dataset and generating the ratio of each genome based on the total number of bases from each dataset.

## Comparing RawHash2 to RawHash, UNCALLED and Sigmap

After generating the `PAF` files by following each step above, run the following command. This will 1) generate the files to evaluate the mapping output of each tool and 2) output the results (i.e., throughput, mean time per read, indexing time, mapping time, relative abundance estimations, precision, recall, and F1 values).

```bash
cd comparison
bash 0_run.sh

#After running 0_run.sh, if you want to just summarize the results again without generating the evaluation files, you can alternatively run the following command:
# bash 2_output_results.sh

cd ..
```

# Sequence Until

You can test the real Sequence Until mechanism and the simulated benefits of Sequence Until using the following commands:

```bash

#Running RawHash2 with the Sequence Until mechanism activated. The mapping will stop as soon as RawHash2 decides that further sequencing is not necessary.
cd sequenceuntil
bash run_rawhash2.sh
cd ..
```

To simulate the behavior of Sequence Until on the fully mapping output from UNCALLED and RawHash2:

```bash

#Running RawHash2 with the Sequence Until mechanism activated. The mapping will stop as soon as RawHash2 decides that further sequencing is not necessary.
cd sequenceuntil/simulate/
bash 1_shuffle.sh
bash 2_output.sh
cd ..
```
