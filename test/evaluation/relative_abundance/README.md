# Relative Abundance Estimation

We assume your current directory points to this directory ([`relative_abundance`](./)).

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

## Comparing RawHash to UNCALLED and Sigmap

After generating the `PAF` files by following each step above, run the following command. This will 1) generate the files to evaluate the mapping output of each tool and 2) output the results (i.e., throughput, mean time per read, indexing time, mapping time, relative abundance estimations, precision, recall, and F1 values).

```bash
cd comparison
bash 0_run.sh

#After running 0_run.sh, if you want to just summarize the results again without generating the evaluation files, you can alternatively run the following command:
# bash 2_output_results.sh

cd ..
```

# Sequence Until
