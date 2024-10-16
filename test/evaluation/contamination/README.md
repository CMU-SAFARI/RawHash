# Contamination Analysis

We assume your current directory points to this directory ([`contamination`](./)).

## RawHash2

Run RawHash2 to map the raw nanopore signals to the corresponding reference genome. This will create the `rawhash2` directory including a `PAF` file that stores the mapping output from RawHash2.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_rawhash2.sh 32
```

## RawHash

Run RawHash to map the raw nanopore signals to the corresponding reference genome. This will create the `rawhash` directory including a `PAF` file that stores the mapping output from RawHash.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_rawhash.sh 32
```

## UNCALLED

Run UNCALLED to map the raw nanopore signals to the corresponding reference genome. This will create the `uncalled` directory including a `PAF` file that stores the mapping output from UNCALLED.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_uncalled.sh 32
```

## Sigmap

Run Sigmap to map the raw nanopore signals to the corresponding reference genome. This will create the `sigmap` directory including a `PAF` file that stores the mapping output from Sigmap.

**Important:** Sigmap will require around 110GB of peak memory during this run.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_sigmap.sh 32
```

## Ground truth mapping (minimap2)

If you have not already done, please follow the `Ground truth mapping (minimap2)` instructions in the [README](../read_mapping/d1_sars-cov-2_r94/README.md) file in [d1_sars-cov-2_r94](../read_mapping/d1_sars-cov-2_r94) to generate the `true_mappings.paf` file that includes the ground truth mapping output from minimap2. Then, run the following command:

```bash
ln -s ../read_mapping/d1_sars-cov-2_r94/true_mappings.paf .
```

## Comparing RawHash2 to RawHash, UNCALLED and Sigmap

After generating the `PAF` files by following each step above, run the following command. This will 1) generate the files to evaluate the mapping output of each tool and 2) output the results (i.e., throughput, mean time per read, indexing time, mapping time, precision, recall, and F1 values).

```bash
cd comparison
bash 0_run.sh

#After running 0_run.sh, if you want to just summarize the results again without generating the evaluation files, you can alternatively run the following command:
# bash 2_output_results.sh

cd ..
```
