# Read Mapping - d4_green_algae_r94

We assume your current directory points to this directory ([`d4_green_algae_r94`](./)).

## RawHash2

Run RawHash2 to map the raw nanopore signals to the corresponding reference genome. This will create the `rawhash2` directory including a `PAF` file that stores the mapping output from RawHash2.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_rawhash2.sh 32
```

## RawHash

Run RawHash to map the raw nanopore signals to the corresponding reference genome. This will create the `rawhash` directory including a `PAF` file that stores the mapping output from RawHash.

**Important:** RawHash will require around 12GB of peak memory during this run.

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

**Important:** Sigmap will require around 30GB of peak memory during this run.

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

## Comparing RawHash2 to RawHash, UNCALLED and Sigmap

After generating the `PAF` files by following each step above, run the following command. This will 1) generate the files to evaluate the mapping output of each tool and 2) output the results (i.e., throughput, mean time per read, indexing time, mapping time, precision, recall, and F1 values).

```bash
cd comparison
bash 0_run.sh

#After running 0_run.sh, if you want to just summarize the results again without generating the evaluation files, you can alternatively run the following command:
# bash 2_output_results.sh

cd ..
```

### Profiling RawHash2

We have created a script to profile the main steps of RawHash2 (i.e., I/O, signal-to-event conversion, sketching, seeding, chaining, and the entire mapping). Before running the following script, RawHash2 needs to be compiled with the profiling mode. Please compile as follows:

```bash
# For profiling
make PROFILE=1
```

After compiling with the profiling mode and adding RawHash2 to your path, you can run the following script to generate the profilied runtimes. The profiling result will be directed to `stderr`.

```bash
bash profile_rawhash2.sh 1
```