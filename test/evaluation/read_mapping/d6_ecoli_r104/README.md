# Read Mapping - d6_ecoli_r104 (R10)

We assume your current directory points to this directory ([`d6_ecoli_r104`](./)).

## RawHash

Run RawHash to map the raw nanopore signals to the corresponding reference genome. This will create the `rawhash` directory including a `PAF` file that stores the mapping output from RawHash.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_rawhash.sh 32
```

<!-- ### Profiling RawHash

We have created a script to profile the main steps of RawHash (i.e., I/O, signal-to-event conversion, sketching, seeding, chaining, and the entire mapping). Before running the following script, RawHash needs to be compiled with the profiling mode. Please modify the corresponding lines of this [Makefile](../../../../src/Makefile) as follows and then recompile RawHash:

```bash
# For profiling
CFLAGS+=-DPROFILERH=1
CPPFLAGS+=-DPROFILERH=1
```

After compiling with the profiling mode and adding RawHash to your path, you can run the following script to generate the profilied runtimes. The profiling result will be directed to `stderr`.

```bash
bash profile_rawhash.sh 1
``` -->

## UNCALLED

UNCALLED currently does not support R10 reads .

## Sigmap

Sigmap currently does not support R10 reads.

## Ground truth mapping (minimap2)

Run minimap2 to map the basecalled sequences to the corresponding reference genome. This will create the `true_mappings.paf` file that includes the ground truth mapping output from minimap2.

The following command will use 32 threads. you can change the maximum threads to use by providing a different value than 32 below.

```bash
bash run_minimap2.sh 32
```

## Comparing RawHash to UNCALLED and Sigmap

Comparison is not possible. However, in order to generate the RawHash results with the existing scripts, we generate the symbolic links of UNCALLED an Sigmap results in the `comparison` directory.

After generating the `PAF` files by following each step above, run the following command. This will 1) generate the files to evaluate the mapping output of RawHash and 2) output the results (i.e., throughput, mean time per read, indexing time, mapping time, precision, recall, and F1 values).

```bash
cd comparison
bash 0_run.sh

#After running 0_run.sh, if you want to just summarize the results again without generating the evaluation files, you can alternatively run the following command:
# bash 2_output_results.sh

cd ..
```
