# Sequence Until

We assume your current directory points to this directory ([`sequenceuntil`](./)).

We provide the scripts to:

* Simulate the benefits of Sequence Until without implementing it in UNCALLED and without activating the Sequence Until mechanism in RawHash. The simulation simply shuffles the mapping output that UNCALLED and RawHash generate after their complete relative abundance estimation runs. Then, it picks the portion of the shuffled mapping output to perform relative abundance estimation using only a portion (i.e., 25%, 10%, 1%, 0.1%, and 0.01%) of the entire mapping output.
* Perform the actual Sequence Until mechanism in RawHash.

## Simulating the benefits of Sequence Until

If you have not already done, please perform the full relative abundance estimations for RawHash and UNCALLED as described in the [README](../README.md) file in [relative_abundance](../) before starting this step.

Shuffle the mapping output of RawHash and UNCALLED to generate a relative abundance estimations using the smaller portions (i.e., 25%, 10%, 1%, 0.1%, and 0.01%) of their mapping outputs. This will create multiple `PAF` files that store the mapping output for RawHash and UNCALLED with various using various portions of the entire mapping output as indicated in the suffix of file names. Running `2_output.sh` will output the relative abundance estimations using smaller portions of the entire mapping output.

```bash
cd simulate
bash 1_shuffle.sh
bash 2_output.sh
cd ..
```

### Ground truth mapping (minimap2)

If you have not already done, please follow the `Ground truth mapping (minimap2)` instructions in the [README](..//README.md) file in [relative_abundance](../) to generate the `true_mappings.paf` file that includes the ground truth mapping output from minimap2.

## Running Sequence Until with RawHash

RawHash provides the opportunity to stop the entire sequencing as soon as it detect that the recently sequenced reads do not significantly change the relative abundance estimations, which we call Sequence Until. To activate the Sequence Until mechanism, we provide the `--sequence-until` that can be passed to `rawhash`.

In the mock [random_community](../../../data/random_community/) dataset, we randomly shuffle reads from each real datasets (datasets D1 to D5) to simulate the behavior of a nanopore sequencer where a strand from any genome may be randomly sequenced at a time. We perform relative abundance estimation with and without the Sequence Until mechanism in RawHash.

The following performs relative abundance estimations for these two modes (i.e., with and without Sequence Until) using RawHash. The `PAF` file with `_sequenceuntil` the substring will include the mapping of raw nanopore signals only until the Sequence Until is activated, while the other `PAF` file will include the mapping of all raw nanopore signals.

```bash
bash run_rawhash.sh
```

The following outputs the relative abundance estimations that each mode generates:

```bash
bash output.sh
```
