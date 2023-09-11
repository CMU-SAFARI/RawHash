# Reproducing the results

## Prerequisites

We compare RawHash2 with [UNCALLED](https://github.com/skovaka/UNCALLED), [Sigmap](https://github.com/haowenz/sigmap), and [RawHash](https://github.com/CMU-SAFARI/RawHash). We specify the versions we use for each tool below.

We list the links to download and compile each tool for comparison below:

* [UNCALLED](https://github.com/skovaka/UNCALLED/tree/74a5d4e5b5d02fb31d6e88926e8a0896dc3475cb)
* [Sigmap](https://github.com/haowenz/sigmap/commit/c9a40483264c9514587a36555b5af48d3f054f6f)
* [RawHash v1](https://github.com/CMU-SAFARI/RawHash/releases/tag/v1.0)

We use minimap2 to generate the ground truth mapping information by mapping basecalled seqeunces to their corresponding reference genomes. We use the following minimap2 version:

* [Minimap2 v2.24](https://github.com/lh3/minimap2/releases/tag/v2.24)

We use various tools to process and analyze the data we generate using each tool. The following tools must also be installed in your machine. We suggest using conda to install these tools with their specified versions as almost all of them are included in the conda repository.

* [Python v3.6.15](https://www.python.org/downloads/release/python-3615/)
* [pip v20.2.3 -- via conda](https://anaconda.org/conda-forge/pip/20.2.3/download/noarch/pip-20.2.3-py_0.tar.bz2)
* [ont_vbz_hdf_plugin v1.0.1 -- via conda](https://anaconda.org/bioconda/ont_vbz_hdf_plugin/files?version=1.0.1)
* [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)

Please make sure that all of these tools are in your `PATH`

## Creating a Virtual Environment (Optional)

You can generate the virtual environment using conda to make sure all the prerequisities are met. Here is an example usage that replicates our virtual environment we use in our evaluations.

```bash
#We will use rawhash2-env directory to download and compile the tools from their repositories
mkdir -p rawhash2-env/bin && cd rawhash2-env

#Important: If you already completed the Step 0 and Step 1 as described, you can skip these steps and add the binaries to your PATH again
#Re-adding the binary to your path is necessary after you start a new shell session.

#Optional Step 0 Creating a conda environment (Note we highly recommend using conda for easy installation of dependencies).
#If not using conda, the packages with the specified versions below (e.g.,  python=3.6.15) must be installed manually in your environment
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda create -n rawhash2-env python=3.6.15 pip=20.2.3 ont_vbz_hdf_plugin=1.0.1
conda activate rawhash2-env

#Installing SRA Toolkit
wget -qO- https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.1/sratoolkit.3.0.1-ubuntu64.tar.gz | tar xzv; cp -r ./sratoolkit.3.0.1-ubuntu64/bin/* bin/; rm -rf sratoolkit.3.0.1-ubuntu64

#Step 1 Compiling the tools
#Cloning and compiling RawHash2
git clone --recursive https://github.com/CMU-SAFARI/RawHash.git rawhash2 && cd rawhash2 && make && cp ./bin/rawhash2 ../bin/ && cd ..

#Cloning and compiling RawHash
wget -qO- https://github.com/CMU-SAFARI/RawHash/releases/download/v1.0/RawHash-1.0.tar.gz | tar -xzv && cd rawhash && make && cp ./bin/rawhash ../bin/ && cd ..

#Cloning and compiling Sigmap
git clone --recursive https://github.com/haowenz/sigmap.git sigmap && cd sigmap && make && cp sigmap ../bin/ && cd ..

#Cloning and compiling UNCALLED
git clone --recursive https://github.com/skovaka/UNCALLED.git uncalled && cd uncalled/submods/bwa && git pull origin master && cd ../../ && pip3 install . && cd ..

#Downloading and compiling minimap2
wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24.tar.bz2; tar -xf minimap2-2.24.tar.bz2; rm minimap2-2.24.tar.bz2; mv minimap2-2.24 minimap2; cd minimap2 && make && cp minimap2 ../bin/ && cd ..

#Step 2 Adding binaries to PATH
#If you are skipping Step 0 and Step 1, uncomment the following line and execute:
# conda activate rawhash2-env
export PATH=$PWD/bin:$PATH
cd rawhash2/test/
```

# Datasets

The scripts and a [README](./data/README.md) can be found in the [data directory](./data/) to download the datasets.

You can start the evaluation step after downloading the datasets as described in [README](./data/README.md).

To start downloading the datasets, enter [data](./data/)

```bash
cd ./data
```

# Evaluation

## Read Mapping

The scripts and a [README](./evaluation/read_mapping/README.md) can be found in the [read mapping directory](./evaluation/read_mapping/) to perform read mapping using UNCALLED, Sigmap, RawHash, and RawHash2.

To start performing the read mapping evaluations, enter [./evaluation/read_mapping](./evaluation/read_mapping/)

```bash
cd ./evaluation/read_mapping
```

## Relative Abundance Estimation

The scripts and a [README](./evaluation/relative_abundance/README.md) can be found in the [relative abundance directory](./evaluation/relative_abundance/) to perform the relative abundance estimation using UNCALLED, Sigmap, RawHash, and RawHash2.

To start performing the relative abundance estimation evaluations, enter [./evaluation/relative_abundance](./evaluation/relative_abundance/)

```bash
cd ./evaluation/relative_abundance
```

## Contamination Analysis

The scripts and a [README](./evaluation/contamination/README.md) can be found in the [contamination directory](./evaluation/contamination/) to perform the contamination analysis using UNCALLED, Sigmap, RawHash, and RawHash2.

To start performing the contamination analysis evaluations, enter [./evaluation/contamination](./evaluation/contamination/)

```bash
cd ./evaluation/contamination/
```

## Sequence Until

The scripts and a [README](./evaluation/relative_abundance/sequenceuntil/README.md) can be found in the [sequence until directory](./evaluation/relative_abundance/sequenceuntil/) to evaluate the benefits of the Sequence Until mechanism using UNCALLED and RawHash2.

To start evaluating Sequence Until, enter [./evaluation/relative_abundance/sequenceuntil](./evaluation/relative_abundance/sequenceuntil/)

```bash
cd ./evaluation/relative_abundance/sequenceuntil
```
