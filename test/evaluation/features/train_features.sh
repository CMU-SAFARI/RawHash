#!/bin/bash

sbatch -o w0_bc0_mc10_occ005_all_logistic.out -J w0_bc0_mc10_occ005_all_logistic -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc0_mc10_occ005_all.features logistic logistic"
sbatch -o w0_bc5_mc10_occ005_all_logistic.out -J w0_bc5_mc10_occ005_all_logistic -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc5_mc10_occ005_all.features logistic logistic"

sbatch -o w0_bc0_mc10_occ01_all_logistic.out -J w0_bc0_mc10_occ01_all_logistic -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc0_mc10_occ01_all.features logistic logistic"
sbatch -o w0_bc5_mc10_occ01_all_logistic.out -J w0_bc5_mc10_occ01_all_logistic -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc5_mc10_occ01_all.features logistic logistic"

sbatch -o w0_bc0_mc10_occ005_all_logistic01.out -J w0_bc0_mc10_occ005_all_logistic01 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc0_mc10_occ005_all.features logistic logistic01 0.01"
sbatch -o w0_bc5_mc10_occ005_all_logistic01.out -J w0_bc5_mc10_occ005_all_logistic01 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc5_mc10_occ005_all.features logistic logistic01 0.01"

sbatch -o w0_bc0_mc10_occ01_all_logistic01.out -J w0_bc0_mc10_occ01_all_logistic01 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc0_mc10_occ01_all.features logistic logistic01 0.01"
sbatch -o w0_bc5_mc10_occ01_all_logistic01.out -J w0_bc5_mc10_occ01_all_logistic01 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc5_mc10_occ01_all.features logistic logistic01 0.01"

sbatch -o w0_bc0_mc10_occ005_all_logistic0001.out -J w0_bc0_mc10_occ005_all_logistic0001 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc0_mc10_occ005_all.features logistic logistic0001 0.0001"
sbatch -o w0_bc5_mc10_occ005_all_logistic0001.out -J w0_bc5_mc10_occ005_all_logistic0001 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc5_mc10_occ005_all.features logistic logistic0001 0.0001"

sbatch -o w0_bc0_mc10_occ01_all_logistic0001.out -J w0_bc0_mc10_occ01_all_logistic0001 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc0_mc10_occ01_all.features logistic logistic0001 0.0001"
sbatch -o w0_bc5_mc10_occ01_all_logistic0001.out -J w0_bc5_mc10_occ01_all_logistic0001 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w0_bc5_mc10_occ01_all.features logistic logistic0001 0.0001"

sbatch -o w3_bc0_mc10_occ005_all_logistic.out -J w3_bc0_mc10_occ005_all_logistic -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc0_mc10_occ005_all.features logistic logistic"
sbatch -o w3_bc5_mc10_occ005_all_logistic.out -J w3_bc5_mc10_occ005_all_logistic -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc5_mc10_occ005_all.features logistic logistic"

sbatch -o w3_bc0_mc10_occ01_all_logistic.out -J w3_bc0_mc10_occ01_all_logistic -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc0_mc10_occ01_all.features logistic logistic"
sbatch -o w3_bc5_mc10_occ01_all_logistic.out -J w3_bc5_mc10_occ01_all_logistic -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc5_mc10_occ01_all.features logistic logistic"

sbatch -o w3_bc0_mc10_occ005_all_logistic01.out -J w3_bc0_mc10_occ005_all_logistic01 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc0_mc10_occ005_all.features logistic logistic01 0.01"
sbatch -o w3_bc5_mc10_occ005_all_logistic01.out -J w3_bc5_mc10_occ005_all_logistic01 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc5_mc10_occ005_all.features logistic logistic01 0.01"

sbatch -o w3_bc0_mc10_occ01_all_logistic01.out -J w3_bc0_mc10_occ01_all_logistic01 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc0_mc10_occ01_all.features logistic logistic01 0.01"
sbatch -o w3_bc5_mc10_occ01_all_logistic01.out -J w3_bc5_mc10_occ01_all_logistic01 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc5_mc10_occ01_all.features logistic logistic01 0.01"


sbatch -o w3_bc0_mc10_occ005_all_logistic0001.out -J w3_bc0_mc10_occ005_all_logistic0001 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc0_mc10_occ005_all.features logistic logistic0001 0.0001"
sbatch -o w3_bc5_mc10_occ005_all_logistic0001.out -J w3_bc5_mc10_occ005_all_logistic0001 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc5_mc10_occ005_all.features logistic logistic0001 0.0001"

sbatch -o w3_bc0_mc10_occ01_all_logistic0001.out -J w3_bc0_mc10_occ01_all_logistic0001 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc0_mc10_occ01_all.features logistic logistic0001 0.0001"
sbatch -o w3_bc5_mc10_occ01_all_logistic0001.out -J w3_bc5_mc10_occ01_all_logistic0001 -p gpu_part --gres gpu:1 --wrap="python train_rh2.py w3_bc5_mc10_occ01_all.features logistic logistic0001 0.0001"
