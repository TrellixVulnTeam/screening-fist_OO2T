#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate sxfst # rdkit 

python pipeline.py ../nb/00-pilot/00.3
