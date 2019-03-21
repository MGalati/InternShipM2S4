#!/bin/bash
#$ -q normal.q
#$ -N mock_ITS
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd

echo "Activation du module conda"
module purge
module load system/conda/5.1.0

echo "Activation de l'environnement conda"
conda activate qiime2-2018.11
