#!/bin/bash

#$ -q short.q
#$ -N sub_16S
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd

mkdir /homedir/galati/data/16S/16S_sub/
SUB=/homedir/galati/data/16S/16S_sub/
TRIM=/homedir/galati/data/16S/16S_primer_trimmed2/
cd ${TRIM}
dir ${TRIM} R*.fastq.gz > filenames.txt

echo "Création d'un fichier avec les noms d'échantillons"

for NAME in `awk '{print $1}' ${TRIM}listfile.txt`

do

echo "Sous-échantillonnage" ${NAME}

vsearch --fastx_subsample ${NAME} --fastqout ${SUB}$NAME --sample_size 1000 --randseed 3

done

