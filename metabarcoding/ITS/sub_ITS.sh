#!/bin/bash

#$ -q short.q
#$ -N sub_ITS
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 12
#$ -l mem_free=6G
#$ -V
#$ -cwd

module purge
module load bioinfo/vsearch/2.9.1

mkdir /homedir/galati/data/metab/ITS/SUB/
SUB=/homedir/galati/data/metab/ITS/SUB/
PRIM=/homedir/galati/data/metab/ITS/PRIM/

cd ${PRIM}
dir ${PRIM} R*.fastq.gz > filenames.txt

echo "Création d'un fichier avec les noms d'échantillons"

for NAME in `awk '{print $1}' ${PRIM}filenames.txt`

do

echo "Sous-échantillonnage" ${NAME}

vsearch --fastx_subsample ${NAME} --fastqout ${SUB}$NAME --sample_size 1000 --randseed 3

done

echo "Suppression des fichiers vides"
find /homedir/galati/data/metab/ITS/SUB/ -empty -type f -delete

# JOB END
date

exit 0
