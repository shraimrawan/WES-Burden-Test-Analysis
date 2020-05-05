#!/bin/sh
#$ -cwd
export LD_LIBRARY_PATH=/home/shraimr/miniconda3/lib
/scr1/users/shraimr/vt/vt decompose veo.remicade.062419.updated.vcf -o veo.remicade.decomposed.101819.vcf
