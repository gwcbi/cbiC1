#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --partition=short
#SBATCH -N 1

[[ -z "$srafile" ]] && echo "Must provide srafile" && exit 1
[[ ! -e "$srafile" ]] && echo "$srafile does not exist" && exit 1

# module load sratoolkit
module load sra

fastq-dump --split-3 -F -R --defline-qual "+"  $srafile

echo "complete" && exit 0

