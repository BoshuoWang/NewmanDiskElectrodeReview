#!/bin/bash
#
#SBATCH --job-name=LegBasMat
#SBATCH --output=LegBasMat.out
#SBATCH --error=LegBasMat.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=20G
#SBATCH -c 1

uname -n 1>&2

/usr/bin/time -v /admin/apps/rhel7/matlabR2019a/bin/matlab -nodisplay -r "LegBasMat; quit;"

