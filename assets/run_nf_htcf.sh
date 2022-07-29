#!/usr/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -J variant_calling_test.out
#SBATCH -o variant_calling_test.out

eval $(spack load --sh openjdk)
# note that the versioning isn't strictly necessary --
# update this to the newest version. Just don't go older, in general
eval $(spack load --sh singularityce@3.8.0)
eval $(spack load --sh nextflow)

tmp=$(mktemp -d /tmp/$USER-singularity-XXXXXX)

mkdir singularity
mkdir local_tmp

export NXF_SINGULARITY_CACHEDIR=singularity
export SINGULARITY_TMPDIR=$tmp
export SINGULARITY_CACHEDIR=$tmp

# note that this isn't bulletproof -- it relies on the assumption that
# the second item in the list of groups of which your $USER on htcf is
# associated is the lab group
lab_group=$(groups | cut -d ' ' -f2)

nextflow run \
    /ref/${lab_group}/software/variant_calling_and_bsa/main.nf \
    -profile test_full,singularity,htcf \
    -resume
