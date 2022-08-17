[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

## Introduction

This pipeline is meant to align genomic reads and perform variant calling against
the reference genome. It can align any number of samples either separately,
together, or in user specified batches.

## Dependencies

If you are in the Brent lab, Nextflow and Singularity are already installed.
If you are in a different lab, would like to use this, and don't know how to
install `git`, `Nextflow >= 20.10.03` and `Singularity`,
please email chasem@wustl.edu.


## Setup

Navigate to your scratch space, create a directory to store pipeline output,
and clone this repo. This look like so:

```
$ cd /scratch/<lab>/<your_scratch>

# load the git module
$ eval $(spack load --sh git)

# it doesn't matter what this is called -- call it anything
$ mkdir variant_calling_pipeline

$ cd variant_calling_pipeline

$ git clone https://github.com/cmatKhan/variant_calling_and_bsa

```


## Test the installation

To test whether the pipeline is appropriately installed, do the following:

```bash
# make sure you are in the directory where you want the output to go. This is
# probably the same directory in which you have the repo. Everything that
# follows assumes this is the case.
$ cd /scratch/<lab>/<user_scratch>/variant_calling_pipeline

# copy the HTCF run script to your directory
$ cp variant_calling_and_bsa/assets/run_nf_htcf.sh .

# launch the test. Note that this script does assume that the code repo is
# in the current directory
$ sbatch run_nf_htcf.sh

```

This will run the pre-configured test data via slurm. You can monitor progress
by using `squeue -u $USER`. Once the jobs start to get scheduled, there will be
a file created in your current directory called `variant_calling_test.out`. You
can watch progress of the pipeline by entering `tail -150 variant_calling_test.out`.
That file is updated as data moves through the pipeline. It will also record
any errors. It should, however, track the data as it moves through. Once all jobs
are complete, that file will display a 'successful completion' message. Throughout
the workflow, intermediate files are deposited in a subdirectory of your current
directory called `work`, and the final results are deposited in a directory called
`results`. Do not delete either of these until you the pipeline has completely
finished (either due to completion or error. If it is finished, there will be
no pipeline jobs in your `squeue -u $USER` output related to the pipeline).

## Running your own data

You will need to create a samplesheet. [See instructions here](docs/usage.md).
Next, move the data from `lts` to `scratch`. I put this into a directory
called `data` in the same directory that I store the pipeline output in. In the
instructions above, this directory was called
`/scratch/<lab>/<user_scratch>/variant_calling_pipeline`.

Create a launch script that is similar to the test script called `run_nf_htcf.sh`.
It might look like so:

```
#!/usr/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -J bsa6.out
#SBATCH -o bsa6.out

eval $(spack load --sh openjdk)
# note that the versioning isn't strictly necessary --
# update this to the newest version. Just don't go older, in general
eval $(spack load --sh singularityce@3.8.0)
eval $(spack load --sh nextflow@22.04.5)

tmp=$(mktemp -d /tmp/$USER-singularity-XXXXXX)

mkdir singularity
mkdir local_tmp

export NXF_SINGULARITY_CACHEDIR=singularity
export SINGULARITY_TMPDIR=$tmp
export SINGULARITY_CACHEDIR=$tmp

nextflow run \
    variant_calling_and_bsa/main.nf \
    -profile singularity,htcf,batch_only \
    --input bsa6_samplesheet.csv \
    --outdir results_bsa6 \
    --aligners "bwamem2" \
    --ploidy 2
    -resume
```

You could alternatively put the input parameters into a json like so:

```json
{
  "fasta": "/ref/mblab/data/KN99/KN99_genome_fungidb.fasta",
  "input": "bsa6_samplesheet.csv",
  "outdir": "results_bsa6",
  "aligners" : "bwamem2",
  "ploidy": 2
}

```

and the launch command would look like:

```
nextflow run \
    variant_calling_and_bsa/main.nf \
    -profile singularity,htcf,batch_only \
    -params-file params.json
    -resume
```



## Citation

Coming soon -- need to cite daniel's paper

## Developer Notes

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool
to run tasks across multiple compute infrastructures in a very portable manner.
It uses Docker/Singularity containers making installation trivial and results
highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
implementation of this pipeline uses one container per process which makes it
much easier to maintain and update software dependencies. Where possible,
these processes have been submitted to and installed from
[nf-core/modules](https://github.com/nf-core/modules) in order to make them
available to all nf-core pipelines, and to everyone within the Nextflow community!

## EVERYTING BELOW IS TEMPLATE

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/mblabcallvariants/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run nf-core/mblabcallvariants -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```console
   nextflow run nf-core/mblabcallvariants --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Documentation

The nf-core/mblabcallvariants pipeline comes with documentation about the pipeline [usage](https://nf-co.re/mblabcallvariants/usage), [parameters](https://nf-co.re/mblabcallvariants/parameters) and [output](https://nf-co.re/mblabcallvariants/output).

## Credits

nf-core/mblabcallvariants was originally written by Chase Mateusiak.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#mblabcallvariants` channel](https://nfcore.slack.com/channels/mblabcallvariants) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/mblabcallvariants for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
