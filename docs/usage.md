# nf-core/mblabcallvariants: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/mblabcallvariants/usage](https://nf-co.re/mblabcallvariants/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet

Not all of these columns will be relevant to your set of sequencing files --
for the time being, you will need to put entries into them. For instance, if
`group`, `pool`, `cond`, `day`, `replicate`, `experiment` don't apply to your
samples, you could put the value `1` into each of those entries. Critically,
if you want to run all samples together as a batch through `freebayes` (more
or less recommended), then make sure that the `group` field is all the same.

__important__: if you have single end reads, simply leave fastq_2 blank (the
column header does need to exist, though)

```console
sample,group,pool,cond,day,replicate,experiment,runNumber,fastq_1,fastq_2
A1-35-8,1,A1,lung,8,1,BSA2,5000,assets/test_data/A1-35-8_R1.fastq.gz,assets/test_data/A1-35-8_R2.fastq.gz
A2-102-5,1,A2,inoculum,15,2,BSA2,5000,assets/test_data/A2-102-5_R1.fastq.gz,assets/test_data/A2-102-5_R2.fastq.gz
A5-35-17,2,A5,brain,0,1,BSA2,5000,assets/test_data/A5-35-17_R1.fastq.gz,assets/test_data/A5-35-17_R2.fastq.gz
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `group` | string or numeral, eg 1 or g1, which groups a sample. For instance, in the example samplesheet, there are two samples in group 1 and only 1 in group 2. Groups with more than 1 will go through additional variant calling steps in which the group is processed together.                                                             |
| `pool` | From Daniel's metadata, the pool of a given sample.                                                             |
| `cond` | The tissue or sample condition. Currently one of inoculum, ypd, lung or brain.                                                             |
| `day` | number of days before sample is extracted and sequenced.                                                             |
| `replicate` | replicate number from Daniel's metadata.                                                             |
| `experiment` | eg BSA2 or BSA6, for instance.                                                             |
| `runNumber` | eg 4869. This should mean that the raw sequencing files may be found in run_4869 somewhere on `lts`.                                                             |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

See these links for example sample sheets:
[example BSA samplesheet](../assets/bsa_samplesheet.csv)  and
[example genotype check samplesheet](../assets/genotype_check_samplesheet.csv).

