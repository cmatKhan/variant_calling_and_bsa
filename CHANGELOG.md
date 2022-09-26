# nf-core/mblabcallvariants: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.2 - [20220925]

### `Fixed`

CNVpytor is now producing tsv output. The module scripts, in particular view,
have been heavily modified from the nf-core versions. Only `tsv` and the pytor
file are produced as output the `tsv` format is hard coded, so the current
output format param for cnvpytor isn't doing anything.

## v1.1.1 - [20220918]
### `Added`
`freebayes_min_base_qual` to the config. Set defaults on mapq, base_qual and
coverage in nextflow.config

## v1.1.0 - [20220914]
### `Added`
Added CNVpytor. However, while the module runs, the result is empty

## v1.0.0 - [20220726]
Initial release. BSA, genotype_check, htcf and kn99 configs working.

Initial release of nf-core/mblabcallvariants, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
