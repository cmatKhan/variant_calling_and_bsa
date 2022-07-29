# snpEff

snpEff requires that an index be created. A number of indicies are maintained
by the software maintainers, but crypto is not one of them. To create a genome
index for KN99, use the biocontainer singularity image which is used in
the nf-core module. Then, do the following:

1. Go to the snpEff documentation on snpEff build and follow instructions on
how to build a snpEff index from a gtf.

2. Download the gtf and genomic fasta from NCBI. Rename the gtf 'genes.gtf' and
the fasta 'sequences.fa' per the snpEff instructions. Place them in a directory
called ASM221672v1 (the accession number).

3. Enter the singularity shell -- note that you'll need to change the bind
path for your machine

```
> singularity shell -B /home/oguzkhan/Desktop/tmp/bsa_tester -B "$PWD" /home/oguzkhan/Desktop/tmp/bsa_tester/work/singularity/depot.galaxyproject.org-singularity-snpeff-5.1--hdfd78af_2.img /bin/bash -c "cd $PWD"
```

4. Run the following command, again changing the path for your machine

```
snpEff build -dataDir /home/oguzkhan/code/mblab_call_variants/assets/snpeff_db/KN99_WITH_NONCODING/data -gtf22 -v ASM221672v1 -noCheckCds -noCheckProtein
```
