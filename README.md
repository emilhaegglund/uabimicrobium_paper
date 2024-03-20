This repository contains code to produce the results published in "Phylogeny and Expansion of Serine/Threonine Kinases in Phagocytotic Bacteria in the Phylum Planctomycetota".

The analysis are run using Snakemake workflows and to ensure that correct version of software are used they are fetched from conda. The workflows can be run using
```
snakemake -s <name_of_workflow.smk> --use-conda --conda-frontend mamba -j <number_of_threads>
```

__Author:__ Emil Hägglund  
__Contact:__ emil.hagglund@icm.uu.se  
__DOI:__ 10.17044/scilifelab.25429441  
__License:__ CC BY 4.0  
This readme file was last updated: 20-03-2024  
