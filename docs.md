---
title:  Configurations
filename: docs.md
--- 


Here are more detailts about conda. 

## Conda

The conda configurations are in `conda.conf`.

Here the process selectors labels:

Labels: 
- `merge_tables`: labels for all the merge metadata processes. Only requires python and pandas
- `r_download`: Label for all the download processes, they more or less all require the same packages. All r packages
  useful for download (TCGAbiolinks, recount3, nsdc, optparse, GenomicDataCommons, tidyverse). 
- `prepare_expression`: prepare recount data (ecount3, recount, data.table, readr (tidyverse))
- `prepare_methylation`: clean and prepare methylation ( data.table, nsdc, stringr(tiyverse), recount3, optparse)
- `netzoopy_panda`: run panda
- `netzoopy_dragon`: run dragon
- `netzoopy_pandalioness`: run pandalioness
- `netzoopy_dragonlionss`: run dragon lioness

For the netzoo functions we are keeping the selectors separated to allow
managing the GPU. Indeed, we recommend using a GPU to process LIONESS networks, 
but in that case the conda environment heavily depends on the cuda version, 
hence we recommend building the environment manually.

Also, for the moment there are some packages that need to be installed manually, 
which is managed by the process themselves. These packages are NetSciDataCompanion, 
...


For further instructions on how to use conda with nextflow 
and to specify your own configurations read the [nextflow
docs](https://www.nextflow.io/docs/latest/config.html#config-conda). 

## Docker


## Tests

## Result folders

Here the structure of all the output folder, as they are created in the tests.
The UUID names in these will likely change for your own analysis, but the structure will 
remain the same. 


### Results: Download


```bash
⚡ tree results/test_download 
results/test_download
├── downloaded_methylation_metadata.csv
├── downloaded_mutation_metadata.csv
├── downloaded_recount_metadata.csv
├── gtex_pancreas
│   └── data_download
│       └── recount3
│           └── gtex_pancreas.rds
└── tcga_paad
    └── data_download
        ├── clinical
        │   ├── clinical_drug_paad.csv
        │   ├── clinical_follow_up_v4.4_nte_paad.csv
        │   ├── clinical_follow_up_v4.4_paad.csv
        │   ├── clinical_nte_paad.csv
        │   ├── clinical_omf_v4.0_paad.csv
        │   ├── clinical_patient_paad.csv
        │   └── clinical_radiation_paad.csv
        ├── methylation
        │   ├── tcga_paad_methylation_manifest.txt
        │   ├── tcga_paad_methylation_metadata.csv
        │   └── tcga_paad_methylations.txt
        ├── mutations
        │   ├── tcga_paad_mutations.txt
        │   ├── tcga_paad_mutations_metadata.csv
        │   └── tcga_paad_mutations_pivot.csv
        └── recount3
            └── tcga_paad.rds


```

### Prepare


This would be the output when you run the `-profile testPrepare` test case

```bash
results/test_prepare
└── tcga_paad
    └── data_prepared
        ├── methylation
        │   ├── tcga_paad_tf_promoter_methylation_clean_all.log
        │   └── tcga_paad_tf_promoter_methylation_raw.csv
        └── recount3
            ├── recount3_pca_tcga_paad_purity01_normtpm_mintpm1_fracsamples00001_tissueall_batchnull_adjnull.png
            ├── recount3_pca_tcga_paad_purity01_normtpm_mintpm1_fracsamples02_tissueall_batchnull_adjnull.png
            ├── recount3_tcga_paad_purity01_normtpm_mintpm1_fracsamples00001_tissueall_batchnull_adjnull.log
            ├── recount3_tcga_paad_purity01_normtpm_mintpm1_fracsamples00001_tissueall_batchnull_adjnull.rds
            ├── recount3_tcga_paad_purity01_normtpm_mintpm1_fracsamples00001_tissueall_batchnull_adjnull.txt
            ├── recount3_tcga_paad_purity01_normtpm_mintpm1_fracsamples02_tissueall_batchnull_adjnull.log
            ├── recount3_tcga_paad_purity01_normtpm_mintpm1_fracsamples02_tissueall_batchnull_adjnull.rds
            └── recount3_tcga_paad_purity01_normtpm_mintpm1_fracsamples02_tissueall_batchnull_adjnull.txt

```

### Analyze

This would be the output when you run the `-profile testAnalyze` test case

```bash
results/test_analyze
└── tcga_paad
    └── analysis
        ├── dragon
        │   ├── tcga_paad_dragon.log
        │   ├── tcga_paad_dragon_filtered_expression.csv
        │   ├── tcga_paad_dragon_input.tsv
        │   └── tcga_paad_dragon_mat.tsv
        ├── lioness_dragon
        │   ├── lioness_dragon
        │   │   ├── lioness-dragon-TCGA-2L-AAQL-01A.csv
        │   │   ├── lioness-dragon-TCGA-HV-A5A3-11A.csv
        │   │   ├── lioness-dragon-TCGA-HV-A7OP-01A.csv
        │   │   ├── lioness-dragon-TCGA-IB-A5ST-01A.csv
        │   │   ├── lioness-dragon-TCGA-US-A779-01A.csv
        │   │   └── lioness-dragon-TCGA-YB-A89D-11A.csv
        │   ├── recount3_tcga_paad_purity01_normtpm_mintpm1_fracsamples02_tissueall_batchnull_adjnull_shuffle.rds
        │   ├── tcga_paad_lioness_dragon.log
        │   ├── tcga_paad_tf_promoter_methylation_clean.csv
        │   └── tcga_paad_tf_promoter_methylation_clean_shuffle.csv
        ├── panda
        │   ├── panda_tcga_paad.log
        │   └── panda_tcga_paad.txt
        └── panda_lioness
            ├── lioness
            │   ├── lioness.TCGA-2L-AAQL-01A-11R-A38C-07.4.h5
            ...
            └── panda.txt
```

### Full 


Here is an example of the output of the full test.
`tcga_luad` is the uuid of the run , hence inside the `tcga_luad` folder you'll find all relevant files
divided in `data_download`, `data_prepare`, `analysis` folders, that correspond to the steps of the pipeline

Here the general out structure. Below we have the full expected output

```bash
results/test_full
└── tcga_luad
    ├── analysis
    │   ├── dragon
    │   ├── lioness_dragon
    │   │   └── lioness_dragon
    │   ├── panda
    │   └── panda_lioness
    │       └── lioness
    ├── data_download
    │   ├── clinical
    │   ├── methylation
    │   ├── mutations
    │   └── recount3
    └── data_prepared
        ├── methylation
        └── recount3
```