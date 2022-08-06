# tcga-data-nf

![](https://github.com/QuackenbushLab/tcga-data-nf/workflows/build/badge.svg)

Workflow to download and prepare TCGA data
This replaces the old recount3-access-nf from Kate.


## Configuration

- param1: this is the parameter description (default: "hello")
- ...
- paramN: this is the parameter description (default: "flow")

## Running the workflow

### Install or update the workflow

```bash
nextflow pull QuackenbushLab/tcga-data-nf
```

### Run the analysis

```bash
nextflow run QuackenbushLab/tcga-data-nf
```

## Results

- `results/analysis.txt`: file with the analysis results.
- `results/tuning.txt`: file with the parameter tuning.
- ...

## Authors

- Kate Hoff Shutta
- Viola Fanfani
- Panagiotis Mandros
