---
title:  Extend workflow
filename: extend_workflow.md
--- 

# Extend Workflow

Being modular, the workflow allows users to add steps. We do not expect or recommend 
to drastically change the structure and length of the workflow, but we would welcome any bug-fix,
"small" addition of GRNs tools, quality control steps.

**BUGS**: please open an issue on Github if you encounter an error. If you know how to solve it, please open a PR with
the bug fix. Refer to the fork-PR-merge cycle section below for guidance on how to do that.

**EXTENSIONS**: first, we would recomment opening an issue or send an email to the maintainers 
to gauge interest in the tool or feature you want to add. If we deem the feature useful for the
broader community we'll ask you to open a PR, following the GENIE3 example below. If, instead, 
you just want to have a personal version of the workflow, we recommend you fork the project 
and appropriately cite it in your work.


## Example: add GENIE3


#### Create a script

You need to have a script ( or a bash command ) that runs the tool of interest. 
This needs to accept data, in the format of the pipeline, as input and save the output results.

For instance, we create the `run_genie3.r` script, in the `bin/r` folder that takes gene expression as input, and runs
the [GENIE3](https://bioconductor.org/packages/release/bioc/html/GENIE3.html) method. 

#### Create (or add) a process to the modules

We have added the `modules/grns.nf` file and added the runGENIE3 process.

You can use the other processes as reference, for instance, we have used the runTCGAPanda process.

Also, add the relevant configuration options. 

For GENIE3, we now need 
- params.genie3.tf_list
- params.genie3.n_cores
- params.genie3.tree_method
- params.genie3.n_trees
- params.genie3.k

#### Add the process to the workflow

You now need to have this step added to the workflow.

We have added GENIE3 to the analyzeExpressionWf, since it only requires
gene expression to be run, and we have also addedd a flag genie3.run_genie3 
that should be set to true if you want the pipeline to include genie3. 

Here is the call inside the workflow:
```{nextflow}
    if (params.run_genie3){
        runGENIE3(data.map{it -> tuple(it[0], it[1])})
    }
```

#### Add testdata and process configurations

Each process needs to be tested. For that, add the relevant test files
inside the `testdata` folder and the configurations into the `test.config` file

For GENIE3 we have addedd a reduced set of TFs that are used to compute the regulatory relationship with the genes and
we select `n_trees=2` such that the inference step is faster. 

#### Modify image/configurations

Add a conda selector and the tool to the docker.

There are many ways that would allow you to test whether the workflow can run, 
one we can recommend for development purposes is the following: 
- Download the latest docker image 
- Run the container from the image and install the required packages (in this case we installed GENIE3 through Bioconductor)
- Commit the container into a new local image [howto](https://stackoverflow.com/questions/63027514/install-package-in-running-docker-container)
- Run the workflow with the updated image
