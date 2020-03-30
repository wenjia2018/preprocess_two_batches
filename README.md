
# Functionality

Here we describe preprocessing steps to generate the basic addhealth RNA expression set. The R functions take a) raw count RNA and b) questionaire ("phenotype") data files and return c) a single standardized [expression set object](http://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) containing both RNA and phenotype data. This object can then interface easily with many downstream bioconductor pipelines. 


# Details of functionality

 - init_rna.R: load raw counts; drop duplicate samples or subjects; create expressionSet object "dat.rds" also containing phenotype data. (Analogously, init_rna_steve.R, does the same but encorportes not the raw counts but, reference-gene normalized counts, c.f. Steve Cole).
 - init_pheno.R: load and merge 5 waves of questionaire data; construct variables of interest to us; write them to disk as "waves.rds".
 - init_data.R: filters genes without ENSG id, haemoglobin genes and genes with insufficient variation for statistical analysis.

# Preprocessing for addhealth RNA and phenotype data

To ease collaboration we provide a portable, reproducible analysis pipeline - using the functionality described below - in a docker image [here](https://hub.docker.com/repository/docker/chumbleycode/preprocess). A more detailed version of the preprocessing is available in pdf [here](https://github.com/chumbleycode/preprocess/blob/master/docs/methods.pdf). A youtube introduction to docker is [![here]()](http://www.youtube.com/watch?v=YFl2mCHdv24).

We assume you have
1. access to the addhealth phenotype and raw RNA counts (tranche 1 of wave 5). This means you have access to the longleaf (UNC) or jccompute (UZH) servers.
2. copied all relevant input data files into one folder: this folder will be called "data_input" inside the docker container 
3. cloned this "preprocess" repository. This will be your working directory: it will be called "workspace" inside the docker container 
4. adjusted docker memory, as per figure below (or if in the terminal using [these commands](https://docs.docker.com/config/containers/resource_constraints/)).


![](screenshot_memory.png)


To open the docker container type the following at the terminal the terminal type:

``` 
docker run -it -v\
/absolute/path/to/data_input:/home/rstudio/data_input \
-v /absolute/path/to/preprocess:/home/rstudio/workspace \
chumbleycode/preprocess:latest bash
```

You should now be in a linux operating system with all the correct packages installed. To reproduce our cleaned datasets use
```
Rscript /home/rstudio/workspace/R/init.R
```

This will write the cleaned files to /home/studio/workspace/data. Alternatively, develop with your own variants of preprocessing, in the workspace folder.

 
