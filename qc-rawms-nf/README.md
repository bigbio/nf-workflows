Identification-free QC  
========================

This workflow downloads full dataset from PRIDE FTP and convert all the files to
 mzML. It computes the a set of metrics using the QuaMeter tool and perform the statistical
 analysis (PCA) and Correlation between each MSRun.


## How to run 

```bash
   
  nextflow main.nf -c nextflow.config -c local        ## If you want to run locally 
  
```
## Prerequisites:

- Docker (For local computer)
- Singularity for cluster setup 

### For local computer: 

- Docker: The crux tool is run using a docker container from [BioContainers](http://biocontainers.pro)
- NextFlow: You need to install [nextflow](https://www.nextflow.io/)

FAQ
==========

- If you are are running in your local computer (PC), please be aware that the conversion step and 
QuaMeter step consumes more than 5GB memory depending of the file size. 

