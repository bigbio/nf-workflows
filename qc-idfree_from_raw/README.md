Identification-free QC with OpenMS
========================

This workflow downloads full dataset from PRIDE FTP and convert all the files to
 mzML. It computes the a set of metrics using the OpenMS QCCalculator tool and plots metrics graphs for each MSRun.


## How to run 

```bash
   
  nextflow main.nf -c nextflow.config -c local        ## If you want to run locally 
  
```
## Prerequisites:

- Docker (For local computer)
- Singularity for cluster setup 

### For local computer: 

- Docker: The OpenMS tool QCCalculator is run using a docker container from [BioContainers](http://biocontainers.pro)
- Docker: QC metric plots are created using an R container with ggplot
- NextFlow: You need to install [nextflow](https://www.nextflow.io/)

FAQ
==========

- If you are are running in your local computer (PC), please be aware that the conversion and QCCalculator step may consume more than 5GB memory depending of the file size. 

