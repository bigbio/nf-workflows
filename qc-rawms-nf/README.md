Quameter IDFree Pipeline
========================

This workflow download a dataset from PRIDE and convert all the files to
 mzML, finally it execute the quameter tool with IDFree option for quality control
 of RAW data previous to Peptide Identification.


Preferquisites:

- Docker: The crux tool is run using a docker container from [BioContainers](http://biocontainers.pro)
- NextFlow: You need to install [nextflow](https://www.nextflow.io/)
