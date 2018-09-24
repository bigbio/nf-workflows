# Proteomics and Multi-omics workflows

This repository store different workflows to analysis proteomics and multi-omics data using [BioContainers](biocontainers.pro) and [Nextflow](nextflow.io). Each workflow resolve one specific proteomics experiment and
should be seen as a tutorial to proteomics, biocontainers and nextflow. We are reciving Issues and PR from contributors. It is important to notice that all pipelines will be run using local setup for nextflow, but this can be easily addapted to other [architectures or executors](https://www.nextflow.io/docs/latest/executor.html)

- [Mass spectrometry Search](ms-crux-id-nf) : This example perform a protein database search unsing an ms2 (mass spectra file) and protein database (fasata database).

- [Converting the raw data from RAW to mzML](thermo-convert-nf) : This example allow you to dowload the raw data from PRIDE ftp and convert the files to mzML and get the json file with the metadata.

