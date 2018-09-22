Simple Mass Spectrometry Identification Workflow
========================

This first workflow take as an input an ms2 spectra file and an fasta file and perform a peptide search using the [Crux search engine](http://crux.ms/). This is a simple workflow to demostrate how to perform peptide search using the popular workflow manager [nextflow](https://www.nextflow.io/).

Preferquisites: 

- Docker: The crux tool is run using a docker container from [BioContainers](http://biocontainers.pro)
- NextFlow: You need to install [nextflow](https://www.nextflow.io/)
