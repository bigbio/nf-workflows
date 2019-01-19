#!/usr/bin/env nextflow

params.px_accession = "PXD000801"
params.file_output = 1

process downloadFiles {
    container 'quay.io/biocontainers/gnu-wget:1.18--3'

    output:
    file '*.raw' into rawFiles

    script:
    """
    wget -v -r -nd -A "*.raw" --no-host-directories --cut-dirs=1 ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/09/PXD010376/
    """
}

process generateMetadata {
    container 'ypriverol/thermorawfileparser:0.1'

    publishDir 'data/', mode:'copy'

    input:
    file rawFile from rawFiles.flatten()

    output:
    file '*.json' into metaResults
    file '*.mzML' into spectraFiles

    script:
    """
    ThermoRawFileParser -i=${rawFile} -m=0 -f=1 -o=./
    """
}