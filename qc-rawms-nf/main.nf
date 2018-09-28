#!/usr/bin/env nextflow
 
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

process qcIdFreeMzMLs{
    container 'quay.io/biocontainers/bumbershoot:3_0_11579--0'

    publishDir 'data/', mode:'copy'

    input:
    file filemzML from spectraFiles.flatten()

    output:
    file '*.tsv' into qcMzMLFiles

    script:
    """
    quameter ${filemzML} -MetricsType idfree -OutputFilepath ${filemzML}-qc.tsv -cpus 5
    """
}
