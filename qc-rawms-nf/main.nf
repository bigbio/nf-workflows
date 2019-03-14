#!/usr/bin/env nextflow

/*
========================================================================================
                 Identification-free QC workflow for proteomics
========================================================================================
 @#### Authors
 Yasset Perez-Riverol <ypriverol@gmail.com>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1:   Download a PRIDE experiment from an FTP URL
 - 2:   Converting the RAW data into mzML files
 - 2.1: Extract the metadata for the mzML files and the RAW files
 - 3:   Extract the QC metrics for each RAW (MSrun) files
 - 4:   Generate the statistics about the Experiment and individual Runs
 - 5:   Output the HTML report for the experiment
 ----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run qc-rawms-nf/main.nf -c qc-rawms-nf/config.nf -profile local

    Mandatory arguments:

    --project_ftp The project ftp folder in PRIDE (e.g. ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2015/11/PXD003133)

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

/*
 * Download files from FTP
 */
process downloadFiles {
    container 'quay.io/biocontainers/gnu-wget:1.18--3'
    
    output:
    file '*.raw' into rawFiles
   
    script: 
    """
    wget -v -r -nd -A "*.raw" --no-host-directories --cut-dirs=1 ${params.project_ftp}
    """ 
}

/*
 * Generate the mzML + metadata for each RAW file
 */
process generateMetadata {

    label 'big_mem'
    container 'ypriverol/thermorawfileparser:0.2'
    memory { 10.GB * task.attempt }
    errorStrategy 'retry'
    
    publishDir 'data/', mode:'copy'
 
    input:
    file rawFile from rawFiles.flatten()
    
    output: 
    file '*.json' into metaResults
    file '*.mzML' into spectraFiles

    script:
    """
    ThermoRawFileParser -i=${rawFile} -m=0 -f=1 -o=./ -v
    """
}

/**
 * Compute QuaMeter metrics for each RAW files
 */
process qcIdFreeMzMLs{

    label 'big_mem'

    container 'quay.io/biocontainers/bumbershoot:3_0_11579--0'
    publishDir 'data/', mode:'copy'

    input:
    file fileMzML from spectraFiles.flatten()

    output:
    file '*.tsv' into qcMzMLFiles

    script:
    """
    quameter ${fileMzML} -MetricsType idfree -OutputFilepath ${fileMzML}-qc.tsv -cpus 5
    """
}

