#!/usr/bin/env nextflow

/*
========================================================================================
            Identification-free QC workflow for proteomics with OpenMS
========================================================================================
 @#### Authors
 Mathias Walzer <walzer@ebi.ac.uk>
 Yasset Perez-Riverol <ypriverol@gmail.com>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1:   Download a PRIDE experiment from an FTP URL
 - 2:   Converting the RAW data into mzML files
 - 2.1: Extract the metadata for the mzML files and the RAW files
 - 3:   Calculate the QC metrics for each MSrun
 - 4:   Generate the QC metric plots
 ----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run qc-idfree_from_raw/main.nf -c qc-idfree_from_raw/config.nf -profile local

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

process IDfreeQualityControl{
   container 'mwalzer/openms-batteries-included:V2.3.0_pepxmlpatch'
   publishDir "${params.result_folder}", mode: 'copy', overwrite: true

   input:
   file mzML_file from spectraFiles.flatten()

   output:
   file "*.qcML" into qcMLs

   script:
   """
   QCCalculator -in ${mzML_file} -out ${mzML}.qcML
   """
}

qc_tool_parameters = ['auctic', 'charge_histogram', 'esiinstability', 
                     'ms1peakcount', 'ms2peakcount', 'rt_events',  
                     'tic', 'ticric', 'topn',
                     'ms1sn', 'ms2sn', 'sn']

process plotQualityControl{

  container 'mwalzer/qc-plotter:latest'
  publishDir "${params.result_folder}", mode: 'copy', overwrite: true

  input:
  file qcML from qcMLs
  each param from qc_tool_parameters


  output:
  file "*.png" into plots

  script:
  """
  qc_plot.sh -$param ${qcML}
  """

}