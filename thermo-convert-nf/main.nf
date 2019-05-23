#!/usr/bin/env nextflow

/*
========================================================================================
                 Metadata Extraction Workflow for PRIDE Data
========================================================================================
 @#### Authors
 Yasset Perez-Riverol <ypriverol@gmail.com>
 Suresh Hewapathirana <sureshhewabi@gmail.com>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1:   Find the PRIDE FTP project path and Download all the raw files from FTP
 - 2:   Converting the RAW data into mzML files and extract the metadata for the mzML files and the RAW files
 - 3:   Format metadata into proper JSON format
 - 4:   Push metadata to the database
 - 5:   Clear Working directories

 ----------------------------------------------------------------------------------------
*/

/*
 * Define the default parameters
 */
params.px_accession = ""
params.pride_username = ""
params.pride_password = ""
params.metadata_path = ""
params.mode = ""
params.files_location = ""
params.ftp_download = false

log.info """\
 ===================================
  M E T A D A T A   P I P E L I N E
 ===================================
 Project Accession  : ${params.px_accession}
 Metadata Files     : ${params.metadata_path}
 Submission Mode    : ${params.mode}
 Raw Files Location : ${params.files_location}
 Enable FTP download: ${params.ftp_download}
 """

process downloadFiles {

    container 'quay.io/pride/pride-py:0.0.8'

    errorStrategy 'retry'
    maxErrors 3

    output:
    file '*.{raw,RAW,Raw}' into rawFiles

    script:
    """
    python3 /pridepy.py download -a $params.px_accession -o ./ -f $params.ftp_download -i $params.files_location
    """
}

process generateMetadata {

    container 'ypriverol/thermorawfileparser:1.1.0'

    memory { 10.GB * task.attempt }
    errorStrategy 'retry'
    queue 'production-rh7'
    publishDir "$params.metadata_path", mode:'copy', overwrite: true

    input:
    file rawFile from rawFiles.flatten()

    output:
    file '*.json' into metaResults mode flatten
    file '*.mzML' into spectraFiles

    script:
    """
    ThermoRawFileParser -i=${rawFile} -m=0 -f=1 -o=./ --ignoreInstrumentErrors
    """
}

process updateMetadata {

     container 'quay.io/pride/pride-py:0.0.8'

     input:
     file metadataFile from metaResults

     when:
     params.mode != 'test'

     script:
     """
     python3 /pridepy.py update-metadata -f $metadataFile -u $params.pride_username -p $params.pride_password
     """
}

workflow.onComplete {
	log.info ( workflow.success ? "\nSuccessful" : "Failed to update metadata in $params.px_accession" )
}
