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

log.info """\
 ===================================
  M E T A D A T A   P I P E L I N E
 ===================================
 Project Accession  : ${params.px_accession}
 Metadata Files     : ${params.metadata_path}
 """

process downloadFiles {

    errorStrategy 'retry'
    maxErrors 3

    output:
        file '*.raw' into rawFiles

    script:
    """
    download_raw_files.py $params.px_accession
    """
}

process generateMetadata {
    container 'quay.io/biocontainers/thermorawfileparser:1.1.0--0'

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
    ThermoRawFileParser.sh -i=${rawFile} -m=0 -f=1 -o=./
    """
}

process updateMetadata {

     input:
     file metadataFile from metaResults

     script:
     """
     update_metadata.py $metadataFile $params.pride_username $params.pride_password
     """
}

workflow.onComplete {
	log.info ( workflow.success ? "\nSuccessful" : "Failed to update metadata in $params.px_accession" )
}
