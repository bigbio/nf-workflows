#!/usr/bin/env nextflow

/*
========================================================================================
                 Proteogenomics Custom database creation
========================================================================================
 Authors
 Yasset Perez-Riverol <ypriverol@gmail.com>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1:   Download the Fasta protein database from ENSEMBL (taxonomy provided).
 - 2:   Add lncRNA, ORFs, sORF   
 - 3:   Add ENSEMBL variant population mutations. 
 - 4:   For cancer specific studies add the mutations from COSMIC or cbioportal studies 
 - 5:   Clear Working directories
 ----------------------------------------------------------------------------------------
*/

/*
 * Define the default parameters
 */

params.taxonomy = "9606"
params.ensembl_downloader_config = "${baseDir}/configs/ensembl_downloader_config.yaml"

ensembl_downloader_config = file(params.ensembl_downloader_config)

process ensembl_protein_fasta_download(){
    
    container 'quay.io/bigbio/pypgatk:0.0.1'
    
    input: 
    file ensembl_downloader_config

    output:
	file "database_ensembl/*.gz" into ensembl_fasta_databases

	script:
	"""
	pypgatk_cli.py ensembl-downloader --config_file "${ensembl_downloader_config}" --taxonomy ${params.taxonomy}
	"""

    
    
}

