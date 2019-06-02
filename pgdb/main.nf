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
 - 1:   Download the Fasta protein database from ENSEMBL (default = 9606).
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

params.lncrna = true 
params.ensembl_lncrna_config = "${baseDir}/configs/ensembl_config.yaml"
ensembl_lncrna_config = file(params.ensembl_lncrna_config)

process ensembl_protein_fasta_download(){
    
    container 'quay.io/bigbio/pypgatk:0.0.1'
    
    input: 
    file ensembl_downloader_config

    output:
	file "database_ensembl/*.gz" into ensembl_fasta_gz_databases

	script:
	"""
	pypgatk_cli.py ensembl-downloader --config_file "${ensembl_downloader_config}" --taxonomy ${params.taxonomy}
	"""
}

process gunzip_ensembl_files{

    container 'quay.io/bigbio/pypgatk:0.0.1'
    publishDir "result", mode: 'copy', overwrite: true

    input: 
    file(fasta_file) from ensembl_fasta_gz_databases

    output: 
    file '*pep.all.fa' into ensembl_protein_database
    file '*cds.all.fa' into ensembl_cds_database
    file '*ncrna.fa' into ensembl_ncrna_database

    script: 
    """
    gunzip -d --force ${fasta_file}
    """ 
}

(lncrna_cds, cds) = ( !params.lncrna 
                 ? [Channel.empty(), ensembl_cds_database] 
                 : [ensembl_cds_database, ensembl_cds_database] ) 

process add_lncrna {
  
  container 'quay.io/bigbio/pypgatk:0.0.1'
  publishDir "result", mode: 'copy', overwrite: true

  input:
  file x from lncrna_cds
  file ensembl_lncrna_config

  output:
  file('*.fa') into optional_lncrna

  script:
  """
  pypgatk_cli.py dnaseq-to-proteindb --config_file "${ensembl_lncrna_config}" --input_fasta ${x} --output_proteindb proteindb_from_lncRNAs_DNAseq.fa --include_biotypes '3prime_overlapping_ncrna, ambiguous_orf, antisense, antisense_RNA, lincRNA, ncrna_host, processed_transcript, sense_intronic, sense_overlapping' --skip_including_all_cds
  """
}




