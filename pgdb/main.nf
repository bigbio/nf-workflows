#!/usr/bin/env nextflow

/*
========================================================================================
                 Proteogenomics Custom database creation
========================================================================================
 Authors
 - Yasset Perez-Riverol <ypriverol@gmail.com>
 - Husen M. Umer <husensofteng@gmail.com>
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

container_path = "path-to-source/py-pgatk/pypgatk/" //temp solution until we put a container on quay.io

params.release = "97"
params.taxonomy = "9606"
params.ensembl_downloader_config = "${baseDir}/configs/ensembl_downloader_config.yaml"
ensembl_downloader_config = file(params.ensembl_downloader_config)

params.lncrna = true 
params.ensembl_lncrna_config = "${baseDir}/configs/ensembl_config.yaml"
ensembl_lncrna_config = file(params.ensembl_lncrna_config)

params.pseudogenes = true 
params.ensembl_pseudogenes_config = "${baseDir}/configs/ensembl_config.yaml"
ensembl_pseudogenes_config = file(params.ensembl_pseudogenes_config)

protein_decoy_config = "${baseDir}/configs/protein_decoy.yaml 

params.final_database_protein = "final_database_protein.fa"
params.final_decoy_database_protein = "final_database_protein_decoy.fa"
params.decoy_prefix = "decoy_"


params.cosmic_config = "${baseDir}/configs/cosmic_config.yaml"
cosmic_config = file(params.cosmic_config)

params.cosmic_user_name = "userName"
params.cosmic_password = "passWord"

/* Pipeline START */

/** 
 * Download data from ensembl for the particular species. 
 */ 
process ensembl_protein_fasta_download(){
    
    container 'quay.io/bigbio/pypgatk:0.0.1'
    
    input: 
    file ensembl_downloader_config

    output:
	file "database_ensembl/*.gz" into ensembl_fasta_gz_databases

	script:
	"""
	pypgatk_cli.py ensembl-downloader --config_file "${ensembl_downloader_config}" --taxonomy ${params.taxonomy} --release ${params.release}
	"""
}

/**
 * Decompress all the data downloaded from ENSEMBL 
 */ 
process gunzip_ensembl_files{

    container 'quay.io/bigbio/pypgatk:0.0.1'
    publishDir "result", mode: 'copy', overwrite: true

    input: 
    file(fasta_file) from ensembl_fasta_gz_databases

    output: 
    file '*pep.all.fa' into ensembl_protein_database
    file '*cdna.all.fa' into ensembl_cdna_database, ensembl_cdna_database_sub
    file '*ncrna.fa' into ensembl_ncrna_database, ensembl_ncrna_database_sub

    script: 
    """
    gunzip -d --force ${fasta_file}
    """ 
}

/* Concatenate cDNA and ncRNA databases */
process merge_cdans{
  
  container 'quay.io/bigbio/pypgatk:0.0.1' 

  input:
  file a from ensembl_cdna_database_sub
  file b from ensembl_ncrna_database_sub
  
  output: 
  file '*.fa' into total_cdans

  script: 
  """
  cat "${a}" "${b}" >> total_cdans.fa
  """ 
}

(lncrna_cdna, cdna_lncrna) = ( !params.lncrna ? [Channel.empty(), total_cdans] : [total_cdans, total_cdans] ) 

process add_lncrna {

  container 'quay.io/bigbio/pypgatk:0.0.1'
  publishDir "result", mode: 'copy', overwrite: true

  input:
  file x from lncrna_cdna
  file ensembl_lncrna_config

  output:
  file('*.fa') into optional_lncrna

  script:
  """
  pypgatk_cli.py dnaseq-to-proteindb --config_file "${ensembl_lncrna_config}" --input_fasta ${x} --output_proteindb proteindb_from_lncRNAs_DNAseq.fa --include_biotypes '3prime_overlapping_ncrna, ambiguous_orf, antisense, antisense_RNA, lincRNA, ncrna_host, processed_transcript, sense_intronic, sense_overlapping' --skip_including_all_cds
  """
}

(pseudogenes_cdna, cdna_pseudogenes) = ( !params.pseudogenes ? [Channel.empty(), cdna_lncrna] : [cdna_lncrna, cdna_lncrna] ) 

/**
 * Creates the pseudogenes protein database. 
 * The default configuration for pypgatk tool is: 
 *  - disrupted_domain 
 *  - IGC_pseudogene
 *  - IGJ_pseudogene 
 *  - IG_pseudogene 
 *  - IGV_pseudogene
 *  - processed_pseudogene
 *  - transcribed_processed_pseudogene
 *  - transcribed_unitary_pseudogene
 *  - transcribed_unprocessed_pseudogene
 *  - translated_processed_pseudogene
 *  - TRJ_pseudogene
 *  - unprocessed_pseudogene  
 */
process add_pseudogenes {

  container 'quay.io/bigbio/pypgatk:0.0.1'
  publishDir "result", mode: 'copy', overwrite: true

  input:
  file x from pseudogenes_cdna
  file ensembl_pseudogenes_config

  output:
  file('*.fa') into optional_pseudogenes

  script:
  """
  pypgatk_cli.py dnaseq-to-proteindb --config_file "${ensembl_pseudogenes_config}" --input_fasta "${x}" --output_proteindb proteindb_from_pseudogenes_DNAseq.fa --include_biotypes 'disrupted_domain, IGC_pseudogene, IGJ_pseudogene, IG_pseudogene, IGV_pseudogene, processed_pseudogene, transcribed_processed_pseudogene, transcribed_unitary_pseudogene, transcribed_unprocessed_pseudogene, translated_processed_pseudogene, TRJ_pseudogene, unprocessed_pseudogene' --skip_including_all_cds
  """
}

/**
 * Merge all the protein databases, for example: ensembl proteins + lncrna + pseudogenes  
 */
process merge_databases_final_proteindb {

    container 'quay.io/bigbio/pypgatk:0.0.1' 
    publishDir "result", mode: 'copy', overwrite: true 

    input:
    file a from ensembl_protein_database
    file b from optional_lncrna
    file c from optional_pseudogenes

    output: 
    file '*.fa' into final_databases

    script: 
    """
    cat "${a}" "${b}" "${c}" >> "${params.final_database_protein}"
    """ 

}


/**
 * Create the decoy database using DecoyPYrat
 * Decoy sequences will have "_DECOY" prefix tag to the protein accession.
 */
process create_decoy_db {
	container 'quay.io/bigbio/pypgatk:0.0.1'
    publishDir "result", mode: 'copy', overwrite: true 

	input:
	file x from final_databases

	output:
	file "*.fasta" into fasta_decoy_db

	script:
	"""
	python pypgatk_cli.py generate-decoy --config_file "${protein_decoy_config}" --input "${x}" --decoy_prefix "${params.decoy_prefix}" --output protein_decoy_database.fa
	"""
}


/* Mutations to proteinDB */

/*** COSMIC Mutations ***/
process cosmic_download {
	
	input:
	file cosmic_config
	
	output:
	file "database_cosmic/*.gz" into cosmic_files
	
	script:
	"""
	python ${container_path}pypgatk_cli.py cosmic-downloader --config_file "${cosmic_config}" --username ${params.cosmic_user_name} --password ${params.cosmic_password}
	"""
}

/** Decompress all the data downloaded from COSMIC ***/ 
process gunzip_cosmic_files{

    publishDir "result", mode: 'copy', overwrite: true

    input: 
    file(data_file) from cosmic_files

    output: 
    file "All_COSMIC_Genes.fasta" into cosmic_genes 
	file "CosmicMutantExport.tsv" into cosmic_mutations
	
	script:
    """
    gunzip -d --force ${data_file}
    """
}

/*** generate proteindb from cosmic mutations ***/
process cosmic_proteindb{
	
	input:
	file g from cosmic_genes
	file m from cosmic_mutations
	
	output:
	file 'cosmic_proteinDB.fa' into cosmic_proteindb
	
	script:
	"""
	python ${container_path}pypgatk_cli.py cosmic-to-proteindb --config_file "${cosmic_config}" --input_mutation ${m} --input_genes ${g} --output_db cosmic_proteinDB.fa --split_by_tissue_type
	"""

}

