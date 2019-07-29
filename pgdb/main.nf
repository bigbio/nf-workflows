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
 	Download the Fasta protein database from ENSEMBL (default species = 9606), latest release
 	Generate ncRNA, psudogenes, ORFs databases   
 	Download ENSEMBL variants (VCFs, default species = 9606)
 	Generate ENSEMBL variant protein database
 	Download gnomAD variants
 	Generate gnomAD variant protein database
 	Download mutations from COSMIC
 	Generate COSMIC mutated protein databases (one database for all mutations, a database per tumor type)
 	Download mutations from cBioPortal (default = all stuides)
 	Generate cBioPortal mutated protein databases (one database for all mutations, a database per tumor type)
 	Concatenate all databases
 	Generate a corresponding decoy database for the concatenated database
 ----------------------------------------------------------------------------------------
*/

/*
 * Define the default parameters
*/

params.ncrna = true 
params.pseudogenes = true
params.altorfs = true
params.ensembl = true
params.gnomad = true
params.cbioportal = true
params.cosmic = true

params.tool_basepath = "./"
container_path = "${params.tool_basepath}/py-pgatk/pypgatk/" //temp solution until we put a container on quay.io

//params.release = "97" //not supported yet
params.taxonomy = "9103"
params.ensembl_downloader_config = "${baseDir}/configs/ensembl_downloader_config.yaml"
ensembl_downloader_config = file(params.ensembl_downloader_config)

params.af_field = "" //set to empty when AF_field does not exist in the INFO filed or filtering on AF is not desired

params.ensembl_config = "${baseDir}/configs/ensembl_config.yaml"
ensembl_config = file(params.ensembl_config)

params.cbioportal_config = "${baseDir}/configs/cbioportal_config.yaml"
cbioportal_config = file(params.cbioportal_config)
	
/* Biotype groups according to: 
 * https://www.ensembl.org/Help/Faq?id=468 and 
 * http://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html
*/
biotypes = [
	'protein_coding': "protein_coding,polymorphic_pseudogene,non_stop_decay,nonsense_mediated_decay,IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene,TEC", 
	'pseudogene': "pseudogene,IG_C_pseudogene,IG_J_pseudogene,IG_V_pseudogene,IG_pseudogene,TR_V_pseudogene,TR_J_pseudogene,processed_pseudogene,rRNA_pseudogene,transcribed_processed_pseudogene,transcribed_unitary_pseudogene,transcribed_unprocessed_pseudogene,translated_unprocessed_pseudogene,unitary_pseudogene,unprocessed_pseudogene,translated_processed_pseudogene", 
	'ncRNA': "lncRNA,Mt_rRNA,Mt_tRNA,miRNA,misc_RNA,rRNA,retained_intron,ribozyme,sRNA,scRNA,scaRNA,snRNA,snoRNA,vaultRNA", 
	]

params.final_database_protein = "final_proteinDB.fa"

protein_decoy_config = "${baseDir}/configs/protein_decoy.yaml"
params.decoy_prefix = "decoy_"

params.cosmic_config = "${baseDir}/configs/cosmic_config.yaml"
cosmic_config = file(params.cosmic_config)

params.gencode_url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19"
//params.gnomad_file_url =  "gs://gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz"
params.gnomad_file_url =  "gs://gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.2*.vcf.bgz" //used only for testing the pipeline

params.cosmic_user_name = ""
params.cosmic_password = ""

ZCAT = (System.properties['os.name'] == 'Mac OS X' ? 'gzcat' : 'zcat')

/* Pipeline START */

/** 
 * Download data from ensembl for the particular species. 
 */ 
process ensembl_protein_fasta_download{
    
    //container 'quay.io/bigbio/pypgatk:0.0.1'
    
    input: 
    file ensembl_downloader_config

    output:
	file "database_ensembl/*.gz" into ensembl_fasta_gz_databases

	script:
	"""
	python ${container_path}pypgatk_cli.py ensembl-downloader --config_file ${ensembl_downloader_config} --taxonomy ${params.taxonomy} -sv -sc
	"""
}

/**
 * Decompress all the data downloaded from ENSEMBL 
 */ 
process gunzip_ensembl_files{
    
    publishDir "result", mode: 'copy', overwrite: true

    input: 
    file(fasta_file) from ensembl_fasta_gz_databases

    output: 
    file '*.pep.all.fa' into ensembl_protein_database
    file '*cdna.all.fa' into ensembl_cdna_database, ensembl_cdna_database_sub
    file '*ncrna.fa' into ensembl_ncrna_database, ensembl_ncrna_database_sub
	file '*.gtf' into gtf
	
    script: 
    """
    gunzip -d --force ${fasta_file}
    """ 
}
/* Concatenate cDNA and ncRNA databases */
process merge_cdnas{
  
  input:
  file a from ensembl_cdna_database_sub
  file b from ensembl_ncrna_database_sub
  
  output: 
  file 'total_cdnas.fa' into total_cdnas

  script: 
  """
  cat "${a}" "${b}" >> total_cdnas.fa
  """ 
}

/**
 * Creates the ncRNA protein database
 */
process add_ncrna{

  //container 'quay.io/bigbio/pypgatk:0.0.1'
  publishDir "result", mode: 'copy', overwrite: true
  
  when:
  params.ncrna
  
  input:
  file x from total_cdnas
  file ensembl_config

  output:
  file 'ncRNAs_proteinDB.fa' into optional_ncrna

  script:
  """
  python ${container_path}pypgatk_cli.py dnaseq-to-proteindb --config_file "${ensembl_config}" --input_fasta ${x} --output_proteindb ncRNAs_proteinDB.fa --include_biotypes "${biotypes['ncRNA']}" --skip_including_all_cds
  """
}

merged_databases = ensembl_protein_database.mix(optional_ncrna)

/**
 * Creates the pseudogenes protein database
 */
process add_pseudogenes {

  //container 'quay.io/bigbio/pypgatk:0.0.1'
  publishDir "result", mode: 'copy', overwrite: true
  
  when:
  params.pseudogenes
  
  input:
  file x from total_cdnas
  file ensembl_config

  output:
  file 'pseudogenes_proteinDB.fa' into optional_pseudogenes

  script:
  """
  python ${container_path}pypgatk_cli.py dnaseq-to-proteindb --config_file "${ensembl_config}" --input_fasta "${x}" --output_proteindb pseudogenes_proteinDB.fa --include_biotypes "${biotypes['pseudogene']}" --skip_including_all_cds
  """
}

merged_databases = merged_databases.mix(optional_pseudogenes)

/**
 * Creates the altORFs protein database
 */
process add_altorfs {

  //container 'quay.io/bigbio/pypgatk:0.0.1'
  publishDir "result", mode: 'copy', overwrite: true
  
  when:
  params.altorfs
  
  input:
  file x from ensembl_cdna_database
  file ensembl_config

  output:
  file('altorfs_proteinDB.fa') into optional_altorfs

  script:
  """
  python ${container_path}pypgatk_cli.py dnaseq-to-proteindb --config_file "${ensembl_config}" --input_fasta "${x}" --output_proteindb altorfs_proteinDB.fa --include_biotypes "${biotypes['protein_coding']}'" --skip_including_all_cds
  """
}

merged_databases = merged_databases.mix(optional_altorfs)

/* Mutations to proteinDB */

/**
 * Download COSMIC Mutations
 */
process cosmic_download {
	
	//container 'quay.io/bigbio/pypgatk:0.0.1'
	
	when:
  	params.cosmic
  		
	input:
	file cosmic_config
	
	output:
	file "database_cosmic/*.gz" into cosmic_files
	
	script:
	"""
	python ${container_path}pypgatk_cli.py cosmic-downloader --config_file "${cosmic_config}" --username ${params.cosmic_user_name} --password ${params.cosmic_password}
	"""
}

/** 
 * Decompress the data downloaded from COSMIC
 */
process gunzip_cosmic_files{

	when:
  	params.cosmic
  	
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

/**
 * Generate proteindb from cosmic mutations
*/
process cosmic_proteindb{
	
	//container 'quay.io/bigbio/pypgatk:0.0.1'
	publishDir "result", mode: 'copy', overwrite: true 
	
	when:
  	params.cosmic
  		
	input:
	file g from cosmic_genes
	file m from cosmic_mutations
	
	output:
	file 'cosmic_proteinDB*.fa' into cosmic_proteindbs
	
	script:
	"""
	python ${container_path}pypgatk_cli.py cosmic-to-proteindb --config_file "${cosmic_config}" --input_mutation ${m} --input_genes ${g} --output_db cosmic_proteinDB.fa
	"""
}

merged_databases = merged_databases.mix(cosmic_proteindbs)

/** 
 * Download VCF files from ensembl for the particular species. 
 */
process ensembl_vcf_download(){
    
    //container 'quay.io/bigbio/pypgatk:0.0.1'
    
    when:
    params.ensembl
    
    input: 
    file ensembl_downloader_config

    output:
	file "database_ensembl/*.vcf.gz" into ensembl_vcf_gz_files
	
	script:
	"""
	python ${container_path}pypgatk_cli.py ensembl-downloader  --config_file ${ensembl_downloader_config} --taxonomy ${params.taxonomy} -sg -sp -sc -sd -sn
	"""
}

/**
 * Decompress vcf files downloaded from ENSEMBL 
 */ 
process gunzip_vcf_ensembl_files{
	
	when:
    params.ensembl
    
    input: 
    file(vcf_file) from ensembl_vcf_gz_files

    output: 
    file "*.vcf" into ensembl_vcf_files
    
    script: 
    """
    gunzip -d --force ${vcf_file}
    """
}

/**
 * Generate protein database(s) from ENSEMBL vcf file(s) 
 */ 

ensembl_vcf_files
	.combine(total_cdnas)
	.combine(gtf)
	.combine(Channel.from(ensembl_config))
	.set { ensembl_files }

process ensembl_vcf_proteinDB {

	//container 'quay.io/bigbio/pypgatk:0.0.1'
	publishDir "result", mode: 'copy', overwrite: true
	
	when:
	params.ensembl
  
  	input:
  	set file(v), file(f), file(g), val(e) from ensembl_files
  	
  	output:
  	file "${v}_proteinDB.fa" into proteinDB_vcf
  	
    script:
  	"""
  	python ${container_path}pypgatk_cli.py vcf-to-proteindb --config_file ${e} --af_field "${params.af_field}" --include_biotypes "${biotypes['protein_coding']}" --input_fasta ${f} --gene_annotations_gtf ${g} --vep_annotated_vcf ${v} --output_proteindb "${v}_proteinDB.fa"
  	"""
}

merged_databases = merged_databases.mix(proteinDB_vcf)

/****** gnomAD variatns *****/

/**
 * Download gencode files (fasta and gtf)
 */
process gencode_download{
	
	when:
	params.gnomad
	
	input:
	val g from params.gencode_url 
	
	output:
	file("gencode.v19.pc_transcripts.fa") into gencode_fasta
	file("gencode.v19.annotation.gtf") into gencode_gtf
	
	script:
	"""
	wget ${g}/gencode.v19.pc_transcripts.fa.gz 
	wget ${g}/gencode.v19.annotation.gtf.gz
	gunzip *.gz
	"""
}

/**
 * Download gnomAD variants (VCF) - requires gsutil
 */
process gnomad_download{
	
	when:
	params.gnomad
	
	input:
	val g from params.gnomad_file_url
	
	output:
	file '*.vcf.bgz' into gnomad_vcf_bgz
	script:
	"""
	gsutil cp ${g} .
	"""
}

/**
 * Extract gnomAD VCF
 */
process extract_gnomad_vcf{
	
	when:
	params.gnomad
	
	input:
	file g from gnomad_vcf_bgz
	
	output:
	file "*.vcf" into gnomad_vcf_files
	
	script:
	"""
	$ZCAT ${g} > ${g}.vcf
	"""
}

/**
 * Generate gmomAD proteinDB
 */

gnomad_vcf_files
	.combine(gencode_fasta)
	.combine(gencode_gtf)
	.combine(Channel.from(ensembl_config))
	.set { gnomad_files }

process gnomad_proteindb{
	
	//container 'quay.io/bigbio/pypgatk:0.0.1'
	publishDir "result", mode: 'copy', overwrite: true
	
	when:
	params.gnomad
	
	input:
	set file(v), file(f), file(g), val(e) from gnomad_files
	
	output:
	file "${v}_proteinDB.fa" into gnomad_vcf_proteindb
	
	script:
	"""
	python ${container_path}pypgatk_cli.py vcf-to-proteindb --config_file ${e} --vep_annotated_vcf ${v} --input_fasta ${f} --gene_annotations_gtf ${g} --output_proteindb "${v}_proteinDB.fa" --af_field controls_AF --transcript_index 6 --biotype_str transcript_type --annotation_field_name vep
	"""
}

merged_databases = merged_databases.mix(gnomad_vcf_proteindb)

/****** cBioPortal mutations *****/
/**
 * Download GRCh37 CDS file from ENSEMBL release 75 
 */
process cds_GRCh37_download{
	
	when:
	params.cbioportal
	
	output:
	file("Homo_sapiens.GRCh37.75.cds.all.fa") into GRCh37_cds
	
	script:
	"""
	wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.75.cds.all.fa.gz 
	gunzip *.gz
	"""
}
/**
 * Download all cBioPortal studies using git-lfs
*/
 process download_all_cbioportal {
 	
 	when:
	params.cbioportal
	
 	output:
 	file('cbioportal_allstudies_data_mutations_mskcc.txt') into cbio_mutations
 	file('cbioportal_allstudies_data_clinical_sample.txt') into cbio_samples
 	
 	script:
 	"""
 	git clone https://github.com/cBioPortal/datahub.git
 	cd datahub
 	git lfs install --local --skip-smudge
 	git lfs pull -I public --include "data*clinical*sample.txt"
 	git lfs pull -I public --include "data_mutations_mskcc.txt"
 	cd ..
 	cat datahub/public/*/data_mutations_mskcc.txt > cbioportal_allstudies_data_mutations_mskcc.txt
 	cat datahub/public/*/*data*clinical*sample.txt > cbioportal_allstudies_data_clinical_sample.txt
 	"""
 }
 
/**
 * Generate proteinDB from cBioPortal mutations
 */
 process cbioportal_proteindb{
	
	publishDir "result", mode: 'copy', overwrite: true 
	
	when:
	params.cbioportal
	
	input:
	file g from GRCh37_cds
	file m from cbio_mutations
	file s from cbio_samples
	
	output:
	file 'cbioPortal_proteinDB*.fa' into cBioportal_proteindb
	
	script:
	"""
	python ${container_path}pypgatk_cli.py cbioportal-to-proteindb --config_file "${cbioportal_config}" --input_mutation ${m} --input_cds ${g} --clinical_sample_file ${s} --output_db cbioPortal_proteinDB.fa
	"""
}

merged_databases = merged_databases.mix(cBioportal_proteindb)

/**
 * Create the decoy database using DecoyPYrat
 * Decoy sequences will have "_DECOY" prefix tag to the protein accession.
 */
 
/**
 * Create the decoy database using DecoyPYrat
 * Decoy sequences will have "_DECOY" prefix tag to the protein accession.
 */

process merge_proteindbs {
	
	//container 'quay.io/bigbio/pypgatk:0.0.1'
    publishDir "result", mode: 'copy', overwrite: true 

	input:
	file("proteindb*") from merged_databases.collect()

	output:
	file "${params.final_database_protein}" into protiendbs
	
	script:
	"""
	cat proteindb* > ${params.final_database_protein}
	"""
}

process decoy {
	
	//container 'quay.io/bigbio/pypgatk:0.0.1'
    publishDir "result", mode: 'copy', overwrite: true 

	input:
	file f from protiendbs

	output:
	file "${params.decoy_prefix}${params.final_database_protein}" into fasta_decoy_db
	
	script:
	"""
	python ${container_path}pypgatk_cli.py generate-decoy --config_file "${protein_decoy_config}" --input $f --decoy_prefix "${params.decoy_prefix}" --output "${params.decoy_prefix}${params.final_database_protein}" 
	"""
}
