#!/usr/bin/env nextflow

/**
 * Input parameters
 */
params.raw_dir = "${baseDir}/test"
params.fasta_file = "${baseDir}/test/crap.fasta"
// X!Tandem template files should not be changed unless for very good reason.
params.xtandem_template = "${baseDir}/config/input.xml"
params.xtandem_taxonomy = "${baseDir}/config/taxonomy.xml"
params.msgf_mods = "${baseDir}/config/Mods.txt"
// precursor tolerance can only be specified in ppm
params.prec_tol = 10
// fragment tolerance can only be specified in Th
params.frag_tol = 0.5
// missed cleavages
params.mc = 1
// TODO: specify a way to define PTMs

// number of threads per search engine
threads = 1

/**
 * Create a channel for all MGF files
 **/
mgf_files = Channel.fromPath("${params.raw_dir}/*.mgf")
// copy channel for MSGF+
mgf_files_msgf = Channel.fromPath("${params.raw_dir}/*.mgf")
fasta_file = file(params.fasta_file)

xtandem_template = file(params.xtandem_template)
xtandem_taxonomy = file(params.xtandem_taxonomy)
msgf_mods = file(params.msgf_mods)

/**
 * Create the decoy database for search with X!Tandem
 * 
 * SearchGUI adds reversed sequences by adding a "_REVERSED" tag to the
 * protein accession.
 */
process createDecoyDb {
	container 'biocontainers/searchgui:v2.8.6_cv2'

	input:
	file "db.fasta" from fasta_file

	output:
	file "db_concatenated_target_decoy.fasta" into fasta_decoy_db

	script:
	"""
	java -cp /home/biodocker/bin/SearchGUI-2.8.6/SearchGUI-2.8.6.jar eu.isas.searchgui.cmd.FastaCLI -decoy -in db.fasta
	"""
}

/**
 * Create the configuration file required to launch X!Tandem
 */
process createTandemConfig {
	input:
	file "settings.xml" from xtandem_template

	output:
	file "adapted_settings.xml" into xtandem_settings

	script:
	"""
	sed -e 's|FRAG_TOL|${params.frag_tol}|' \
	    -e 's|PREC_TOL|${params.prec_tol}|' \
	    -e 's|MISSED_CLEAV|${params.mc}|' \
	    -e 's|THREADS|$threads|' \
	    ${xtandem_template} > adapted_settings.xml
	"""
}

/**
 * Search every MGF file using X!Tandem
 */
process searchTandem {
	container 'jgriss/tandem:v17-02-01-4'

	input:
	file xtandem_settings
	file xtandem_taxonomy
	file fasta_decoy_db
	file mgf_file from mgf_files

	output:
        file "${mgf_file}.xml.mzid" into xtandem_result

	script:
	"""
	sed -e 's|ORG_NAME|${mgf_file}|' ${xtandem_settings} > ${mgf_file}.settings.xml && \
	tandem ${mgf_file}.settings.xml && \
	sed -i 's|value="^XXX"|value="_REVERSED"|' ${mgf_file}.xml.mzid
	"""	
}

/**
 * Create the MSGF+ database index
 */
process createMsgfDbIndex {
	container 'biocontainers/msgfp:v9949_cv3'
	// MSGF+ will raise an exception since the MGF file is empty
	validExitStatus 0,1
	
	input:
	file "user.fasta" from fasta_decoy_db

	output:
	file "user.*" into msgf_fasta_index

	script:
	"""
	touch /tmp/test.mgf
	java -jar /home/biodocker/bin/MSGFPlus_9949/MSGFPlus.jar -s /tmp/test.mgf -d user.fasta -tda 0
	"""
}

/**
 * Search every MGF file using MSGF+
 * 
 * Notes on parameters:
 *   * -inst 3 = QExactive, 1 = Orbitrap
 *   * -e 1 = Trypsin
 *   * -ntt 2 = Termini
 */
process searchMsgf {
	container 'biocontainers/msgfp:v9949_cv3'
	publishDir "result"

	input:
	file "user.fasta" from fasta_decoy_db
	file msgf_fasta_index
	file mgf_file_msgf from mgf_files_msgf
	file msgf_mods

	output:
	file "*.mzid" into msgf_result
	
	script:
	"""
	java -jar /home/biodocker/bin/MSGFPlus_9949/MSGFPlus.jar \
	-d user.fasta -s ${mgf_file_msgf} -t ${params.prec_tol}ppm -ti 0,1 -thread ${threads} \
	-tda 0 -inst 3 -e 1 -ntt ${params.mc} -mod ${msgf_mods} -minCharge 2 -maxCharge 4 
	"""
}
