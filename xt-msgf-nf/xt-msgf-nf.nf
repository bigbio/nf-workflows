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
params.pia_config = "${baseDir}/config/pia_config.xml"

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
(mgf_files, mgf_files_msgf) = Channel.fromPath("${params.raw_dir}/*.mgf").into(2)
fasta_file = file(params.fasta_file)

xtandem_template = file(params.xtandem_template)
xtandem_taxonomy = file(params.xtandem_taxonomy)
msgf_mods = file(params.msgf_mods)
pia_config = file(params.pia_config)

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
	container 'biocontainers/tandem:v17-02-01-4_cv4'

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
	container 'quay.io/biocontainers/msgf_plus:2017.07.21--3'
	// MSGF+ will raise an exception since the MGF file is empty
	validExitStatus 0,1
	
	input:
	file "user.fasta" from fasta_decoy_db

	output:
	file "user.*" into msgf_fasta_index

	script:
	"""
	touch /tmp/test.mgf
	msgf_plus -s /tmp/test.mgf -d user.fasta -tda 0
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
	container 'quay.io/biocontainers/msgf_plus:2017.07.21--3'
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
	msgf_plus -d user.fasta -s ${mgf_file_msgf} -t ${params.prec_tol}ppm -ti 0,1 -thread ${threads} \
	-tda 0 -inst 3 -e 1 -ntt ${params.mc} -mod ${msgf_mods} -minCharge 2 -maxCharge 4 
	"""
}

/**
 * Combine the search results for every MGF files
 */
// Change each result channel into a set with the base MGF name as the key
xtandem_key = xtandem_result.map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mgf\.xml\.mzid/)[0]
		[index, file]
	}
msgf_key = msgf_result.map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mzid/)[0]
		[index, file]
	}

// merge the results based on the MGF filename as key
combined_results = xtandem_key.combine(msgf_key, by: 0)

process mergeSearchResults {
	container 'biocontainers/pia:v1.3.8_cv1'

	input:
	set val(mgf_name), file(xtandem_mzid), file(msgf_mzid) from combined_results

	output:
	file "${mgf_name}.xml" into pia_compilation

	script:
	"""
	java -cp /home/biodocker/pia/pia-1.3.8.jar de.mpc.pia.intermediate.compiler.PIACompiler \
	-infile ${xtandem_mzid} -infile ${msgf_mzid} -name "${mgf_name}" -outfile ${mgf_name}.xml
	"""
}

process filterPiaResuls {
	container 'biocontainers/pia:v1.3.8_cv1'
	publishDir "result"

	input:
	file pia_xml from pia_compilation
	file pia_config

	output:
	file "*.mzTab" into final_result

	// TODO: Failed to download unimod.xml - use cached version

	script:
	"""
	java -jar /home/biodocker/pia/pia-1.3.8.jar -infile ${pia_xml} -paramFile ${pia_config} \
	-psmExport ${pia_xml}.mzTab mzTab
	"""
}
