#!/usr/bin/env nextflow

/**
 * Input parameters
 */
params.mgf_dir = "${baseDir}/test"
params.fasta_file = "${baseDir}/test/uniprot-human-reviewed.fasta"

// X!Tandem template files should not be changed unless for very good reason.
params.xtandem_template = "${baseDir}/config/input.xml"
params.xtandem_taxonomy = "${baseDir}/config/taxonomy.xml"

// MGF parameters
params.msgf_mods = "${baseDir}/config/Mods.txt"

// PIA Parameters
params.pia_config = "${baseDir}/config/pia_config.xml"

/**
 * START General Parameters Settings
 */

// precursor tolerance can only be specified in ppm
params.prec_tol = 10

// fragment tolerance can only be specified in Th
params.frag_tol = 0.1

// missed cleavages
params.mc = 2

// number of threads per search engine
threads = 1

// minimun charge to be search
params.min_charge = 1

// maximum charge to be search
params.max_charge = 4

/**
 * Create a channel for all MGF files
 **/
(mgf_files, mgf_files_msgf) = Channel.fromPath("${params.mgf_dir}/*.mgf").into(2)
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
	publishDir "result"

	memory '5 GB'

	input:
	file xtandem_settings
	file xtandem_taxonomy
	file fasta_decoy_db
	file mgf_file from mgf_files

	output:
    file "${mgf_file}.xml.mzid" into xtandem_result
	file "${mgf_file}.xml" into xtandem_xml_result
	file "${mgf_file}.xml" into all_xtandem_result

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

	memory '4 GB'
	
	input:
	file "user.fasta" from fasta_decoy_db

	output:
	file "user.*" into msgf_fasta_index

	script:
	"""
	touch ./temp.mgf
	msgf_plus -s ./temp.mgf -d user.fasta -tda 0
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

	memory { 16.GB * task.attempt }
    errorStrategy 'retry'

	input:
	file "user.fasta" from fasta_decoy_db
	file msgf_fasta_index
	file mgf_file_msgf from mgf_files_msgf
	file msgf_mods

	output:
	file "*.mzid" into msgf_result
	file "*.mzid" into all_msgf_result

	script:
	"""
	msgf_plus -Xmx${10 * task.attempt}g -d user.fasta -s ${mgf_file_msgf} -t ${params.prec_tol}ppm -ti 0,1 -thread ${threads} \
	-tda 0 -inst 3 -e 1 -ntt ${params.mc} -mod ${msgf_mods} -minCharge ${params.min_charge} -maxCharge ${params.max_charge}
	"""
}

/**
 * Combine the search results for every MGF files
 */
// Change each result channel into a set with the base MGF name as the key
xtandem_key = xtandem_xml_result.map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mgf\.xml/)[0]
		[index, file]
	}
msgf_key = msgf_result.map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mzid/)[0]
		[index, file]
	}

// merge the results based on the MGF filename as key
combined_results = xtandem_key.combine(msgf_key, by: 0)

process mergeSearchResults {
	container 'ypriverol/pia:1.3.10'
	publishDir "result"

	memory { 16.GB * task.attempt }
    errorStrategy 'retry'

	input:
	set val(mgf_name), file(xtandem_mzid), file(msgf_mzid) from combined_results

	output:
	file "${mgf_name}.xml" into pia_compilation

	script:
	"""
	pia -Xmx${10 * task.attempt}g compiler -infile ${xtandem_mzid} -infile ${msgf_mzid} -name "${mgf_name}" -outfile ${mgf_name}.xml
	"""
}

/**
 * Inference and filtering of the results
 */
process filterPiaResults {
	container 'ypriverol/pia:1.3.10'
	publishDir "result", mode: 'copy', overwrite: true

	memory { 16.GB * task.attempt }
    errorStrategy 'retry'

	input:
	file pia_xml from pia_compilation
	file pia_config

	output:
	file "*.mztab" into final_result

	script:
	"""
	pia -Xmx${10 * task.attempt}g inference -infile ${pia_xml} -paramFile ${pia_config} -psmExport ${pia_xml}.mztab mzTab
	"""
}

/**
 * Collect all the files from search engines.
 */

all_search_files = all_msgf_result.concat(all_xtandem_result)

/**
 * Compile all the results in to one xml file.
 */
 process mergeCompleteSearchResults {
 	container 'ypriverol/pia:1.3.10'
 	publishDir "result"

 	memory { 16.GB * task.attempt }
    errorStrategy 'retry'

 	input:
 	file(input_files) from all_search_files.collect()

 	output:
 	file "*.xml" into merged_pia_compilation

 	script:
 	"""
 	pia -Xmx${10 * task.attempt}g compiler -infile ${(input_files as List).join(" -infile ")} -name pia-fina-complete -outfile pia-final-merged.xml
 	"""
 }

 /**
  * Inference and filtering the final mztab results
  */
 process filterPiaCompleteResults {
 	container 'ypriverol/pia:1.3.10'
 	publishDir "result", mode: 'copy', overwrite: true

 	memory { 16.GB * task.attempt }
     errorStrategy 'retry'

 	input:
 	file pia_xml from merged_pia_compilation
 	file pia_config

 	output:
 	file "*.mztab" into complete_final_result

 	script:
 	"""
 	pia -Xmx${10 * task.attempt}g inference -infile ${pia_xml} -paramFile ${pia_config} -psmExport ${pia_xml}.mztab mzTab
 	"""
 }