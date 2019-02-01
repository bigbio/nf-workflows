#!/usr/bin/env nextflow

/**
 * Input parameters
 */
params.mgf_dir = "${baseDir}/test"
params.fasta_file = "${baseDir}/test/uniprot-human-reviewed.fasta"

params.tide_config = "${baseDir}/config/tide_config.txt"

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
(mgf_files, mgf_files_2) = Channel.fromPath("${params.mgf_dir}/*.mgf").into(2)
fasta_file = file(params.fasta_file)

tide_config_file = file(params.tide_config)
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
 * Create the MSGF+ database index
 */
process createMsgfDbIndex {
	container 'quay.io/biocontainers/msgf_plus:2017.07.21--3'
	memory '4 GB'
	
	input:
	file "user.fasta" from fasta_decoy_db

	output:
	file "user.*" into msgf_fasta_index

	script:
	"""
	msgf_plus edu.ucsd.msjava.msdbsearch.BuildSA -d user.fasta -tda 0
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

	memory { 4.GB * task.attempt }
    	errorStrategy 'retry'

	input:
	file "user.fasta" from fasta_decoy_db
	file msgf_fasta_index
	file mgf_file_msgf from mgf_files
	file msgf_mods

	output:
	file "*.mzid" into msgf_result, all_msgf_result

	script:
	"""
	msgf_plus -Xmx${4 * task.attempt}g -d user.fasta -s ${mgf_file_msgf} -t ${params.prec_tol}ppm -ti 0,1 -thread ${threads} \
	-tda 0 -inst 3 -e 1 -ntt ${params.mc} -mod ${msgf_mods} -minCharge ${params.min_charge} -maxCharge ${params.max_charge}
	"""
}

process createTideIndex {
	container 'biocontainers/crux:v3.2_cv3'
	memory { 4.GB }

	input:
	file fasta_decoy_db
	file tide_config_file

	output:
	file "fasta-index/*" into fasta_tide_index

	script:
	"""
	crux tide-index --parameter-file "${tide_config_file}" "${fasta_decoy_db}" fasta-index
	"""
}

process searchTide {
	container 'biocontainers/crux:v3.2_cv3'
	memory { 4.GB }
	publishDir "result"

	input:
	file fasta_tide_index
	file mgf_file from mgf_files_2

	output:
	file "*.txt" into tide_result, all_tide_result

	script:
	"""
	mkdir fasta-index
	mv ${fasta_tide_index} fasta-index/
	crux tide-search --parameter-file "${tide_config_file}" "${mgf_file}" fasta-index

	# move and rename the result file
	mv crux-output/tide-search.txt ${mgf_file}.tide.txt
	"""
}

/**
 * Combine the search results for every MGF files
 */
// Change each result channel into a set with the base MGF name as the key
tide_key = tide_result.map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mgf\.tide.txt/)[0]
		[index, file]
	}
msgf_key = msgf_result.map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mzid/)[0]
		[index, file]
	}

// merge the results based on the MGF filename as key
combined_results = tide_key.combine(msgf_key, by: 0)

process mergeSearchResults {
	container 'ypriverol/pia:1.3.10'
	publishDir "result"

	memory { 16.GB * task.attempt }
    errorStrategy 'retry'

	input:
	set val(mgf_name), file(tide_mzid), file(msgf_mzid) from combined_results

	output:
	file "${mgf_name}.xml" into pia_compilation

	script:
	"""
	pia -Xmx${10 * task.attempt}g compiler -infile ${tide_mzid} -infile ${msgf_mzid} -name "${mgf_name}" -outfile ${mgf_name}.xml
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

all_search_files = all_msgf_result.concat(all_tide_result)

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
