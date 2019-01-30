#!/usr/bin/env nextflow

/**
 * Input parameters
 */
params.mgf_dir = "${baseDir}/test"
params.fasta_file = "${baseDir}/test/uniprot-human-reviewed.fasta"
// SearchGUI parameter file
params.search_params = "${baseDir}/test/search_gui_params.par"

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

search_gui_config = file(params.search_params)
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

process searchMSGF {
	container 'biocontainers/searchgui:v2.8.6_cv2'
	publishDir "result"

	memory '4 GB'

	input:
	file search_gui_config
	file fasta_decoy_db
	file mgf_file from mgf_files_msgf

	output:
	file "${mgf_file}.mzid" into msgf_result, all_msgf_result

	script:
	"""
	# fix the fasta file path
	sed -i 's|"/[^"]*\\.fasta|"${fasta_decoy_db}|' ${search_gui_config}

	# run the search
	java -Xmx3500M -cp /home/biodocker/bin/SearchGUI-2.8.6/SearchGUI-2.8.6.jar eu.isas.searchgui.cmd.SearchCLI \
	-spectrum_files ${mgf_file} \
	-output_folder . \
	-id_params ${search_gui_config} \
	-xtandem 0 -msgf 1 -comet 0 -ms_amanda 0 -myrimatch 0 -andromeda 0 -omssa 0 -tide 0

	# extract the mzid file
	unzip -o searchgui_out.zip

	# simply fix the filename to fit our pattern
	rename 's/msgf/mgf/' *.mzid
	"""
}

process searchComet {
	container 'biocontainers/searchgui:v2.8.6_cv2'
	publishDir "result"

	memory '4 GB'

	input:
	file search_gui_config
	file fasta_decoy_db
	file mgf_file from mgf_files

	output:
	file "${mgf_file}.pep.xml" into comet_xml_result, all_comet_result
	
	script:
	"""
	# fix the fasta file path
	sed -i 's|"/[^"]*\\.fasta|"${fasta_decoy_db}|' ${search_gui_config}
	
	# run the search
	java -Xmx3500M -cp /home/biodocker/bin/SearchGUI-2.8.6/SearchGUI-2.8.6.jar eu.isas.searchgui.cmd.SearchCLI \
	-spectrum_files ${mgf_file} \
	-output_folder . \
	-id_params ${search_gui_config} \
	-xtandem 0 -msgf 0 -comet 1 -ms_amanda 0 -myrimatch 0 -andromeda 0 -omssa 0 -tide 0

	# extract the result file
	unzip -o searchgui_out.zip

	# rename the result file
	rename 's|\\.comet|.mgf|' *.comet.pep.xml
	"""
}

/**
 * Combine the search results for every MGF files
 */
// Change each result channel into a set with the base MGF name as the key
comet_key = comet_xml_result.map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mgf\.pep\.xml/)[0]
		[index, file]
	}
msgf_key = msgf_result.map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mgf\.mzid/)[0]
		[index, file]
	}

// merge the results based on the MGF filename as key
combined_results = comet_key.combine(msgf_key, by: 0)

process mergeSearchResults {
	container 'ypriverol/pia:1.3.10'
	publishDir "result"

	memory { 4.GB * task.attempt }
	errorStrategy 'retry'

	input:
	set val(mgf_name), file(comet_xml), file(msgf_mzid) from combined_results

	output:
	file "${mgf_name}.xml" into pia_compilation

	script:
	"""
	pia -Xmx${4 * task.attempt}g compiler -infile ${comet_xml} -infile ${msgf_mzid} -name "${mgf_name}" -outfile ${mgf_name}.xml
	"""
}

/**
 * Inference and filtering of the results
 */
process filterPiaResults {
	container 'ypriverol/pia:1.3.10'
	publishDir "result", mode: 'copy', overwrite: true

	memory { 4.GB * task.attempt }
    errorStrategy 'retry'

	input:
	file pia_xml from pia_compilation
	file pia_config

	output:
	file "*.mztab" into final_result

	script:
	"""
	pia -Xmx${4 * task.attempt}g inference -infile ${pia_xml} -paramFile ${pia_config} -psmExport ${pia_xml}.mztab mzTab
	"""
}

/**
 * Collect all the files from search engines.
 */

all_search_files = all_msgf_result.concat(all_comet_result)

/**
 * Compile all the results in to one xml file.
 */
 process mergeCompleteSearchResults {
 	container 'ypriverol/pia:1.3.10'
 	publishDir "result"

 	memory { 4.GB * task.attempt }
    	errorStrategy 'retry'

 	input:
 	file(input_files) from all_search_files.collect()

 	output:
 	file "*.xml" into merged_pia_compilation

 	script:
 	"""
 	pia -Xmx${4 * task.attempt}g compiler -infile ${(input_files as List).join(" -infile ")} -name pia-fina-complete -outfile pia-final-merged.xml
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
