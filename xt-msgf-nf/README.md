# Tide / MSGF+ Reprocessing pipeline

This pipeline reprocesses a set of MS/MS RAW files and processes them using MSGF+ and X!Tandem through SearchGUI. Search results are combined using PIA.

## How to run the pipeline

> nextflow



/**
 * Combine the search results for every MGF files
 */
// Change each result channel into a set with the base MGF name as the key
comet_key = comet_xml_result.map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mgf\.target\.text/)[0]
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

 	memory { 60.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3

 	input:
 	file(input_files) from all_search_files.collect()

 	output:
 	file "*.xml" into merged_pia_compilation

 	script:
 	"""
 	pia -Xmx${40 * task.attempt}g compiler -infile ${(input_files as List).join(" -infile ")} -name pia-fina-complete -outfile pia-final-merged.xml
 	"""
 }

 /**
  * Inference and filtering the final mztab results
  */
 process filterPiaCompleteResults {
 	container 'ypriverol/pia:1.3.10'
 	publishDir "result", mode: 'copy', overwrite: true

 	memory { 60.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3

 	input:
 	file pia_xml from merged_pia_compilation
 	file pia_config

 	output:
 	file "*.mztab" into complete_final_result

 	script:
 	"""
 	pia -Xmx${40 * task.attempt}g inference -infile ${pia_xml} -paramFile ${pia_config} -psmExport ${pia_xml}.mztab mzTab
 	"""
 }

process peptideShakerConvert {
	container 'quay.io/biocontainers/peptide-shaker:1.16.36--0'
	publishDir "result"

	memory '4 GB'

	input:
	file "${mgf_file}.cpsx" from peptideshaker_results

	output:
	file "${mgf_file}.mzid" into combined_search_results

	script:
	"""
	peptide-shaker eu.isas.peptideshaker.cmd.MzidCLI -Xms512m -Xmx1g \
	-in ${mgf_file}.cpsx -output_file ${mgf_file}.mzid \
	-contact_first_name  "BigBio team" -contact_last_name "Bioinformatics" \
    -contact_email       "ypriverol@gmail.com" -contact_address "Cambridge, UK" \
    -organization_name   "BigBio Team" -organization_email "ypriverol@gmail.com" \
    -organization_address "Cambridge, UK"
	"""
}
