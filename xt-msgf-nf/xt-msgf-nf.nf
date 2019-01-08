#!/usr/bin/env nextflow

/**
 * Input parameters
 */
params.raw_dir = "${baseDir}/test"
params.fasta_file = "${baseDir}/test/crap.fasta"
// precursor tolerance can only be specified in ppm
params.prec_tol = 10
// fragment tolerance can only be specified in Th
params.frag_tol = 0.5
// missed cleavages
params.mc = 1
// TODO: specify a way to define PTMs

/**
 * Create a channel for all MGF files
 **/
mgf_files = Channel.fromPath("${params.raw_dir}/*.mgf")
fasta_file = file(params.fasta_file)

/**
 * Use SearchGUI to create a decoy database
 */
process createDecoyDb {
	container 'biocontainers/searchgui:v2.8.6_cv2'

	input:
	file "db.fasta" from fasta_file

	output:
	file "*concatenated_target_decoy*" into fasta_db

	script:
	"""
	java -cp /home/biodocker/bin/SearchGUI-2.8.6/SearchGUI-2.8.6.jar eu.isas.searchgui.cmd.FastaCLI -decoy -in db.fasta
	"""
}

/**
 * Create the parameter file for all SearchGUI runs
 */
process createParamsFile {
	container 'biocontainers/searchgui:v2.8.6_cv2'

	input:
	val prec_tol from params.prec_tol
	val frag_tol from params.frag_tol
        val mc from params.mc
	file fasta_db

	output:
	file "parameters.par" into params_file	

	script:
	"""
	java -cp /home/biodocker/bin/SearchGUI-2.8.6/SearchGUI-2.8.6.jar eu.isas.searchgui.cmd.IdentificationParametersCLI \
	-db db_concatenated_target_decoy.fasta \
	-prec_tol ${prec_tol} \
	-prec_ppm 1 \
	-frag_tol ${frag_tol} \
	-enzyme "Trypsin" \
	-fixed_mods "Carbamidomethylation of C" \
	-variable_mods "Oxidation of M,Acetylation of protein N-term" \
	-min_charge 2 \
	-max_charge 4 \
	-mc ${mc} \
	-out 'parameters.par'
	"""
}

/**
 * Search the mgf files
 */

process search {
	container 'biocontainers/searchgui:v2.8.6_cv2'
	publishDir 'result'

	input:
	file params_file
	file fasta_db
	file mgf_file from mgf_files

	output:
	file "*.zip" into search_result

	script:
	"""
	java -cp /home/biodocker/bin/SearchGUI-2.8.6/SearchGUI-2.8.6.jar eu.isas.searchgui.cmd.SearchCLI \
	-spectrum_files ${mgf_file} \
	-output_folder . \
	-id_params ${params_file} \
	-xtandem_refine 0 \
	-xtandem 1 -msgf 1 -comet 0 -ms_amanda 0 -myrimatch 0 -andromeda 0 -omssa 0 -tide 0
	"""
}
