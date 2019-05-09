#!/usr/bin/env nextflow

/**
 * Input parameters
 */
params.raw_dir = "${baseDir}/test"
params.fasta_file = "${baseDir}/test/sp_human_19-01.fasta"

/**
 * Search engine + clustering parameters
 */
// precursor tolerance can only be specified in ppm
params.prec_tol = 10
// fragment tolerance can only be specified in Th
params.frag_tol = 0.5
// missed cleavages
params.mc = 1

/**
 * Id transferer parameters
 */
params.min_ident = 2
params.min_ratio = 0.7

// X!Tandem template files should not be changed unless for very good reason.
params.xtandem_template = "${baseDir}/config/input.xml"
params.xtandem_taxonomy = "${baseDir}/config/taxonomy.xml"
params.msgf_mods = "${baseDir}/config/Mods.txt"

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
 * Create the MSGF+ database index adding decoys that are
 * automatically added by MSGF+
 */
process createMsgfDbIndex {
	container 'biocontainers/msgfp:v9949_cv3'
	// MSGF+ will raise an exception since the MGF file is empty
	validExitStatus 0,1
	
	input:
	file "user.fasta" from fasta_file

	output:
	file "user.revCat*" into msgf_fasta_index

	script:
	"""
	touch /tmp/test.mgf
	java -jar /home/biodocker/bin/MSGFPlus_9949/MSGFPlus.jar -s /tmp/test.mgf -d user.fasta -tda 1
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

	input:
	file "user.fasta" from fasta_file
	file msgf_fasta_index
	file mgf_file_msgf from mgf_files_msgf
	file msgf_mods

	output:
	file "*.mzid" into msgf_result
	
	script:
	"""
	java -jar /home/biodocker/bin/MSGFPlus_9949/MSGFPlus.jar \
	-d user.fasta -s ${mgf_file_msgf} -t ${params.prec_tol}ppm -ti 0,1 -thread ${threads} \
	-tda 1 -inst 3 -e 1 -ntt ${params.mc} -mod ${msgf_mods} -minCharge 2 -maxCharge 4 \
	-addFeatures 1
	"""
}

/**
 * Create channels containing the MGF files and the search result files
 */
// Change each result channel into a set with the base MGF name as the key
(mgf_keys1, mgf_keys2) = Channel.fromPath("${params.raw_dir}/*.mgf").map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mgf/)[0]
		[index, file]
	}.into(2)

xtandem_combined = xtandem_result.map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mgf\.xml\.mzid/)[0]
		[index, file]
	}.combine(mgf_keys1, by: 0)

msgf_combined = msgf_result.map { file ->
		(whole_name, index) = (file =~ /.*\/(.*)\.mzid/)[0]
		[index, file]
	}.combine(mgf_keys2, by: 0)

/**
 * Create annotated MGF files
 */

process annotateMsgf {
	container 'jgriss/spectra-cluster-py:latest'

	input:
	set val(mgf_name), file(msgf_file), file(mgf_file) from msgf_combined

	output:
	file "${mgf_name}.msgf.mgf" into msgf_annotated

	script:
	"""
	mgf_search_result_annotator --input ${mgf_file} --search ${msgf_file} --output ${mgf_name}.msgf.mgf --format "MSGF_ident" --fdr 0.01
	"""
}

process annotateTandem {
	container 'jgriss/spectra-cluster-py:latest'

	input:
	set val(mgf_name), file(xtandem_file), file(mgf_file) from xtandem_combined

	output:
	file "${mgf_name}.xt.mgf" into xtandem_annotated

	script:
	"""
	mgf_search_result_annotator --input ${mgf_file} --search ${xtandem_file} --output "${mgf_name}.xt.mgf" --format xtandem --fdr 0.01
	"""
}

/**
 * Perform the clustering on all files coming from one search engine
 */
annotated_mgf = Channel.create()
annotated_mgf << msgf_annotated.toList()
annotated_mgf << xtandem_annotated.toList()
annotated_mgf.close()

process runClustering {
	container 'biocontainers/spectra-cluster-cli:vv1.1.2_cv2'
	publishDir "result"

	input:
	file ("*") from annotated_mgf

	output:
	file "*.clustering" into clustering_result

	script:
	"""
	if [ `ls -1 *.mgf | grep -c ".xt.mgf"` -gt 0 ]; then
		ENGINE="xtandem"
	else
		ENGINE="msgf"
	fi

	spectra-cluster-cli -major_peak_jobs ${threads} \
	-threshold_start 1 -threshold_end 0.99 -rounds 5 \
	-precursor_tolerance ${params.prec_tol} -precursor_tolerance_unit ppm \
	-fragment_tolerance ${params.frag_tol} -filter mz_150 -output_path \${ENGINE}.clustering \
	*.mgf
	"""
}

/**
 * Transfer the identifications based on the clustering results
 */
process transferIds {
	container 'jgriss/spectra-cluster-py:latest'
	publishDir "result"

	input:
	file clustering_file from clustering_result
	file fasta_file

	output:
	file '*.spec_counts.tsv' into spectral_counts

	script:
	"""
	id_transferer_cli --input "${clustering_file}" --output psm_quant.tsv \
		--min_identified ${params.min_ident} --min_ratio ${params.min_ratio} \
		--return_all_identified --only_unidentified
	protein_annotator --input psm_quant.tsv --output ${clustering_file}.spec_counts.tsv \
		--fasta "${fasta_file}" --ignore_il
	"""
}
