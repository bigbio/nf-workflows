#!/usr/bin/env nextflow

/*
 *  This pipeline is designed to analyze a Label-free Quant experiment using the OpenMS framework. It has been use to analyse the
 *  benchmark dataset iPRG2015.
 *
 *   Params:
 *
 *   --mzml_folder The folder containing the mzML folder
 *   --exp_design The Experimental design file
 *   --fasta      Fasta database
 *   --id_config  Identification step configuration file (see example configs/msgf.ini)
 *   --index_config Peptide indexer config file.
 *   --fdr_config The FDR config provides the parameters to filter peptides/proteins
 *   --idfilter_config The Identification filter config
 *
 */

/*
 *
 * files = Channel.fromPath("${params.raw_folder}/*.mzML").map { file -> [ name:file.baseName, file:file]}
 *
 * group_info = Channel
 *	.fromPath(params.exp_design)
 *	.splitCsv(header: true)
 *	.map { row -> [ name:row[1], group:row[] ]}
 *
 * comb = files.phase(groupInfo) {it -> it.name}
 *	     .map { it -> tuple( it.group[1], it.file[0])}
 *	     .groupTuple(by:0)
 *
 */

(mz_files, mz_result_files) = Channel.fromPath("${params.raw_folder}/*.mzML").into(2)
fasta_file   = file(params.fasta)

/*
 * Config files for each processing step
 */
id_config    = file(params.id_config)
index_config = file(params.index_config)
fdr_config   = file(params.fdr_config)
idfilter_config = file(params.idfilter_config)

/**
 * Identification step using MSGF+
 */
process peptideIdentification {

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq'
   publishDir "results", mode: 'copy', overwrite: true

   memory { 10.GB * task.attempt }

   input:
   file mz_ml from mz_files
   file "database.fasta" from fasta_file
   file id_config

   output:
   file "*.idXML" into id_xmls

   script:
   """
   MSGFPlusAdapter -ini "${id_config}" -database database.fasta -in ${mz_ml} -out ${mz_ml}.idXML
   """
}

/**
 * PeptideIndexer is an step that is used to map the identified peptides to protein ids.
 */
process peptideIndexer {

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq'
   publishDir "results", mode: 'copy', overwrite: true

   input:
   file id_xml from id_xmls
   file "database.fasta" from fasta_file
   file index_config

   output:
   file "*.idXML" into index_xmls

   script:
   """
   PeptideIndexer -ini "${index_config}" -fasta database.fasta -in ${id_xml} -out ${id_xml.baseName}-index.idXML
   """
}


/**
 * FalseDiscoveryRate This step compute the FDR for the identified peptides.
 */
process peptideFDRCompute {

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq'
   publishDir "results", mode: 'copy', overwrite: true

   input:
   file index_xml from index_xmls
   file fdr_config

   output:
   file "*.idXML" into fdr_xmls

   script:
   """
   FalseDiscoveryRate -ini "${fdr_config}" -in ${index_xml} -out ${index_xml.baseName}-fdr.idXML
   """
}


/**
 * IDFilter This step filter the peptides using the FDR computation
 */
process peptideFDRFilter {

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq'
   publishDir "results", mode: 'copy', overwrite: true

   input:
   file fdr_xml from fdr_xmls
   file idfilter_config

   output:
   file "*.idXML" into peptide_xmls

   script:
   """
   IDFilter -ini "${idfilter_config}" -in ${fdr_xml} -out ${fdr_xml.baseName}-filter.idXML
   """
}

(mzMLs, mzML_print) = mz_result_files.map { file -> tuple(file.name, file)}.into(2)
(peptide_collection, peptide_print) = peptide_xmls.map { file -> tuple(file.baseName.replaceFirst("-index-fdr-filter",""), file)}.into(2)
(combined_results, print_combined) = mzMLs.combine(peptide_collection, by: 0).into(2)


mzML_print.subscribe{ println "value: $it"}
peptide_print.subscribe{ println "value: $it"}
print_combined.subscribe{ println "value: $it"}

process idQualityControl{
   container 'mwalzer/openms-batteries-included:V2.3.0_pepxmlpatch'
   publishDir "results", mode: 'copy', overwrite: true

   input:
   set val(mzML), file(mzML_file), file(idx_file) from combined_results

   output:
   file "*.qcML" into qc_MLs

   script:
   """
   QCCalculator -in ${mzML_file} -id ${idx_file} -out ${mzML}.qcML
   """
}