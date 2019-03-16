#!/usr/bin/env nextflow

/*
 *  This pipeline benchmark two different tools for mzML data conversion ThermoRawParser and ms-convert. A prerequisite of the pipeline
 *  is a folder where the data of each tool is deposited with the name defined in the Experimental Design file.
 *
 *   - msconvert
 *   - thermorawparser
 *
 *   Note: Please be aware that files should keep the same name in both experiments for the statistical analysis.
 *
 *   Params:
 *
 *   --raw_folder The folder containing the raw data
 *   --exp_design The Experimental design file
 *   --fasta      Fasta database
 *   --id_config  Identification step configuration file
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

mz_files = Channel.fromPath("${params.raw_folder}/*.mzML")
fasta_file = file(params.fasta)
id_config = file(params.id_config)

process peptideIdentification {

   container 'mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq'
   publishDir "results", mode: 'copy', overwrite: true

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