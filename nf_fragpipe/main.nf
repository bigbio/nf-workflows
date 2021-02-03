#!/usr/bin/env nextflow

/* ---------------------------------------------
   A workflow to use MSFragger and IonQuant to
   perform a closed search and LFQ
   --------------------------------------------- */

params.fasta_file = "swissprot_iRT.fasta"
params.msfragger_config = "closed_fragger.params"
params.raw_file_dir = "raw_files"
params.ptm_shepherd = "ptmshepherd-0.3.5.jar"
params.thermo_ext_dir = "/smb/tools/MSFragger-3.0/MSFragger-3.0/ext/thermo"
params.file_pattern = "*.mzML"
params.mode = "narrow_search"
params.create_spec_lib = "True"

irt_file = Channel.fromPath("../GasPhase/fragpipe_14_DIA_mode/20-10-01_DIA_GasPhase/irt.tsv")

// find all raw files
Channel
    .fromPath("${params.raw_file_dir}/${params.file_pattern}")
    .into { raw_files; raw_files_II; raw_files_III }

// convert the parameters to file objects
fasta_in_file = file(params.fasta_file)
msfragger_config_file = file(params.msfragger_config)
ptm_shepherd_jar = file(params.ptm_shepherd)

process createFasta {
    cpus 1
    memory '10 GB'
    container 'prvst/philosopher:3.2.9'

    input:
    file fasta_in_file

    output:
    file("*-decoys-contam-*.fasta") into fasta_file

    script:
    """
    philosopher workspace --init

    # add decoys + contaminants
    philosopher database --custom ${fasta_in_file} --contam
    """
}

// run the search
process runSearch {
    cpus 5
    memory '50 GB'
    stageInMode 'link'
    errorStrategy 'retry'
    maxRetries 3

    input:
    file("*") from raw_files.collect()
    file("db.fasta") from fasta_file
    file msfragger_config_file

    output:
    file("*.pepXML") into search_result_out

    script:
    """
    java -Xmx48G -jar /smb/tools/MSFragger-3.0/MSFragger-3.0/MSFragger-3.0.jar ${msfragger_config_file} ${params.file_pattern}
    """
}

(search_result, search_result_II) = search_result_out.into(2)


process philosopher {
    cpus 10
    memory '10 GB'
    container 'prvst/philosopher:3.2.9'
    publishDir "${params.mode}", mode: 'copy'

    input:
    file("*") from search_result
    file fasta_file

    output:
    file("*.tsv") into final_result
    file("psm.tsv") into psm_result
    file("combined.prot.xml") into combined_prot_result
    file("*.pep.xml") into pep_xml_results

    script:
    """
    # init the workspace
    philosopher workspace --init

    # annotate the database
    philosopher database --annotate ${fasta_file}

    # peptide prophet
    for PEPXML in *.pepXML; do
        philosopher peptideprophet \\
            --database ${fasta_file} \\
            --nonparam \\
            --expectscore \\
            --decoyprobs \\
            --masswidth 1000.0 --clevel -2 \${PEPXML}
    done

    # protein prophet
    philosopher proteinprophet --maxppmdiff 2000000 --output combined *.pep.xml

    # create the combined peptide result
    philosopher filter --sequential --prot 0.01 --tag rev_ \
        --pepxml ./ \
        --protxml combined.prot.xml

    # create the reports
    philosopher report
    """
}


process combinedPeptideResult {
    cpus 10
    memory '10 GB'
    container 'prvst/philosopher:3.2.9'
    publishDir "${params.mode}", mode: 'copy'

    input:
    file fasta_file
    file prot_xml from combined_prot_result
    file ("*") from search_result_II

    output:
    file("*.pep.xml") into combined_pep_xml
    file("protein.fas") into protein_fasta

    script:
    """
    # init the workspace
    philosopher workspace --init

    # annotate the database
    philosopher database --annotate ${fasta_file}

    # run peptide prophet to combine all search results
    philosopher peptideprophet --database ${fasta_file} \
        --combine \
        --decoyprobs \
        --accmass \
        --nonparam \
        --output combined \
        *.pepXML

    # filter on the peptide and protein level
    philosopher filter --razor --pepxml *.pep.xml --protxml ${prot_xml}
    """
}

process ptmShepherd {
    cpus 10
    memory '10 GB'
    stageInMode 'link'
    publishDir "${params.mode}", mode: 'copy'
    when params.mode != "narrow_search"

    input:
    file("psm.tsv") from psm_result
    file("*") from raw_files_II.collect()
    file ptm_shepherd_jar
    
    output:
    file("*.tsv")

    script:
    """
    # use the default config
    echo -e "dataset = DATASET1 psm.tsv .\\n" > ptm_shepherd.txt

    java -Xmx10G -Dbatmass.io.libs.thermo.dir=${params.thermo_ext_dir} -jar ${ptm_shepherd_jar} ptm_shepherd.txt
    """
}

// merge raw files and pep.xml files
combined_result = pep_xml_results
    .flatten()
    .map { fileobj ->
        (whole_path, filename) = (fileobj =~ /.*interact-(.*)\.pep\.xml/)[0]
        [filename, fileobj]
    }
    .join(raw_files_III
        .map {
            (whole_path, filename) = (it =~ /.*\/(.*)\.mzML/)[0]
            [filename, it]
        }
    )

process convertSpeclibFiles {
    cpus 2
    container 'grosenberger/easypqp'
    memory '5 GB'
    stageInMode 'link'
    when params.create_spec_lib == "True"

    input:
    tuple val(name), file(pep_xml_file), file(mzml_file) from combined_result

    output:
    file("*pkl") into pkl_files

    script:
    """
    easypqp convert \
        --pepxml ${pep_xml_file} \
        --spectra ${mzml_file} \
        --exclude-range -1.5,3.5 \
        --psms ${name}.psmpkl \
        --peaks ${name}.peakpkl
    """
}//*/

process createLibrary {
    cpus 1
    memory '10 GB'
    container 'grosenberger/easypqp'
    publishDir "${params.mode}", mode: 'copy'
    when params.create_spec_lib == "True"

    input:
    file("*") from pkl_files.collect()
    file("*") from final_result
    file irt_file
    
    output:
    file("${params.mode}_easypqp_lib.tsv") into spec_lib

    script:
    """
    # build the library
    easypqp library \
        --psmtsv psm.tsv \
        --peptidetsv peptide.tsv \
        --rt_reference ${irt_file} \
        --out ${params.mode}_easypqp_lib.tsv \
        --rt_lowess_fraction 0.01 \
        *.psmpkl *.peakpkl
    """
}//*/