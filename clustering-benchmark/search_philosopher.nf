// ---------------------------------------------
// Parameters
// ---------------------------------------------

// TODO: A) change to input file (TSV with filenames + batch ID)
//       B) use IsoProt structure (better)
params.input_directory = "./testfiles"
params.input_filetype = ".mzML"

params.analysis_name = "DefaultSearch"
params.fasta_file = "testfiles/human_sp_21-02-08.fasta"
params.msgf_param_file = "search_params/msgfplus_param.txt"
params.comet_param_file = "search_params/comet_param.txt"
params.fragger_param_file = "search_params/fragger.params"
params.fragger_jar = "/nfs/tools/MSFragger-3.2/MSFragger-3.2.jar"

// add contaminants to FASTA database
params.add_contam = "True"
params.search_engine = "comet" // alternative 'msgf' or 'fragger'

// containers to use
philosopher_container = "prvst/philosopher:3.2.9"
msgf_container = "quay.io/biocontainers/msgf_plus:2021.01.08--0"
pwiz_container = "proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses"
// TODO: not sure how to better support MSFragger as it cannot be put into
// a container due to license issues
fragger_jar_path = "${params.fragger_jar}"

/**
 *  -------------------- Get the input files + structure ---------------------
 */
print("Discovering files in ${params.input_directory}/*...\n")

// split between files and directories in the input directory
(input_files_msgf, input_files_comet, input_files_fragger, input_files_mgf) = Channel
    .fromPath("${params.input_directory}/*${params.input_filetype}", type: "file")
    .into(4)

// convert files
user_fasta_file = file(params.fasta_file)
msgf_param_file = file(params.msgf_param_file)
comet_param_file = file(params.comet_param_file)
fragger_param_file = file(params.fragger_param_file)

// make sure the file exist
if (!user_fasta_file.exists()) {
    log.error("Failed to find FASTA file")
    exit(1)
} else {
    log.info("Using FASTA file ${user_fasta_file}")
}

if (params.search_engine == "comet" && !comet_param_file.exists()) {
    logger.error("Failed to find comet parameter file '${comet_param_file}'")
    exit(1)
}
if (params.search_engine == "msgf" && !msgf_param_file.exists()) {
    logger.error("Failed to find MSGF parameter file '${msgf_param_file}'")
    exit(1)
}
if (params.search_engine == "fragger" && !fragger_param_file.exists()) {
    logger.error("Failed to find fragger parameter file '${fragger_param_file}'")
    exit(1)
}

process createFasta {
    cpus 1
    memory '10 GB'
    container "${philosopher_container}"

    input:
    file user_fasta_file

    output:
    file("*.fasta") into fasta_file,fasta_file_II,fasta_file_III,fasta_file_IV,fasta_file_comet,fasta_file_msgf

    script:
    def contam = ""

    if (params.add_contam == "True") {
        contam = "--contam"
    }

    """
    philosopher workspace --init

    # add decoys
    philosopher database --custom ${user_fasta_file} ${contam}
    """
}

process createMsgfDb {
    container "${msgf_container}"
    memory "8GB"
    validExitStatus 1 // this is expected to fail

    when:
    params.search_engine == "msgf"

    input:
    file fasta_file

    output:
    file("*") into msgf_fasta_db

    script:
    """
    # create a mock input file
    touch test.mgf

    java -Xmx7900M -jar /usr/local/share/msgf_plus-2021.01.08-0/MSGFPlus.jar -s test.mgf -d ${fasta_file}
    """
}

process runMsgfSearch {
    container "${msgf_container}"
    memory "8GB"
    cpus 4

    when:
    params.search_engine == "msgf"

    input:
    file fasta_file_II
    file("*") from msgf_fasta_db.collect()
    file(input_file) from input_files_msgf
    file msgf_param_file

    output:
    file("*.mzid") into search_result_mzIdent

    script:
    """
    INSTRUMENT="3" # Q-Exactive
    PROTOCOL="5" # Standard
    MODE="3" # HCD

    java -Xmx7900M -jar /usr/local/share/msgf_plus-2021.01.08-0/MSGFPlus.jar -conf ${msgf_param_file} -s ${input_file} -d ${fasta_file_II} -thread ${task.cpus}
    """
}

search_result_pepxml = Channel.empty()

process convertMzId {
    container "${pwiz_container}"
    cpus 1
    memory "8GB"

    input:
    file(mzid_file) from search_result_mzIdent

    output:
    file("*.pepXML") into search_result_msgf_pepXML

    script:
    """
    wine idconvert --pepXML ${mzid_file}
    sed -i 's|max_num_internal_cleavages="-1"|max_num_internal_cleavages="1"|' *.pepXML
    """
}

// add to global search results
search_result_pepxml = search_result_pepxml.mix(search_result_msgf_pepXML)

process runCometSearch {
    container "${philosopher_container}"
    cpus 10
    memory "16GB"

    when:
    params.search_engine == "comet"

    input:
    file fasta_file_comet
    file(input_file) from input_files_comet
    file comet_param_file

    output:
    file("*.pepXML") into search_result_comet_pepXML

    script:
    """
    # adapt parameter file
    cp ${comet_param_file} comet.params

    sed -i \\
        -e 's|PATH_TO_FASTA|${fasta_file_comet}|' \\
        -e 's|NUM_THREADS|${task.cpus}|' \\
        comet.params

    # initialize the workspace
    philosopher workspace --init

    # run the search
    philosopher comet --param comet.params ${input_file}

    # keep file ending consistent with MSGF+
    rename 's|pep\\.xml|pepXML|' *.pep.xml
    """
}

process runMsFraggerSearch {
    cpus 10
    memory "16GB"
    // TODO: this currently only works if java + MSFragger are available locally

    when:
    params.search_engine == "fragger"

    input:
    file "fasta.fas" from fasta_file_msgf
    file(input_file) from input_files_fragger
    fragger_param_file

    output:
    file("*.pepXML") into search_result_fragger_pepXML

    script:
    """
    # adapt the config file
    cp ${fragger_param_file} local_fragger.params
    sed -i 's|THE_THREADS|${task.cpus}|' local_fragger.params

    # run the search
    java -jar ${fragger_jar_path} local_fragger.params ${input_file}
    """
}
//*/

// add to global search results
search_result_pepxml = search_result_pepxml
    .mix(search_result_comet_pepXML)
    .mix(search_result_fragger_pepXML)

process peptideProphet {
    cpus 10
    memory '10 GB'
    container "${philosopher_container}"

    input:
    file(pepxml_file) from search_result_pepxml
    file fasta_file_III

    output:
    file("*.pep.xml") into pep_xml_result,pep_xml_result_II

    script:
    """
    # init the workspace
    philosopher workspace --init

    # annotate the database
    philosopher database --annotate ${fasta_file_III}

    # peptide prophet
    philosopher peptideprophet \\
        --database ${fasta_file_III} \\
        --ppm \\
        --accmass \\
        --nonparam \\
        --expectscore \\
        --decoyprobs \\
        ${pepxml_file}
    """
}
// quay.io/biocontainers/proteowizard:3.0.9992

process proteinProphet {
    cpus 10
    memory '40 GB'
    container "${philosopher_container}"

    input:
    file("*") from pep_xml_result.collect()

    output:
    file("*.prot.xml") into protein_prophet_result

    script:
    """
    philosopher workspace --init
    philosopher proteinprophet --maxppmdiff 2000000 --minprob 0.9 --output combined interact-*.pep.xml
    """
}

process fdrFiltering {
    cpus 2
    memory "8GB"
    container "${philosopher_container}"
    publishDir "${params.analysis_name}_search_result", mode: 'link'

    input:
    file ("*") from pep_xml_result_II.collect()
    file protein_prophet_result
    file fasta_file_IV

    output:
    file("psm.tsv") into psm_tsv
    file("protein.tsv") into protein_tsv
    file("peptide.tsv") into peptide_tsv
    file("psm.tsv") into experiment_psm_tsv

    script:
    """
    # philosopher.exe filter
    philosopher workspace --init

    philosopher database --annotate ${fasta_file_IV}

    philosopher filter \
        --sequential \
        --prot 0.01 \
        --tag rev_ \
        --razor \
        --picked \
        --pepxml ./ \
        --protxml ${protein_prophet_result}

    philosopher report
    """
}