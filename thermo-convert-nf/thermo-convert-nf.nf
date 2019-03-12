#!/usr/bin/env nextflow

params.px_accession = ""

process findProjectPath {

    output:
    file 'project_path.txt' into public_dataset_path

    """
    python ../../../scripts/Project.py $params.px_accession
    """
}

process downloadFiles {
    container 'quay.io/biocontainers/gnu-wget:1.18--3'

    input:
    file directory from public_dataset_path.flatten()

    output:
    file '*.raw' into rawFiles

    script:
    """
    wget -v -r -nd -A "Campy_X2_Fr01.raw" --no-host-directories --cut-dirs=1 `cat $directory`
    """
}

process generateMetadata {
    container 'ypriverol/thermorawfileparser:0.1'

    publishDir 'data/', mode:'copy'

    input:
    file rawFile from rawFiles.flatten()

    output:
    file '*.json' into metaResults
    file '*.mzML' into spectraFiles

    script:
    """
    ThermoRawFileParser -i=${rawFile} -m=0 -f=1 -o=./
    """
}

public_dataset_path.subscribe {
    println it
}

metaResults.subscribe {
    println it
}