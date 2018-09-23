#!/usr/bin/env nextflow
 
process downloadFiles {
    container 'quay.io/biocontainers/gnu-wget:1.18--3'
    
    output:
    file '*.raw' into rawFiles
   
    script: 
    """
    wget -v -r -nd -A "*.raw" --no-host-directories --cut-dirs=1 ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/09/PXD010376/
    """ 
}


process generateMetadata {
    container 'quay.io/biocontainers/thermorawfileparser:0.0.2018.09.07--0' 

    input:
    file '*.raw' from rawFiles
    
    output: 
    file '*.json' into metaResults
    
    script:
    """
    bash -c ThermoRawFileParser.sh -m 1 -i ${rawFile} -o $params.metaFolder
    """

}

metaResults.subscribe { results ->
    results.copyTo('./results.txt')
    println "Final results at: results.txt"
    
}

