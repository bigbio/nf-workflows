#!/usr/bin/env nextflow

params.rawFolder = 'data/'
params.metaFolder = 'meta/'
 
process downloadFiles {
    container 'quay.io/biocontainers/gnu-wget:1.18--3'
   
    script: 
    """
    wget -v -r -nd -A "*.raw" --no-host-directories --cut-dirs=1 ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/09/PXD010376/ -P $params.rawFolder
    """ 
}

rawFiles = Channel.fromPath( 'data/*.raw' )

process generateMetadata {
    container 'quay.io/biocontainers/thermorawfileparser:0.0.2018.09.07--0' 

    input:
    file 'query.raw' from rawFiles
    
    script:
    """
    thermorawparser -m 1 -i query.raw -o $params.metaFolder
    """

}


