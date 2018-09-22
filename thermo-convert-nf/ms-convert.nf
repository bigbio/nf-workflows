#!/usr/bin/env nextflow

params.spectra = 'data/*.RAW'

spectra = file(params.spectra)


process indexPeptides {
    container 'biocontainers/thermoraw'
    
    input:
    file 'small-yeast.fasta' from peptides
    file 'demo.ms2' from spectra

    output:
    file 'crux-output/tide-search.target.txt' into searchResults
    file 'crux-output/tide-search.decoy.txt' into decoyResults

    script:
    """
    crux tide-index small-yeast.fasta yeast-index
    crux tide-search --compute-sp T demo.ms2 yeast-index
    """
}

process postProcess {
    container 'biocontainers/crux'

    input:
    file 'search.target.txt' from searchResults        
    file 'search.decoy.txt' from decoyResults

    output:
    file 'crux-output/percolator.target.psms.txt' into percResults

    script:
    """
    crux percolator search.target.txt
    """
}

percResults.subscribe { results ->
    results.copyTo('./results.txt')
    println "Final results at: results.txt"
    
}
