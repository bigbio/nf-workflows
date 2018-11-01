#!/usr/bin/env nextflow

params.spectra_dir = "./nf-test/"

input_dir = params.spectra_dir.replaceFirst(/\/$/, "")
input_dir = "${input_dir}/"

println "Processing MGF files in ${input_dir}..."

mgf_files = Channel.fromPath("${input_dir}*.mgf")

process runClustering {
    container 'biocontainers/spectra-cluster-cli:1.1.2'

    input:
    file "spectra*.mgf" from mgf_files

    output:
    file '*.clustering' into clusteringResults1, clusteringResults2

    script:
    """
    spectra-cluster-cli -output_path ./result_95.clustering -rounds 5 -threshold_start 1 -threshold_end 0.95 *.mgf
    spectra-cluster-cli -output_path ./result_97.clustering -rounds 5 -threshold_start 1 -threshold_end 0.97 *.mgf
    spectra-cluster-cli -output_path ./result_99.clustering -rounds 5 -threshold_start 1 -threshold_end 0.99 *.mgf
    """
}

process getStatistics {
    container 'jgriss/spectra-cluster-py:latest'

    input:
    file "result*.clustering" from clusteringResults1

    output:
    file '*.tsv' into clusteringStats

    script:
    """
    clustering_stats --output stats.tsv --min_size 3 *.clustering
    """
}

clusteringResults2.collect {
    println "Copying clustering results to ${params.spectra_dir}"
    it.collect { it.copyTo(input_dir) }
}

clusteringStats.collect {
    println "Copying statistics results to ${params.spectra_dir}"
    it.copyTo(input_dir)
}

