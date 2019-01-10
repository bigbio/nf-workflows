#!/usr/bin/env nextflow

params.spectra_dir = "./nf-test/"

input_dir = params.spectra_dir.replaceFirst(/\/$/, "")
input_dir = "${input_dir}/"

println "Processing MGF files in ${input_dir}..."

mgf_files = Channel.fromPath("${input_dir}*.mgf")
mgf_files_2 = mgf_files

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

process runNewClustering {
    container 'jgriss/spectra-cluster-2:develop'

    input:
    file  "spectra*.mgf" from mgf_files2

    output:
    file '*.clustering' into newClusteringResults, newClusteringResult2

    script:
    """
    spectra-cluster-2 -output.path ./result_95.clustering -rounds 5 -threshold.start 1 -threshold.end 0.95 *.mgf
    spectra-cluster-2 -output.path ./result_97.clustering -rounds 5 -threshold.start 1 -threshold.end 0.97 *.mgf
    spectra-cluster-2 -output.path ./result_99.clustering -rounds 5 -threshold.start 1 -threshold.end 0.99 *.mgf
    """
}

process getStatistics {
    container 'jgriss/spectra-cluster-py:latest'

    input:
    file "spectra-cluster_*.clustering" from clusteringResults1
    file "spectra-cluster-2_*.clustering" from newClusteringResults1

    output:
    file '*.tsv' into clusteringStats, clusteringStatsForR

    script:
    """
    clustering_stats --output stats.tsv --min_size 3 *.clustering
    """
}

process createPlot {
    input:
    file "stats.tsv" from clusteringStatsForR

    output:
    file "*.png" into statPlots

    script:
    '''
    #!/usr/bin/env Rscript
    library(ggplot2)
    stats <- read.csv("stats.tsv", sep = "\\t")
    stats[, "rel_clustered"] <- stats[, "clustered_spectra"] / stats[, "total_spectra"]
    stats[, "rel_incorr"] <- stats[, "incorrect_spectra"] / stats[, "clustered_identified_spectra"]

    stats[, "filename"] <- gsub(".clustering", "", stats[, "filename"])

    stats[, "threshold"] <- gsub(".*_([0-9]*)\$", "\\\\1", stats[, "filename"])
    stats[, "algorithm"] <- gsub("_[^_]*\$", "", stats[, "filename"])

    plot.obj <- ggplot(stats, aes(x = rel_incorr, y = rel_clustered, group = algorithm, color = threshold)) +
        geom_line(color = "black", alpha = "0.5") +
        geom_point(aes(color = threshold)) +
        labs(x = "Rel. incorr. spectra", y = "Rel. clustered spectra")

    png("stats.png", width = 2000, height = 2000, res = 300)
    print(plot.obj)
    dev.off()
    '''
}

clusteringResults2.collect {
    println "Copying clustering results to ${params.spectra_dir}"
    it.collect { it.copyTo(input_dir) }
}

newClusteringResult2.collect {
    println "Copying new clustering results to ${params.spectra_dir}"
    it.collect { it.copyTo(input_dir) }
}

clusteringStats.collect {
    println "Copying statistics results to ${params.spectra_dir}"
    it.copyTo(input_dir)
}

statPlots.collect {
    println "Copying plots to ${params.spectra_dir}"
    it.copyTo(input_dir)
}
