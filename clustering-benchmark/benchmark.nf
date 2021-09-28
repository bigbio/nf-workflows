// ---------------------------------------------
// Parameters
// ---------------------------------------------

params.input_directory = "./testfiles"
params.input_filetype = ".mzML"

params.clustering_parameters = "cluster_parameters.json"
params.psm_result = "DefaultSearch_search_result/psm.tsv"

// ---------------------------------------------
// Basic config parameters
python_container = "python:3.9.5"
pwiz_container = "proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses"
falcon_container = "jgriss/falcon:latest" // TODO: This container is not public
spectra_cluster_container = "quay.io/biocontainers/spectra-cluster-cli:1.1.2--0"
spectra_cluster_py_container = "jgriss/spectra-cluster-py"

// ---------------------------------------------
// change to files
clustering_parameters = file(params.clustering_parameters)

if (!clustering_parameters.exists()) {
    log.error("Failed to find clustering_parameters file '${clustering_parameters}'")
    exit(1)
}

psm_result = file(params.psm_result)

if (!psm_result.exists()) {
  log.error("Failed to find PSM file '${psm_result}'")
  exit(1)
}

// load all files
log.info("Discovering files in ${params.input_directory}...")

(input_files_spectra_cluster, input_files_falcon) = Channel
    .fromPath("${params.input_directory}/*${params.input_filetype}", type: "file")
    .into(2)

process convertToMGF {
    container "${pwiz_container}"
    cpus 1
    memory "8GB"

    input:
    file input_file from input_files_spectra_cluster

    output:
    file("*.mgf") into mgf_input_files

    script:
    """
    # fix the mzML file
    cp ${input_file} input.mzML

    sed -i '/1003145/d' input.mzML

    wine msconvert --mgf input.mzML
    """
}

// run the clustering
spectra_cluster_threshold = [0.99, 0.96, 0.93, 0.9, 0.8, 0.75, 0.5]

process runSpectraCluster {
    container "${spectra_cluster_container}"
    memory '12 GB'
    cpus 8
    containerOptions = "--user 1000"

    input:
    file("*") from mgf_input_files.collect()
    each threshold from spectra_cluster_threshold

    output:
    tuple val("spectra-cluster_${threshold}"), file("*.clustering") into spectra_cluster_raw_result

    script:
    """
    java -Xmx11G -jar /usr/local/share/spectra-cluster-cli-1.1.2-0/spectra-cluster-cli-1.1.2.jar \
        -threshold_start 1 \
        -threshold_end ${threshold} \
        -rounds 5 \
        -precursor_tolerance 20 \
        -precursor_tolerance_unit "ppm" \
        -fragment_tolerance 0.5 \
        -output_path result.clustering \
        *.mgf
    """
}

process convertSpectraCluster {
    container "${spectra_cluster_py_container}"
    memory '1GB'
    cpus 1
    publishDir "clustering_results", mode: 'link'

    input:
    tuple val(settings), file(clustering_file) from spectra_cluster_raw_result

    output:
    tuple val(settings), file("*.tsv") into spectra_cluster_result

    script:
    """
#!/usr/bin/python3

import spectra_cluster.clustering_parser as clustering_parser
import re
import sys


scan_pattern = re.compile(r'.*\\.([0-9]*)\\.[0-9]\$')


input_file = "${clustering_file}"
parser = clustering_parser.ClusteringParser(input_file)

# process all clusters
with open("spectra-cluster.tsv", "w") as writer:
  # write the header
  writer.write("filename\\tspectrum\\tprecursor_mz\\tcluster\\n")
  current_cluster_index = 0

  for cluster in parser:
    spectra = cluster.get_spectra()

    for s in spectra:
      match = scan_pattern.match(s.get_title())

      if not match:
        print("Error: Failed to process title: " + s.get_title())
        sys.exit(1)

      writer.write("{}\\t{}\\t{}\\t{}\\n".format(s.get_filename()[:-4], match.group(1), str(s.precursor_mz), str(current_cluster_index)))

    current_cluster_index += 1

"""
}

falcon_eps = ["0.1", "0.15", "0.2", "0.25", "0.3", "0.35"]

process runFalcon {
    container "${falcon_container}" // TODO: use "proper" falcon container
    memory '16GB'
    cpus 32 // TODO: Missing option in Falcon to set the humber of threads
    
    input:
    file("*") from input_files_falcon.collect()
    each eps from falcon_eps

    output:
    tuple val("falcon_${eps}"), file("*.csv") into falcon_raw_result

    script:
    """
    falcon --precursor_tol 20 ppm --fragment_tol 0.05 --eps ${eps} *.mzML falcon_esult
    """
}

process convertFalcon {
    container "${python_container}"
    memory '1GB'
    cpus 1
    containerOptions = "--user 1000"
    publishDir "clustering_results", mode: "link"

    input:
    tuple val(tool), file(falcon_result) from falcon_raw_result

    output:
    tuple tool, file("falcon.tsv") into falcon_result

    script:
    """
#!/usr/bin/python3

import os
import sys
import csv
import re

identifier_pattern = re.compile(r'mzspec:USI000000:([^:]*):scan:(.*)')

# open the result file
with open("falcon.tsv", "w") as writer:
  # write the header
  writer.write("filename\\tspectrum\\tprecursor_mz\\tcluster\\n")

  with open("${falcon_result}", "r") as input_file:
    reader = csv.DictReader(filter(lambda row: row[0]!='#', input_file))

    for entry in reader:
      # get the id
      match = identifier_pattern.match(entry["identifier"])

      if not match:
        print("Error: Failed to decode identifier " + entry["identifier"])
        sys.exit(1)

      filename = match.group(1)
      scan = match.group(2)

      writer.write("{}\\t{}\\t{}\\t{}\\n".format(filename, scan, entry["precursor_mz"], entry["cluster"]))
    """
}


clustering_result_files = spectra_cluster_result.mix(falcon_result)


process addFragpipePsmData {
  container "${python_container}"
  cpus 1
  memory "1GB"
  publishDir "clustering_results", mode: "link"

  input:
  tuple val(tool), file(result_file) from clustering_result_files
  file psm_result

  output:
  tuple val(tool), file("${tool}_psm.tsv") into clustering_psm_results

  script:
  """
#!/usr/bin/env python3


import os
import sys
import csv
import re


scan_pattern = re.compile(r'.*\\.([0-9]*)\\.[0-9]\$')


# load the PSMs first
psms = dict()

with open("${psm_result}", "r") as psm_file:
  reader = csv.DictReader(psm_file, delimiter="\\t")

  for psm_entry in reader:
    # filename of the actual input file
    org_file = psm_entry["Spectrum File"][9:-8]

    # extract the scan number
    scan_match = scan_pattern.match(psm_entry["Spectrum"])

    if not scan_match:
      print("Error: Failed to parse spectrum id: " + psm_entry["Spectrum"])
      sys.exit(1)

    spectrum_number = int(scan_match.group(1))

    psm = {
      "filename": org_file,
      "spectrum": str(spectrum_number),
      "peptide": psm_entry["Peptide"],
      "ptms": psm_entry["Assigned Modifications"],
      "retention": psm_entry["Retention"]
    }

    # create the key and save the PSM
    key = psm["filename"] + "." + psm["spectrum"]

    if key in psms:
      print("Error: Duplicated PSMs for " + key)
      sys.exit(1)

    psms[key] = psm

# add the PSM info to the clustering result file
with open("${result_file}", "r") as reader:
  with open("${tool}_psm.tsv", "w") as writer:
    # update the header
    header = reader.__next__().strip()
    new_header = header + "\\tpeptide\\tptms\\tretention\\n"
    writer.write(new_header)

    # process the results
    for line in reader:
      fields = line.strip().split("\\t")
      
      # create the PSM key
      filename = fields[0]
      spectrum = int(fields[1])
      psm_key = filename + "." + str(spectrum)

      if psm_key in psms:
        psm = psms[psm_key]
        fields += [psm["peptide"], psm["ptms"], psm["retention"]]
      else:
        fields += ["NA", "NA", "NA"]

      writer.write("\\t".join(fields) + "\\n")

  """
}

process analyseClustering {
  cpus 1
  memory "1GB"
  // TODO: run in R container
  publishDir "clustering_results", mode: "link"

  input:
  file("*") from clustering_psm_results.collect()

  output:
  file("tool_stats.tsv") into tool_stats
  file("*.png") into clustering_plots

  script:
  template 'clustering_stats.R'
}
