# spectra-cluster LFQ pipeline

This pipeline uses spectral-clustering to infer additional identifications from
MS/MS search engine results and improve label-free quantitation accuracy.

## Usage

The pipeline uses [nextflow](https://nextflow.io) and [docker](https://docker.com) and therefore
only runs on Linux systems.

### Test dataset

Once nextflow and docker are installed, the example workflow can simply be run using:

```bash
# The pipeline expects that the current working directory is the directory
# containing the "config" and the "test" folder
nextflow run -resume lfq-clustering.nf
```

This will process one included MGF file of the CPTAC Study 6 dataset and search it against the human SwissProt database (**without contaminants**). Using a database without contaminants is generally not recommended and only included for testing purposes.

### Processing own data

```bash
# The pipeline expects that the current working directory is the directory
# containing the "config" folder
nextflow run -resume lfq-clustering.nf --prec_rol 10 --frag_tol 0.5 --mc 1 \
	--min_ident 2 --min_ratio 0.7 \
	--raw_dir /my/mgf/file_directory --fasta_file /my/species.fasta
```

The above command will perform the following steps:

  1) Append reversed decoy sequences to the FASTA database
  2) Search all MGF files using X!Tandem and MSGF+ - both result sets are processed separately
  3) Cluster the MGF files
  4) Infer additional PSMs based on the clustering results
  5) Report the clustering results and a table with the complete list of PSMs in the `result` directory.

### Parameters

  * prec_tol: The precursor tolerance to use (in ppm)
  * frag_tol: The fragment tolerance to use (in Da)
  * mc: Number of missed cleavages
  * min_ident: Minimum number of identified spectra within a cluster to tranfer identifications
  * min_ratio: Minimum proportion of spectra that have to be identified as the same peptide within a cluster to transfer identifications.
  * raw_dir: Path to the directory containing the MGF files (no trailing "/")
  * fasta_file: Path to the FASTA file to use (without decoy sequences).

### Specifying PTMs

By default, the workflow uses Carbamidomethylation (C) as fixed modification, and Oxidation (M) and N-term Acetylation as variable modifications. To change this, you need to manually edit the `Mods.txt` and the `input.xml` file in the `config` directory for MSGF+ and X!Tandem respectively.

Alterantively, you can create a copy of these file, edit this copy and set the `--xtandem_template` and `--msgf_mods` command line parameters accordingly.

### Multi-threading

By default, all tasks that can be run in parallel will be run in parallel by nextflow. Additionally, you can change the `threads` variable in the workflow script. This will then increase the number of threads used by the search engines and the clustering process. Note that the searches and the clustering might run at the same time.
