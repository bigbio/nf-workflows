# Tide / MSGF+ Reprocessing pipeline

This pipeline reprocesses a set of MS/MS RAW files and processes them using MSGF+ and X!Tandem through SearchGUI. Search results are combined using PIA.

## How to run the pipeline

nextflow run main.nf -c nextflow.config -profile local,standard
