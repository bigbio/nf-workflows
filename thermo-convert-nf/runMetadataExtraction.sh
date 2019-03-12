#!/usr/bin/env bash

PX_ACCESSION=$1
nextflow run thermo-convert-nf.nf --px_accession $PX_ACCESSION -name $PX_ACCESSION -c nextflow.config
#&& nextflow clean -f $PX_ACCESSION