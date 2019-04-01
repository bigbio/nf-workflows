#!/usr/bin/env bash

PX_ACCESSION=$1
nextflow run thermo-convert-nf.nf -name $PX_ACCESSION \
                                  -c nextflow.config \
                                  --px_accession $PX_ACCESSION \
                                  --pride_username "sureshhewabi@gmail.com" \
                                  --pride_password "********"
nextflow clean -f $PX_ACCESSION