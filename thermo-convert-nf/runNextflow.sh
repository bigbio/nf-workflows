#!/usr/bin/env bash

ACCESSION=$1
USERNAME=$2
PASSWORD=$3
NEXTFLOW_DIR=$4
${NEXTFLOW_DIR}nextflow run thermo-convert-nf.nf -name $ACCESSION \
                                  -c nextflow.config \
                                  --px_accession $ACCESSION \
                                  --pride_username $USERNAME \
                                  --pride_password $PASSWORD \
                                  --metadata_path $NEXTFLOW_DIR + "data/" + $ACCESSION
${NEXTFLOW_DIR}nextflow clean -f $ACCESSION