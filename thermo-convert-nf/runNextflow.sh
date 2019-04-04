#!/usr/bin/env bash

ACCESSION=$1
USERNAME=$2
PASSWORD=$3
METADATA_DIR=$4

${NEXTFLOW}nextflow run thermo-convert-nf.nf -name $ACCESSION \
                                  -c nextflow.config \
                                  --px_accession $ACCESSION \
                                  --pride_username $USERNAME \
                                  --pride_password $PASSWORD \
                                  --metadata_path $METADATA_DIR + "/data/" + $ACCESSION
${NEXTFLOW}nextflow clean -f $ACCESSION