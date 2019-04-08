#!/usr/bin/env bash

##### VARIABLES
ACCESSION=$1
USERNAME=$2
PASSWORD=$3
METADATA_DIR=$4
SCRIPT_DIR="$( cd "$(dirname "$0")" ; pwd -P )"

# run Nextflow pipeline
${NEXTFLOW}nextflow run ${SCRIPT_DIR}/thermo-convert-nf.nf -name $ACCESSION \
                                  -c ${SCRIPT_DIR}/nextflow.config \
                                  -profile cluster \
                                  --px_accession $ACCESSION \
                                  --pride_username $USERNAME \
                                  --pride_password $PASSWORD \
                                  --metadata_path $METADATA_DIR + "/data/" + $ACCESSION
# Clean working directories
echo "/nCleanning the working directory..."
${NEXTFLOW}nextflow clean -f $ACCESSION