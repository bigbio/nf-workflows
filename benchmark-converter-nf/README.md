Label-Free Protein Quantitation Pipeline
==================================================

This pipeline uses the popular proteomics framework OpenMS to perform protein quantitative analysis. The current pipeline 
containers the following steps: 

- Peptide Identification using the MSGF+ search engine.
   - The database containing the decoys proteins should be provided. 
- Peptide FDR calculation and filtering 1% PSM FDR. 
- LFQ with OpenMS ProteomicsLFQ for iPRG    

## Peptide identification in LSF

### ID
Note: the identification db is defined in the `.ini` file in configs folder.
```
cd /hps/nobackup/proteomics/conversion_comparison/<peakfile input>/
mkdir ../msgf_newconverts
for f in *.mzML; do if [[ ! -r ../msgf_newconverts/${f%.*}.idXML  ]] ; then  bsub -e ${f%.*}.msgf.err -o ${f%.*}.msgf.out -R "select[singularity]" -M 8192 -R "rusage[mem=8192]" "singularity exec docker://mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq MSGFPlusAdapter -ini ../msgf.ini -in ${f%.*}.mzML -out ../msgf_newconverts/${f%.*}.idXML"; fi; done
```

### FDR calc.
Note: this will remove non-scoring PSM
```
cd /hps/nobackup/proteomics/conversion_comparison/msgf_newconverts/
mkdir ../ids_newconverts
for f in *.idXML; do if [[ ! -r ../ids_newconverts/${f%.*}.idXML  ]] ; then  bsub -e ${f%.*}.msgf.err -o ${f%.*}.msgf.out -R "select[singularity]" -M 8192 -R "rusage[mem=8192]" "singularity exec docker://mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq FalseDiscoveryRate -ini ../fdr.ini -in $f -out ../ids_newconverts/${f%.*}.idXML"; fi; done
```

## LFQ
Note: fixed relative paths in msstats config file 
```
cd /hps/nobackup/proteomics/conversion_comparison/
bsub -R "select[singularity]" -M 18192 -R "rusage[mem=18192]" -n 10 "singularity exec docker://mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq bash run_lfq.sh"
```
