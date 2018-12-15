#!bin/bash
###########################
#

#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample1_A.mzML -out ./iPRG2015/JD_06232014_sample1_A.idXML &
#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample2_A.mzML -out ./iPRG2015/JD_06232014_sample2_A.idXML &
#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample3_A.mzML -out ./iPRG2015/JD_06232014_sample3_A.idXML &
#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample4_A.mzML -out ./iPRG2015/JD_06232014_sample4_A.idXML &
#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample1_B.mzML -out ./iPRG2015/JD_06232014_sample1_B.idXML &
#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample2_B.mzML -out ./iPRG2015/JD_06232014_sample2_B.idXML &
#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample3_B.mzML -out ./iPRG2015/JD_06232014_sample3_B.idXML &
#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample4_B.mzML -out ./iPRG2015/JD_06232014_sample4_B.idXML &
#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample1_C.mzML -out ./iPRG2015/JD_06232014_sample1_C.idXML &
#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample2_C.mzML -out ./iPRG2015/JD_06232014_sample2_C.idXML &
#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample3_C.mzML -out ./iPRG2015/JD_06232014_sample3_C.idXML &
#MSGFPlusAdapter -ini msgf.ini -in ./iPRG2015/JD_06232014_sample4_C.mzML -out ./iPRG2015/JD_06232014_sample4_C.idXML &
#mkdir ../msgf
#cd iPRG2015
#for f in *.mzML; do if [[ ! -r ../msgf/${f%.*}.idXML  ]] ; then  bsub -e ${f%.*}.msgf.err -o ${f%.*}.msgf.out -R "select[singularity]" -M 8192 -R "rusage[mem=8192]" "singularity exec docker://mwalzer/openms-batteries-included:V2.3.0_pepxmlpatch MSGFPlusAdapter -ini ../msgf.ini -in ${f%.*}.mzML -out ../msgf/${f%.*}.idXML"; fi; done
mkdir new_msgf
cd new_iPRG2015 
for f in *.mzML; do if [[ ! -r ../new_msgf/${f%.*}.idXML  ]] ; then  bsub -e ${f%.*}.msgf.err -o ${f%.*}.msgf.out -P rh74 -M 8192 -R "rusage[mem=8192]" "singularity exec docker://mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq MSGFPlusAdapter -ini ../msgf.ini -in ${f%.*}.mzML -out ../new_msgf/${f%.*}.idXML"; fi; done
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample1_A.idXML -out ./iPRG2015/JD_06232014_sample1_A.idXML
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample2_A.idXML -out ./iPRG2015/JD_06232014_sample2_A.idXML
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample3_A.idXML -out ./iPRG2015/JD_06232014_sample3_A.idXML
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample4_A.idXML -out ./iPRG2015/JD_06232014_sample4_A.idXML
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample1_B.idXML -out ./iPRG2015/JD_06232014_sample1_B.idXML
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample2_B.idXML -out ./iPRG2015/JD_06232014_sample2_B.idXML
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample3_B.idXML -out ./iPRG2015/JD_06232014_sample3_B.idXML
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample4_B.idXML -out ./iPRG2015/JD_06232014_sample4_B.idXML
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample1_C.idXML -out ./iPRG2015/JD_06232014_sample1_C.idXML
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample2_C.idXML -out ./iPRG2015/JD_06232014_sample2_C.idXML
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample3_C.idXML -out ./iPRG2015/JD_06232014_sample3_C.idXML
#PeptideIndexer -ini pi.ini -in ./iPRG2015/JD_06232014_sample4_C.idXML -out ./iPRG2015/JD_06232014_sample4_C.idXML
#mkdir ../ids
#cd ../msgf
#for f in *.idXML; do if [[ ! -r ../ids/${f%.*}.idXML ]] ; then  bsub -e ${f%.*}.pi.err -o ${f%.*}.pi.out -R "select[singularity]" -M 8192 -R "rusage[mem=8192]" "singularity exec docker://mwalzer/openms-batteries-included:V2.3.0_pepxmlpatch PeptideIndexer -ini pi.ini  -in $f -out ../ids/$f"; fi; done
cd ../new_msgf
for f in *.idXML; do if [[ ! -r ../new_ids/${f%.*}.idXML ]] ; then  bsub -e ${f%.*}.pi.err -o ${f%.*}.pi.out -P rh74 -M 8192 -R "rusage[mem=8192]" "singularity exec docker://mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq PeptideIndexer -ini ../pi.ini  -in $f -out ../new_ids/$f"; fi; done
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample1_A.idXML -out ./iPRG2015/JD_06232014_sample1_A.idXML
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample2_A.idXML -out ./iPRG2015/JD_06232014_sample2_A.idXML
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample3_A.idXML -out ./iPRG2015/JD_06232014_sample3_A.idXML
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample4_A.idXML -out ./iPRG2015/JD_06232014_sample4_A.idXML
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample1_B.idXML -out ./iPRG2015/JD_06232014_sample1_B.idXML
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample2_B.idXML -out ./iPRG2015/JD_06232014_sample2_B.idXML
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample3_B.idXML -out ./iPRG2015/JD_06232014_sample3_B.idXML
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample4_B.idXML -out ./iPRG2015/JD_06232014_sample4_B.idXML
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample1_C.idXML -out ./iPRG2015/JD_06232014_sample1_C.idXML
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample2_C.idXML -out ./iPRG2015/JD_06232014_sample2_C.idXML
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample3_C.idXML -out ./iPRG2015/JD_06232014_sample3_C.idXML
#FalseDiscoveryRate -ini fdr.ini -in ./iPRG2015/JD_06232014_sample4_C.idXML -out ./iPRG2015/JD_06232014_sample4_C.idXML
#cd ../ids
#for f in *.idXML; do if [[ -r ../ids/${f%.*}.idXML ]] ; then  bsub -e ${f%.*}.fdr.err -o ${f%.*}.fdr.out -R "select[singularity]" -M 8192 -R "rusage[mem=8192]" "singularity exec docker://mwalzer/openms-batteries-included:V2.3.0_pepxmlpatch FalseDiscoveryRate -ini fdr.ini  -in $f -out $f"; fi; done
cd ../new_ids
for f in *.idXML; do if [[ -r ../new_ids/${f%.*}.idXML ]] ; then  bsub -e ${f%.*}.fdr.err -o ${f%.*}.fdr.out -P rh74 -M 8192 -R "rusage[mem=8192]" "singularity exec docker://mwalzer/openms-batteries-included:V2.4.0_proteomic_lfq FalseDiscoveryRate -ini ../fdr.ini  -in $f -out $f"; fi; done
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample1_A.idXML -out ./iPRG2015/JD_06232014_sample1_A.idXML
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample2_A.idXML -out ./iPRG2015/JD_06232014_sample2_A.idXML
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample3_A.idXML -out ./iPRG2015/JD_06232014_sample3_A.idXML
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample4_A.idXML -out ./iPRG2015/JD_06232014_sample4_A.idXML
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample1_B.idXML -out ./iPRG2015/JD_06232014_sample1_B.idXML
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample2_B.idXML -out ./iPRG2015/JD_06232014_sample2_B.idXML
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample3_B.idXML -out ./iPRG2015/JD_06232014_sample3_B.idXML
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample4_B.idXML -out ./iPRG2015/JD_06232014_sample4_B.idXML
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample1_C.idXML -out ./iPRG2015/JD_06232014_sample1_C.idXML
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample2_C.idXML -out ./iPRG2015/JD_06232014_sample2_C.idXML
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample3_C.idXML -out ./iPRG2015/JD_06232014_sample3_C.idXML
#IDFilter -ini idf.ini -in ./iPRG2015/JD_06232014_sample4_C.idXML -out ./iPRG2015/JD_06232014_sample4_C.idXML

cd new_ids
for f in *.idXML; do if [[ -r ../new_ids/${f%.*}.idXML ]] ; then  bsub -e ${f%.*}.qc.err -o ${f%.*}.qc.out -P rh74 -M 8192 -R "rusage[mem=8192]" "singularity exec docker://mwalzer/openms-batteries-included:V2.3.0_pepxmlpatch QCCalculator -in  ../new_iPRG2015/${f%.*}.mzML  -id $f -out ../new_qc/${f%.*}.qcML"; fi; done

cd ids
for f in *.idXML; do if [[ -r ../ids/${f%.*}.idXML ]] ; then  bsub -e ${f%.*}.qc.err -o ${f%.*}.qc.out -P rh74 -M 8192 -R "rusage[mem=8192]" "singularity exec docker://mwalzer/openms-batteries-included:V2.3.0_pepxmlpatch QCCalculator -in  ../iPRG2015/${f%.*}.mzML  -id $f -out ../qc/${f%.*}.qcML"; fi; done

#
###########################

ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./ids/JD_06232014_sample1_A.idXML  ./ids/JD_06232014_sample1_B.idXML  ./ids/JD_06232014_sample1_C.idXML \
./ids/JD_06232014_sample2_A.idXML  ./ids/JD_06232014_sample2_B.idXML  ./ids/JD_06232014_sample2_C.idXML \
./ids/JD_06232014_sample3_A.idXML  ./ids/JD_06232014_sample3_B.idXML  ./ids/JD_06232014_sample3_C.idXML \
./ids/JD_06232014_sample4_A.idXML  ./ids/JD_06232014_sample4_B.idXML  ./ids/JD_06232014_sample4_C.idXML \
-design ./experimental_design.tsv \
-Alignment:max_rt_shift 0.1 \
-fasta ./database/iPRG2015_target_decoy_nocontaminants.fasta -targeted_only "true" \
-transfer_ids "false" \
-out_msstats "iPRG2015_targeted_only.csv" \
-out_cxml "iPRG2015_targeted_only.consensusXML" \
-out iPRG2015_targeted_only.mzTab -threads 10 > iPRG2015_targeted_only.log 

FileInfo -in iPRG2015_targeted_only.consensusXML > iPRG2015_targeted_only.fileinfo

exit

ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./ids/JD_06232014_sample1_A.idXML  ./ids/JD_06232014_sample1_B.idXML  ./ids/JD_06232014_sample1_C.idXML \
./ids/JD_06232014_sample2_A.idXML  ./ids/JD_06232014_sample2_B.idXML  ./ids/JD_06232014_sample2_C.idXML \
./ids/JD_06232014_sample3_A.idXML  ./ids/JD_06232014_sample3_B.idXML  ./ids/JD_06232014_sample3_C.idXML \
./ids/JD_06232014_sample4_A.idXML  ./ids/JD_06232014_sample4_B.idXML  ./ids/JD_06232014_sample4_C.idXML \
-design ./experimental_design.tsv \
-Alignment:max_rt_shift 0.1 \
-fasta ./database/iPRG2015_target_decoy_nocontaminants.fasta -targeted_only "true" \
-transfer_ids "merged" \
-out_msstats "iPRG2015_targeted_only_merged.csv" \
-out_cxml "iPRG2015_targeted_only_merged.consensusXML" \
-out iPRG2015_targeted_only_merged.mzTab -threads 10 > iPRG2015_targeted_only_merged.log 

FileInfo -in iPRG2015_targeted_only_merged.consensusXML > iPRG2015_targeted_only_merged.fileinfo

ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./ids/JD_06232014_sample1_A.idXML  ./ids/JD_06232014_sample1_B.idXML  ./ids/JD_06232014_sample1_C.idXML \
./ids/JD_06232014_sample2_A.idXML  ./ids/JD_06232014_sample2_B.idXML  ./ids/JD_06232014_sample2_C.idXML \
./ids/JD_06232014_sample3_A.idXML  ./ids/JD_06232014_sample3_B.idXML  ./ids/JD_06232014_sample3_C.idXML \
./ids/JD_06232014_sample4_A.idXML  ./ids/JD_06232014_sample4_B.idXML  ./ids/JD_06232014_sample4_C.idXML \
-design ./experimental_design.tsv \
-Alignment:max_rt_shift 0.1 \
-fasta ./database/iPRG2015_target_decoy_nocontaminants.fasta -targeted_only "true" \
-transfer_ids "SVM" \
-out_msstats "iPRG2015_targeted_only_SVM.csv" \
-out_cxml "iPRG2015_targeted_only_SVM.consensusXML" \
-out iPRG2015_targeted_only_SVM.mzTab -threads 10 > iPRG2015_targeted_only_SVM.log 

FileInfo -in iPRG2015_targeted_only_SVM.consensusXML > iPRG2015_targeted_only_SVM.fileinfo
exit 1

ProteomicsLFQ -in \
./new_iPRG2015/JD_06232014_sample1_A.mzML  ./new_iPRG2015/JD_06232014_sample1_B.mzML  ./new_iPRG2015/JD_06232014_sample1_C.mzML \
./new_iPRG2015/JD_06232014_sample2_A.mzML  ./new_iPRG2015/JD_06232014_sample2_B.mzML  ./new_iPRG2015/JD_06232014_sample2_C.mzML \
./new_iPRG2015/JD_06232014_sample3_A.mzML  ./new_iPRG2015/JD_06232014_sample3_B.mzML  ./new_iPRG2015/JD_06232014_sample3_C.mzML \
./new_iPRG2015/JD_06232014_sample4_A.mzML  ./new_iPRG2015/JD_06232014_sample4_B.mzML  ./new_iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./new_ids/JD_06232014_sample1_A.idXML  ./new_ids/JD_06232014_sample1_B.idXML  ./new_ids/JD_06232014_sample1_C.idXML \
./new_ids/JD_06232014_sample2_A.idXML  ./new_ids/JD_06232014_sample2_B.idXML  ./new_ids/JD_06232014_sample2_C.idXML \
./new_ids/JD_06232014_sample3_A.idXML  ./new_ids/JD_06232014_sample3_B.idXML  ./new_ids/JD_06232014_sample3_C.idXML \
./new_ids/JD_06232014_sample4_A.idXML  ./new_ids/JD_06232014_sample4_B.idXML  ./new_ids/JD_06232014_sample4_C.idXML \
-design ./new_experimental_design.tsv \
-Alignment:max_rt_shift 0.1 \
-fasta ./database/iPRG2015_target_decoy_nocontaminants.fasta -targeted_only "true" \
-transfer_ids "false" \
-out_msstats "new_iPRG2015_targeted_only.csv" \
-out_cxml "new_iPRG2015_targeted_only.consensusXML" \
-out new_iPRG2015_targeted_only.mzTab -threads 10 > new_iPRG2015_targeted_only.log

FileInfo -in new_iPRG2015_targeted_only.consensusXML > new_iPRG2015_targeted_only.fileinfo

exit

ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./ids/JD_06232014_sample1_A.idXML  ./ids/JD_06232014_sample1_B.idXML  ./ids/JD_06232014_sample1_C.idXML \
./ids/JD_06232014_sample2_A.idXML  ./ids/JD_06232014_sample2_B.idXML  ./ids/JD_06232014_sample2_C.idXML \
./ids/JD_06232014_sample3_A.idXML  ./ids/JD_06232014_sample3_B.idXML  ./ids/JD_06232014_sample3_C.idXML \
./ids/JD_06232014_sample4_A.idXML  ./ids/JD_06232014_sample4_B.idXML  ./ids/JD_06232014_sample4_C.idXML \
-design ./experimental_design.tsv \
-Alignment:max_rt_shift 0 \
-fasta ./database/iPRG2015_target_decoy_nocontaminants.fasta \
-out iPRG2015_2.mzTab -debug 667 -threads 10 > iPRG2015.log 

FileInfo -in debug_consensus.consensusXML > iPRG2015.fileinfo

ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./ids/JD_06232014_sample1_A.idXML  ./ids/JD_06232014_sample1_B.idXML  ./ids/JD_06232014_sample1_C.idXML \
./ids/JD_06232014_sample2_A.idXML  ./ids/JD_06232014_sample2_B.idXML  ./ids/JD_06232014_sample2_C.idXML \
./ids/JD_06232014_sample3_A.idXML  ./ids/JD_06232014_sample3_B.idXML  ./ids/JD_06232014_sample3_C.idXML \
./ids/JD_06232014_sample4_A.idXML  ./ids/JD_06232014_sample4_B.idXML  ./ids/JD_06232014_sample4_C.idXML \
-design ./experimental_design.tsv \
-Alignment:min_run_occur 9 \
-Alignment:max_rt_shift 0 \
-fasta ./database/iPRG2015_target_decoy_nocontaminants.fasta \
-out iPRG2015_mrc9.mzTab -debug 667 -threads 10 > iPRG2015mrc9.log 

FileInfo -in debug_consensus.consensusXML > iPRG2015mrc9.fileinfo



#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/1A.idXML -out ./iPRG2015/PeptideAtlas/idXML/1A_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/2A.idXML -out ./iPRG2015/PeptideAtlas/idXML/2A_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/3A.idXML -out ./iPRG2015/PeptideAtlas/idXML/3A_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/4A.idXML -out ./iPRG2015/PeptideAtlas/idXML/4A_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/1B.idXML -out ./iPRG2015/PeptideAtlas/idXML/1B_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/2B.idXML -out ./iPRG2015/PeptideAtlas/idXML/2B_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/3B.idXML -out ./iPRG2015/PeptideAtlas/idXML/3B_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/4B.idXML -out ./iPRG2015/PeptideAtlas/idXML/4B_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/1C.idXML -out ./iPRG2015/PeptideAtlas/idXML/1C_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/2C.idXML -out ./iPRG2015/PeptideAtlas/idXML/2C_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/3C.idXML -out ./iPRG2015/PeptideAtlas/idXML/3C_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/4C.idXML -out ./iPRG2015/PeptideAtlas/idXML/4C_IDF.idXML
#
ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./PeptideAtlas/idXML/1A_IDF.idXML  ./PeptideAtlas/idXML/1B_IDF.idXML  ./PeptideAtlas/idXML/1C_IDF.idXML \
./PeptideAtlas/idXML/2A_IDF.idXML  ./PeptideAtlas/idXML/2B_IDF.idXML  ./PeptideAtlas/idXML/2C_IDF.idXML \
./PeptideAtlas/idXML/3A_IDF.idXML  ./PeptideAtlas/idXML/3B_IDF.idXML  ./PeptideAtlas/idXML/3C_IDF.idXML \
./PeptideAtlas/idXML/4A_IDF.idXML  ./PeptideAtlas/idXML/4B_IDF.idXML  ./PeptideAtlas/idXML/4C_IDF.idXML \
-design ./iPRG2015/experimental_design.tsv \
-Alignment:max_rt_shift 0 \
-fasta ./iPRG2015/iPRG2015_decoy.fasta \
-out iPRG2015_PeptideAtlas.mzTab -debug 667 -threads 10 > iPRG2015_PeptideAtlas.log 

FileInfo -in debug_consensus.consensusXML > iPRG2015PeptideAtlas.fileinfo


