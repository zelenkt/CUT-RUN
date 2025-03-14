ml BEDTools/2.30.0-GCC-11.2.0

# combine peaks from different replicates
cd /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k27ac/raw_files 
cat *_001_non_stringent.stringent.bed | sortBed | mergeBed -c 4,5,6 -o sum,max,collapse | sortBed > final_h3k27ac_cat_merge_001_non_stringent.stringent.bed
cat raw/seacr*ko_001* | sortBed | mergeBed -c 4,5,6 -o sum,max,collapse | sortBed > final_ko_h3k27ac_cat_merge_001_non_stringent.stringent.bed
cat raw/seacr*wt_001* | sortBed | mergeBed -c 4,5,6 -o sum,max,collapse | sortBed > final_wt_h3k27ac_cat_merge_001_non_stringent.stringent.bed

cd /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k4me1/raw_files 
cat *_001_non_stringent.stringent.bed | sortBed | mergeBed -c 4,5,6 -o sum,max,collapse | sortBed > final_h3k4me1_cat_merge_001_non_stringent.stringent.bed
cat raw/seacr*ko_001* | sortBed | mergeBed -c 4,5,6 -o sum,max,collapse | sortBed > final_ko_h3k4me1_cat_merge_001_non_stringent.stringent.bed
cat raw/seacr*wt_001* | sortBed | mergeBed -c 4,5,6 -o sum,max,collapse | sortBed > final_wt_h3k4me1_cat_merge_001_non_stringent.stringent.bed


# for the purpose of regulatory elements analysis, take all regions with H3K4me1 signal
awk -v OFS="\t" '{print $1,"wtTE"NR,"wtTE"NR,$2,$3,$4,".","wtTE"NR,"wtTE"NR}' final_wt_h3k4me1_cat_merge_001_non_stringent.stringent.bed > wt_TE_K4m1.gff
awk -v OFS="\t" '{print $1,"koTE"NR,"koTE"NR,$2,$3,$4,".","koTE"NR,"koTE"NR}' final_ko_h3k4me1_cat_merge_001_non_stringent.stringent.bed > ko_TE_K4m1.gff


cd /share/lab_avram/HPC_Cluster/user/tomas/rose-docker
singularity run rose_latest.sif ROSE_main.py -g mm10 -i ./tz_data/wt_TE_K4m1.gff -r ./tz_data/merged_wt-60-91_h3k27ac.bam -o output_wt_TE-with27acBAM -s 12500 -t 2500
singularity run rose_latest.sif ROSE_main.py -g mm10 -i ./tz_data/wt_TE_K4m1.gff -r ./tz_data/merged_wt-17-60_h3k27me3.bam -o output_wt_TE-withK27me3BAM_PE -s 12500 -t 2500
singularity run rose_latest.sif ROSE_main.py -g mm10 -i ./tz_data/ko_TE_K4m1.gff -r ./tz_data/merged_ko-47-48_h3k27ac.bam -o output_ko_TE-with27acBAM -s 12500 -t 2500
singularity run rose_latest.sif ROSE_main.py -g mm10 -i ./tz_data/ko_TE_K4m1.gff -r ./tz_data/merged_ko-38-61_h3k27me3.bam -o output_ko_TE-withK27me3BAM_PE -s 12500 -t 2500


# FOR ROSE PLOTS with the ranking, add order to the ROSE ranked SEs and plot them in R by myself
cd output_wt_TE-with27acBAM 
awk -F'\t' -v OFS="\t" 'NR>1{print $0}' wt_TE_K4m1_AllStitched.table_withGENES.txt | awk -F'\t' -v OFS="\t" '{print $0, NR}'> wt_TE_K4m1-with27acBAM_AllStitched.table_withGENES_with_order.txt
cd ../output_ko_TE-with27acBAM 
awk -F'\t' -v OFS="\t" 'NR>1{print $0}' ko_TE_K4m1_AllStitched.table_withGENES.txt | awk -F'\t' -v OFS="\t" '{print $0, NR}'> ko_TE_K4m1-with27acBAM_AllStitched.table_withGENES_with_order.txt
cd ../output_wt_TE-withK27me3BAM_PE
awk -F'\t' -v OFS="\t" 'NR>1{print $0}' wt_TE_K4m1_AllStitched.table_withGENES.txt | awk -F'\t' -v OFS="\t" '{print $0, NR}'> wt_TE_K4m1-withK27me3BAM_PE_AllStitched.table_withGENES_with_order.txt
cd ../output_ko_TE-withK27me3BAM_PE
awk -F'\t' -v OFS="\t" 'NR>1{print $0}' ko_TE_K4m1_AllStitched.table_withGENES.txt | awk -F'\t' -v OFS="\t" '{print $0, NR}'> ko_TE_K4m1-withK27me3BAM_PE_AllStitched.table_withGENES_with_order.txt




# # to SUMMARIZE THE REGULATORY ELEMENTS FILES GENERATED:
# # these are all TE elements, i.e. non-filtered H3K4me1 peaks (non-stitched)
# /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/rose-docker/all_bed_super-elements/h3k4me1_merged_all.bed

# # these are the filtered SE (i.e. they have H3K4me1 signal and strong H3K27ac). To make it comparable between WT and KO, I have to combine WT and KO set (used for deeptools):
# # output_wt_TE-with27acBAM.bed
# # output_ko_TE-with27acBAM.bed
# cd /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/rose-docker/all_bed_super-elements
# cat output_wt_TE-with27acBAM.bed output_ko_TE-with27acBAM.bed | sortBed | mergeBed -c 4,5,6,7 -o collapse,mean,collapse,collapse | sortBed  | uniq > SE-Filt_k4me1_with_k27acBAM_merged_all.bed

# # these are all the stitched TE which includes SE but also regions that were not classified as SE
# # output_wt_TE-with27acBAM_nonFilt.bed
# # output_ko_TE-with27acBAM_nonFilt.bed
# cd /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/rose-docker/all_bed_super-elements
# cat output_wt_TE-with27acBAM_nonFilt.bed output_ko_TE-with27acBAM_nonFilt.bed | sortBed | mergeBed -c 4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse | sortBed  | uniq > SE-nonFilt_k4me1_with_k27acBAM_merged_all.bed

# # these are the filtered PE (i.e. they have H3K4me1 signal and strong H3K27me3). To make it comparable between WT and KO, I have to combine WT and KO set (used for deeptools):
# # output_wt_TE-withK27me3BAM_PE.bed
# # output_ko_TE-withK27me3BAM_PE.bed
# cd /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/rose-docker/all_bed_super-elements
# cat output_wt_TE-withK27me3BAM_PE.bed output_ko_TE-withK27me3BAM_PE.bed | sortBed | mergeBed -c 4,5,6,7 -o collapse,mean,collapse,collapse | sortBed  | uniq > PE-Filt_k4me1_with_k27me3BAM_merged_all.bed

# # these are all the stitched TE which includes PE but also regions that were not classified as PE, coordinates-wise, this should be the same as the stitched-TE file which includes SE but also regions that were not classified as SE, i.e. SE-nonFilt_k4me1_with_k27acBAM_merged_all.bed
# # output_wt_PE-with27me3BAM_nonFilt.bed
# # output_ko_PE-with27me3BAM_nonFilt.bed
# cd /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/rose-docker/all_bed_super-elements
# cat output_wt_PE-with27me3BAM_nonFilt.bed output_ko_PE-with27me3BAM_nonFilt.bed | sortBed | mergeBed -c 4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse | sortBed  | uniq > PE-nonFilt_k4me1_with_k27me3BAM_merged_all.bed











# in the meantime, get RNAseq as pseudobulk extracted from scRNAseq which will be better for the plot
# modify scRNAseq output so that it looks like the bulk data
# the original scRNA version is:
# gene geneID geneSymbol_with_dupl  pval    avg_log2FC  pct.1   pct.2   adjPval
# on 12-17-2024 I updated these commands using the most recent scRNAseq data after some additional filtering we performed
# cp /Users/4476125/Documents/BERLIN_Shaw_lab_2023/BERLIN/2-Single_Cell_RNAseq_Pipeline/Output/ovarT_integrated_noRefUsed/reclustered/markers/allKOvsWT_ANNOTATED_ovarT_assayRNA.txt ./
# or run it on shared drive directly and distribute the scRNAseq file in all SE, TE, PE directories such as DESeq2_poisedenhancer-noFilt_01-16-2025
# cd /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27me3/cd8_ovar_tumor_mm10/DESeq2_poisedenhancer-noFilt_01-16-2025
cat allKOvsWT_ANNOTATED_ovarT_assayRNA.txt | tr -d '"'| awk  -F'\t' -v OFS="\t" 'NR>1{print $2,$1,$6-$7,$5,$6,$7,$4,$8}' > tmp1315; awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' /Volumes/lab_avram/HPC_Cluster/user/tomas/lib/geneNames_refgenie_ensembl_mm10_description_longestIsoLength.gtf tmp1315 | awk -F'\t' -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$11}' | awk -F'\t' -v OFS="\t" '{for(i=1; i<=9; i++) if ($i=="") $i="NA"; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > scRNAseq_adapted_markers.seurat_obj_reclust2_12-17-2024.rna_allClustTogether.txt; rm tmp1315
# new columns are geneID    gene  pct.1-pct.2  avg_log2FC  pct.1   pct.2  pval   adjPval    description


rna_dataset=scRNAseq_adapted_markers.seurat_obj_reclust2_12-17-2024.rna_allClustTogether.txt
# peaks=all_diff_peaks_noFilt_transcriptID.bed
peaks1=UPDATE_12-24-2024_TSS0bp___all_diff_peaks_noFilt_transcriptID.bed
peaks2=UPDATE_12-24-2024_TSSall___all_diff_peaks_noFilt_transcriptID.bed








# for TE (traditional enhancers)
awk -v OFS="\t" '{print $1,$2,$3,"h3k4me1_"NR,$4,"."}' final_h3k4me1_cat_merge_001_non_stringent.stringent.bed | sortBed | uniq > h3k4me1_merged_all.bed
cd /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10
sample1="CRK27AM8TW60"
sample2="CRK27AM8TW91"
sample3="M8CTK47A"
sample4="M8CTK48A"
myArray=(${sample1} ${sample2} ${sample3} ${sample4})
peakset="h3k4me1"
for sample in ${myArray[@]}; do
  bedtools coverage -a /share/lab_avram/HPC_Cluster/user/tomas/rose-docker/all_bed_super-elements/${peakset}_merged_all.bed -b ${sample}/${sample}.bam | uniq > ${sample}_${peakset}_merged_peaks_coverage.bed
done
# extract raw reads and calculate differential signal at TE 
# change line 140-141 to column $7 in ~/bin/CUTRUN_DESeq2_for_seacr.sh
sbatch --chdir=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10 --export=s=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10/samples.csv,f=_h3k4me1_merged_peaks_coverage.bed ~/bin/CUTRUN_DESeq2_for_seacr.sh
mv DESeq2_seacr DESeq2_enhancer_01-16-2025
cd DESeq2_enhancer_01-16-2025
# paste all_differential_peaks_noFilt_for_annotation.bed ../../../../rose-docker/all_bed_super-elements/${peakset}_merged_all.bed | awk -F'\t' -v OFS="\t" '{print $1,$2,$3,$4,$5,$7,$8,$9,$10,$14,$15,$20,$19}' > all_differential_peaks_noFilt_for_annotation_TE-details.bed
paste all_differential_peaks_noFilt_for_annotation.bed ../../../../rose-docker/all_bed_super-elements/${peakset}_merged_all.bed | awk -F'\t' -v OFS="\t" '{print $1,$2,$3,$4,$5,$7,$8,$9,$10,$14,$14,$14,$15}' > all_differential_peaks_noFilt_for_annotation_TE-details.bed
bedtools closest -a all_differential_peaks_noFilt_for_annotation_TE-details.bed -b /share/lab_avram/HPC_Cluster/user/tomas/lib/coordinates_TSS0bp_transcriptID_refgenie_ensembl_mm10.bed -d -t all -k 1000 | awk -F'\t' -v OFS="\t" '$21<50000{print $17,$20,$21,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | awk '{gsub(",","",$(NF))}1' OFS='\t' |  awk -F'\t' -v OFS="\t" '{print $0, "1"}' > $peaks1
bedtools closest -a all_differential_peaks_noFilt_for_annotation_TE-details.bed -b /share/lab_avram/HPC_Cluster/user/tomas/lib/coordinates_transcriptID_refgenie_ensembl_mm10.bed -d -t all -k 1000 | awk -F'\t' -v OFS="\t" '$21<50000{print $17,$20,$21,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | awk '{gsub(",","",$(NF))}1' OFS='\t' |  awk -F'\t' -v OFS="\t" '{print $0, "1"}' > $peaks2

##### next SCRNA rerun the upper commands for enhnacers (not stitched SE)
awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $rna_dataset $peaks1 | awk -F'\t' -v OFS="\t" '!($18=="#NA") {print $0}' > 01b_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt
awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $rna_dataset $peaks2 | awk -F'\t' -v OFS="\t" '!($18=="#NA") {print $0}' > 01b_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt
# now intersect the merged all-diff k27ac peaks with the SEs...first put the coordinates in front
awk -F'\t' -v OFS="\t" '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' 01b_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt| sortBed | uniq > 02b_coord-first_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt
# now intersect them
bedtools intersect -a 02b_coord-first_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt -b all-diff_h3k27ac_logFC02_Pval03__Pval005.bed -wo  > 03b_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp_with_diff_h3k27ac_logFC02_Pval03__Pval005_peaks.txt
# now intersect the merged all-diff k27ac peaks with the SEs...first put the coordinates in front
awk -F'\t' -v OFS="\t" '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' 01b_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt | sortBed | uniq > 02b_coord-first_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt
# now intersect them
bedtools intersect -a 02b_coord-first_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt -b all-diff_h3k27ac_logFC02_Pval03__Pval005.bed -wo  > 03b_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall_with_diff_h3k27ac_logFC02_Pval03__Pval005_peaks.txt









# for SE (super-enhancers)
# make adjustments in filtering and take all stitched peaks not just defined SE..and rerun the bedtools coverage and DESEQ pipeline
cd /share/lab_avram/HPC_Cluster/user/tomas/rose-docker/output_wt_TE-with27acBAM
awk -F'\t' -v OFS="\t" 'NR>1{print $2, $3, $4, $1, $2":"$3"-"$4,$7,$8,$9,$11,"h3k4me1_wt"NR}' *AllStitched_REGION_TO_GENE.txt |  sort -k1,1 -k2,2n | uniq > output_wt_TE-with27acBAM_nonFilt.bed
cd /share/lab_avram/HPC_Cluster/user/tomas/rose-docker/output_ko_TE-with27acBAM
awk -F':|-|\t' -v OFS="\t" 'NR>1{print $2, $3, $4, $1, $2":"$3"-"$4,$7,$8,$9,$11,"h3k4me1_ko"NR}' *AllStitched_REGION_TO_GENE.txt |  sort -k1,1 -k2,2n | uniq > output_ko_TE-with27acBAM_nonFilt.bed
cd /share/lab_avram/HPC_Cluster/user/tomas/rose-docker/all_bed_super-elements
cat output_wt_TE-with27acBAM_nonFilt.bed output_ko_TE-with27acBAM_nonFilt.bed | sortBed | mergeBed -c 4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse | sortBed  | uniq > SE-nonFilt_k4me1_with_k27acBAM_merged_all.bed

cd /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10
sample1="CRK27AM8TW60"
sample2="CRK27AM8TW91"
sample3="M8CTK47A"
sample4="M8CTK48A"
myArray=(${sample1} ${sample2} ${sample3} ${sample4})
peakset="SE-nonFilt_k4me1_with_k27acBAM"
for sample in ${myArray[@]}; do
  bedtools coverage -a /share/lab_avram/HPC_Cluster/user/tomas/rose-docker/all_bed_super-elements/${peakset}_merged_all.bed -b ${sample}/${sample}.bam | uniq > ${sample}_${peakset}_merged_peaks_coverage.bed; 
done
# extract raw reads and calculate differential signal at SE 
# change line 170 to column $11 in ~/bin/cut_and_run_pipeline/02_CUTRUN_DESeq2_for_seacr.sh
sbatch --chdir=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10 --export=s=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10/samples.csv,f=_SE-nonFilt_k4me1_with_k27acBAM_merged_peaks_coverage.bed ~/bin/CUTRUN_DESeq2_for_seacr.sh
mv DESeq2_seacr DESeq2_superenhancer-noFilt_01-16-2025
cd DESeq2_superenhancer-noFilt_01-16-2025
paste all_differential_peaks_noFilt_for_annotation.bed ../../../../rose-docker/all_bed_super-elements/${peakset}_merged_all.bed | awk -F'\t' -v OFS="\t" '{print $1,$2,$3,$4,$5,$7,$8,$9,$10,$14,$15,$20,$19}' > all_differential_peaks_noFilt_for_annotation_SE-details.bed
bedtools closest -a all_differential_peaks_noFilt_for_annotation_SE-details.bed -b /share/lab_avram/HPC_Cluster/user/tomas/lib/coordinates_TSS0bp_transcriptID_refgenie_ensembl_mm10.bed -d -t all -k 1000 | awk -F'\t' -v OFS="\t" '$21<50000{print $17,$20,$21,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | awk '{gsub(",","",$(NF))}1' OFS='\t' |  awk -F'\t' 'BEGIN{OFS="\t"}{if($(NF) ~ /1/){$17="1"}else{$17="0"};print $0}' > $peaks1
bedtools closest -a all_differential_peaks_noFilt_for_annotation_SE-details.bed -b /share/lab_avram/HPC_Cluster/user/tomas/lib/coordinates_transcriptID_refgenie_ensembl_mm10.bed -d -t all -k 1000 | awk -F'\t' -v OFS="\t" '$21<50000{print $17,$20,$21,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | awk '{gsub(",","",$(NF))}1' OFS='\t' |  awk -F'\t' 'BEGIN{OFS="\t"}{if($(NF) ~ /1/){$17="1"}else{$17="0"};print $0}' > $peaks2

# scRNA with K27ac SE
# use the vlookup function  and further filter out elements without any expression data
awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $rna_dataset $peaks1 | awk -F'\t' -v OFS="\t" '!($18=="#NA") {print $0}' > 01b_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt
awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $rna_dataset $peaks2 | awk -F'\t' -v OFS="\t" '!($18=="#NA") {print $0}' > 01b_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt
# now intersect the merged all-diff k27ac peaks with the SEs...first put the coordinates in front
awk -F'\t' -v OFS="\t" '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' 01b_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt | sortBed | uniq > 02b_coord-first_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt
# now intersect them
bedtools intersect -a 02b_coord-first_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt -b all-diff_h3k27ac_logFC02_Pval03__Pval005.bed -wo  > 03b_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp_with_diff_h3k27ac_logFC02_Pval03__Pval005_peaks.txt
# now intersect the merged all-diff k27ac peaks with the SEs...first put the coordinates in front
awk -F'\t' -v OFS="\t" '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' 01b_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt| sortBed | uniq > 02b_coord-first_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt
# now intersect them
bedtools intersect -a 02b_coord-first_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt -b all-diff_h3k27ac_logFC02_Pval03__Pval005.bed -wo  > 03b_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall_with_diff_h3k27ac_logFC02_Pval03__Pval005_peaks.txt











# poised-SE = PE
cd /share/lab_avram/HPC_Cluster/user/tomas/rose-docker/output_wt_TE-withK27me3BAM_PE
awk -F'\t' -v OFS="\t" 'NR>1{print $2, $3, $4, $1, $2":"$3"-"$4,$7,$8,$9,$11,"h3k4me1_wt"NR}' *AllStitched_REGION_TO_GENE.txt |  sort -k1,1 -k2,2n | uniq > output_wt_PE-with27me3BAM_nonFilt.bed
cd /share/lab_avram/HPC_Cluster/user/tomas/rose-docker/output_ko_TE-withK27me3BAM_PE
awk -F':|-|\t' -v OFS="\t" 'NR>1{print $2, $3, $4, $1, $2":"$3"-"$4,$7,$8,$9,$11,"h3k4me1_ko"NR}' *AllStitched_REGION_TO_GENE.txt |  sort -k1,1 -k2,2n | uniq > output_ko_PE-with27me3BAM_nonFilt.bed
cd /share/lab_avram/HPC_Cluster/user/tomas/rose-docker/all_bed_super-elements
cat output_wt_PE-with27me3BAM_nonFilt.bed output_ko_PE-with27me3BAM_nonFilt.bed  | sortBed | mergeBed -c 4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse | sortBed  | uniq > PE-nonFilt_k4me1_with_k27me3BAM_merged_all.bed

cd /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27me3/cd8_ovar_tumor_mm10
sample1="M8CTK38M"
sample2="M8CTK61M"
sample3="M8CTW17M"
sample4="M8CTW60M"
myArray=(${sample1} ${sample2} ${sample3} ${sample4})
peakset="PE-nonFilt_k4me1_with_k27me3BAM"
for sample in ${myArray[@]}; do
  bedtools coverage -a /share/lab_avram/HPC_Cluster/user/tomas/rose-docker/all_bed_super-elements/${peakset}_merged_all.bed -b ${sample}/${sample}.bam | uniq > ${sample}_${peakset}_merged_peaks_coverage.bed
done
# change line 170 to column $11 in ~/bin/cut_and_run_pipeline/02_CUTRUN_DESeq2_for_seacr.sh
sbatch --chdir=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27me3/cd8_ovar_tumor_mm10 --export=s=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27me3/cd8_ovar_tumor_mm10/samples.csv,f=_PE-nonFilt_k4me1_with_k27me3BAM_merged_peaks_coverage.bed ~/bin/CUTRUN_DESeq2_for_seacr.sh
mv DESeq2_seacr DESeq2_poisedenhancer-noFilt_01-16-2025
cd DESeq2_poisedenhancer-noFilt_01-16-2025
paste all_differential_peaks_noFilt_for_annotation.bed ../../../../rose-docker/all_bed_super-elements/${peakset}_merged_all.bed | awk -F'\t' -v OFS="\t" '{print $1,$2,$3,$4,$5,$7,$8,$9,$10,$14,$15,$20,$19}' > all_differential_peaks_noFilt_for_annotation_PE-details.bed
bedtools closest -a all_differential_peaks_noFilt_for_annotation_PE-details.bed -b /share/lab_avram/HPC_Cluster/user/tomas/lib/coordinates_TSS0bp_transcriptID_refgenie_ensembl_mm10.bed -d -t all -k 1000 | awk -F'\t' -v OFS="\t" '$21<50000{print $17,$20,$21,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | awk '{gsub(",","",$(NF))}1' OFS='\t' |  awk -F'\t' 'BEGIN{OFS="\t"}{if($(NF) ~ /1/){$17="1"}else{$17="0"};print $0}' > $peaks1
bedtools closest -a all_differential_peaks_noFilt_for_annotation_PE-details.bed -b /share/lab_avram/HPC_Cluster/user/tomas/lib/coordinates_transcriptID_refgenie_ensembl_mm10.bed -d -t all -k 1000 | awk -F'\t' -v OFS="\t" '$21<50000{print $17,$20,$21,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | awk '{gsub(",","",$(NF))}1' OFS='\t' |  awk -F'\t' 'BEGIN{OFS="\t"}{if($(NF) ~ /1/){$17="1"}else{$17="0"};print $0}' > $peaks2

awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $rna_dataset $peaks1 | awk -F'\t' -v OFS="\t" '!($18=="#NA") {print $0}' > 01b_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt
awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $rna_dataset $peaks2 | awk -F'\t' -v OFS="\t" '!($18=="#NA") {print $0}' > 01b_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt

awk -F'\t' -v OFS="\t" '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' 01b_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt | sortBed | uniq > 02b_coord-first_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt
bedtools intersect -a 02b_coord-first_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt -b all-diff_h3k27ac_logFC02_Pval03__Pval005.bed -wo  > 03b_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp_with_diff_h3k27ac_logFC02_Pval03__Pval005_peaks.txt
bedtools intersect -a 02b_coord-first_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.txt -b all-diff_h3k27me3_logFC02_Pval03__Pval005.bed -wo  > 04b_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp_with_diff_h3k27me3_logFC02_Pval03__Pval005_peaks.txt

awk -F'\t' -v OFS="\t" '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' 01b_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt | sortBed | uniq > 02b_coord-first_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt
bedtools intersect -a 02b_coord-first_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt -b all-diff_h3k27ac_logFC02_Pval03__Pval005.bed -wo  > 03b_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall_with_diff_h3k27ac_logFC02_Pval03__Pval005_peaks.txt
bedtools intersect -a 02b_coord-first_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall.txt -b all-diff_h3k27me3_logFC02_Pval03__Pval005.bed -wo  > 04b_final_PE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSSall_with_diff_h3k27me3_logFC02_Pval03__Pval005_peaks.txt