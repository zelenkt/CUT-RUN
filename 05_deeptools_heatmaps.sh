#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=30GB
#SBATCH --qos=medium
#SBATCH -t 24:00:00
#SBATCH --job-name=Deeptools_heatmaps
#SBATCH --output=Deeptools_heatmaps.%j.out
#SBATCH --error=Deeptools_heatmaps.%j.err

##### LAST UPDATED 02-26-2023 Tomas Zelenka  #####

# to run, uncomment the section of interest and run the following command:
# sbatch /share/lab_avram/HPC_Cluster/user/tomas/bin/cut_and_run_pipeline/05_deeptools_heatmaps.sh



module load deepTools/3.5.0-foss-2021a
ml BEDTools/2.30.0-GCC-11.2.0
cd /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/deeptools_plots


######## FIRST PREPARE INPUT FILES
# # these are scRNAseq data differentially expressed genes from  /Users/4476125/Documents/BERLIN_Shaw_lab_2023/BERLIN/2-Single_Cell_RNAseq_Pipeline/Output/ovarT_integrated_noRefUsed/reclustered/markers/allKOvsWT_ovarT_assayRNA.txt
# # these were generated as part of the scRNAseq pipeline and annotated also within R as described in 10_berlin_v5.0.1_combined.R
# # 01-31-2025_allKOvsWT_ANNOTATED_melanT_assayRNA.txt
# # 01-31-2025_allKOvsWT_ANNOTATED_ovarT_assayRNA.txt

# # they were further modified by this to add coordinates from lib/geneID_coordinates_refgenie_ensembl_mm10.txt
# cat 01_01-31-2025_allKOvsWT_ANNOTATED_ovarT_assayRNA.txt | tr -d '"'| awk  -F'\t' -v OFS="\t" 'NR>1{print $2,$1,$6-$7,$5,$6,$7,$4,$8}' > tmp1315; awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' /share/lab_avram/HPC_Cluster/user/tomas/lib/geneNames_refgenie_ensembl_mm10_description_longestIsoLength.gtf tmp1315 | awk -F'\t' -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$11}' | awk -F'\t' -v OFS="\t" '{for(i=1; i<=9; i++) if ($i=="") $i="NA"; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > 02_ovarT_scRNAseq_adapted_markers.seurat_obj_reclust2_01-31-2025.rna_allClustTogether.txt; rm tmp1315
# cat 01_01-31-2025_allKOvsWT_ANNOTATED_melanT_assayRNA.txt | tr -d '"'| awk  -F'\t' -v OFS="\t" 'NR>1{print $2,$1,$6-$7,$5,$6,$7,$4,$8}' > tmp1315; awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' /share/lab_avram/HPC_Cluster/user/tomas/lib/geneNames_refgenie_ensembl_mm10_description_longestIsoLength.gtf tmp1315 | awk -F'\t' -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$11}' | awk -F'\t' -v OFS="\t" '{for(i=1; i<=9; i++) if ($i=="") $i="NA"; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > 02_melanT_scRNAseq_adapted_markers.seurat_obj_reclust2_01-31-2025.rna_allClustTogether.txt; rm tmp1315
# # new columns are geneID    gene  pct.1-pct.2  avg_log2FC  pct.1   pct.2  pval   adjPval    description

# cat 02_ovarT_scRNAseq_adapted_markers.seurat_obj_reclust2_01-31-2025.rna_allClustTogether.txt > tmp1315; awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' /share/lab_avram/HPC_Cluster/user/tomas/lib/geneID_coordinates_refgenie_ensembl_mm10.txt tmp1315 | awk -F'\t' -v OFS="\t" '{print $11,$12,$13,$1,$14,$2,$3,$4,$5,$6,$7,$8,$9}' | awk -F'\t' -v OFS="\t" '{for(i=1; i<=13; i++) if ($i=="") $i="NA"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | awk '!($1 ~ /^chrUn/ || $1 ~ /random$/ || $1~/^chrMT/ || $1 ~ /^NA/ || $1 ~ /^chrG/ || $1 ~ /^chrJ/ || $1 ~ /^chrY/) {print $0}' | sort -k1,1 -k2,2n | uniq > 03_ovarT_coord_scRNAseq_adapted_markers.seurat_obj_reclust2_01-31-2025.rna_allClustTogether.txt; rm tmp1315
# cat 02_melanT_scRNAseq_adapted_markers.seurat_obj_reclust2_01-31-2025.rna_allClustTogether.txt > tmp1315; awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' /share/lab_avram/HPC_Cluster/user/tomas/lib/geneID_coordinates_refgenie_ensembl_mm10.txt tmp1315 | awk -F'\t' -v OFS="\t" '{print $11,$12,$13,$1,$14,$2,$3,$4,$5,$6,$7,$8,$9}' | awk -F'\t' -v OFS="\t" '{for(i=1; i<=13; i++) if ($i=="") $i="NA"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | awk '!($1 ~ /^chrUn/ || $1 ~ /random$/ || $1~/^chrMT/ || $1 ~ /^NA/ || $1 ~ /^chrG/ || $1 ~ /^chrJ/ || $1 ~ /^chrY/) {print $0}' | sort -k1,1 -k2,2n | uniq > 03_melanT_coord_scRNAseq_adapted_markers.seurat_obj_reclust2_01-31-2025.rna_allClustTogether.txt; rm tmp1315
# # new columns are chr start end geneID strand geneSymbol  pct.1-pct.2  avg_log2FC  pct.1   pct.2  pval   adjPval    description

# # filter up or downregulated genes
# # $8 = log2FC
# # $11 = pval
# awk '($8 > 0 && $11<0.05) {print $0}' 03_melanT_coord_scRNAseq_adapted_markers.seurat_obj_reclust2_01-31-2025.rna_allClustTogether.txt > 04_up_melanT_scRNAseq_Pval005.bed
# awk '($8 < 0 && $11<0.05) {print $0}' 03_melanT_coord_scRNAseq_adapted_markers.seurat_obj_reclust2_01-31-2025.rna_allClustTogether.txt > 04_down_melanT_scRNAseq_Pval005.bed
# awk '($8 > 0 && $11<0.05) {print $0}' 03_ovarT_coord_scRNAseq_adapted_markers.seurat_obj_reclust2_01-31-2025.rna_allClustTogether.txt > 04_up_ovarT_scRNAseq_Pval005.bed
# awk '($8 < 0 && $11<0.05) {print $0}' 03_ovarT_coord_scRNAseq_adapted_markers.seurat_obj_reclust2_01-31-2025.rna_allClustTogether.txt > 04_down_ovarT_scRNAseq_Pval005.bed

# # FROM THIS POINT ON FOCUS ONLY ON OVARIAN DATA
# # create promoter file with -1000 bp to TSS to TSS+200bp
# cd ~/lib/
# awk -F'\t' -v OFS="\t" '{print $0, $1,$2,$3}' coordinates_geneID_refgenie_ensembl_mm10.bed | awk 'BEGIN {OFS="\t"} 
#     {
#         if ($6 == "+") {
#             start = ($2 - 1000 > 0) ? $2 - 1000 : 0;  # Ensure start is not negative (0-based)
#             end = $2 + 200;
#         } else if ($6 == "-") {
#             start = ($3 - 200 > 0) ? $3 - 200 : 0;  # Ensure start is not negative
#             end = $3 + 1000;
#         }
#         print $1, start, end, $4, $5, $6, $7, $8, $9, $10;
#     }' > coordinates_promoters1000bpTSS200bp_geneID_refgenie_ensembl_mm10_with_origCoord.bed
# cd /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/deeptools_plots





path_bw="/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/pygenometracks/input/mouse_mm10"
promoters="/share/lab_avram/HPC_Cluster/user/tomas/lib/coordinates_promoters1000bpTSS200bp_geneID_refgenie_ensembl_mm10_with_origCoord.bed"
bcl11b_peaks="/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/pygenometracks/input/mouse_mm10/wt_bcl11b_cat_merge_0005_non_stringent.stringent.bed"

# Heatmap for H3K4me3 CUT&RUN signal visualized for TSS (±3 kb) of genes with Bcl11b binding at their promoter (TSS, −1000 to +200 bp).
# this gives the promoters intersecting with Bcl11b, but I also had to move the original coordinates of these genes to the first position with the awk command
# bedtools intersect -a $promoters -b $bcl11b_peaks -wa | awk -F'\t' -v OFS="\t" '{print $8, $9, $10, $4, $5, $6, $7, $1, $2, $3}' | sortBed | uniq  > 05_inter_bcl11b-bound_promoters_all_genes_origCoord.bed
# bedtools intersect -a $promoters -b $bcl11b_peaks -v -wa | awk -F'\t' -v OFS="\t" '{print $8, $9, $10, $4, $5, $6, $7, $1, $2, $3}' | sortBed | uniq  > 05_inter_bcl11b-unbound_promoters_all_genes_origCoord.bed



# prepare TE, SE and PE files to compare them for CUT&RUN signal
# these are all TE elements, i.e. non-filtered H3K4me1 peaks
TE=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/rose-docker/all_bed_super-elements/h3k4me1_merged_all.bed

# these are the filtered SE (i.e. they have H3K4me1 signal and strong H3K27ac). To make it comparable between WT and KO, I have to combine WT and KO set (used for deeptools):
# output_wt_TE-with27acBAM.bed
# output_ko_TE-with27acBAM.bed
SE=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/rose-docker/all_bed_super-elements/SE-Filt_k4me1_with_k27acBAM_merged_all.bed

# these are the filtered PE (i.e. they have H3K4me1 signal and strong H3K27me3). To make it comparable between WT and KO, I have to combine WT and KO set (used for deeptools):
# output_wt_TE-withK27me3BAM_PE.bed
# output_ko_TE-withK27me3BAM_PE.bed
PE=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/rose-docker/all_bed_super-elements/PE-Filt_k4me1_with_k27me3BAM_merged_all.bed

# these are all the stitched TE which includes SE but also regions that were not classified as SE
# output_wt_TE-with27acBAM_nonFilt.bed
# output_ko_TE-with27acBAM_nonFilt.bed
# these are already combined into a single file here
stitchedTE=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/rose-docker/all_bed_super-elements/SE-nonFilt_k4me1_with_k27acBAM_merged_all.bed 

# these are all the stitched TE which includes PE but also regions that were not classified as PE, coordinates-wise, this should be the same as the stitched-TE file which includes SE but also regions that were not classified as SE, i.e. SE-nonFilt_k4me1_with_k27acBAM_merged_all.bed
# output_wt_PE-with27me3BAM_nonFilt.bed
# output_ko_PE-with27me3BAM_nonFilt.bed
# these are already combined into a single file here
# /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/rose-docker/all_bed_super-elements/PE-nonFilt_k4me1_with_k27me3BAM_merged_all.bed

# intersect the regulatory elements with Bcl11b peaks
# bedtools intersect -a $TE -b $bcl11b_peaks -wa | sortBed | uniq  > 06_inter_TE_bcl11b_bound.bed
# bedtools intersect -a $SE -b $bcl11b_peaks -wa | sortBed | uniq  > 07_inter_SE_bcl11b_bound.bed
# bedtools intersect -a $PE -b $bcl11b_peaks -wa | sortBed | uniq  > 08_inter_PE_bcl11b_bound.bed
# bedtools intersect -a $stitchedTE -b $bcl11b_peaks -wa | sortBed | uniq  > 09_inter_stitchedTE_bcl11b_bound.bed

# bedtools intersect -a $TE -b $bcl11b_peaks -v -wa | sortBed | uniq  > 06_inter_TE_bcl11b_unbound.bed
# bedtools intersect -a $SE -b $bcl11b_peaks -v -wa | sortBed | uniq  > 07_inter_SE_bcl11b_unbound.bed
# bedtools intersect -a $PE -b $bcl11b_peaks -v -wa | sortBed | uniq  > 08_inter_PE_bcl11b_unbound.bed
# bedtools intersect -a $stitchedTE -b $bcl11b_peaks -v -wa | sortBed | uniq  > 09_inter_stitchedTE_bcl11b_unbound.bed



# consider checking signal on a list of genes - on stemness genes, effector genes, NK genes..would it make sense?

# list of bw files to be checked that are present in $path_bw
# atac_mouse_CD8tumor_KO_mC8ATK_smooth.bw
# atac_mouse_CD8tumor_WT_mC8ATW_smooth.bw
# bcl11b_wt_mean_rep_mm10_treat_pileup.bw
# h3k27ac_ko_mean_rep_mm10_treat_pileup.bw
# h3k27ac_wt_mean_rep_mm10_treat_pileup.bw
# h3k27me3_ko_mean_rep_mm10_treat_pileup.bw
# h3k27me3_wt_mean_rep_mm10_treat_pileup.bw
# h3k4me1_ko_mean_rep_mm10_treat_pileup.bw
# h3k4me1_wt_mean_rep_mm10_treat_pileup.bw


################# MAIN PART TO GENERATE HEATMAPS USED FOR PUBLICATION

# computeMatrix scale-regions --regionsFileName 04_up_ovarT_scRNAseq_Pval005.bed 04_down_ovarT_scRNAseq_Pval005.bed --scoreFileName "$path_bw"/h3k27ac_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27ac_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 04_UDEG-DDEG_vs_K27ac-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 04_UDEG-DDEG_vs_K27ac-signal_3kbp.gz --colorMap Reds --outFileName 04_UDEG-DDEG_vs_K27ac-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 04_up_ovarT_scRNAseq_Pval005.bed 04_down_ovarT_scRNAseq_Pval005.bed --scoreFileName "$path_bw"/h3k27me3_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27me3_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 04_UDEG-DDEG_vs_K27me3-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 04_UDEG-DDEG_vs_K27me3-signal_3kbp.gz --colorMap Reds --outFileName 04_UDEG-DDEG_vs_K27me3-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 04_up_ovarT_scRNAseq_Pval005.bed 04_down_ovarT_scRNAseq_Pval005.bed --scoreFileName "$path_bw"/atac_mouse_CD8tumor_WT_mC8ATW_smooth.bw "$path_bw"/atac_mouse_CD8tumor_KO_mC8ATK_smooth.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 04_UDEG-DDEG_vs_ATAC-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 04_UDEG-DDEG_vs_ATAC-signal_3kbp.gz --colorMap Reds --outFileName 04_UDEG-DDEG_vs_ATAC-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 04_up_ovarT_scRNAseq_Pval005.bed 04_down_ovarT_scRNAseq_Pval005.bed --scoreFileName "$path_bw"/h3k4me1_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k4me1_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 04_UDEG-DDEG_vs_K4me1-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 04_UDEG-DDEG_vs_K4me1-signal_3kbp.gz --colorMap Reds --outFileName 04_UDEG-DDEG_vs_K4me1-signal_3kbp.svg




# computeMatrix reference-point --regionsFileName 05_inter_bcl11b-bound_promoters_all_genes_origCoord.bed 05_inter_bcl11b-unbound_promoters_all_genes_origCoord.bed --scoreFileName "$path_bw"/h3k27me3_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27me3_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 05_inter_bcl11b-promoters-all_vs_K27me3-signal_refPointTSS.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 05_inter_bcl11b-promoters-all_vs_K27me3-signal_refPointTSS.gz --refPointLabel TSS --colorMap Reds --outFileName 05_inter_bcl11b-promoters-all_vs_K27me3-signal_refPointTSS.svg

# computeMatrix reference-point --regionsFileName 05_inter_bcl11b-bound_promoters_all_genes_origCoord.bed 05_inter_bcl11b-unbound_promoters_all_genes_origCoord.bed --scoreFileName "$path_bw"/h3k27ac_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27ac_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 05_inter_bcl11b-promoters-all_vs_K27ac-signal_refPointTSS.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 05_inter_bcl11b-promoters-all_vs_K27ac-signal_refPointTSS.gz --refPointLabel TSS --colorMap Reds --outFileName 05_inter_bcl11b-promoters-all_vs_K27ac-signal_refPointTSS.svg

# computeMatrix reference-point --regionsFileName 05_inter_bcl11b-bound_promoters_all_genes_origCoord.bed 05_inter_bcl11b-unbound_promoters_all_genes_origCoord.bed --scoreFileName "$path_bw"/h3k4me1_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k4me1_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 05_inter_bcl11b-promoters-all_vs_K4me1-signal_refPointTSS.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 05_inter_bcl11b-promoters-all_vs_K4me1-signal_refPointTSS.gz --refPointLabel TSS --colorMap Reds --outFileName 05_inter_bcl11b-promoters-all_vs_K4me1-signal_refPointTSS.svg

# computeMatrix reference-point --regionsFileName 05_inter_bcl11b-bound_promoters_all_genes_origCoord.bed 05_inter_bcl11b-unbound_promoters_all_genes_origCoord.bed --scoreFileName "$path_bw"/atac_mouse_CD8tumor_WT_mC8ATW_smooth.bw "$path_bw"/atac_mouse_CD8tumor_KO_mC8ATK_smooth.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 05_inter_bcl11b-promoters-all_vs_ATAC-signal_refPointTSS.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 05_inter_bcl11b-promoters-all_vs_ATAC-signal_refPointTSS.gz --refPointLabel TSS --colorMap Reds --outFileName 05_inter_bcl11b-promoters-all_vs_ATAC-signal_refPointTSS.svg





# computeMatrix scale-regions --regionsFileName 06_inter_TE_bcl11b_bound.bed 06_inter_TE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k27ac_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27ac_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 06_TE_bcl11b_vs_K27ac-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 06_TE_bcl11b_vs_K27ac-signal_3kbp.gz --colorMap Reds --outFileName 06_TE_bcl11b_vs_K27ac-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 06_inter_TE_bcl11b_bound.bed 06_inter_TE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k27me3_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27me3_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 06_TE_bcl11b_vs_K27me3-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 06_TE_bcl11b_vs_K27me3-signal_3kbp.gz --colorMap Reds --outFileName 06_TE_bcl11b_vs_K27me3-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 06_inter_TE_bcl11b_bound.bed 06_inter_TE_bcl11b_unbound.bed --scoreFileName "$path_bw"/atac_mouse_CD8tumor_WT_mC8ATW_smooth.bw "$path_bw"/atac_mouse_CD8tumor_KO_mC8ATK_smooth.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 06_TE_bcl11b_vs_ATAC-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 06_TE_bcl11b_vs_ATAC-signal_3kbp.gz --colorMap Reds --outFileName 06_TE_bcl11b_vs_ATAC-signal_3kbp.svg


# computeMatrix scale-regions --regionsFileName 06_inter_TE_bcl11b_bound.bed 06_inter_TE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k4me1_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k4me1_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 06_TE_bcl11b_vs_K4me1-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 06_TE_bcl11b_vs_K4me1-signal_3kbp.gz --colorMap Reds --outFileName 06_TE_bcl11b_vs_K4me1-signal_3kbp.svg







# computeMatrix scale-regions --regionsFileName 07_inter_SE_bcl11b_bound.bed 07_inter_SE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k27ac_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27ac_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 07_SE_bcl11b_vs_K27ac-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 07_SE_bcl11b_vs_K27ac-signal_3kbp.gz --colorMap Reds --outFileName 07_SE_bcl11b_vs_K27ac-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 07_inter_SE_bcl11b_bound.bed 07_inter_SE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k27me3_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27me3_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 07_SE_bcl11b_vs_K27me3-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 07_SE_bcl11b_vs_K27me3-signal_3kbp.gz --colorMap Reds --outFileName 07_SE_bcl11b_vs_K27me3-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 07_inter_SE_bcl11b_bound.bed 07_inter_SE_bcl11b_unbound.bed --scoreFileName "$path_bw"/atac_mouse_CD8tumor_WT_mC8ATW_smooth.bw "$path_bw"/atac_mouse_CD8tumor_KO_mC8ATK_smooth.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 07_SE_bcl11b_vs_ATAC-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 07_SE_bcl11b_vs_ATAC-signal_3kbp.gz --colorMap Reds --outFileName 07_SE_bcl11b_vs_ATAC-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 07_inter_SE_bcl11b_bound.bed 07_inter_SE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k4me1_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k4me1_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 07_SE_bcl11b_vs_K4me1-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 07_SE_bcl11b_vs_K4me1-signal_3kbp.gz --colorMap Reds --outFileName 07_SE_bcl11b_vs_K4me1-signal_3kbp.svg






# computeMatrix scale-regions --regionsFileName 08_inter_PE_bcl11b_bound.bed 08_inter_PE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k27ac_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27ac_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 08_PE_bcl11b_vs_K27ac-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 08_PE_bcl11b_vs_K27ac-signal_3kbp.gz --colorMap Reds --outFileName 08_PE_bcl11b_vs_K27ac-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 08_inter_PE_bcl11b_bound.bed 08_inter_PE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k27me3_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27me3_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 08_PE_bcl11b_vs_K27me3-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 08_PE_bcl11b_vs_K27me3-signal_3kbp.gz --colorMap Reds --outFileName 08_PE_bcl11b_vs_K27me3-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 08_inter_PE_bcl11b_bound.bed 08_inter_PE_bcl11b_unbound.bed --scoreFileName "$path_bw"/atac_mouse_CD8tumor_WT_mC8ATW_smooth.bw "$path_bw"/atac_mouse_CD8tumor_KO_mC8ATK_smooth.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 08_PE_bcl11b_vs_ATAC-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 08_PE_bcl11b_vs_ATAC-signal_3kbp.gz --colorMap Reds --outFileName 08_PE_bcl11b_vs_ATAC-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 08_inter_PE_bcl11b_bound.bed 08_inter_PE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k4me1_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k4me1_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 08_PE_bcl11b_vs_K4me1-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 08_PE_bcl11b_vs_K4me1-signal_3kbp.gz --colorMap Reds --outFileName 08_PE_bcl11b_vs_K4me1-signal_3kbp.svg







# computeMatrix scale-regions --regionsFileName 09_inter_stitchedTE_bcl11b_bound.bed 09_inter_stitchedTE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k27ac_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27ac_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 09_stitchedTE_bcl11b_vs_K27ac-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 09_stitchedTE_bcl11b_vs_K27ac-signal_3kbp.gz --colorMap Reds --outFileName 09_stitchedTE_bcl11b_vs_K27ac-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 09_inter_stitchedTE_bcl11b_bound.bed 09_inter_stitchedTE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k27me3_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27me3_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 09_stitchedTE_bcl11b_vs_K27me3-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 09_stitchedTE_bcl11b_vs_K27me3-signal_3kbp.gz --colorMap Reds --outFileName 09_stitchedTE_bcl11b_vs_K27me3-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 09_inter_stitchedTE_bcl11b_bound.bed 09_inter_stitchedTE_bcl11b_unbound.bed --scoreFileName "$path_bw"/atac_mouse_CD8tumor_WT_mC8ATW_smooth.bw "$path_bw"/atac_mouse_CD8tumor_KO_mC8ATK_smooth.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 09_stitchedTE_bcl11b_vs_ATAC-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 09_stitchedTE_bcl11b_vs_ATAC-signal_3kbp.gz --colorMap Reds --outFileName 09_stitchedTE_bcl11b_vs_ATAC-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName 09_inter_stitchedTE_bcl11b_bound.bed 09_inter_stitchedTE_bcl11b_unbound.bed --scoreFileName "$path_bw"/h3k4me1_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k4me1_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 09_stitchedTE_bcl11b_vs_K4me1-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 09_stitchedTE_bcl11b_vs_K4me1-signal_3kbp.gz --colorMap Reds --outFileName 09_stitchedTE_bcl11b_vs_K4me1-signal_3kbp.svg



# computeMatrix scale-regions --regionsFileName $bcl11b_peaks --scoreFileName "$path_bw"/h3k27ac_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27ac_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 10_bcl11b-peaks_vs_K27ac-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 10_bcl11b-peaks_vs_K27ac-signal_3kbp.gz --colorMap Reds --outFileName 10_bcl11b-peaks_vs_K27ac-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName $bcl11b_peaks --scoreFileName "$path_bw"/h3k27me3_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k27me3_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 10_bcl11b-peaks_vs_K27me3-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 10_bcl11b-peaks_vs_K27me3-signal_3kbp.gz --colorMap Reds --outFileName 10_bcl11b-peaks_vs_K27me3-signal_3kbp.svg
 
# computeMatrix scale-regions --regionsFileName $bcl11b_peaks --scoreFileName "$path_bw"/atac_mouse_CD8tumor_WT_mC8ATW_smooth.bw "$path_bw"/atac_mouse_CD8tumor_KO_mC8ATK_smooth.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 10_bcl11b-peaks_vs_ATAC-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 10_bcl11b-peaks_vs_ATAC-signal_3kbp.gz --colorMap Reds --outFileName 10_bcl11b-peaks_vs_ATAC-signal_3kbp.svg

# computeMatrix scale-regions --regionsFileName $bcl11b_peaks --scoreFileName "$path_bw"/h3k4me1_wt_mean_rep_mm10_treat_pileup.bw "$path_bw"/h3k4me1_ko_mean_rep_mm10_treat_pileup.bw --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros --missingDataAsZero --outFileName 10_bcl11b-peaks_vs_K4me1-signal_3kbp.gz --numberOfProcessors max --blackListFileName /share/lab_avram/HPC_Cluster/annotations_genomes/alias/hg38/blacklist/default/hg38_blacklist.bed.gz
# plotHeatmap --matrixFile 10_bcl11b-peaks_vs_K4me1-signal_3kbp.gz --colorMap Reds --outFileName 10_bcl11b-peaks_vs_K4me1-signal_3kbp.svg













mv *.err logs/
mv *.out logs/

date
