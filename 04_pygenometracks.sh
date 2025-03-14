
# This is a summary of commands to be run to generate genomics tracks figures using pygenometrakcs software
# the figures and config files are saved in:
cd /share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/pygenometracks 

conda activate pygenometracks

# # if needed to modify the paths in .ini files, just run this in the ini_files_modify_pathways directory
# chmod +x replace_paths.sh 
#  ./replace_paths.sh

# RUN to create trakcs for HUMAN data
pyGenomeTracks --tracks BACH2.ini --region chr6:89,900,001-91,000,962 --outFileName BACH2.png
pyGenomeTracks --tracks BACH2.ini --region chr6:89,900,001-91,000,962 --outFileName BACH2.svg

pyGenomeTracks --tracks CCR7.ini --region chr17:40,549,342-40,597,718 --outFileName CCR7.png
pyGenomeTracks --tracks CCR7.ini --region chr17:40,549,342-40,597,718 --outFileName CCR7.svg

pyGenomeTracks --tracks CTLA4.ini --region chr2:203,850,090-203,990,817 --outFileName CTLA4.png    # some decreased in ATAC but not in both patients significant
pyGenomeTracks --tracks CTLA4.ini --region chr2:203,850,090-203,990,817 --outFileName CTLA4.svg

pyGenomeTracks --tracks FCER1G.ini --region chr1:161,193,349-161,265,677 --outFileName FCER1G.png  # no nice Bcl11b peaks
pyGenomeTracks --tracks FCER1G.ini --region chr1:161,193,349-161,265,677 --outFileName FCER1G.svg 

pyGenomeTracks --tracks IL2RA.ini --region chr10:5,987,437-6,298,410 --outFileName IL2RA.png
pyGenomeTracks --tracks IL2RA.ini --region chr10:5,987,437-6,298,410 --outFileName IL2RA.svg

pyGenomeTracks --tracks KLRK1.ini --region chr12:10,351,634-10,506,513 --outFileName KLRK1.png
pyGenomeTracks --tracks KLRK1.ini --region chr12:10,351,634-10,506,513 --outFileName KLRK1.svg

pyGenomeTracks --tracks NCR1.ini --region chr19:54,893,436-55,029,675 --outFileName NCR1.png
pyGenomeTracks --tracks NCR1.ini --region chr19:54,893,436-55,029,675 --outFileName NCR1.svg

pyGenomeTracks --tracks NCR3.ini --region chr6:31,570,765-31,620,000 --outFileName NCR3.png   # also increase ATAC in KO
pyGenomeTracks --tracks NCR3.ini --region chr6:31,570,765-31,620,000 --outFileName NCR3.svg

pyGenomeTracks --tracks PDCD1.ini --region chr2:241,840,000-241,906,301 --outFileName PDCD1.png
pyGenomeTracks --tracks PDCD1.ini --region chr2:241,840,000-241,906,301 --outFileName PDCD1.svg

pyGenomeTracks --tracks PRF1.ini --region chr10:70,584,809-70,662,552 --outFileName PRF1.png
pyGenomeTracks --tracks PRF1.ini --region chr10:70,584,809-70,662,552 --outFileName PRF1.svg

pyGenomeTracks --tracks XCL1.ini --region chr1:168,546,250-168,690,906 --outFileName XCL1.png    # no nice Bcl11b peaks
pyGenomeTracks --tracks XCL1.ini --region chr1:168,546,250-168,690,906 --outFileName XCL1.svg 




# RUN to create trakcs for MOUSE data
pyGenomeTracks --tracks Fcer1g.ini --region chr1:171,213,440-171,324,814 --outFileName Fcer1g.png
pyGenomeTracks --tracks Tcf7.ini --region chr11:52,209,926-52,669,437 --outFileName Tcf7.png
pyGenomeTracks --tracks Tyrobp.ini --region chr7:30,385,644-30,561,256 --outFileName Tyrobp.png
pyGenomeTracks --tracks Slamf6.ini --region chr1:171,902,126-172,138,148 --outFileName Slamf6.png
pyGenomeTracks --tracks Xcl1.ini --region chr1:164,914,243-165,150,265 --outFileName Xcl1.png
pyGenomeTracks --tracks Ccr7.ini --region chr11:99,084,436-99,379,252 --outFileName Ccr7.png
pyGenomeTracks --tracks NK.ini --region chr6:128,587,253-131,480,042 --outFileName NK.png
pyGenomeTracks --tracks Tox.ini --region chr4:6,642,850-7,531,941 --outFileName Tox.png
pyGenomeTracks --tracks Ctla4.ini --region chr1:60,720,001-61,192,045 --outFileName Ctla4.png
pyGenomeTracks --tracks Pdcd1.ini --region chr1:93993232-94,254,122 --outFileName Pdcd1.png
pyGenomeTracks --tracks Entpd1.ini --region chr19:40,535,109-41,128,512 --outFileName Entpd1.png
pyGenomeTracks --tracks Bach2.ini --region chr4:32,032,372-33,021,028 --outFileName Bach2.png

pyGenomeTracks --tracks Fcer1g.ini --region chr1:171,213,440-171,324,814 --outFileName Fcer1g.svg
pyGenomeTracks --tracks Tcf7.ini --region chr11:52,209,926-52,669,437 --outFileName Tcf7.svg
pyGenomeTracks --tracks Tyrobp.ini --region chr7:30,385,644-30,561,256 --outFileName Tyrobp.svg
pyGenomeTracks --tracks Slamf6.ini --region chr1:171,902,126-172,138,148 --outFileName Slamf6.svg
pyGenomeTracks --tracks Xcl1.ini --region chr1:164,914,243-165,150,265 --outFileName Xcl1.svg
pyGenomeTracks --tracks Ccr7.ini --region chr11:99,084,436-99,379,252 --outFileName Ccr7.svg
pyGenomeTracks --tracks NK.ini --region chr6:128,587,253-131,480,042 --outFileName NK.svg
pyGenomeTracks --tracks Tox.ini --region chr4:6,642,850-7,531,941 --outFileName Tox.svg
pyGenomeTracks --tracks Ctla4.ini --region chr1:60,720,001-61,192,045 --outFileName Ctla4.svg
pyGenomeTracks --tracks Pdcd1.ini --region chr1:93993232-94,254,122 --outFileName Pdcd1.svg
pyGenomeTracks --tracks Entpd1.ini --region chr19:40,535,109-41,128,512 --outFileName Entpd1.svg
pyGenomeTracks --tracks Bach2.ini --region chr4:32,032,372-33,021,028 --outFileName Bach2.svg















######################## everything below are details on how to setup the files etc





# # I prepared tracks in the Bach2 ini file...add also gene gtf files and test the regions I prepared for the selected genes
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/human/bcl11b/hg38/HTC19_mean_rep_hg38_subtract_notSkipped.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/human/bcl11b/hg38/HTC19_mean_rep_hg38_subtract_SkippedNCR.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/human/HTC19B_bcl11b-cutrun_cat_merge_0005_non_stringent.stringent.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/human/HTC19B_bcl11b-cutrun_cat_merge_001_non_stringent.stringent.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/human/bcl11b/hg38/HTC41_mean_rep_hg38_subtract_notSkipped.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/human/bcl11b/hg38/HTC41_mean_rep_hg38_subtract_SkippedNCR.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/human/HTC41B_bcl11b-cutrun_cat_merge_0005_non_stringent.stringent.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/human/HTC41B_bcl11b-cutrun_cat_merge_001_non_stringent.stringent.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/human/bcl11b/hg38/atac/atac_HTA19C_mean_BR1-BR2.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/human/p19_hg38_atac_down_peaks_P005_for_annotation.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/human/bcl11b/hg38/atac/atac_HTA19B_mean_BR1-BR2.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/human/p19_hg38_atac_up_peaks_P005_for_annotation.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/human/bcl11b/hg38/atac/atac_HTA41C_mean_BR1-BR2.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/human/p41_hg38_atac_down_peaks_P005_for_annotation.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/human/bcl11b/hg38/atac/atac_HTA41B_mean_BR1-BR2.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/human/p41_hg38_atac_up_peaks_P005_for_annotation.bed


# # regions for HUMAN that were not used
# chr5:134,090,325-134,270,791	TCF7  decreased atac in both KO and no Bcl11b peaks
# chr19:35,887,070-35,918,717	TYROBP no atac changes
# chr1:160,444,061-160,733,375	SLAMF6 no atac changes
# chr8:58,755,017-59,596,529	TOX  .rather UP atac in KO
# chr10:95,673,199-96,295,147	ENTPD1 rather UP atac in KO
# chr14:24,627,868-24,658,966	GZMB nice Bcl11b peaks but no ATAC changes
# pyGenomeTracks --tracks 03_template_human.ini --region chr6:89,699,001-91,286,962 --outFileName BACH2.png




# # default installation installed very old version 3.0 so I need to specify the version I want
# # I was getting this error: warning  libmamba Problem type not implemented SOLVER_RULE_STRICT_REPO_PRIORITY which was solved by making chanel flexible
# conda config --set channel_priority flexible
# conda create -n pygenometracks -c bioconda -c conda-forge pygenometracks=3.9

# conda activate pygenometracks
# cd /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks


# # this uses list of regions to plot them all - use this with the template ini file to find the optimal scaling and then put WT and KO on the same scale for each gene individually
# pyGenomeTracks --tracks 03_template.ini --BED regions.bed --outFileName SETUP.png
# pyGenomeTracks --tracks 03_template_human.ini --BED regions_human.bed --outFileName SETUP.png



# # original before I made some adjustments
# # pyGenomeTracks --tracks Fcer1g.ini --region chr1:171,224,562-171,283,158 --outFileName Fcer1g.png
# # pyGenomeTracks --tracks Pdcd1.ini --region chr1:94,019,737-94,254,122 --outFileName Pdcd1.png



# # regions for HUMAN
# chr1:161,193,349-161,265,677	FCER1G
# chr5:134,090,325-134,270,791	TCF7
# chr19:35,887,070-35,918,717	TYROBP
# chr1:160,444,061-160,733,375	SLAMF6
# chr1:168,546,250-168,690,906	XCL1
# chr17:40,549,342-40,597,718	CCR7
# chr8:58,755,017-59,596,529	TOX
# chr2:203,850,090-203,990,817	CTLA4
# chr2:241,835,938-241,906,301	PDCD1
# chr10:95,673,199-96,295,147	ENTPD1
# chr6:89,699,001-91,286,962	BACH2
# chr12:10,351,634-10,506,513	KLRK1
# chr19:54,893,436-55,029,675	NCR1
# chr6:31,570,765-31,620,388	NCR3
# chr10:5,987,437-6,298,410	IL2RA
# chr10:70,584,809-70,662,552	PRF1
# chr14:24,627,868-24,658,966	GZMB


# # regions to be used:
# chr1:171,224,562-171,283,158	Fcer1g
# chr11:52,209,926-52,669,437	Tcf7
# chr7:30,385,644-30,561,256	Tyrobp
# chr1:171,902,126-172,138,148	Slamf6
# chr1:164,914,243-165,150,265	Xcl1
# chr11:99,084,436-99,379,252	Ccr7
# chr6:128,587,253-131,480,042	NK
# chr4:6,642,850-7,531,941	Tox
# chr1:60,720,001-61,192,045	Ctla4
# chr1:94,019,737-94,254,122	Pdcd1
# chr19:40,535,109-41,128,512	Entpd1
# chr4:32,032,372-33,021,028	Bach2

# # regions not to be used:
# chr10:61,219,580-61,535,195	Prf1
# chr15:9,477,719-9,728,973	Il7r
# chr6:124,895,817-124,985,589	Lag3
# chr16:43,592,753-43,829,914	Tigit
# chr11:46,411,506-46,706,322	Havcr2
# chr14:56,112,121-56,413,747	Granzymes
# chr1:106,469,241-107,413,331	Bcl2




# # note that bed files require to have only 6 fields so I had to modify them and place them in the input folder c('chr','start','end','id','score','strand')
# cd /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/overlaps_intersections/differential/
# for i in *;do awk -F'\t' -v OFS="\t" '{print $1, $2, $3, $4, $5, $6 }' $i > /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/$i; done

# cd /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/ROSE
# for i in output*;do awk -F'\t' -v OFS="\t" '{print $1, $2, $3, $4, "1000", "." }' $i > /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/$i; done

# cd /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/ROSE/SE_withRNAseq/
# for i in 02b_coord-first_final_*;do awk -F'\t' -v OFS="\t" '{print $1, $2, $3, $4, "1000", "." }' $i > /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/$i; done

# $ modify bcl11b peaks separately
# awk -F'\t' -v OFS="\t" '{print $1, $2, $3, $6, $4, "." }' /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/wt_bcl11b_cat_merge_0005_non_stringent.stringent.bed > wt_bcl11b_cat_merge_0005_non_stringent.stringent2.bed

# # FOR HUMANS
# cd /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/human/before_modif
# for i in p*;do awk -F'\t' -v OFS="\t" '{print $1, $2, $3, $4, $5, $6 }' $i > ../$i; done

# # modify bcl11b peaks separately
# cd /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/human/before_modif
# for i in HTC*;do awk  -F'\t' -v OFS="\t" '{print $1, $2, $3, $6, $4, "." }' $i > ../$i; done




# # this autogenerates config file that can then be adapted
# make_tracks_file --trackFiles <file1.bed> <file2.bw> ... -o tracks.ini


# make_tracks_file --trackFiles /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/bcl11b/bcl11b_wt_mean_rep_mm10_treat_pileup.bw /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/overlaps_intersections/wt_bcl11b_cat_merge_0005_non_stringent.stringent.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/atacseq/mouse/mC8ATW_smooth.bw /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/down_atac_CD8T_peaks_P005_for_annotation_genrich.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/atacseq/mouse/mC8ATK_smooth.bw /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/up_atac_CD8T_peaks_P005_for_annotation_genrich.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k27me3/h3k27me3_wt_mean_rep_mm10_treat_pileup.bw /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/down_h3k27me3_logFC02_Pval03__Pval005.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k27me3/h3k27me3_ko_mean_rep_mm10_treat_pileup.bw /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/up_h3k27me3_logFC02_Pval03__Pval005.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k27ac/h3k27ac_wt_mean_rep_mm10_treat_pileup.bw /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/down_h3k27ac_logFC02_Pval03__Pval005.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k27ac/h3k27ac_ko_mean_rep_mm10_treat_pileup.bw /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/up_h3k27ac_logFC02_Pval03__Pval005.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k4me1/h3k4me1_wt_mean_rep_mm10_treat_pileup.bw /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/down_h3k4me1_logFC02_Pval03__Pval005.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k4me1/h3k4me1_ko_mean_rep_mm10_treat_pileup.bw /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/up_h3k4me1_logFC02_Pval03__Pval005.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/output_wt_TE-with27acBAM.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/output_ko_TE-with27acBAM.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/output_wt_TE-withK27me3BAM_PE.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/output_ko_TE-withK27me3BAM_PE.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/02b_coord-first_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.bed /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/02b_coord-first_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.bed -o template_tracks.ini


# make_tracks_file --trackFiles /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/pygenometracks/input/gencode.vM25.annotation.gtf.gz -o template_gtf.ini


# # there are many predicted or non-coding genes which I don't want to show so I filter them
# gunzip gencode.vM25.basic.annotation.gtf.gz
# grep -Ev '\b\w*Rik\b|\bGm[0-9]+' gencode.vM25.basic.annotation.gtf > gencode.vM25.basic.filt-Rik-Gm_genes.gtf

# # for humans
# # generate a grep command that will filter out rows that contain words that start with any MIR, RPL, RN7SK or any of the following, followed by numbers: "RF, AL, AC, LINC"
# grep -Ev '\b(MIR|RPL|DNAJC19P6|RN7SK|RF[0-9]|AL[0-9]|AC[0-9]|LINC[0-9])' filt_refgenie_hg38.gtf > refgenie_hg38_filtRiboNonCodeGenes.gtf
# gzip refgenie_hg38_filtRiboNonCodeGenes.gtf



# ########### backup links and names


# path_bw <- '/Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/'
# bigwig_files <- paste0(path_bw,rep(c('bcl11b/bcl11b_wt_mean_rep_mm10_treat_pileup.bw','h3k4me1/h3k4me1_wt_mean_rep_mm10_treat_pileup.bw','h3k4me1/h3k4me1_ko_mean_rep_mm10_treat_pileup.bw','h3k27ac/h3k27ac_wt_mean_rep_mm10_treat_pileup.bw', 'h3k27ac/h3k27ac_ko_mean_rep_mm10_treat_pileup.bw','h3k27me3/h3k27me3_wt_mean_rep_mm10_treat_pileup.bw','h3k27me3/h3k27me3_ko_mean_rep_mm10_treat_pileup.bw')))
# track_names <-  c('bcl11b','h3k4me1_wt','h3k4me1_ko', 'h3k27ac_wt', 'h3k27ac_ko', 'h3k27me3_wt','h3k27me3_ko')








# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/bcl11b/bcl11b_wt_mean_rep_mm10_treat_pileup.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/overlaps_intersections/wt_bcl11b_cat_merge_0005_non_stringent.stringent.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/atacseq/mouse/mC8ATW_smooth.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/overlaps_intersections/differential/down_atac_CD8T_peaks_P005_for_annotation_genrich.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/atacseq/mouse/mC8ATK_smooth.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/overlaps_intersections/differential/up_atac_CD8T_peaks_P005_for_annotation_genrich.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k27me3/h3k27me3_wt_mean_rep_mm10_treat_pileup.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/overlaps_intersections/differential/down_h3k27me3_logFC02_Pval03__Pval005.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k27me3/h3k27me3_ko_mean_rep_mm10_treat_pileup.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/overlaps_intersections/differential/up_h3k27me3_logFC02_Pval03__Pval005.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k27ac/h3k27ac_wt_mean_rep_mm10_treat_pileup.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/overlaps_intersections/differential/down_h3k27ac_logFC02_Pval03__Pval005.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k27ac/h3k27ac_ko_mean_rep_mm10_treat_pileup.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/overlaps_intersections/differential/up_h3k27ac_logFC02_Pval03__Pval005.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k4me1/h3k4me1_wt_mean_rep_mm10_treat_pileup.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/overlaps_intersections/differential/down_h3k4me1_logFC02_Pval03__Pval005.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/cutrun/mouse/h3k4me1/h3k4me1_ko_mean_rep_mm10_treat_pileup.bw
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/overlaps_intersections/differential/up_h3k4me1_logFC02_Pval03__Pval005.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/ROSE/output_wt_TE-with27acBAM.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/ROSE/output_ko_TE-with27acBAM.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/ROSE/output_wt_TE-withK27me3BAM_PE.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/ROSE/output_ko_TE-withK27me3BAM_PE.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/ROSE/SE_withRNAseq/02b_coord-first_final_TE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.bed
# /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/ROSE/SE_withRNAseq/02b_coord-first_final_SE_nonFilt_with_scRNAseq_ovarTumCD8_transcripID_TSS0bp.bed
