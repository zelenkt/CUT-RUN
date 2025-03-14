#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=8:00:00
#SBATCH --job-name=DESeq2_SEACR
#SBATCH --output=DESEq2_SEACR.%j.out
#SBATCH --error=DESEq2_SEACR.%j.err

##### LAST UPDATED 01/16/2025 Tomas Zelenka #####


##### ARGUMENTS AND USAGE #####
usage="$(basename $0) [-h] [-d dir] [-r file] [-t string] [-c string]


Optional Arguments:
    -h  Show this help text
    -d  Set target directory (directory to write analysis files) (default: $PWD/DESeq2)
    -s  Set ATAC-seq / ChIP-seq sample list file (default: /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/tools/test-samples.csv)
    -t  Set test condition (column from ATAC-seq / ChIP-seq sample list file) (i.e. genotype or cell_subset) (default: genotype)
    -c  Set control group (i.e. WT if test condition is genotype, or Th0 if test condition is cell_subset (default: WT)
    -r  Set Raw reads counts as an output from CUT&RUN pipeline, from feature counts or from other sources

Example usage:
    sbatch --chdir=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10 --export=s=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10/samples.csv,f=_h3k4me1_merged_peaks_coverage.bed ~/bin/CUTRUN_DESeq2_for_seacr.sh
    
    sbatch --chdir=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10 --export=s=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10/samples.csv,f=_SE-nonFilt_k4me1_with_k27acBAM_merged_peaks_coverage.bed ~/bin/CUTRUN_DESeq2_for_seacr.sh

    sbatch --chdir=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10 --export=s=/share/lab_avram/HPC_Cluster/user/tomas/08_cut_run_5-30-2023/h3k27ac/cd8_ovar_tumor_mm10/samples.csv,f=_PE-nonFilt_k4me1_with_k27me3BAM_merged_peaks_coverage.bed ~/bin/CUTRUN_DESeq2_for_seacr.sh
    "



while getopts ':h:d:s:t:c:' option; do
    case "$option" in
    h)  echo "$usage"
        exit 0
        ;;
	d)  if [[ -d "$OPTARG" ]]; then
            d="$OPTARG"
            printf "Target directory manually set to %s\n" "$OPTARG"
        else
            printf "%s doesn't exist or is not a valid directory. Please specify a valid directory and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"
            exit 1
        fi
        ;;
    
    r)  if [[ -r "$OPTARG" ]]; then
            r="$OPTARG"
            printf "Raw read counts manually set to: %s\n" "$OPTARG"
        else
            printf "%s is empty or does not exist. Please specify a valid file with Raw read counts and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"
        fi
        ;;
    s)  if [[ -s "$OPTARG" ]]; then
            s="$OPTARG"
            printf "ATAC-seq / ChIP-seq sample list file manually set to: %s\n" "$OPTARG"
        else
            printf "%s is empty or does not exist. Please specify valid ATAC-seq / ChIP-seq sample list file and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"
        fi
        ;;
    t)  t="$OPTARG"
        printf "Test condition manually set to: %s\n" "$OPTARG"
        ;;
    f)  f="$OPTARG"
        printf "Test condition manually set to: %s\n" "$OPTARG"
        ;;
    c)  c="$OPTARG"
        printf "Control group manually set to: %s\n" "$OPTARG"
        ;;
   \?)  printf "Invalid option: -%s.\n" "$OPTARG"
        echo "$usage"
        exit 1
        ;;
    :)  printf "Invalid option: -%s requires an argument.\n" "$OPTARG"
        echo "$usage"
        exit 1
        ;;
    esac
done
shift $((OPTIND -1))

date;pwd

##### LOAD MODULES #####
ml purge
module load R/4.1.0-foss-2020b

##### Use Avram Lab's shared R library #####
export R_LIBS="/share/lab_avram/HPC_Cluster/share/R/x86_64-pc-linux-gnu-library/4.1"


# suffix="_SE-nonFilt_k4me1_with_k27acBAM_merged_peaks_coverage.bed"
suffix="${f:-_merged_peaks_coverage.bed}"

species="mouse"

if [[ $species == "mouse" ]]; then
  printf "%sSelected species is mouse\n"
deseq2_genome="mm10"

else
if [[ $species == "human" ]]; then
  printf "%sSelected species is human\n"
  deseq2_genome="hg38"
  else
  printf "%sMake sure to enter the correct species\n"
fi
fi

##### SET VARIABLES #####
target_dir="${d:-${PWD}/DESeq2_seacr}"
if [[ "$target_dir" == "$PWD/DESeq2_seacr" && ! -d "$target_dir" ]]; then # IF TARGET DIR DOES NOT EXIST
    mkdir -p "$target_dir"
elif [[ ! -d "$target_dir" ]]; then # THIS CAN ONLY BE TRUE IF SOMEONE CALLED THIS FUNCTION THROUGH SBATCH AND DIRECTORY EXPORTED AN INVALID "d" variable. (e.g. sbatch --export=d=/some/invalid/dir DESeq2.sh)
    printf "%s doesn't exist or is not a valid directory. Please specify a valid directory and rerun %s. Exiting...\n" "$target_dir" "$(basename $0)"
    exit 1
fi

samples="${s:-${PWD}/samples.csv}"  # samples="${r:-${PWD}/DESeq2/test-samples.csv}"
[[ -s "$samples" ]] || { printf "%s is empty or does not exist. Please specify valid ATAC-seq / ChIP-seq sample list file and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"; exit 1; } # VERIFY SAMPLE LIST FILE IS VALID

condition="${t:-genotype}"
read -ra condition_list <<<$(head -n1 "$samples") # GENERATE ARRAY OF POSSIBLE CONDITIONS
[[ "${condition_list[@]}" =~ "$condition" ]] || { printf "%s is not a valid condition.  Please select a condition from the following list: %s. Exiting...\n" "$condition" "${condition_list[*]}"; exit 1; } #VERIFY SELECTED CONDITION IS VALID

baseline="${c:-WT}"
baseline_check=$(awk -F$',' -v condition=$condition -v baseline=$baseline 'BEGIN {s=0} 
    NR==1{for(i=1; i<=NF; i++) if($i==condition) {a[i]++; break} }
    {for (i in a) if ($i==baseline) {s++; exit} }
    END {print s}' "$samples")
[[ $baseline_check > 0 ]] || { printf "There are no ATAC-seq / ChIP-seq samples of type \"%s\" under condition \"%s\". Please select a valid control group and rerun %s. Exiting...\n" "$baseline" "$condition" "$(basename $0)"; exit 1; } #VERIFY THAT AT LEAST ONE SAMPLE UNDER SELECTED CONDITION MATCHES BASELINE TYPE


## Extract the read_count values from individual ref_coverage files - coverage of the WT+KO consensus peaks generated by project level pepatac or from CUT&RUN pipeline
#samples="/share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/tools/samples_C4P.csv"; 
if [[ ! -f ${target_dir}/raw_read_counts_cutrun_consensus_peaks.txt ]]; then
    printf "Extracting raw read counts for consensus peaks from individual samples...\n" "$base"

    # note that for MACS2 output I should print awk '{print $11} but for seacr it is awk '{print $7}
    sampleCount=$(awk -F',' 'NR>1{print $1}' $samples | wc -l); awk -F',' 'NR>1{print $1}' $samples | sort > tmp_4245277sampleNames.txt; for ((sampleOrder=1; sampleOrder<=$sampleCount; sampleOrder++)); do sampleName=$(sed -n ${sampleOrder}p tmp_4245277sampleNames.txt); awk '{print $7}' ${sampleName}${suffix} > tmp_3541643543534154_${sampleName}; echo -e '' > tmp_9878798798798_consensus_${sampleOrder}; awk -v OFS="_" '{print $1,$2,$3}' ${sampleName}${suffix} >> tmp_9878798798798_consensus_${sampleOrder}; done; 
    awk '
    { 
        for (i=1; i<=NF; i++)  {
            a[NR,i] = $i
        }
    }
    NF>p { p = NF }
    END {    
        for(j=1; j<=p; j++) {
            str=a[1,j]
            for(i=2; i<=NR; i++){
                str=str" "a[i,j];
            }
            print str
        }
    }' tmp_4245277sampleNames.txt | sed -e 's/ /\t/g' > tmp_13213213213132sampleNames.txt; paste tmp_3541643543534154* > tmp999898018194245242; cat tmp_13213213213132sampleNames.txt tmp999898018194245242 > tmp8184414654654read_counts.txt; rm tmp_13213213213132sampleNames.txt; rm tmp_3541643543534154*; rm tmp_4245277sampleNames.txt; rm tmp999898018194245242; paste tmp_9878798798798_consensus_1 tmp8184414654654read_counts.txt | awk '!($1 ~ /^chrUn/ || $1 ~ /random/ || $1~/^chrM/ || $1~/^chrY/) {print $0}' > raw_read_counts_cutrun_consensus_peaks.txt; rm tmp_9878798798798_consensus*; rm tmp8184414654654read_counts.txt; mv raw_read_counts_cutrun_consensus_peaks.txt "$target_dir"/
    
    # this one keeps the tmp files
    # tmp_4245277sampleNames.txt | sed -e 's/ /\t/g' > tmp_13213213213132sampleNames.txt; paste tmp_3541643543534154* > tmp999898018194245242; cat tmp_13213213213132sampleNames.txt tmp999898018194245242 > tmp8184414654654read_counts.txt; paste tmp_9878798798798_consensus_1 tmp8184414654654read_counts.txt | awk '!($1 ~ /^chrUn/ || $1 ~ /random/ || $1~/^chrM/ || $1~/^chrY/) {print $0}' > raw_read_counts_cutrun_consensus_peaks.txt; mv raw_read_counts_cutrun_consensus_peaks.txt "$target_dir"/

# modify to line below to remove tmp files
# tmp_4245277sampleNames.txt | sed -e 's/ /\t/g' > tmp_13213213213132sampleNames.txt; paste tmp_3541643543534154* > tmp999898018194245242; cat tmp_13213213213132sampleNames.txt tmp999898018194245242 > tmp8184414654654read_counts.txt; rm tmp_13213213213132sampleNames.txt; rm tmp_3541643543534154*; rm tmp_4245277sampleNames.txt; rm tmp999898018194245242; paste tmp_9878798798798_consensus_1 tmp8184414654654read_counts.txt | awk '!($1 ~ /^chrUn/ || $1 ~ /random/ || $1~/^chrM/ || $1~/^chrY/) {print $0}' > raw_read_counts_cutrun_consensus_peaks.txt; rm tmp_9878798798798_consensus*; rm tmp8184414654654read_counts.txt; mv raw_read_counts_cutrun_consensus_peaks.txt "$target_dir"/

    printf "Done.\n"
else
    printf "Raw read counts file already exists\n" "$base"
fi


readcounts="${r:-${target_dir}/raw_read_counts_cutrun_consensus_peaks.txt}"


## here are some backup commands that can be used to perform differential analysis with the simplified consensus peaks from pepatac
####readcounts="${r:-${PWD}/summary/PEPATAC_tomas_mm10_ref_peaks_coverage.tsv}"  # samples="${r:-${PWD}/DESeq2/test-samples.csv}"
####[[ -r "$readcounts" ]] || { printf "%s is empty or does not exist. Please specify  a valid file with Raw read counts and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"; exit 1; } # VERIFY RAW READ COUTNS FILE IS VALID

# prepare the raw read count output from pepatac (the original TSV file needs to be converted like this otherwise it creates problems when preparing matrix)
####colCount=$(awk '{print NF; exit}' ${readcounts}); for ((col=1; col<=$colCount; col++)); do awk -v col="$col" '{print $(col)}' ${readcounts} > tmp_columns_for_readcountcoverage_$col; done; paste tmp_columns_for_readcountcoverage* > $target_dir/pepatac_raw_read_counts.txt; rm tmp_columns_for_readcountcoverage*


##### BEGIN SCRIPT #####
printf "Condition = %s\n" "$condition"
printf "Baseline = %s\n" "$baseline"
printf "ATAC-seq / ChIP-seq Sample list file = %s\n" "$samples"
printf "Raw read counts file = %s\n" "$readcounts"

( cd "$target_dir" && Rscript "$HOME"/bin/cutrun_DESeq2.R -c "$condition" -b "$baseline" -s "$samples" -r "$readcounts" ) # RUN DESEQ2.R WITHIN SUBSHELL TO CHANGE DIRECTORY TO ANALYSIS DIR

# ( cd "$target_dir" && Rscript "$HOME"/bin/annotate_genes.R -i DESeq2_analysis.txt -o DESeq2_analysis_annotated.txt ) # RUN ANNOTATE_GENES_DESE2.R WITHIN SUBSHELL TO CHANGE DIRECTORY TO ANALYSIS DIR



# prepare filtered DESeq2 output using p < 0.05 and separately padj < 0.05 - these significantly differential peaks can be then used for peak annotation
# the files contain the following columns
# chr start end peakID -log10P strand baseMean log2FC pvalue padj

cd $target_dir
ml BEDTools/2.30.0-GCC-11.2.0


# the order of columns in the all_differential_peaks is
# chr start end name -logPval strand baseMean($4) log2FC($5) pval($8) padj($9)
cat DESeq2_analysis.txt | tr -d '"' | awk -F'_|\t' -v OFS="\t" 'NR>1 {print $1, $2, $3, "peak_"NR,-log($8)/log(10), ".", $4, $5, $8, $9}' | sortBed | uniq > all_differential_peaks_noFilt_for_annotation.bed

cat DESeq2_analysis.txt | tr -d '"' | awk -F'_|\t' -v OFS="\t" 'NR>1 && $8<0.05 {print $1, $2, $3, "peak_"NR,-log($8)/log(10), ".", $4, $5, $8, $9}' | sortBed | uniq > all_differential_peaks_P005_for_annotation.bed

awk -v OFS="\t" '$8<0 {print $0}' all_differential_peaks_P005_for_annotation.bed | sortBed | uniq > down_peaks_P005_for_annotation.bed
awk -v OFS="\t" '$8>0 {print $0}' all_differential_peaks_P005_for_annotation.bed | sortBed | uniq > up_peaks_P005_for_annotation.bed


cat DESeq2_analysis.txt | tr -d '"' | awk -F'_|\t' -v OFS="\t" 'NR>1 && $9<0.05 {print $1, $2, $3, "peak_"NR,-log($8)/log(10), ".", $4, $5, $8, $9}' | sortBed | uniq > all_differential_peaks_FDR005_for_annotation.bed

awk -v OFS="\t" '$8<0 {print $0}' all_differential_peaks_FDR005_for_annotation.bed | sortBed | uniq > down_peaks_FDR005_for_annotation.bed
awk -v OFS="\t" '$8>0 {print $0}' all_differential_peaks_FDR005_for_annotation.bed | sortBed | uniq > up_peaks_FDR005_for_annotation.bed





mv ../*.err $target_dir/
mv ../*.out $target_dir/



date

