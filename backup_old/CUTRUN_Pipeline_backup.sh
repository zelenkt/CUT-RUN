#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=16GB
#SBATCH -t 24:00:00
#SBATCH --job-name=CUTRUN_Pipeline
#SBATCH --output=CUTRUN_Pipeline.%j.out
#SBATCH --error=CUTRUN_Pipeline.%j.err

##### LAST UPDATED 2022.11.09 Tomas Zelenka #####
##### LAST UPDATED 2021.11.30 Brian Connolly #####

##### OPTIONAL ARGUMENTS #####
usage="$(basename $0) [-h] [-p] [-d dir] [-n string] [-p PE] -- Program to analyze fastq files from ATAC-seq and ChIP-seq.

Where:
    -h  Show this help text
    -m  Active 'pre-macs' mode. Halts program after final BAM file is produced. (default: off)
    -d  Set target directory (directory containing ATAC-seq fastq files) (default: Current directory)
    -n  Set sample base name (i.e. ATAC_Th2_WT_1452) (default: Name of directory containing files)
    -p  Set sample as PE (paired end) or SE (single end) (default: auto-detect)

Example usage:
    sbatch --chdir=/directory/containing/fastq/files --export=pre_macs=true,is_paired_end=false $HOME/bin/$(basename $0)

    or

    $(basename $0) -d /directory/containing/fastq/files -n ATAC_Th2_WT_1452 -m -p SE"

while getopts ':hmd:n:p:' option; do
    case "$option" in
	h) echo "$usage"
	   exit 0
	   ;;
	m) printf "Pre-macs mode selected. Program will terminate before calling peaks.\n."
	   pre_macs=true
	   ;;
	d) if [[ -d "$OPTARG" ]]; then
	       cd "$OPTARG"
	       printf "Working directory set to %s\n" "$OPTARG"
	   else
	       printf "%s doesn't exist or isn't a valid directory.\n" "$OPTARG"
	       exit 1
	   fi
	   ;;
	n) base="$OPTARG"
	   printf "Sample base name set to %s.\n" "$OPTARG"
	   ;;
	p) if [[ "$OPTARG" == "PE" ]]; then
	       is_paired_end=true
	       printf "Paired-end mode has been selected.\n"
	   elif [[ "$OPTARG" == "SE" ]]; then
	       is_paired_end=false
	       printf "Single-end mode has been selected.\n"
	   else
	       printf "Invalid mode selected \"%s\". -p option must be either PE or SE.\n" "$OPTARG"
	       exit 1
	   fi
	   ;;
       \?) printf "Invalid option: -%s.\n" "$OPTARG"
	   echo "$usage"
	   exit 1
	   ;;
	:) printf "Invalid option: -%s requires an argument.\n" "$OPTARG"
	   echo "$usage"
	   exit 1
	   ;;
    esac
done
shift $((OPTIND -1))

##### LOAD MODULES #####
# Modules are loaded in-line below as needed 

# At UF, these modules were loaded to ACCESS A BETTER SORT COMMAND. 
# commenting out at Moffitt. Will check performance.  
# Used in sort command on line 247
# module load gcc/4.7.2 coreutils 


##### SET VARIABLES #####
base="${base-${PWD##*/}}"
pre_macs="${pre_macs-false}"
trim_file=$HOME/lib/trim_file.fa
index=/share/lab_avram/HPC_Cluster/reference/bowtie2/mm10/mm10
chr_sizes=/share/lab_avram/HPC_Cluster/reference/genomes/mm10/ucsc/mm10.chrom.sizes

#CREATE CLEANED UP CHROMOSOME LIST FOR MM10 IN USER LIB IF IT DOESN'T ALREADY EXIST
[[ -f $HOME/lib/mouse_chr_list.txt ]] || awk '{if ($1 !~ /[M,Y,_]/) print $1}' $chr_sizes > $HOME/lib/mouse_chr_list.txt
chr_list=$HOME/lib/mouse_chr_list.txt

##### SET JAVA&TMPDIR OPTIONS #####
export _JAVA_OPTIONS="-Xmx3g"
#CREATE TMP DIR IF IT DOESN'T ALREADY EXIST
[[ -d tmp ]] || mkdir tmp
export TMPDIR="$PWD"/tmp

date; pwd

##### BEGIN SCRIPT #####
#CREATES "Logs" DIRECTORY IF IT DOESN'T EXIST
[[ -d "Logs" ]] || mkdir Logs

#DETERMINE IF SAMPLE IS PAIRED END OR SINGLE END
if [[ -z $is_paired_end ]]; then
    if test -n "$(find -L "$PWD" -maxdepth 1 -name '*_R2_001.fastq.gz' -o -iname '*_Reverse.fastq.gz')"; then
	is_paired_end=true
	printf "%s detected to be paired end.\n" "$base"
    else
	is_paired_end=false
	printf "%s detected to be single end.\n" "$base"
    fi
fi

#CONCATENATE READ FILES
#FORWARD READS
if ! [[ -f "$base"_Forward.fastq.gz || -f "$base".sam || -f "$base"_nodups.bam || -f "$base"_full.bam || -f "$base".bam ]]; then
    printf "Concatenating forward reads..."
    cat *R1_001.fastq.gz > "$base"_Forward.fastq.gz
    printf "Done.\n"
fi

#REVERSE READS
if $is_paired_end && ! [[ -f "$base"_Reverse.fastq.gz || -f "$base".sam || -f "$base"_nodups.bam || -f "$base"_full.bam || -f "$base".bam ]]; then
    printf "Concatenating reverse reads..."
    cat *R2_001.fastq.gz > "$base"_Reverse.fastq.gz
    printf "Done.\n"
fi

#BEGIN TRIMMING
if [[ ! -f "$base"_1P.fastq.gz ]]; then
    printf "Trimming read via Trimmomatic..."
    module load Trimmomatic/0.39-Java-1.8
    if $is_paired_end; then
	trimmomatic PE -threads 7 -phred33 -baseout "$base".fastq.gz "$base"_Forward.fastq.gz "$base"_Reverse.fastq.gz ILLUMINACLIP:"$trim_file":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> Logs/"$base".trim.log
    else
	trimmomatic SE -threads 7 -phred33 "$base"_Forward.fastq.gz "$base"_1P.fastq.gz ILLUMINACLIP:"$trim_file":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> Logs/"$base".trim.log
    fi
    printf "Done.\n"
else
    printf "Trimmed fastq files already detected...skipping trimmomatic.\n"
fi

#BEGIN BOWTIE2 ALIGNING
printf "Reloading bowtie2"
if ! [[ -f "$base".sam || -f "$base"_nodups.bam || -f "$base".bam || -f "$base"_full.bam ]]; then
    module load Bowtie2
    printf "Aligning read to mm10 genome via bowtie2..."
    if $is_paired_end; then
	bowtie2  --dovetail --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 7 --met-file Logs/"$base".b2.metrics -x "$index" -1 "$base"_1P.fastq.gz -2 "$base"_2P.fastq.gz 2> Logs/"$base".b2.log | awk 'NR==FNR {chr[$1]; next}; {if ($1 ~ /@/ || ($3 in chr && $5>30)) print $0}' "$chr_list" - > "$base".sam
    else
	bowtie2 --dovetail --local --very-sensitive-local --no-unal --no-mixed  --no-discordant --phred33 -I 10 -X 700 -p 7 --met-file Logs/"$base".b2.metrics -x "$index" -U "$base"_1P.fastq.gz 2> Logs/"$base".b2.log | awk 'NR==FNR {chr[$1]; next}; {if ($1 ~ /@/ || ($3 in chr && $5>30)) print $0}' "$chr_list" - > "$base".sam
    fi # tomas added parameters to the bowtie mapping to be consistent with the original protocol
    module unload Bowtie2
    printf "Done.\n"
else
    printf "Alignment file already detected...skipping bowtie2.\n"
fi

#SORT SAM FILE INTO BAM FILE
if ! [[ -f "$base"_nodups.bam || -f "$base".bam || -f "$base"_full.bam ]]; then
    module load SAMtools
    printf "Sorting and compressing %s.sam into bam file via samtools..." "$base"
    samtools sort -@ 7 -m 3G -l 0 -O bam -o "$base"_full.bam -T $SLURM_JOB_ID "$base".sam
    printf "Done.\n"
    module unload SAMtools
else
    printf "Sorted alignment file already detected...skipping samtools sort.\n"
fi

#REMOVE DUPLICATES
if ! [[ -f "$base"_nodups.bam || -f "$base".bam ]]; then
    module load picard
    # Touch output files due to bug in picard
    touch "$base"_nodups.bam
    touch Logs/"$base".picard.metrics
    printf "Removing duplicates from %s_full.bam via Picardtools...\n" "$base"
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates I="$base"_full.bam O="$base"_nodups.bam REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate m=Logs/"$base".picard.metrics TMP_DIR="$PWD"/tmp VALIDATION_STRINGENCY=SILENT
    printf "Done.\n"
else
    printf "Deduplicated alignment file already detected...skipping Picardtools.\n"
fi

#RESORT BAM FILE
if [[ ! -f "$base".bam ]]; then
    module load SAMtools
    printf "Resorting %s_nodups.bam via samtools..." "$base"
    samtools sort -@ 7 -m 3G -l 9 -O bam -o "$base".bam -T $SLURM_JOBID "$base"_nodups.bam
    printf "Done.\n"
else
    printf "Final sorted bam file already detected...skipping samtools sort.\n"
fi

#REMOVE INTERMEDIATE FILES
rm "$base".sam
#rm "$base"_full.bam
#rm "$base"_nodups.bam
rm "$base"_Forward.fastq.gz
rm "$base"_Reverse.fastq.gz

#BEGIN FASTQC
if [[ ! -f FastQC/"$base"_fastqc.html ]]; then
    printf "Analyzing %s.bam via FastQC...\n" "$base"
    module load FastQC
    [[ -d "FastQC" ]] || mkdir FastQC #Create FastQC directory if it doesn't already exist
    # fastqc "$base".bam -o FastQC
    fastqc "$base".bam
    mv -v "$base"_fastqc.* FastQC/
    printf "Done.\n"
else
    printf "%s_fastqc.html already exists...skipping FastQC analysis.\n" "$base"
fi

#EXIT IF PRE-MACS SETTING SELECTED
if $pre_macs; then
    printf "Pre-macs mode selected. Exiting...\n"
    rm -r tmp
    date
    exit 0
fi

#CALL PEAKS
printf "Calling peaks section...\n"
module load MACS2
module load parallel/20180822-foss-2018b
if ! [[ -f Peaks/"$base"_treat_pileup.bdg || -f Peaks/"$base"_treat_pileup.bdg.gz ]]; then
    printf "Calling peaks with macs2...\n"
    macs2 callpeak -t "$base".bam --call-summits -g mm -p 0.00001 -f BAMPE --SPMR --outdir Peaks -B -n "$base"
    # version 1 eric:  macs2 callpeak -t "$base".bam --call-summits --nomodel -g mm -q 0.05 -f BAM --SPMR --outdir Peaks -B -n "$base"
    # version 2 nomodel: macs2 callpeak -t "$base".bam --call-summits --nomodel -g mm -p 1e-5 -f BAMPE --SPMR --outdir Peaks -B -n "$base"
    # version 3 model bampe: macs2 callpeak -t "$base".bam --call-summits -g mm -p 1e-5 -f BAMPE --SPMR --outdir Peaks -B -n "$base" 
    printf "Done.\n"
else
    printf "Pileup files already detected...skipping macs2 peak calling.\n"
fi

#CLEAN PILEUP FILE (VIA SUBTRACT), CLEAN, SORT, THEN CONVERT TO BIGWIG
if ! [[ -f Peaks/"$base"_subtract.bdg || -f Peaks/"$base"_subtract.bw ]]; then
    printf "Cleaning pileup files via subtract...\n"
    macs2 bdgcmp -t Peaks/"$base"_treat_pileup.bdg -c Peaks/"$base"_control_lambda.bdg -m subtract -o Peaks/"$base"_subtract.bdg
    printf "Done.\n"

    temp_file="$(mktemp --tmpdir="$PWD")"
    #REMOVE NEGATIVE VALUES FROM SUBTRACT FILE
    printf "Removing negative values from %s_subtract.bdg..." "$base"
    parallel --pipepart -q --block 5M -a Peaks/"$base"_subtract.bdg awk '$4>0' > "$temp_file" && mv "$temp_file" Peaks/"$base"_subtract.bdg
    printf "Done.\n"
else
    printf "%s_subtract.bdg already exists...skipping cleanup with macs2 bdgcmp.\n" "$base"
fi

#CONVERT SUBTRACT.BDG FILE TO BIGWIG
if [[ ! -f Peaks/"$base"_subtract.bw ]]; then
    temp_file=$(mktemp --tmpdir="$PWD")
    #SORT SUBTRACT FILE
    printf "Sorting %s_subtract.bdg for conversion to bigWig..." "$base"
    sort --parallel=$(nproc) -k1,1 -k2,2n Peaks/"$base"_subtract.bdg > "$temp_file" && mv "$temp_file" Peaks/"$base"_subtract.bdg
    printf "Done.\n"
	
    #CONVERT BEDGRAPH TO BIGWIG
    module load Kent_tools
    printf "Converting %s_subtract.bdg to bigWig via bedGraphToBigWig..." "$base"
    bedGraphToBigWig Peaks/"$base"_subtract.bdg "$chr_sizes" Peaks/"$base"_subtract.bw
    printf "Done.\n"
else
    printf "s_subtract.bw already exists...skipping conversion.\n" "$base"
fi

#COMPRESS REMAINING BEDGRAPH FILES FOR STORAGE
gzip Peaks/*.bdg

#REMOVE tmp DIRECTORY
rm -rf tmp

date

