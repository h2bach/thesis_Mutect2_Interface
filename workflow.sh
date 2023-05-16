dir = "."
ref = "references"
input_fastq = ""
reads = "data/reads"
aligned_reads = "data/aligned_reads"

# -------------------
# STEP 0: sra-tools
# -------------------


# -------------------
# STEP 0: Seqtk
# -------------------
seqtk sample ${reads}/run_forward.fastq $rate > ${reads}/run_forward.fastq
seqtk sample ${reads}/run_backward.fastq $rate > ${reads}/run_backward.fastq

# -------------------
# STEP 1: QC - fastqc
# -------------------
fastqc {seqfile1} {seqfile2} {seqfileN}

# -------------------
# STEP 1.1: Trimmomatic (optional)
# Using PE param for paired-end
# Using phred33 as default for usual PE-seq
# -------------------

trimmomatic PE -threads ${num_thread} -phred33 ${reads}/run_forward.fastq ${reads}/run_backward.fastq ${reads}/output_forward.fastq ${reads}/output_backward.fastq 


# -------------------
# STEP 2: BWA - MEM
# -------------------
bwa mem -t ${num_thread} -R ${read_group_info} ${ref} ${reads}/run_forward.fastq.fastq ${reads}/run_backward.fastq ${aligned_reads}/run.paired.bam
# e.g.: bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

############################
# Docker for GATK from here
############################
# docker run -v ~/data:/gatk/data -it broadinstitute/gatk:4.4.0.0

# -------------------
# STEP 3: Mark Duplicates and Sort - GATK4 
# In GATK 3, usually used Picard in this step
# -------------------
gatk MarkDuplicatesSpark -I ${aligned_reads}/run.paired.bam -O ${aligned_reads}/run.paired_sorted_dedup_reads.bam

# -------------------
# STEP 4: BAMSurgeon 
# -------------------
python3 -O ../scripts/check_dependencies.py
python3 -O ../scripts/randomsites.py 
python3 -O ../bin/addsnv.py -v ../test_data/random_snvs.txt -f ${aligned_reads}/run.paired_sorted_dedup_reads.bam -r $REF -o ${aligned_reads}/run.paired_sorted_dedup_reads_bs.bam -n 5 --picardjar $1 --aligner mem --seed 1234


# -------------------
# STEP 5: BQSR
# -------------------
# 1. build the model
gatk BaseRecalibrator -I ${aligned_reads}/run.paired_sorted_dedup_reads_bs.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/run.paired_sorted_dedup_reads_bs.bam -R ${ref} --bqsr-recal-file {$data}/recal_data.table -O ${aligned_reads}/run.paired_sorted_dedup_bqsr_reads_bs.bam


# -------------------
# STEP 6.1: Create PON (optional)
# -------------------
# # Tumor-only Mutect2 calls on normals
# gatk Mutect2 -R reference.fasta -I normal1.bam -max-mnp-distance 0 -O normal1.vcf.gz
# gatk Mutect2 -R reference.fasta -I normal2.bam -max-mnp-distance 0 -O normal2.vcf.gz
# gatk Mutect2 -R reference.fasta -I normal3.bam -max-mnp-distance 0 -O normal3.vcf.gz
# gatk Mutect2 -R reference.fasta -I normal4.bam -max-mnp-distance 0 -O normal4.vcf.gz
# # ............
# gatk Mutect2 -R reference.fasta -I normal40.bam -max-mnp-distance 0 -O normal40.vcf.gz
# # GenomicsDBImport
# gatk GenomicsDBImport -R reference.fasta -L intervals.interval_list \
# --genomicsdb-workspace-path pon_db \
# -V normal1.vcf.gz \
# -V normal2.vcf.gz \
# # ...
# -V normal40.vcf.gz
# # CreateSomaticPanelOfNormals
# gatk CreateSomaticPanelOfNormals -R reference.fasta -V gendb://pon_db -O pon.vcf.gz

# -------------------
# STEP 6: Mutect2
# -------------------

gatk Mutect2 /path/to/reference_genome.fa \
    -L /path/to/supporting_files/intervals.interval_list \
    -I /path/to/tumor.bam \
    -germline-resource /path/to/supporting_files/af-only-gnomad.vcf \
    -pon /path/to/supporting_files/panel_of_normals.vcf \
    --f1r2-tar-gz /path/to/data/f1r2.tar.gz \
    -O /path/to/results/unfiltered_variants.vcf


# -------------------
# STEP 7.1: LearnReadOrientationModel
# -------------------
gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz

# -------------------
# STEP 7.2: CalculateContamination
# -------------------
gatk GetPileupSummaries \
    -I tumor.bam \
    -V chr17_small_exac_common_3_grch38.vcf.gz \
    -L chr17_small_exac_common_3_grch38.vcf.gz \
    -O getpileupsummaries.table

gatk CalculateContamination \
    -I getpileupsummaries.table \
    -tumor-segmentation segments.table \ 
    -O calculatecontamination.table

# -------------------
# STEP 7: FilterMutectCalls
# -------------------

gatk FilterMutectCalls -V unfiltered.vcf \
    --tumor-segmentation segments.table] \
    --contamination-table contamination.table] \
    --ob-priors read-orientation-model.tar.gz \
    -O filtered.vcf

# -------------------
# STEP 8: Funcotator & ANNOVAR
# -------------------

############# FUNCOTATOR #############

gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download

gatk Funcotator \
    --variant /path/to/variants.vcf \
    --reference /path/to/reference_genome.fa \
    --ref-version hg38 \
    --data-sources-path /path/to/funcotator_dataSources \
    --output /path/to/result/variants.funcotated.vcf \
    --output-file-format VCF

######################################

############# ANNOVAR #############

# Tai ve du lieu chu thich
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ${protocol} humandb/

# Chu thich dua tren du lieu
./table_annovar.pl /path/to/filtered_variants.vcf humandb/ -buildver ref-buidver \
    -out /path/to/results/var_anno_folder/anno_prefix \
    --thread any-number-of-available-threads \
    -remove -protocol {names_of_protocols} \
    -operation g|r|f -nastring . -vcfinput -polish
    
######################################

!zgrep chr 17 /home/jupyter-user/sandbox/12_somatic_oncefiltered_funcotate.vcf.gz | zgrep 7674220