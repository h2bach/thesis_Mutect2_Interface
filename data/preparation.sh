###########################
# references
wget -P ${ref}/ $ref_link
gunzip ${ref}/ref.fa.gz
samtools faidx ${ref}/ref.fa
gatk CreateSequenceDictionary -R ${ref}/ref.fa -O ${ref}/ref.dict

# known sites data for BQSR
wget -P ${ref} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P ${ref} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx


###########################
# Picard
wget https://github.com/broadinstitute/picard/releases/download/2.27.3/picard.jar

#sra-toolkit
prefetch {sra_codes}
fasterq-dump {sra_codes}.sra
