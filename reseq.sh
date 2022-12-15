####################################################
# Create a directory, and prepare some path variables for subsequent calls
####################################################
workdir=/work/my_reseq  #Setting a Working Path
refdir=$workdir/ref
datadir=$workdir/data
scriptdir=$workdir/scripts
tmpdir=$workdir/tmp 

###################################################################
# Download reference genome data and build index on the genome
###################################################################
cd $refdir 
wget -c ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa.gz
wget -c ftp://ftp.ensemblgenomes.org/pub/plants/release-47/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.47.chromosome.4.gff3.gz
gunzip *gz
sh $scriptdir/index.sh Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa Arabidopsis_thaliana.TAIR10.47.chromosome.4.gff3

#Reference gene related file variables are set to facilitate subsequent use
REF=$refdir/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa
BWA_INDEX=$refdir/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa
GFF=$refdir/Arabidopsis_thaliana.TAIR10.47.chromosome.4.gff3
GTF=$refdir/Arabidopsis_thaliana.TAIR10.47.chromosome.4.gtf
INDEX_FAI=$refdir/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa.fai
PICARD_DICT=$refdir/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.dict

####################################################################
#Perform fastqc quality control on raw data
#####################################################################
cd $workdir
mkdir 1.fastqc

fastqc $datadir/*.gz  -o $workdir/1.fastqc

############################################################################
# Data quality control: de-coupling of original sequences, deletion of low-quality reads, etc
############################################################################
cd $workdir 
mkdir 2.data_qc
cd  2.data_qc
#Remove the adapter using the fastp tool
#--qualified_quality_phred the quality value that a base is qualified. 
#            Default 15 means phred quality >=Q15 is qualified. (int [=15])
#--unqualified_percent_limit how many percents of bases are allowed to be unqualified
#--n_base_limit if one read's number of N base is >n_base_limit, 
#            then this read/pair is discarded 
#--detect_adapter_for_pe 

for i in p1 p2 pool1 pool2; do 
echo "RUN CMD: fastp --thread 1 --qualified_quality_phred 10 \
    --unqualified_percent_limit 50 \
    --n_base_limit 10 \
    -i $datadir/${i}_r1.fq.gz \
    -I $datadir/${i}_r2.fq.gz \
    -o ${i}_1.clean.fq.gz \
    -O ${i}_2.clean.fq.gz \
    --adapter_fasta $workdir/data/illumina_multiplex.fa \
    -h ${i}.html -j ${i}.json"

fastp --thread 1 --qualified_quality_phred 10 \
    --unqualified_percent_limit 50 \
    --n_base_limit 10 \
    -i $datadir/${i}_r1.fq.gz \
    -I $datadir/${i}_r2.fq.gz \
    -o ${i}_1.clean.fq.gz \
    -O ${i}_2.clean.fq.gz \
    --adapter_fasta $workdir/data/illumina_multiplex.fa \
    -h ${i}.html -j ${i}.json    
done

#Quality control data statistical summary：
python $scriptdir/qc_stat.py -d $workdir/2.data_qc/ -o $workdir/2.data_qc/ -p all_sample_qc

###################################################################
# reads were mapped with the genome
##################################################################
cd $workdir 
mkdir -p $workdir/3.map/bwa
cd $workdir/3.map/bwa

#map bwa 
for i in p1 p2 pool1 pool2; do
echo  "RUN CMD: bwa mem  $BWA_INDEX  $workdir/2.data_qc/${i}_1.clean.fq.gz \
        $workdir/2.data_qc/${i}_2.clean.fq.gz -t 1 -M \
        -R '@RG\tID:${i}\tLB:${i}\tPL:ILLUMINA\tSM:${i}' \
        |samtools view -bS -h - > $workdir/3.map/bwa/${i}.bam"

bwa mem  $BWA_INDEX  $workdir/2.data_qc/${i}_1.clean.fq.gz \
        $workdir/2.data_qc/${i}_2.clean.fq.gz -t 1 -M \
        -R "@RG\tID:${i}\tLB:${i}\tPL:ILLUMINA\tSM:${i}" \
        |samtools view -bS -h - > $workdir/3.map/bwa/${i}.bam
done

#bam file sorting
for i in p1 p2 pool1 pool2; do
echo "RUN CMD: picard SortSam -Xmx1g VALIDATION_STRINGENCY=LENIENT I=$workdir/3.map/bwa/${i}.bam \
    O=$workdir/3.map/bwa/${i}.sorted.bam SORT_ORDER=coordinate \
    TMP_DIR=$tmpdir"

picard SortSam -Xmx1g VALIDATION_STRINGENCY=LENIENT I=$workdir/3.map/bwa/${i}.bam \
    O=$workdir/3.map/bwa/${i}.sorted.bam SORT_ORDER=coordinate TMP_DIR=$tmpdir
done


#Remove Duplicates
cd $workdir  
mkdir -p $workdir/3.map/result
cd $workdir/3.map/result
for i in p1 p2 pool1 pool2; do
echo "RUN CMD: picard MarkDuplicates -Xmx4g MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=512 VALIDATION_STRINGENCY=LENIENT  \
    INPUT=$workdir/3.map/bwa/${i}.sorted.bam \
    OUTPUT=$workdir/3.map/result/${i}.sorted.dedup.bam \
    METRICS_FILE=$workdir/3.map/result/${i}.sorted.dedup.metrics"

picard MarkDuplicates -Xmx4g MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=512 VALIDATION_STRINGENCY=LENIENT  \
    INPUT=$workdir/3.map/bwa/${i}.sorted.bam \
    OUTPUT=$workdir/3.map/result/${i}.sorted.dedup.bam \
    METRICS_FILE=$workdir/3.map/result/${i}.sorted.dedup.metrics
done

# bam index
for i in p1 p2 pool1 pool2; do
    echo "RUN CMD: samtools index ${i}.sorted.dedup.bam"
    samtools index ${i}.sorted.dedup.bam
done

###################################################3
#GATK call SNP
####################################################
cd $workdir  
mkdir -p $workdir/4.snp_indel/GATK
cd $workdir/4.snp_indel/GATK

for i in `cat sample_list.txt`
do
	echo "gatk --java-options "-Xmx4g" HaplotypeCaller -R $REF   \
  -I $workdir/3.map/result/${i}.sorted.dedup.bam \
  -O ${i}.g.vcf.gz --max-alternate-alleles 4  --sample-ploidy 2 \
  -ERC GVCF --tmp-dir $tmpdir"
done > SNPcalling.sh

# ############Merge GVCF files##########################

# gatk --java-options "-Xmx4g" CombineGVCFs -R $REF \
    # --variant p1.g.vcf.gz \
    # --variant p2.g.vcf.gz \
    # --variant pool1.g.vcf.gz \
    # --variant pool2.g.vcf.gz \
  # -O all.g.vcf.gz --tmp-dir $tmpdir

# #将gvcf转换成过VCF
# gatk --java-options "-Xmx4g" GenotypeGVCFs  -R $REF \
  # -V all.g.vcf.gz \
  # -O all.raw.vcf.gz --tmp-dir $tmpdir


#############################################################
#Variation result quality control filter, remove low quality variation result
#############################################################
cd $workdir 
mkdir -p $workdir/4.snp_indel/var_qc
cd $workdir/4.snp_indel/var_qc


#filteration：vcftools
#--max-missing Exclude sites on the basis of the proportion of missing data 
#(defined to be between 0 and 1, where 0 allows sites that are completely missing 
#and 1 indicates no missing data allowed).

vcftools --gzvcf all.varFilter.vcf.gz --recode --recode-INFO-all --stdout \
    --maf 0.05  --max-missing 0.7  --minDP 4  --maxDP 1000  \
    --minQ 30 --minGQ 80 --min-alleles 2  --max-alleles 2 |gzip - > all.clean.vcf.gz

#Separate indel and SNP to different files
vcftools --remove-indels --recode --recode-INFO-all --gzvcf all.clean.vcf.gz --stdout |gzip - >all.clean.snp.vcf.gz
vcftools --keep-only-indels  --recode --recode-INFO-all --gzvcf all.clean.vcf.gz --stdout |gzip - >all.clean.indel.vcf.gz

#################################################################
#Annotation of variation results: ANNOVAR
#################################################################
cd $workdir 
mkdir -p $workdir/5.var_ann/
cd $workdir/5.var_ann/

table_annovar.pl $workdir/4.snp_indel/var_qc/all.clean.snp.vcf.gz $refdir  \
    -buildver unknown -out $workdir/5.var_ann/snp  \
    -remove -protocol refGene -operation g -nastring . -vcfinput

table_annovar.pl $workdir/4.snp_indel/var_qc/all.clean.indel.vcf.gz $refdir  \
    -buildver unknown -out $workdir/5.var_ann/indel  \
    -remove -protocol refGene -operation g -nastring . -vcfinput
