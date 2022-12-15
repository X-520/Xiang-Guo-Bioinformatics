####################################################
# Prepare some path variables for subsequent calls
####################################################

workdir=/work/rnaseq_demo  #Setting a Working Path
refdir=$workdir/ref
datadir=$workdir/data
scriptdir=$workdir/scripts

###################################################################
# Download reference genome data and build on the genome: HISAT index
###################################################################
mkdir $refdir
cd $refdir  

#hisat2 Provides common species genome index files：

#wget -c ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
#wget -c ftp://ftp.ensembl.org/pub/release-99/gff3/homo_sapiens/Homo_sapiens.GRCh38.99.chromosome.22.gff3.gz
#gunzip *gz

#Build index common species index download： http://daehwankimlab.github.io/hisat2/download/

sh $scriptdir/index.sh Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.99.chromosome.22.gff3

#Reference gene related file variables are set to facilitate subsequent use

REF_INDEX=$refdir/Homo_sapiens.GRCh38.dna.chromosome.22
GFF=$refdir/Homo_sapiens.GRCh38.99.chromosome.22.gff3
GTF=$refdir/Homo_sapiens.GRCh38.99.chromosome.22.gtf
GENE_BED=$refdir/gene.bed
GENE_LENGTH=$refdir/gene_length.txt



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
# If there are many samples, in practice the for loop is given a list of sample names

for i in normal_rep1 normal_rep2 normal_rep3 tumor_rep1 tumor_rep2 tumor_rep3; do 
echo "RUN CMD: fastp --thread 1 --qualified_quality_phred 10 \
--unqualified_percent_limit 50 \
--n_base_limit 10 \
-i $datadir/${i}_r1.fastq.gz \
-I $datadir/${i}_r2.fastq.gz \
-o ${i}_1.clean.fq.gz \
-O ${i}_2.clean.fq.gz \
--adapter_fasta $workdir/data/illumina_multiplex.fa \
-h ${i}.html -j ${i}.json"

fastp --thread 1 --qualified_quality_phred 10 \
--unqualified_percent_limit 50 \
--n_base_limit 10 \
-i $datadir/${i}_r1.fastq.gz \
-I $datadir/${i}_r2.fastq.gz \
-o ${i}_1.clean.fq.gz \
-O ${i}_2.clean.fq.gz \
--adapter_fasta $workdir/data/illumina_multiplex.fa \
-h ${i}.html -j ${i}.json	
done

#Quality control data statistical summary：
python $scriptdir/qc_stat.py -d $workdir/2.data_qc/ -o $workdir/2.data_qc/ -p all_sample_qc

###################################################################
# Match the reads to the genome
##################################################################
cd $workdir  
mkdir -p $workdir/3.map/hisat2
cd $workdir/3.map/hisat2


for i in normal_rep1 normal_rep2 normal_rep3 tumor_rep1 tumor_rep2 tumor_rep3; do
echo "RUN CMD: hisat2 -p 1 --rg-id=${i} --rg SM:${i} --rg LB:${i} --rg PL:ILLUMINA \
-x $REF_INDEX --dta --rna-strandness RF \
-1 $workdir/2.data_qc/${i}_1.clean.fq.gz \
-2 $workdir/2.data_qc/${i}_2.clean.fq.gz \
-S ${i}.sam 2>${i}.summary"

hisat2 -p 1 --rg-id=${i} --rg SM:${i} --rg LB:${i} --rg PL:ILLUMINA \
-x $REF_INDEX --dta --rna-strandness RF \
-1 $workdir/2.data_qc/${i}_1.clean.fq.gz \
-2 $workdir/2.data_qc/${i}_2.clean.fq.gz \
-S ${i}.sam 2>${i}.summary 
done

# sam2bam : Format and sort
for i in normal_rep1 normal_rep2 normal_rep3 tumor_rep1 tumor_rep2 tumor_rep3; do
echo "RUN CMD: samtools sort --threads 1 -m 3G -o ${i}.bam ${i}.sam"
samtools sort  --threads 1 -m 3G -o ${i}.bam ${i}.sam
done

# bam index
for i in normal_rep1 normal_rep2 normal_rep3 tumor_rep1 tumor_rep2 tumor_rep3; do
echo "RUN CMD: samtools index ${i}.bam"
samtools index ${i}.bam
done

######################################################################
# The map results were analyzed by QC
######################################################################
## RSeQC was used to compare the result files for quality control analysis
cd $workdir  
mkdir -p $workdir/3.map/map_QC
cd $workdir/3.map/map_QC

#1.Fragment inner size, fragment selection is abnormal

for i in normal_rep1 normal_rep2 normal_rep3 tumor_rep1 tumor_rep2 tumor_rep3; do
echo "RUN CMD: inner_distance.py -i $workdir/3.map/hisat2/${i}.bam  -r $GENE_BED  -o ${i}_inner_size"
inner_distance.py -i $workdir/3.map/hisat2/${i}.bam  -r $GENE_BED  -o ${i}_inner_size
done

#2.Gene coverage, whether the RNA is degraded

for i in normal_rep1 normal_rep2 normal_rep3 tumor_rep1 tumor_rep2 tumor_rep3; do
echo "RUN CMD: geneBody_coverage.py -r $GENE_BED -i $workdir/3.map/hisat2/${i}.bam  -o ${i}.genebody"	
geneBody_coverage.py -r $GENE_BED -i $workdir/3.map/hisat2/${i}.bam  -o ${i}.genebody
done


###############################################################
# Gene expression quantification and result presentation
## The known gene expression was quantified by Htseq-count
##################################################################
cd $workdir/
mkdir 4.expression
cd 4.expression

#--stranded yes or no   Library type setting, whether it is chain specific library

for i in normal_rep1 normal_rep2 normal_rep3 tumor_rep1 tumor_rep2 tumor_rep3; do
echo "RUN CMD: htseq-count --format bam --order pos --mode intersection-strict \
--stranded reverse --minaqual 1 --type exon \
--idattr gene_id $workdir/3.map/hisat2/${i}.bam $GTF > ${i}_gene.tsv"
		
htseq-count --format bam --order pos --mode intersection-strict \
--stranded yes --minaqual 1 --type exon \
--idattr gene_id $workdir/3.map/hisat2/${i}.bam $GTF > ${i}_gene.tsv
done

### The quantitative results of different samples were incorporated to facilitate subsequent gene expression analysis
python $scriptdir/merge_gene_count.py -p all_gene_count \
-f  normal_rep1_gene.tsv -l normal_rep1 \
-f  normal_rep2_gene.tsv -l normal_rep2 \
-f  normal_rep3_gene.tsv -l normal_rep3 \
-f  tumor_rep1_gene.tsv -l tumor_rep1 \
-f  tumor_rep2_gene.tsv -l tumor_rep2 \
-f  tumor_rep3_gene.tsv -l tumor_rep3


###Presentation of quantitative results of gene expression
#1.The expression density of each sample
#2.box distribution of expression quantity of each sample
#4 Each sample expression correlation analysis heat map and cluster map

/usr/bin/Rscript $scriptdir/fpkm_and_plot.R -i all_gene_count.tsv  -l $GENE_LENGTH  -o ./


##################################################
## Gene differential expression analysis (DESeq2), and mapping volcano map and MA map, as well as differential gene expression heat map
###################################################
cd $workdir/
mkdir 5.deg
cd 5.deg
#The grouping file normal_vs_tumor.pare. txt was created to search for differentially expressed genes among different groups
#ID	group
#normal_rep1	normal
#normal_rep2	normal
#normal_rep3	normal
#tumor_rep1	tumor
#tumor_rep2	tumor
#tumor_rep3	tumor

/usr/bin/Rscript $scriptdir/deseq_analysis.r -i $workdir/4.expression/all_gene_count.tsv -g normal_vs_tumor.compare.txt  -k $workdir/4.expression/all_gene_fpkm.tsv -r normal --fdr 1 --fc 1.1 -p normal_vs_tumor

#Heat map of differentially expressed genes
/usr/bin/Rscript $scriptdir/heatmap.r -i $workdir/4.expression/all_gene_fpkm.tsv -l normal_vs_tumor.DEG.final.tsv -p normal_vs_tumor.deg_gene_heatmap -o ./

##################################################
#Functional enrichment analysis of differentially expressed genes
##################################################
cd $workdir/
mkdir 6.enrich
cd 6.enrich


/usr/bin/Rscript $scriptdir/enrichGO_pip.R --deg.file $workdir/5.deg/normal_vs_tumor.DEG.final.tsv -o GO/ -n normal_vs_tumor --pvalueCutoff 0.5 --ann.db org.Hs.eg.db  --idtype ENSEMBL --totype ENTREZID 

/usr/bin/Rscript $scriptdir/enrichKEGG_pip.R --deg.file $workdir/5.deg/normal_vs_tumor.DEG.final.tsv -o KEGG -n normal_vs_tumor --pvalueCutoff 1 --ann.db org.Hs.eg.db --organism hsa --idtype ENSEMBL --totype ENTREZID 

#GSEA analysis
# download GSEA software:http://software.broadinstitute.org/gsea/downloads.jsp
GSEA software