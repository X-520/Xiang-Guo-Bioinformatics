#ChIP-Seq
conda create -n chip seqkit bowtie2 macs2 deeptools homer samtools  python=3 
conda activate chip

cd /mnt/d/online_course/ChIP

mkdir 1.mapping
mkdir 2.callPeak
mkdir 3.PeakAnno
mkdir 4.motif

#1, Alignment: bowtie2, Ensembl
#Build a databas
bowtie2-build --threads 4 data/genome.fa  data/genome

#Alignment
bowtie2 -p 4 -x data/genome -1 data/treat.1.fq   -2  data/treat.2.fq|samtools view -bS |samtools sort  >1.mapping/treat.bam
bowtie2 -p 4 -x data/genome -1 data/control.1.fq   -2   data/control.2.fq |samtools view -bS |samtools sort  >1.mapping/control.bam

#2, PeakCalling
#The effective size of genome was calculated and the non-N base length of genome fasta was removed
Base=`cat  data/genome.fa|grep -v ">"|sed '1i\>seq'|seqkit fx2tab -B ATCG  -n -l|awk '{print $2*$3/100}'`

#MACS2
macs2 callpeak -f  BAMPE -c  1.mapping/control.bam -t  1.mapping/treat.bam -q 0.05 -g $Base -n 2.callPeak/zmays --nomodel 2> macs2.log


#3, peak annotation to obtain the annotated genes in chip enrichment results
#Rstudio:chipAnno.R analysis


#4, motif analysis
cd /mnt/d/online_course/ChIP
#Homer
findMotifsGenome.pl 2.callPeak/zmays_summits.bed  data/genome.fa   4.motif

#5, visualization（IGV）
