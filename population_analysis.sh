####################################################
# Create a directory, and prepare some path variables for subsequent calls
####################################################
workdir=/work/pop_cucu_demo  
refdir=$workdir/ref
datadir=$workdir/data
scriptdir=$workdir/scripts
export PATH=$scriptdir:$PATH   

GFF=$refdir/Cucumber_v2.gff3
REF=$refdir/Cucumber_v2.chr.fa
FAI=$refdir/Cucumber_v2.chr.fa.fai
GROUP=$datadir/pop_group.txt

####################################################################
#data filtering
#####################################################################
cd $workdir  
mkdir 00.filter
cd 00.filter

#filtration：vcftools
#--max-missing Exclude sites on the basis of the proportion of missing data 
#(defined to be between 0 and 1, where 0 allows sites that are completely missing 
#and 1 indicates no missing data allowed).

vcftools --gzvcf all.varFilter.vcf.gz --recode --recode-INFO-all --stdout     --maf 0.05  --max-missing 0.8  --minDP 2  --maxDP 1000      \
    --minQ 30 --minGQ 0 --min-alleles 2  --max-alleles 2 --remove-indels |gzip - > clean.vcf.gz  

####################################################################
#phylogenetic analysis
#####################################################################

cd $workdir  
mkdir 01.phylo_tree
cd 01.phylo_tree
#file format conversion
run_pipeline.pl  -Xmx5G -importGuess  $workdir/00.filter/clean.vcf.gz  \
    -ExportPlugin -saveAs supergene.phy -format Phylip_Inter

#Maximum likelihood method to construct evolutionary tree
#method 1：fasttree
fasttree -nt -gtr  supergene.phy   >  fasttree.nwk
#method 2：iqtree
iqtree2 -s supergene.phy -st DNA -T 2  -mem 8G \
    -m  GTR  -redo \
    -B 1000 -bnni \
    --prefix iqtree 
    

#Evolution tree drawing
ggtree.r -t fasttree.nwk -f $GROUP -g Group -n fasttree
ggtree.r -t fasttree.nwk -f $GROUP -g Group -n fasttree.unrooted -l unrooted

ggtree.r -t iqtree.treefile -f $GROUP -g Group -n iqtree
ggtree.r -t iqtree.treefile -f $GROUP -g Group -n iqtree.unrooted -l unrooted

########################################################
#PCA
#########################################################
cd $workdir  
mkdir 02.PCA
cd 02.PCA

## plink: PCA analysis
plink --vcf  $workdir/00.filter/clean.vcf.gz --pca 10 --out  plink_pca   \
    --allow-extra-chr --set-missing-var-ids @:#    --vcf-half-call missing

#visualize
pca_plink_plot.r -i plink_pca.eigenvec -f $GROUP -g Group --name plink_pca

#########################################################
#STRUCTURE analysis
######################################################
cd $workdir 
mkdir 03.STRUCTURE
cd 03.STRUCTURE

## filter LD 
#50 10 0.2   50 SNP Windows step 10 r2 0.2
plink --vcf  $workdir/00.filter/clean.vcf.gz  --indep-pairwise 50 10 0.2 --out ld   \
    --allow-extra-chr --set-missing-var-ids @:# 
plink --vcf  $workdir/00.filter/clean.vcf.gz  --make-bed --extract ld.prune.in  \
    --out LDfiltered --recode vcf-iid  --keep-allele-order  --allow-extra-chr --set-missing-var-ids @:#  

#Convert to the format of plink
vcftools --vcf LDfiltered.vcf --plink \
    --out plink
#Convert to the bed format required by admixture
plink --noweb --file plink  --recode12 --out admixture \
     --allow-extra-chr  --keep-allele-order

#admixture Population structure analysis
for k in {2..10};do
    admixture -j2 -C 0.01 --cv admixture.ped $k >admixture.log$k.out
done

#visualize 
structure_plot.r  -d ./ -s admixture.nosex 
structure_plot.r  -d ./ -s admixture.nosex -f $GROUP -g Group  #按照分组顺序显示

#Determine the optimal K, and the K corresponding to the smallest CV value is the optimal K
grep "CV error" *out


##########################################################
#linkage disequilibrium analysis LDdecay
########################################################

cd $workdir 
mkdir 04.LDdecay
cd 04.LDdecay

#Group: Split according to its group name
cat $GROUP |grep EastAsian |cut -f 1 >EastAsian_popid.txt
cat $GROUP |grep Eurasian |cut -f 1 >Eurasian_popid.txt
cat $GROUP |grep Xishuangbanna |cut -f 1 >Xishuangbanna_popid.txt
cat $GROUP |grep Indian |cut -f 1 >Indian_popid.txt


for i in EastAsian Eurasian Xishuangbanna  Indian;do 

    PopLDdecay -InVCF  $workdir/00.filter/clean.vcf.gz \
        -SubPop  ${i}_popid.txt -MaxDist 500 -OutStat ${i}.stat
    Plot_OnePop.pl -inFile ${i}.stat.gz -output ${i}.ld
done 

#Multiple groups are placed in a single graph
echo EastAsian.stat.gz EastAsian >ld_stat.list
echo Eurasian.stat.gz Eurasian >>ld_stat.list
echo Xishuangbanna.stat.gz Xishuangbanna >>ld_stat.list
echo Indian.stat.gz Indian >>ld_stat.list
Plot_MultiPop.pl -inList ld_stat.list -output ld_stat.multi


######################################################
#selection analysis
########################################################
cd $workdir 
mkdir 05.select_sweep
cd 05.select_sweep

#The samples from different groups were separated into different files
cat $GROUP |grep wild |cut -f 1 >wild_popid.txt
cat $GROUP |grep cultivated|grep -v "Indian"|grep -v "Xishuangbanna" |cut -f 1 >cultivated_popid.txt

#fst  pi  tajimaD analysis
mkdir fst_pi_tajimaD
cd fst_pi_tajimaD

#Set the input vcf file with the calculated window and step size
gzvcf=$workdir/00.filter/clean.vcf.gz
window=100000
step=10000

#pi
vcftools  --gzvcf $gzvcf \
    --window-pi $window --window-pi-step  $step  \
    --keep ../wild_popid.txt   --out pi.wild
vcftools  --gzvcf $gzvcf \
    --window-pi $window --window-pi-step  $step  \
    --keep ../cultivated_popid.txt  --out pi.cultivated

#pi:visualization
pi_manhattan_plot.r -i pi.wild.windowed.pi -F $FAI -f 19226500 -n pi.wild
pi_smooth_line_plot.r -i pi.wild.windowed.pi -F $FAI -f 19226500  -n pi.wild.smoothline

#Fst
vcftools  --gzvcf $gzvcf --fst-window-size $window --fst-window-step $step  \
    --weir-fst-pop  ../wild_popid.txt --weir-fst-pop ../cultivated_popid.txt --out  Fst.wild.cultivated

#fst:visualization
fst_manhattan_plot.r -i Fst.wild.cultivated.windowed.weir.fst -F $FAI -f 19226500 -n Fst.wild.cultivated
fst_manhattan_plot.r -i Fst.wild.cultivated.windowed.weir.fst -F $FAI -f 19226500  -n Fst.wild.cultivated_vline --vline
fst_smooth_line_plot.r -i Fst.wild.cultivated.windowed.weir.fst -F $FAI -f 19226500  -n Fst.wild.cultivated_smoothline


#Tajima's D 
vcftools --gzvcf $gzvcf --TajimaD  100000  --keep  ../wild_popid.txt  --out wild
vcftools --gzvcf $gzvcf --TajimaD  100000  --keep  ../cultivated_popid.txt  --out cultivated

tajimaD_manhattan_plot.r -i cultivated.Tajima.D -F $FAI -f 19226500  -n TajimaD.cultivated 
tajimaD_smooth_line_plot.r -i cultivated.Tajima.D -F $FAI -f 19226500  -n TajimaD.cultivated.smoothline

tajimaD_manhattan_plot.r -i wild.Tajima.D -F $FAI -f 19226500  -n TajimaD.wild 
tajimaD_smooth_line_plot.r -i wild.Tajima.D -F $FAI -f 19226500  -n TajimaD.wild.smoothline
