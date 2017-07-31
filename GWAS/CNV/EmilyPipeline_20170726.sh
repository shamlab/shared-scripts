###	Date : Jul/26/2017
###	Purpose : To have a record of all commands
###	This is not a auto script. Please run the command line by line and modify it accordingly.  

hmm=/psychipc01/disk2/projects/merce_project/20141219_exome_chip/output/HumanCoreExome_v12-A_from_PennCNV/From_PennCNV/exome.hmm
pfb=/psychipc01/disk2/projects/merce_project/20141219_exome_chip/output/HumanCoreExome_v12-A_from_PennCNV/From_PennCNV/humancoreexome-12v1-1_a.pfb
sample=/home2/groups/paultam/Jacob/CNV/output_20170726/GenomeStudio_export
dir=/home2/groups/paultam/Jacob/CNV/output_20170726/
sample_list=$sample/sample_list.txt
plink=/home/jacobhsu/software/plink-1.07-x86_64/plink
#plink=/home/jacobhsu/software/plink-1.9-beta-x86_64_20150613/plink

####	Split by sample
####	Preparing input signal intensity files from Illumina intensity file (if there are more :than 1 sample**)
cd $sample
kcolumn.pl IS-141125_BalleleFreq_LogRRatio.txt split 3 -heading 3 -tab -out IS-141125 -name_by_header --beforestring .GType
kcolumn.pl IS-151208_BalleleFreq_LogRRatio.txt split 3 -heading 3 -tab -out IS-151208 -name_by_header --beforestring .GType
kcolumn.pl IS-160329_BalleleFreq_LogRRatio.txt split 3 -heading 3 -tab -out IS-160329 -name_by_header --beforestring .GType
kcolumn.pl IS-161202_BalleleFreq_LogRRatio.txt split 3 -heading 3 -tab -out IS-161202 -name_by_header --beforestring .GType

####	CNV Calling by PennCNV, call all CNV together
cd $dir
mkdir -p results
cd results
/software/cnvpipeline/penncnv/detect_cnv.pl --test -hmm $hmm -pfb $pfb --list $sample_list -log $dir/log/CLOACA_vs_70trios.log -out $dir/results/CLOACA_vs_70trios.rawcnv --coordinate_from_input  --confidence 2>&1

####    CNV filtering
cd $dir/results/AfterQC
/software/cnvpipeline/penncnv/filter_cnv.pl $dir/results/CLOACA_vs_70trios.rawcnv -qclogfile $dir/log/CLOACA_vs_70trios.log -qclrrsd 0.35 --qcwf 0.1 -qcnumcnv 1000 -length 1k -qcpassout $dir/results/AfterQC/CLOACA_vs_70trios.qcpass -qcsumout $dir/results/AfterQC/CLOACA_vs_70trios.qcsum -out $dir/results/AfterQC/CLOACA_vs_70trios.goodcnv

####	Convert the Penncnv cnv summary into tab delimited file for the input of PLINK
cd $dir/results/AfterQC
wget http://www.openbioinformatics.org/penncnv/download/penncnv_to_plink.pl
chmod 744 penncnv_to_plink.pl
perl $dir/results/AfterQC/penncnv_to_plink.pl -i $dir/results/AfterQC/CLOACA_vs_70trios.goodcnv -o $dir/results/AfterQC/CLOACA_vs_70trios.goodcnv.plink.cnv
                                               

####	Modify the file CLOACA_vs_70trios.goodcnv.plink.cnv
cd $dir/results/AfterQC
sed 's/\/home2\/groups\/paultam\/Jacob\/CNV\/output_20170726\/GenomeStudio_export\///' CLOACA_vs_70trios.goodcnv.plink.cnv > tmp.cnv 
sed 's/\/home2\/groups\/paultam\/Jacob\/CNV\/output_20170726\/GenomeStudio_export\///' tmp.cnv > tmp1.cnv
sed  's/A\tIS/\tIS/' tmp1.cnv > tmp2.cnv
sed  's/B\tIS/\tIS/' tmp2.cnv > tmp3.cnv
sed  's/C\tIS/\tIS/' tmp3.cnv > tmp4.cnv
sed  's/C_blood\t/\t/' tmp4.cnv > tmp5.cnv
sed  's/C_tissue\t/\t/' tmp5.cnv > tmp6.cnv
sed  's/lood//' tmp6.cnv > tmp7.cnv
sed  's/issue//' tmp7.cnv > tmp8.cnv
sed  's/B_new\t/\t/' tmp8.cnv > tmp9.cnv
sed  's/C_new\t/\t/' tmp9.cnv > tmp10.cnv
sed 's/IS-141125\.//' tmp10.cnv > tmp11.cnv
sed 's/IS-151208\.//' tmp11.cnv > tmp12.cnv
sed 's/IS-160329\.//' tmp12.cnv > tmp13.cnv
sed 's/IS-161202\.//' tmp13.cnv > tmp14.cnv
sed  's/B_new\t/B\t/' tmp14.cnv > tmp15.cnv
sed  's/C_new\t/C\t/' tmp15.cnv > tmp16.cnv
sed 's/IS-141125\.//' tmp16.cnv > tmp17.cnv
sed 's/IS-151208\.//' tmp17.cnv > tmp18.cnv
sed 's/IS-160329\.//' tmp18.cnv > tmp19.cnv
sed 's/IS-161202\.//' tmp19.cnv > tmp20.cnv
cp tmp20.cnv CLOACA_vs_70trios.goodcnv.plink.cnv
rm tmp*

####	Create or copy a fam fille
cd $dir/results/
mkdir -p Plink
cd Plink
cp $dir/Phenotype/CNV_calling_70trios.ped  ./CLOACA_vs_70trios.goodcnv.plink.fam
cp ../AfterQC/CLOACA_vs_70trios.goodcnv.plink.cnv ./


#--------------------------------------------------------
####	Modify the disease status of non-CLOACA patient |
#--------------------------------------------------------


####	make .map for the subsequent analysis
$plink --noweb --cfile CLOACA_vs_70trios.goodcnv.plink --cnv-make-map --out CLOACA_vs_70trios.goodcnv.plink
####	load the CNV data (see PLINK website)
$plink --noweb --cfile CLOACA_vs_70trios.goodcnv.plink --allow-no-sex --out CLOACA_vs_70trios.goodcnv.plink
####	checking for the overlapping CNV within an individual for sanity check of cnv file
$plink --noweb --cfile CLOACA_vs_70trios.goodcnv.plink --cnv-check-no-overlap --allow-no-sex --out CLOACA_vs_70trios.goodcnv.plink_QC
####	To perform a set of global test of CNV burden in cases versus controls ( To drop individuals from the file who do not have at least one segment after filtering)
$plink --noweb --cfile CLOACA_vs_70trios.goodcnv.plink --cnv-drop-no-segment --cnv-indiv-perm --mperm 1000000 --allow-no-sex --out CLOACA_vs_70trios.goodcnv.plink_globalCNVburden


##########  Unique CNV
####    To identify the CNV that are unique to case or controls
$plink --noweb --cfile CLOACA_vs_70trios.goodcnv.plink --cnv-unique --cnv-write --out unique_cnv
####    make .map for the unique CNV
$plink --noweb --cnv-list unique_cnv.cnv --cnv-make-map --out unique_cnv
####    write summary files for the unique CNV
$plink --noweb --map unique_cnv.cnv.map --fam unique_cnv.fam --cnv-list unique_cnv.cnv --allow-no-sex --out unique_cnv
####    Unique segment: To perform a set of global test of CNV burden in cases versus controls
$plink --noweb --map unique_cnv.cnv.map --fam unique_cnv.fam --cnv-list unique_cnv.cnv --cnv-drop-no-segment --cnv-indiv-perm --mperm 1000000 --allow-no-sex --out unique_cnv_globalCNVburden
####    write summary files for unique CNV (only in Cases)
$plink --noweb --map unique_cnv.cnv.map --fam unique_cnv.fam --cnv-list unique_cnv.cnv --allow-no-sex --cnv-drop-no-segment --filter-cases --cnv-write --out unique_cnv_cases
$plink --noweb --cnv-list unique_cnv_cases.cnv --cnv-make-map --out unique_cnv_cases
$plink --noweb --map unique_cnv_cases.cnv.map --fam unique_cnv_cases.fam --cnv-list unique_cnv_cases.cnv --allow-no-sex --out unique_cnv_cases
####    check if it intersects with any gene regions
wget https://www.cog-genomics.org/static/bin/plink/glist-hg19
mv glist-hg19 glist_hg19.txt
$plink --noweb --cfile unique_cnv_cases --cnv-verbose-report-regions --cnv-intersect glist_hg19.txt --out uniqueCase_gene &
$plink --noweb --cfile unique_cnv_cases --cnv-verbose-report-regions --cnv-intersect glist_hg19.txt --cnv-border 20 --out uniqueCase_gene_20kbaround &
$plink --noweb --map unique_cnv_cases.cnv.map --fam unique_cnv_cases.fam --cnv-list unique_cnv_cases.cnv --allow-no-sex --cnv-drop-no-segment --cnv-freq-exclude-below 2  --cnv-write --out uniqueCase_gene_20kbaround_moreThan1seg
$plink --noweb --cnv-list uniqueCase_gene_20kbaround_moreThan1seg.cnv --cnv-make-map --out uniqueCase_gene_20kbaround_moreThan1seg
$plink --noweb --cfile uniqueCase_gene_20kbaround_moreThan1seg --cnv-verbose-report-regions --cnv-intersect glist_hg19.txt --cnv-border 20 --out uniqueCase_gene_20kbaround_moreThan1seg &

