# 2015.4.16 # 
#Here is the demo analysis pipeline applying Asmvar to the 37 human de novo assemblies.
#Because the individual genotypes for the 10 trio de novo assemblies from the GenomeDK consortium requires further access authorization and to reduce the size for the software packge , we only present the result for chromosome 20 of a widely used cell line NA12878 in the output directory. 
#For more information, access to the full result of the 37 human de novo assemblies, technical support, invitation to project analysis, suggestions etc , please feel free to contact liusiyang@genomics.cn and huangshujia@genomics.cn 
#In the following text, the exectuable programes and the directories start with captials while the other individual files are in small letters.


##############################################################################
# 0. Last Alignment- Preparation for AsmVar
##############################################################################
#We suggest you apply LAST aligner for assembly-vs-assembly alignment http://last.cbrc.jp/ 
#Below is the best practice kindly advised by the author to align human assemblies
#Assembly-vs-reference alignment
ref= `wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz`
#Human reference Fasta File eg GRCh37, GRCh38 or hs37d5
DataDir=../../Data
OutDir=../../GigaDB
for sam in `cat $DataDir/DemoSamples.list`
do
	lastal -e25 -v -q3 -j4 $ref.lastdb $asm.fa | last-split -s35 -v |gzip >$DataDir/$sam.maf.gz
done

##############################################################################
# 1. AsmVar Detection
##############################################################################
AsmVar=../../../src/AsmvarDetect/ASV_VariantDetector 
echo "Step1: Running Detection\n";
for sam in `cat $DataDir/DemoSamples.list`
do
	for chr in 1..22,X,Y
	do
		#We suggest that you use parameter "-r" to do the analysis by separate chromosome.
		$AsmVar -s $sam -i $DataDir/$sam.maf.gz -t $ref -q $DataDir/$sam.fa.gz -o $OutDir/$sam -r $chr > $OutDir/$sam.age 2> $OutDir/$sam.AsmVarDetection.log	
	done
done

#After you finish the AsmVar Detection, you will obtain the raw variants vcf file for each assembly. For example, "NA12878.APLG.20.vcf" below indicate the raw vcf for chr20 of individual assembly NA12878.APLG. $OutDir/$sam in the above command is the prefix of the raw vcf file. 
#Should you have more than one sample, you would need to apply GATK_CombineVariants to merge vcf from multiple samples into one multi-sample vcf file as the input file for the next Altalignment step


##############################################################################
# 2. AsmVar Altalignment
##############################################################################
echo "Step2: AsmVar Altalignment\n";
ProgAltAlign=../../../src/AsmvarAlterAlign/ASV_AlterAlign.py
bam= `wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.high_coverage.20100311.bam` 
#This is a bam file example for one of the individual NA12878, you will have to generate your own bam file for your samples. Bam files of each samples will be applied in the alternative alignment procedue to generate the genotype likelihood
rawvcf= $OutDir/NA12878.APLG.20.vcf 

#To maximize efficiency, AsmVar does alternative alignment by individual for difference alleles.
for sam in `cat $DataDir/DemoSamples.list`
do
	for chr in 1..22,X,Y
	do
		#Here, we also suggest that you use "-c" to specify the chromosome
		python $ProgAltAlign -r $ref -v $rawvcf -b $bam -q 30 -s $sam -o $OutDir/$sam.AltAlign -c $chr 2> $OutDir/$sam.AltAlign.log
	done
done

#After you finish the Altalignment, you will obtain the vcf file with the read alignment information, e.g. "$OutDir/NA12878.APLG.AltAlign.20.vcf.gz". The following field will be added to the vcf file
# ##FORMAT=<ID=AA,Number=4,Type=Integer,Description="Information of Alternate Align. Format: Ref_perfect,Alt_Perfect,Both_Perfect,Both_Imperfect">
#Should you have more than one sample, you would need to merge all the altalign vcf files as well

##############################################################################
# 3. AsmVar Genotyping
##############################################################################
echo "Step3: AsmVar Genotyping\n";
altvcf=$OutDir/NA12878.APLG.AltAlign.20.vcf #One merge alt vcf file
Genoyping=../../../src/AsmvarGenotype/GMM/SVGenotyping.py
ped= #A ped-format file including your family information
for sam in `cat $DataDir/DemoSamples.list`
do
	#Also, we suggest you run by each chromosomes using "-c"
	python $Genoyping -p $ped -f $OutDir/$sam.Genotyping.Figure -o $OutDir/$sam.Genotyping -c 20 $altvcf 2> $OutDir/$sam.Genotyping.log
done

#After you finish the genotyping, you will obtain the vcf file with the genotype information including the genotype, genotype likelihood, genotype quality e.g. "$OutDir/NA12878.APLG.Genotyping.20.vcf". The following fields will be added to the vcf
# ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality. -10*log10(1-p), p=w*N/(sigma(W*N)) N is gaussian density (p: The posterior probability of Genotype call is correct)">
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
#Before you run the Genotype recalibration, you have to "tabix" your vcf file.(http://sourceforge.net/projects/samtools/files/tabix/)
bgzip -c $altvcf.vcf > $altvcf.vcf.gz
tabix -p vcf $altvcf.vcf.gz

##############################################################################
# 4. AsmVar RecalibrationVQ
##############################################################################
echo "Step4: AsmVar RecalibrationVQ\n";
Vqprog=../../../src/AsmvarVarScore/ASV_VarScore.py
Roc=../../../src/AsmvarVarScore/ROC.py
FeatureToScore=../../../src/AsmvarVarScore/FeatureToScore.py
allInOneTabixvcf=$OutDir/$altvcf.vcf.gz #One tabix genotyping vcf file
traingSet=$DataDir/DemoDoublehit.txt # Positive training dataset which indicates the highly confident variants such as the double-hit variants from the 37 human de novo assemblies
queryLenList=$DataDir/DemoLength.list # The list file of each query length, generated by GATK_CreateSequenceDictionary; this information is used for calculating the relative distance between the variant towards the end of the scaffold; Generally speaking, variants on the edges of the assemblies are of less confidence.
name= "NA12878.APLG" #Name of your output vcf
python $Vqprog -T $traingSet -q $queryLenList -i $allInOneTabixvcf -f $OutDir/$name.Recal.Figure 2> $OutDir/$name.Recal.log | $bgzip -f > $OutDir/$name.Recal.vcf.gz && $tabix -f -p vcf $OutDir/$name.Recal.vcf.gz && echo "** Recal-VQ DONE **" 
python $Roc -i $OutDir/$name.Recal.vcf.gz -f $OutDir/$name.Recal.ROC.Figure > $OutDir/$name.Recal.ROC.txt && echo "** ROC done **" 
python $FeatureToScore $OutDir/$name.Recal.vcf.gz $queryLenList $OutDir/$name.recal.fs.Figure > $OutDir/$name.recal.fs.txt && echo "** Featrue Score Done **"
#You will obtain one vcf and some demo figures in the Recalibraion step.


##############################################################################
# 5. AsmVar VQSR
##############################################################################
echo "Step5: AsmVar VQSR\n";
VQSR=../../../src/AsmvarFinalQC/AsmvarSVQC.pl
allInOneRecalvcf=$OutDir/$name.Recal.vcf.gz #One recalibration vcf
VQ= 3 #Variants quality specified by the ROC figure
perl $VQSR -v $allInOneRecalvcf -q $VQ 2> $OutDir/$name.svqc.log | $bgzip -f > $OutDir/$name.svqc.vcf.gz && echo "** SVQC Done **" 

echo "** All Done **"
