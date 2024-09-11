AsmVar
==========


AsmVar is a software for discovery, genotyping and characterization of structural variants and novel sequence at single nucleotide resolution from de novo assemblies on a population scale

__Contributors:__ Shujia Huang, Siyang Liu, Junhua Rao & Weijian Ye (Contribute equally)  <br/>
__Contact              :__ liusiyang@genomics.cn & huangshujia@genomics.cn                <br/>
__Institute            :__ BGI & KU in GenomeDK consortium                 <br/>
__Latest Version         :__ 2015-04-16                                      <br/>


Dependency
---------------------------------------
* **[last](http://last.cbrc.jp/)**
* **[scipy](http://www.scipy.org/)**
* **[numpy](http://www.numpy.org/)**
* **[scikit-learn](http://scikit-learn.org/stable/)**
* **[matplotlib](http://matplotlib.org/)**
* **[tabix](http://sourceforge.net/projects/samtools/files/tabix/)**

Installation
---------------------------------------
	$ cd src/AsmvarDetect
	$ make

One short shell ~/demo/AsmVarDetector/Test.sh can be run to check the correct installation of AsmVar. It can be used to verify whether the installation is working correctly. This shell will run automatically and the output files are t.age,t4.1025.vcf,t4.1025.svd,t4.1025.summary,t4.1025.gap.bed.

Running AsmVar workflow
---------------------------------------
The demo shell  ./demo/Cookbook/scripts/DemoPipelineGuideline.sh. provides the example for application of Asmvar to detect, genotype and recalibrate the variation calls <br/>
Due to the github storage limitations, the full demo including the output files can be found through this link https://www.dropbox.com/sh/fapzzvvo1whizmn/AAAeBsTSRmTyOn8zVJc24fgca?dl=0 <br/>
We are currently improving the characterization modules and therefore will close the codes for those modules for a while before the 1st, July, 2015 <br/>

There are 5 main steps about Asmvar workflow, <br/>
- Step1: AsmVar Detection <br/>
- Step2: AsmVar Altalignment <br/>
- Step3: AsmVar Genotyping <br/>
- Step4: AsmVar RecalibrationVQ <br/>
- Step5: AsmVar VQSR <br/>


Instant plans for further developments    
---------------------------------------  
1. Open codes again for the characterization modules
2. Improve user experience for AsmVar based on the feedback including possibility of using Bam instead MAF as input files; Level of user friendliness.
3. Streamline evaluations of de novo assemblies with regards to continuity, completeness and accuracy. 
4. Exploration of novel statistical approaches especially for the alternative align, genotyping and recalibration modules to resolve multi-allelic structural variations and to improve time and memory efficiency
5. Application of AsmVar in haplotype-resolved de novo assemblies

Contribution
------------
The AsmVar is initially developed for the DanishPanGenome project in the GenomeDK platform. The AsmVar framework construction, applications, statistical methods and the evaluation protocols are established by [Siyang Liu](https://github.com/SiyangLiu) and [Shujia Huang](https://github.com/ShujiaHuang), Junhua Rao and [Weijian Ye](https://github.com/WeijianYe) led by [Siyang Liu](https://github.com/SiyangLiu) and [Shujia Huang](https://github.com/ShujiaHuang).  The coding contributions for different modules are written in the source code title. Most of the initial codes are written in C++, python and perl by [Shujia Huang](https://github.com/ShujiaHuang) and perl, python and R by [Siyang Liu](https://github.com/SiyangLiu)(All modules especially SV discovery and genotyping), Junhua Rao (especially SV mechanism), [Weijian Ye](https://github.com/WeijianYe) (especially Ancestral state and Novel sequence) . [Shujia Huang](https://github.com/ShujiaHuang) and [Siyang Liu](https://github.com/SiyangLiu) take charge of the software architecture and the efficiency of the algorithms. The work is supervised by [Anders Krogh](http://www.binf.ku.dk/staff/?pure=en/persons/8330) and Jun Wang.

Please cite the paper:  

- Siyang Liu, Shujia Huang, Junhua Rao, Weijian Ye, The Genome Denmark Consortium, Anders Krogh, Jun Wang, Discovery, genotyping and characterization of structural variation and novel sequence at single nucleotide resolution from de novo genome assemblies on a population scale, *GigaScience*, Volume 4, Issue 1, December 2015, s13742–015–0103–4, https://doi.org/10.1186/s13742-015-0103-4


## LICENSE 
Released under the [MIT License](http://opensource.org/licenses/MIT)
Copyright &copy; 2014-2015
