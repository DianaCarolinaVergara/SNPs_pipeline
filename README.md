# Phylogeny Construction with SNPs 

Better explanaitions in our **wiki**:
[Wiki](https://github.com/DianaCarolinaVergara/SNPs_pipeline/wiki)

Or the R Markdown with all the commands >[file](github.com/DianaCarolinaVergara/SNPs_pipeline/blob/master/Muricea_combined_aligned_to_reference_PIPELINE.Rmd)<


Or a full population genetics R pipeline in [GrunLab](http://grunwaldlab.github.io/Population_Genetics_in_R/qc.html)

#### Pipeline

We received a `vcf file` where SNPsaurus converted genomic DNA into nextRAD genotyping-by-sequencing libraries (SNPsaurus, LLC).
The **nextRAD libraries** were sequenced on a HiSeq 4000 on two lanes of 150bp reads.

![](https://www.rbcbioscienceusa.com/wp-content/uploads/2019/03/SNP-961x480.png)

With this pipeline you would be able to: 

> Modified from population genetics http://grunwaldlab.github.io/Population_Genetics_in_R/qc.html

* Using RStudio 3.4.4 and packages: 

> genepop, parallel, poppr, dartR, devtools, phytools, seqinr, phylotools, adegenet, pegas, hierfstat


   1. Remove `indels` (insertions-deletions)

   2. SNPs upper and lower 20% of depth distribution

   3. Delete: 

      3.1 Samples (missingness >70%)

      3.2 SNPs (>90%) with high degree of missingness information

   4. Rewrite `vcf file`


![](https://d33wubrfki0l68.cloudfront.net/62bcc8535a06077094ca3c29c383e37ad7334311/a263f/assets/img/logo.svg)

* Using **RStudio** or VCFTOOLS 0.1.17

   1. Run Minor Allele Frequency (MAF)

* Convert vcf file to `PHYLIP format`

> Suitable for RAxML

![](http://www.phylo.org/img/logo_cipres.gif)

![](https://cyverseuk.org/wp-content/uploads/2016/11/raxml_banner.png)

Or you can use another tool from CIPRESS Phylogenetic Collection (BEAST2, MRBAYES, RAXML) [CyVerse](https://cyverseuk.org/applications/cipress-phylogenetic-collection-beast-mrbayes-raxml/)

![](https://cyverseuk.org/wp-content/uploads/2016/11/beast2_banner.png) ![](https://cyverseuk.org/wp-content/uploads/2016/11/mrbayes_banner.png)

Need:

- [X] Desktop RStudio or R in a HPC. 
- [X] Load Packages 
- [X] VCF file 
- [ ] PGDSpider (or any format conversor)
- [X] CIPRES account 

__________________________________________________________


### References

*	RStudio: A Platform‐Independent IDE for R and Sweave - Racine - 2012 - Journal of Applied Econometrics - Wiley Online Library. https://onlinelibrary.wiley.com/doi/abs/10.1002/jae.1278.

*	variant call format and VCFtools | Bioinformatics | Oxford Academic. https://academic.oup.com/bioinformatics/article/27/15/2156/402296.

*	Lischer, H. E. L. & Excoffier, L. PGDSpider: an automated data conversion tool for connecting population genetics and genomics programs. Bioinformatics 28, 298–299 (2012).

*	Stamatakis, A. RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics 30, 1312–1313 (2014).

*	Miller, M. A., Pfeiffer, W. & Schwartz, T. The CIPRES Science Gateway: Enabling High-impact Science for Phylogenetics Researchers with Limited Resources. in Proceedings of the 1st Conference of the Extreme Science and Engineering Discovery Environment: Bridging from the eXtreme to the Campus and Beyond 39:1–39:8 (ACM, 2012). doi:10.1145/2335755.2335836.

*	Letunic, I. & Bork, P. Interactive Tree Of Life (iTOL) v4: recent updates and new developments. Nucleic Acids Res. 47, W256–W259 (2019).

* https://www.onworks.net/programs/vcftools-online?amp=0


![](https://s3.amazonaws.com/user-media.venngage.com/434503-b8f72eedfbceddd1c2371f8f115ce244.png)
_______________________________________________

.

.

.


#### Future pipeline: 

1. Learn the SNPs format

2. SNPs calling

3. SNPs filtering
   3.1. Samples quality
   3.2. Variants quality

