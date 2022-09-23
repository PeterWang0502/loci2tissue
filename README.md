# loci2tissue

Annotating a set of genomic loci is a critical bioinformatics problem with many use cases. The emergence of biomedical big data created unprecedented opportunities for annotation. Although tools are available, options are still limited. Among all types of emerging biomedical data, expression quantitative trait loci (eQTLs) become increasingly popular resources for researchers. Biologically speaking, eQTL’s nature of affecting gene expression level makes it an attractive resource for analyzing links between genetic variants and complex diseases and traits. 

`loci2tissue` is a bioinformatic tool that aims to utilize eQTLs to annotate genomic loci to obtain tissue enrichment ranking. It uses the input genomic loci and query size to form query regions, extracts tissue-specific eQTLs within these regions, calculates the average normalized expression level of these eQTLs’ corresponding eGenes, and then perform statistical tests to rank tissues by their gene enrichment. Loci2tissue generates more accurate results compared to other tools given its two advantages: first, it utilizes tissue-specific eQTLs instead of physical proximity to link loci to target genes, which generates more accurate mapping; second, loci2tissue uses percentiles instead of the original gene expression data and thus decreases the bias caused by extreme gene expression data and allow cross-tissue comparison.

# Protocols to use loci2tissue
## input query regions
loci2tissue takes the input query regions using genomic variants’ chromosomal locations. The variants should be filtered by a p-value threshold of 10^-6, which has proven to generate most enriched gene associations by multiple test trials. The input dataset should have two columns: one being chromosome number, the other being the position. Such data could be obtained by downloading disease or tissue’s association results from Phenotype-Genotype Integrator (PheGenI): https://www.ncbi.nlm.nih.gov/gap/phegeni. You could use the following code to process the data if you download input query region data from PheGenI.

```
yfile <- "..."    # this should be the local address of your association file
list <- read.delim(yfile, header=TRUE)
list <- list[list[,11]<=10^-6,]   # filter out variants whose p-values are higher than 10^-6
list <- list[,c(9,10)]
list[,1] = paste("chr", list[,1], sep="")
list
```
## prepare eQTL datasets
Tissue-specific eQTLs are used in loci2tissue to obtain significant eQTLs within the input query regions. The eQTLs will be filtered with a p-value threshold of 10^-4. An example for the eQTL set’s format id below.

```
snp.id	snp.chr	snp.pos	gene.name	gene.entrez.id
chr1_769577_G_A_b38,ENSG00000177757.2	chr1	769577	FAM87B	400728
chr1_769577_G_A_b38,ENSG00000225880.5	chr1	769577	LINC00115	79854
chr1_769577_G_A_b38,ENSG00000230092.7	chr1	769577	AL669831.4	NA
chr1_769577_G_A_b38,ENSG00000237491.8	chr1	769577	LINC01409	105378580
...
```

Multi-tissue eQTL data could be downloaded from GTEx: https://www.gtexportal.org/home/. In the example codes below, we only used three tissue-specific eQTL set to form the eQTL datasets. When actually using loci2tissue, you will likely use a lot more eQTL sets.

```
# read in eQTL set
aa.fileb <- ".../eQTLsets/Artery_Aorta.txt"
ac.fileb <- ".../eQTLsets/Artery_Coronary.txt"
ag.fileb <- ".../eQTLsets/Adrenal_Gland.txt"
eqtl.aa <- read.delim(aa.fileb, header=TRUE)
eqtl.ac <- read.delim(ac.fileb, header=TRUE)
eqtl.ag <- read.delim(ag.fileb, header=TRUE)

# construct eQTL dataset
eqtllist <- list(Artery_Aorta=eqtl.aa, Artery_Coronary=eqtl.ac, Adrenal_Gland=eqtl.ag)
eqtllist
```
## prepare normalized gene expression datasets

For each tissue-specific eQTL set used in the eQTL datasets, we will need a corresponding tissue-specific normalized gene expression set. The normalized gene expression matrices could be obtained from GTEx as well. 

The expression matrices will need to be processed first. We first sort the normalized gene expression data within each sample and then replace actual expression levels with their percentiles within each sample. The average percentile for each gene across samples is then calculated within the tissue. The average percentiles for n genes within the tissue are ranked into percentiles again and treated as the gene expression levels. Here is the format of the gene expression set.

```
X.chr	start	end	gene_id	gene.name	mean	percentile
chr17	19190244	19190245	ENSG00000264940.4	SNORD3C	0.0709829515369367	1
chr3	27802761	27802762	ENSG00000225548.5	LINC01980	0.0703849176627305	0.999953609203934
chr13	75876885	75876886	ENSG00000223458.2	LMO7DN-IT1	0.0684152453116919	0.999907218407868
chr4	188400735	188400736	ENSG00000249378.5	LINC01060	0.0662115890963626	0.999860827611802
chr14	36521643	36521644	ENSG00000253563.2	NKX2-1-AS1	0.0650326505390475	0.999814436815736
...
```

You will also need to obtain the means of gene expression levels for each tissue as a vector, which will later be used in the statistical tests. Below is the example code to construct the gene expression datasets.

```
# read in gene expression set
aa.filea <- "/Users/peter/Desktop/Loci2tissue/normData/aa_rank.txt"
ac.filea <- "/Users/peter/Desktop/Loci2tissue/normData/ac_rank.txt"
ag.filea <- "/Users/peter/Desktop/Loci2tissue/normData/ag_rank.txt"
rank.aa <- read.delim(aa.filea, header=TRUE)
rank.ac <- read.delim(ac.filea, header=TRUE)
rank.ag <- read.delim(ag.filea, header=TRUE)

# construct normalized gene expression dataset
ranklist <- list(Artery_Aorta=rank.aa, Artery_Coronary=rank.ac, Adrenal_Gland=rank.ag)
ranklist

# generate mean of gene expression levels vector
meanpercentile <- c()
for(i in 1:length(ranklist)){
  meanpercentile <- c(meanpercentile, mean(ranklist[[i]][,7]))
}
meanpercentile
```
## prepare tissue enrichment analysis
Below are the codes for function `tissueQuery` that performs tissue enrichment analysis.

```
tissueQuery <- function(list, size, eqtllist, ranklist, meanper){
  average <- c()
  egenesize <- c()
  for(i in 1:length(eqtllist)){
    output <- NULL
    newout <- NULL
    eqtl <- eqtllist[[i]]
    egene <- ranklist[[i]]
    egene.name <- egene[,5]
    for(n in 1:length(list[,1])){
      selection <- eqtl[list[n,1]==eqtl[,2],]
      eqtl.range <- selection[,3]
      selection <- selection[eqtl.range>=list[n,2]-size & eqtl.range<=list[n,2]+size,]
      output <- rbind(output, selection)
      output <- output[!duplicated(output[,4]),]
    }
    print(paste(names(eqtllist)[i], "eQTLs have been found", sep=" "))
    for(n in 1:length(output[,1])){
      b <- output[n,4]
      if(sum(is.na(egene[b==egene.name,]))==0){
        newout <- rbind(newout, egene[b==egene.name,])
      }
      newout <- newout[!duplicated(newout[,5]),]
    }
    print(paste(names(ranklist)[i], "eGenes have been located", sep=" "))
    average <- c(average, mean(newout[,7]))
    egenesize <- c(egenesize, length(newout[,1]))
  }
  tab <- data.frame(names(eqtllist), average, egenesize)
  tab <- tab[!is.na(tab[,2]),]
  p.val = rep(1, length(tab[,1]))
  for(n in 1:length(tab[,1])){
    p.val[n] = pnorm(meanper[n], mean=tab[n,2], sd=1/sqrt(12*tab[n,3]), lower.tail=TRUE)
  }
  log.p.val = rep(0, length(p.val))
  for(n in 1:length(p.val)){
    log.p.val[n] = -log(p.val[n])
  }
  tab <- data.frame(tab, p.val, log.p.val)
  tab <- tab[rev(order(tab[,5])),]
  return(tab)
}
```

now you can perform loci2tissue.
```
tab <- tissueQuery(list, 10000, eqtllist, ranklist, meanpercentile)
tab
```
