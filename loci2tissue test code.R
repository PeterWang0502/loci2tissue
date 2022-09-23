# define function tissueQuery
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

# input query region
yfile <- "/Users/peter/Downloads/Crohn Disease.txt"
list <- read.delim(yfile, header=TRUE)
list <- list[list[,11]<=10^-6,] #filter out variants whose p-values are higher than 10^-6
list <- list[,c(9,10)]
list[,1] = paste("chr", list[,1], sep="")
list

# import eQTL datasets
aa.fileb <- "/Users/peter/Desktop/Loci2path/eQTLsets/Artery_Aorta.txt"
ac.fileb <- "/Users/peter/Desktop/Loci2path/eQTLsets/Artery_Coronary.txt"
ag.fileb <- "/Users/peter/Desktop/Loci2path/eQTLsets/Adrenal_Gland.txt"
as.fileb <- "/Users/peter/Desktop/Loci2path/eQTLsets/Adipose_Subcutaneous.txt"
at.fileb <- "/Users/peter/Desktop/Loci2path/eQTLsets/Artery_Tibial.txt"
eqtl.aa <- read.delim(aa.fileb, header=TRUE)
eqtl.ac <- read.delim(ac.fileb, header=TRUE)
eqtl.ag <- read.delim(ag.fileb, header=TRUE)
eqtl.as <- read.delim(as.fileb, header=TRUE)
eqtl.at <- read.delim(at.fileb, header=TRUE)

eqtllist <- list(Artery_Aorta=eqtl.aa, Artery_Coronary=eqtl.ac, Adrenal_Gland=eqtl.ag, 
                 Adipose_Subcutaneous=eqtl.as, Artery_Tibial=eqtl.at)
eqtllist

# import normalized gene expression datasets
aa.filea <- "/Users/peter/Desktop/Loci2tissue/normData/aa_rank.txt"
ac.filea <- "/Users/peter/Desktop/Loci2tissue/normData/ac_rank.txt"
ag.filea <- "/Users/peter/Desktop/Loci2tissue/normData/ag_rank.txt"
as.filea <- "/Users/peter/Desktop/Loci2tissue/normData/as_rank.txt"
at.filea <- "/Users/peter/Desktop/Loci2tissue/normData/at_rank.txt"
rank.aa <- read.delim(aa.filea, header=TRUE)
rank.ac <- read.delim(ac.filea, header=TRUE)
rank.ag <- read.delim(ag.filea, header=TRUE)
rank.as <- read.delim(as.filea, header=TRUE)
rank.at <- read.delim(at.filea, header=TRUE)

ranklist <- list(Artery_Aorta=rank.aa, Artery_Coronary=rank.ac, Adrenal_Gland=rank.ag, 
                 Adipose_Subcutaneous=rank.as, Artery_Tibial=rank.at)
ranklist

meanpercentile <- c()
for(i in 1:length(ranklist)){
  meanpercentile <- c(meanpercentile, mean(ranklist[[i]][,7]))
}
meanpercentile

# perform tissue enrichment analysis
tab <- tissueQuery(list, 10000, eqtllist, ranklist, meanpercentile)
tab
