library(haplo.stats)

# Command line parameter for capture prefix of plink files (.map .ped)
argv = commandArgs()
prefix = argv[6]

# Load in name and SNPs of MHs
load("MHsInfo.Rdata")

# Load in plink file, remove unneeded columns and add colnames
dataInput <- read.csv(paste0(prefix, ".ped"), sep=" ", header=FALSE)
rownames(dataInput) <- dataInput[,2]
dataInput <- dataInput[, -1:-6]
dataMap <- read.csv(paste0(prefix, ".map"), sep="\t", header=FALSE)
colName = c()
for (rsID in dataMap$V2){
  colName <- c(colName, paste0(rsID, ".a1"))
  colName <- c(colName, paste0(rsID, ".a2"))
}
colnames(dataInput) <- colName

# Function for phasing
GTbyHS <- function(label, geno){
  save.em <- haplo.em(geno=geno, locus.label= label, miss.val=c(0,NA))
  hap1codeAll <- sapply(1:nrow(geno), function(i) {idx = save.em$subj.id == i; save.em$hap1code[idx][which.max(save.em$post[idx])];})
  hap2codeAll <- sapply(1:nrow(geno), function(i) {idx = save.em$subj.id == i; save.em$hap2code[idx][which.max(save.em$post[idx])];})
  MH_label <- sapply(1:nrow(save.em$haplotype), function(x) paste(save.em$haplotype[x,], collapse = '-'))
  Genotype1 <- factor(hap1codeAll, levels=1:length(MH_label),labels=MH_label)
  Genotype2 <- factor(hap2codeAll, levels=1:length(MH_label),labels=MH_label)
  GenotypeAll <- data.frame(Genotype1, Genotype2)
  rownames(GenotypeAll) <- rownames(geno)
  return(GenotypeAll)
}

# Phase available MHs
geno_phased <- data.frame(1:nrow(dataInput))
geno_phased <- geno_phased[,-1] # Empty dataframe for recoding results
colName <- c("IDs")
for (MH in names(MHsInfo)){
  # Find MHs having all SNPs that make it up
  flag = TRUE
  for (SNP in MHsInfo[[MH]]){
    if (!SNP %in% dataMap$V2){
      flag = FALSE
    }
  }
  # Phase
  if (flag){
    colName <- c(colName, MH)
    selectcol <- c()
    for (SNP in MHsInfo[[MH]]){
      selectcol <- c(selectcol,paste0(SNP, ".a1"))
      selectcol <- c(selectcol,paste0(SNP, ".a2"))
    }
    geno_sel <- subset(dataInput, select = selectcol)
    geno_tem <- GTbyHS(label=MHsInfo[[MH]], geno=geno_sel)
    geno_phased <- cbind(geno_phased, geno_tem)
  }
}

# Reshape and output
geno_phased <- cbind(rownames(geno_phased), rownames(geno_phased), geno_phased)
rownames(geno_phased) <- NULL
colnames(geno_phased) <- NULL
geno_phased_part1 <- geno_phased[,seq(1,ncol(geno_phased),2)]
colnames(geno_phased_part1) <- colName
geno_phased_part2 <- geno_phased[,seq(0,ncol(geno_phased),2)]
colnames(geno_phased_part2) <- colName
dataOutput <- rbind(geno_phased_part1, geno_phased_part2)
dataOutput <- dataOutput[order(dataOutput$IDs), ]
write.table(dataOutput, file=paste0(prefix, "_phased.txt"), sep=' ', quote=FALSE, row.names=FALSE, col.names=TRUE)
