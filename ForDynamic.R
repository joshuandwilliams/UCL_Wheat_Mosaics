#arina vs chinese missing 7B
#sy mattis v chinese missing 7B
getwd()
setwd("E://Prof Mott Files//snps")

library(tidyverse)

#----
# PREPARING DATA
# Load in the datasets
filenames <- list.files(pattern="*.delta.snps")
datalist <- lapply(filenames, read_tsv, col_names=c("POSR", "BASER", "BASEQ", "TAGR", "TAGQ"), col_types="icc_______cc")

# Select only those dataframes where the reference is chinese spring
refs <- c()
for(i in 1:length(filenames)){
  name <- filenames[[i]]
  if(substr(name[[1]], 1, 7) == "chinese"){
    refs <- c(refs, i)
  }
}
datalist <- datalist[refs]

# Modify rows for more useful information
for(j in 1:length(datalist)){
  chromname <- rep(substr(datalist[[j]]$TAGR[1], 1, 5), length(rownames(datalist[[j]])))
  queryname <- rep(substr(datalist[[j]]$TAGQ[1], 8, 10), length(rownames(datalist[[j]])))
  datalist[[j]] <- datalist[[j]] %>% add_column(CHROM = chromname, QUERY = queryname, .before = "POSR")
  datalist[[j]] <- datalist[[j]][,1:5]
}

# Group dataframes in list by chromosome and output as separate lists
CHRNAMES <- c("chr1A", "chr1B", "chr1D", "chr2A", "chr2B", "chr2D", "chr3A", "chr3B", "chr3D", "chr4A", "chr4B", "chr4D", "chr5A", "chr5B", "chr5D", "chr6A", "chr6B", "chr6D", "chr7A", "chr7B", "chr7D")
for(Chromname in CHRNAMES){
  templist <- list()
  for(n in 1:length(datalist)){
    if(datalist[[n]]$CHROM[1] == Chromname){
      templist[[length(templist)+1]] <- datalist[[n]]
    }
  }
  assign(paste0("list_", Chromname, sep = ""), templist)
}
rm(templist)
rm(datalist)


#----
# SUBSETTING DATA

# Load in founders data
founders <- read.table("C:\\Users\\joshu\\OneDrive\\Josh\\Summer 2021\\Richard Mott\\0 Original Data\\FOUNDERS\\Founders.traw", header = T)
founders$CHROM <- paste("chr", substr(founders$SNP, 1, 2), sep = "")

# Create unique row identifier
founders$CombinedPos <- paste(founders$CHROM, founders$POS)

# Remove indels, non-biallelic SNPs, and SNPs that aren't present in the founders data
SNPstotal <- 0
subsetlist <- list()
founderslist <- list()

subsetchrom <- function(founderdata, chromdata){
  full <- Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2), chromdata)
  full <- full[order(full$POSR),]
  
  # Delete indels
  full <- full[!(full$BASER == "." | full$BASEQ == "."),] # These are insertions/deletions
  
  # Delete non-biallelic snps - table needs to be sorted by POSR
  toremove <- c()
  for(row in 2:length(rownames(full))){
    if(full$POSR[row] == full$POSR[row-1]){
      if(full$BASER[row] != full$BASER[row-1] | full$BASEQ[row] != full$BASEQ[row-1]){
        toremove <- c(toremove, row, row-1)
      }
    }
  }
  full <- full[-toremove,]
  
  fullcopy <- full
  
  # Remove rows where the pos is not in both datasets
  full$QUERY <- NULL
  full$CombinedPos <- paste(full$CHROM, full$POSR)
  full <- full[!duplicated(full$CombinedPos),]
  
  full <- full[full$CombinedPos %in% founders$CombinedPos,]
  founderssubset <- founders[founders$CombinedPos %in% full$CombinedPos,]
  
  # Remove rows where the REF and ALT alleles arent the same
  full <- subset(full, founderssubset$COUNTED == full$BASEQ & founderssubset$ALT == full$BASER)
  founderssubset <- founders[founders$CombinedPos %in% full$CombinedPos,]
  
  fullcopy$CombinedPos <- paste(fullcopy$CHROM, fullcopy$POSR)
  fullcopy <- fullcopy[fullcopy$CombinedPos %in% full$CombinedPos,]
  
  #SNPstotal <<- SNPstotal + length(rownames(founderssubset))
  subsetlist[[length(subsetlist)+1]] <<- fullcopy
  founderslist[[length(founderslist)+1]] <<- founderssubset
}

#subsetchrom(founderdata=founders, chromdata=list_chr6D)

applylist <- list(list_chr1A, list_chr1B, list_chr1D, list_chr2A, list_chr2B, list_chr2D, list_chr3A, list_chr3B, list_chr3D, list_chr4A, list_chr4B, list_chr4D, list_chr5A, list_chr5B, list_chr5D, list_chr6A, list_chr6B, list_chr6D, list_chr7A, list_chr7B, list_chr7D)
mapply(subsetchrom, applylist, MoreArgs=list(founderdata=founders))

subsetlist

#----
# EMISSIONS
finalfunc <- function(datalist){
  for(data in datalist){
    chromname <- data$CHROM[1]
    rows <- sort(unique(data$QUERY))
    cols <- sort(unique(data$POSR))
    em <- matrix(0, nrow=length(rows), ncol=length(cols))
    rownames(em) <- rows
    colnames(em) <- cols
    po <- data.frame(CHROM=data$CHROM[which(!duplicated(data$POSR))], POS=unique(data$POSR))
    
    for(index in 1:length(rownames(data))){
      row <- match(data$QUERY[index], rows)
      col <- match(data$POSR[index], cols)
      em[row, col] <- 1
    }
    save(em, po, file = paste("E:\\Prof Mott Files\\Mosaics Data\\emissiondir\\", chromname, ".RData", sep = ""))
  }
}

finalfunc(subsetlist)

# DOSAGES
for(chr in founderslist){
  chromname <- chr$CHROM[1]
  snp.data <- data.frame(CHROM=chr$CHROM, POS=chr$POS)
  chr.split <- chr[,c(7:22)]
  chr.split <- t(chr.split)
  
  save(chr.split, snp.data, file = paste("E:\\Prof Mott Files\\Mosaics Data\\dosagedir\\", chromname, ".RData", sep = ""))
}
