file.cer <-"Scaled_Winsorized_MayoRNAseq_CER.csv"
file.tcx <- "Scaled_Winsorized_MayoRNAseq_TCX.csv"


tbl.counts <- read.table(file.cer, header=TRUE,
                         stringsAsFactors=FALSE, quote="", sep = ",", check.names=FALSE)
names(tbl.counts)[1] <- "ENSMBL_ID"
rownames(tbl.counts) <- tbl.counts$ENSMBL_ID
tbl.counts$ENSMBL_ID <- NULL
tbl.counts <- t(tbl.counts)
library(org.Hs.eg.db)

# get names for the genes
ensembl.geneIDs <- rownames(tbl.counts)
tbl.map <- select(org.Hs.eg.db, keys=ensembl.geneIDs, keytype="ENSEMBL", columns=c("ENSEMBL", "SYMBOL"))
dups <- which(duplicated(tbl.map$ENSEMBL))
if(length(dups) > 0)
   tbl.map <- tbl.map[-dups,]
tbl <- cbind(tbl.map, tbl.counts)
dim(tbl)
unmapped <- which(is.na(tbl$SYMBOL))
length(unmapped)
if(length(unmapped) > 0)
  tbl <- tbl[-unmapped,]
dups <- which(duplicated(tbl$SYMBOL))
if(length(dups) > 0)
    tbl <- tbl[-dups,]

rownames(tbl) <- tbl$SYMBOL
delete.these.columns <- match(c("ENSEMBL", "SYMBOL"), colnames(tbl))
if(length(delete.these.columns) > 0)
    tbl <- tbl[, -delete.these.columns]

dim(tbl)
mtx.cer <- tbl
mtx.tcx
dim(mtx.cer)
dim(mtx.tcx)
mean(mtx.cer)
mean(mtx.tcx)


