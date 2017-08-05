#------------------------------------------------------------------------------------------------------------------------
chain.19to38 <- import.chain(system.file(package="trenaUtilities", "data", "hg19ToHg38.over.chain"))
chain.38to19 <- import.chain(system.file(package="trenaUtilities", "data", "hg38ToHg19.over.chain"))
#------------------------------------------------------------------------------------------------------------------------
liftoverBedTable.hg19.hg38 <- function(tbl)
{
     # some minimal sanity checks:  chrX or chromX or  X values in 1st column
   incoming.column.count <- ncol(tbl)
   incoming.column.names <- colnames(tbl)

   stopifnot(incoming.column.count >= 3)

   tbl.otherColumns <- data.frame()
   if(incoming.column.count > 3){
      tbl.otherColumns <- tbl[, 4:incoming.column.count]
      tbl <- tbl[, 1:3]
      colnames(tbl) <- c("chrom", "start", "end")
      }

   chromValues <- sort(unique(sub("^chr", "", tbl[,1])))
   chromValues <- sub("^chrom", "", chromValues)
   chromValues <- toupper(chromValues)
   expected.chrom.values <- all(chromValues %in% as.character(c(1:22, "X", "Y", "M")))

   if(!expected.chrom.values)
      stop(sprintf("unexpected chromosome names in input bed file: %s", paste(chromValues, collapse=",")))

     # all numeric in second and third columsn

   stopifnot(all(is.numeric(tbl[,2])))
   stopifnot(all(is.numeric(tbl[,3])))


   if(incoming.column.count == 3)
       colnames(tbl) <- c("chrom", "start", "end")

   gr <- GRanges(Rle(tbl$chrom), IRanges(tbl$start, tbl$end))
   seqlevelsStyle(gr) <- "UCSC"
   seqinfo(gr) <- SeqinfoForUCSCGenome("hg19")[seqlevels(gr)]
   gr.38 <- unlist(liftOver(gr, chain.19to38))
   seqinfo(gr.38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr)]
   tbl.out <- data.frame(chrom=as.character(seqnames(gr.38)), start=start(gr.38), end=end(gr.38), stringsAsFactors=FALSE)
   if(ncol(tbl.otherColumns) > 0)
      tbl.out <- cbind(tbl.out, tbl.otherColumns)

   colnames(tbl.out) <- incoming.column.names

   tbl.out

} # liftoverBedTable.hg19.hg38
#------------------------------------------------------------------------------------------------------------------------
liftoverBedFile.hg19.hg38 <- function(filename.hg19, filename.hg38)
{
   tbl.19 <- read.table(filename.hg19, sep="\t", as.is=TRUE)

   tbl.38 <- liftoverBedTable.hg19.hg38(tbl.19)
   printf("writing %d rows to %s", nrow(tbl.out), filename.hg38)

   write.table(tbl.out, file=filename.hg38, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

} # liftoverBedFile.hg19.hg38
#------------------------------------------------------------------------------------------------------------------------
