#------------------------------------------------------------------------------------------------------------------------
chain.19to38 <- import.chain(system.file(package="trenaUtilities", "data", "hg19ToHg38.over.chain"))
chain.38to19 <- import.chain(system.file(package="trenaUtilities", "data", "hg38ToHg19.over.chain"))
#------------------------------------------------------------------------------------------------------------------------
liftoverBedTable.hg19.hg38 <- function(tbl)
{
     # some minimal sanity checks:  chrX or chromX or  X values in 1st column

   chromValues <- sort(unique(sub("^chr", "", tbl[,1])))
   chromValues <- sub("^chrom", "", chromValues)
   chromValues <- toupper(chromValues)
   expected.chrom.values <- all(chromValues %in% as.character(c(1:22, "X", "Y", "M")))

   if(!expected.chrom.values)
      stop(sprintf("unexpected chromosome names in input bed file: %s", paste(chromValues, collapse=",")))

     # all numeric in second and third columsn

   stopifnot(all(is.numeric(tbl[,2])))
   stopifnot(all(is.numeric(tbl[,3])))

   if(ncol(tbl) == 3)
       colnames(tbl) <- c("chrom", "start", "end")

   if(ncol(tbl) == 4){
      colnames(tbl) <- c("chrom", "start", "end", "name")
      }

   gr <- GRanges(Rle(tbl$chrom), IRanges(tbl$start, tbl$end))
   seqlevelsStyle(gr) <- "UCSC"
   seqinfo(gr) <- SeqinfoForUCSCGenome("hg19")[seqlevels(gr)]
   gr.38 <- unlist(liftOver(gr, chain.19to38))
   seqinfo(gr.38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr)]
   tbl.out <- data.frame(chrom=as.character(seqnames(gr.38)), start=start(gr.38), end=end(gr.38), stringsAsFactors=FALSE)
   if("name" %in% colnames(tbl))
      tbl.out <- cbind(tbl.out, name=tbl$name)

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
