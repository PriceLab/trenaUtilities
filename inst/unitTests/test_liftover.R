# test_liftover.R
#------------------------------------------------------------------------------------------------------------------------
library(trenaUtilities)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   printf("--- runTests")
   test_liftoverBedTable.19.38()
   test_liftoverNarrowPeakBedTable.19.38()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_liftoverBedTable.19.38 <- function()
{
   printf("--- test_liftoverBedTable.19.38")

     # submit a 1-row bed table with the coordinates of the AGL gene
     # this provides a simple check on hg19 -> hg38 coordinates, which
     # are easily available in the UCSC genome browser:

   tbl.agl.19 <- data.frame(chrom="chr1", start=99850084, end=99924020, stringsAsFactors=FALSE)
   tbl.agl.38 <- liftoverBedTable.hg19.hg38(tbl.agl.19)
   checkEquals(as.list(tbl.agl.38), list(chrom="chr1", start=99384528, end=99458464))

     # now a ginned up 3-region bed table, arbitrarily moved far down chromosome 1

   tbl.19 <- data.frame(chrom=rep("chr1", 3),
                        start=c(10155, 10175, 10435),
                        end=c(10305, 10325, 10585),
                        stringsAsFactors=FALSE)

     # move those locations far down chr1, into regions which were remapped between hg19 and hg38
   tbl.19$start <- 100000000 + tbl.19$start
   tbl.19$end   <- 100000000 + tbl.19$end

   tbl.38 <- liftoverBedTable.hg19.hg38(tbl.19)

   checkEquals(dim(tbl.19), dim(tbl.38))

   hg.19.lengths <- 1 + tbl.19$end - tbl.19$start
   hg.38.lengths <- 1 + tbl.38$end - tbl.38$start
   checkTrue(all(hg.19.lengths == hg.38.lengths))

} # test_liftoverBedTable.19.38
#------------------------------------------------------------------------------------------------------------------------
test_liftoverNarrowPeakBedTable.19.38 <- function()
{
   printf("--- test_liftoverNarrowPeakBedTable.19.38")

   tbl.narrowPeak.19 <- data.frame(chrom=rep("chr1", 4),
                                  chromStart=c(99862695, 99872315, 99907355, 99923255),
                                    chromEnd=c(99862845, 99872465, 99907505, 99923405),
                                        name=rep(".", 4),
                                       score=rep(0, 4),
                                      strand=rep(".", 4),
                                 signalValue=c(67, 43, 21, 38),
                                      pValue=rep(-1, 4),
                                      qValue=rep(-1, 4),
                                        peak=rep(75, 4))

   tbl.narrowPeak.38 <- liftoverBedTable.hg19.hg38(tbl.narrowPeak.19)
   checkEquals(colnames(tbl.narrowPeak.19), colnames(tbl.narrowPeak.38))
   widths.19 <- 1 + tbl.narrowPeak.19$chromEnd - tbl.narrowPeak.19$chromStart;
   widths.38 <- 1 + tbl.narrowPeak.38$chromEnd - tbl.narrowPeak.38$chromStart;
   checkTrue(all(widths.19 == widths.38))
   checkTrue(all(tbl.narrowPeak.19$chromStart != tbl.narrowPeak.38$chromStart))
   checkTrue(all(tbl.narrowPeak.19$chromEnd != tbl.narrowPeak.38$chromEnd))

} # test_liftoverNarrowPeakBedTable.19.38
#------------------------------------------------------------------------------------------------------------------------
