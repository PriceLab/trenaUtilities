library(TrenaHelpers)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   aqp4.tss <- 26865884
   prep <- TrenaPrep("AQP4", "chr18", aqp4.tss-100, aqp4.tss+100, regulatoryRegionSources=list())
   checkTrue(is(prep, "TrenaPrep"))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------

