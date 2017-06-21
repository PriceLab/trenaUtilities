library(TrenaHelpers)
library(RPostgreSQL)   # TODO: part of the nasty hack.  remove
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx")){
   load("~/github/projects/examples/microservices/trenaGeneModel/datasets/coryAD/rosmap_counts_matrix_normalized_geneSymbols_25031x638.RData")
   # load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   mtx <- asinh(mtx)
   mtx.var <- apply(mtx, 1, var)
   deleters <- which(mtx.var < 0.01)
   if(length(deleters) > 0)   # 15838 x 638
      mtx <- mtx[-deleters,]
   }

#------------------------------------------------------------------------------------------------------------------------
# temporary hack.  the database-accessing classes should clean up after themselves
closeAllPostgresConnections <- function()
{
   connections <- RPostgreSQL::dbListConnections(RPostgreSQL::PostgreSQL())
   printf(" TODO, nasty hack!  found %d Postgres connections open, now closing...", length(connections))
   if(length(connections) > 0) for(i in 1:length(connections)){
      dbDisconnect(connections[[i]])
      } # for i

} # closeAllPostgresConnections
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_canonicalTableColumnNames()

   test_.callFootprintFilterAndTFexpander()
   test_getRegulatoryRegions_oneFootprintSource()
   test_getRegulatoryRegions_twoFootprintSources()

   test_.callHumanDHSFilterAndTFexpander()
   test_getRegulatoryRegions_encodeDHS
   test_getRegulatoryRegions_twoFootprintSources_oneDHS()

   test_createGeneModel()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   aqp4.tss <- 26865884
   prep <- TrenaPrep("AQP4", targetGeneTSS=aqp4.tss, "chr18", aqp4.tss-100, aqp4.tss+100, regulatoryRegionSources=list())
   checkTrue(is(prep, "TrenaPrep"))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_canonicalTableColumnNames <- function()
{
   printf("--- test_canonicalTableColumnNames")
   aqp4.tss <- 26865884
   prep <- TrenaPrep("AQP4", targetGeneTSS=aqp4.tss, "chr18", aqp4.tss-100, aqp4.tss+100, regulatoryRegionSources=list())
   checkTrue(is(prep, "TrenaPrep"))
   checkTrue(all(getRegulatoryTableColumnNames(prep) ==
                 c("chrom", "motifStart", "motifEnd", "motifName", "strand", "score", "length", "distance.from.tss", "id", "tf")))
   checkTrue(all(getGeneModelTableColumnNames(prep) ==
                 c("tf", "randomForest", "pearson", "spearman", "betaLasso", "pcaMax", "concordance")))


} # test_canonicalTableColumnNames
#------------------------------------------------------------------------------------------------------------------------
test_.callFootprintFilterAndTFexpander <- function()
{
   printf("--- test_.callFootprintFilterAndTFexpander")
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   chromosome <- "chr18"
   chromStart <- aqp4.tss - 100
   chromEnd   <- aqp4.tss + 100
   source <- "postgres://whovian/brain_hint_20"

   prep <- TrenaPrep("AQP4", aqp4.tss, chromosome, chromStart, chromEnd, regulatoryRegionSources=list(source))

   tbl.fp <- TrenaHelpers:::.callFootprintFilter(prep, source, chromosome, chromStart, chromEnd,
                                                 targetGene, targetGeneTSS=aqp4.tss)
   checkTrue(all(colnames(tbl.fp) == getRegulatoryTableColumnNames(prep)))
   checkTrue(nrow(tbl.fp) > 8 & nrow(tbl.fp) <20)
   tbl.fpExpanded <- expandRegulatoryRegionsTableByTF(prep, tbl.fp)
   checkTrue(nrow(tbl.fpExpanded) > 100)
   checkTrue(min(nchar(tbl.fpExpanded$tf)) > 0)
   checkTrue(max(nchar(tbl.fpExpanded$tf)) < 20)
   checkEquals(length(grep(";", tbl.fpExpanded$tf)), 0)
   checkTrue(length(grep(";", tbl.fp$tf)) > 0)

} # test_.callFootprintFilterAndTFexpander
#------------------------------------------------------------------------------------------------------------------------
test_.callHumanDHSFilterAndTFexpander <- function()
{
   printf("--- test_.callHumanDHSFilterAndTFexpander")
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   chromosome <- "chr18"
   chromStart <- aqp4.tss - 100
   chromEnd   <- aqp4.tss + 100
   source <- "encodeHumanDHS"

   prep <- TrenaPrep("AQP4", aqp4.tss, chromosome, chromStart, chromEnd, regulatoryRegionSources=list(source))

   tbl.fp <- TrenaHelpers:::.callHumanDHSFilter(prep, source, chromosome, chromStart, chromEnd,
                                                targetGene, targetGeneTSS=aqp4.tss)

   checkTrue(all(colnames(tbl.fp) == getRegulatoryTableColumnNames(prep)))
   checkTrue(nrow(tbl.fp) > 10 & nrow(tbl.fp) < 30)
   tbl.fpExpanded <- expandRegulatoryRegionsTableByTF(prep, tbl.fp)
   checkTrue(nrow(tbl.fpExpanded) > 100)
   checkTrue(min(nchar(tbl.fpExpanded$tf)) > 0)
   checkTrue(max(nchar(tbl.fpExpanded$tf)) < 20)
   checkEquals(length(grep(";", tbl.fpExpanded$tf)), 0)
   checkTrue(length(grep(";", tbl.fp$tf)) > 0)

} # test_.callFootprintFilterAndTFexpander
#------------------------------------------------------------------------------------------------------------------------
test_getRegulatoryRegions_oneFootprintSource <- function()
{
   printf("--- test_getRegulatoryRegions_oneFootprintSource")
   aqp4.tss <- 26865884
   sources <- list("postgres://whovian/brain_hint_20")
   prep <- TrenaPrep("AQP4", aqp4.tss, "chr18", aqp4.tss-100, aqp4.tss+100, regulatoryRegionSources=sources)
   x <- getRegulatoryRegions(prep)
   checkTrue(is(x, "list"))
   checkEquals(names(x), as.character(sources))
   checkTrue(all(unlist(lapply(x, function(element) is(element, "data.frame")), use.names=FALSE)))

   tbl.reg <- x[[sources[[1]]]]
   checkTrue(all(colnames(tbl.reg) == getRegulatoryTableColumnNames(prep)))
   closeAllPostgresConnections()

} # test_getRegulatoryRegions_oneFootprintSource
#------------------------------------------------------------------------------------------------------------------------
test_getRegulatoryRegions_twoFootprintSources <- function()
{
   printf("--- test_getRegulatoryRegions_twoFootprintSources")
   aqp4.tss <- 26865884
   sources <- list("postgres://whovian/brain_hint_20",
                   "postgres://whovian/brain_wellington_16")
   prep <- TrenaPrep("AQP4", aqp4.tss, "chr18", aqp4.tss-100, aqp4.tss+100, regulatoryRegionSources=sources)

      # get, but don't combined, footprints from both sources
   x <- getRegulatoryRegions(prep, combine=FALSE)  # the default
   checkTrue(is(x, "list"))
   checkEquals(length(x), 2)
   checkEquals(names(x), as.character(sources))
   checkTrue(all(unlist(lapply(x, function(element) is(element, "data.frame")), use.names=FALSE)))

   tbl.reg <- x[[sources[[1]]]]
   checkTrue(all(colnames(tbl.reg) == getRegulatoryTableColumnNames(prep)))

      # get AND combine  footprints from both sources
   x <- getRegulatoryRegions(prep, combine=TRUE)  # the default
   checkTrue(is(x, "list"))
   checkEquals(length(x), 3)
   checkEquals(names(x), c(as.character(sources), "all"))
   checkTrue(all(unlist(lapply(x, function(element) is(element, "data.frame")), use.names=FALSE)))

       # the "all" table should have rowcount the sum of the other two
   summary <- lapply(x, nrow)
   checkEquals(summary[["postgres://whovian/brain_hint_20"]], 13)
   checkEquals(summary[["postgres://whovian/brain_wellington_16"]], 11)
   checkEquals(summary[["all"]], 24)

   checkTrue(all(colnames(x[[1]]) == getRegulatoryTableColumnNames(prep)))
   checkTrue(all(colnames(x[[2]]) == getRegulatoryTableColumnNames(prep)))
   checkTrue(all(colnames(x[[3]]) == getRegulatoryTableColumnNames(prep)))

   closeAllPostgresConnections()

} # test_getRegulatoryRegions_twoFootprintSources
#------------------------------------------------------------------------------------------------------------------------
test_getRegulatoryRegions_encodeDHS <- function()
{
   printf("--- test_getRegulatoryRegions_encodeDHS")

   aqp4.tss <- 26865884
   sources <- list("encodeHumanDHS")
   prep <- TrenaPrep("AQP4", aqp4.tss, "chr18", aqp4.tss-100, aqp4.tss+100, regulatoryRegionSources=sources)
   x <- getRegulatoryRegions(prep)
   checkTrue(is(x, "list"))
   checkEquals(names(x), as.character(sources))
   checkTrue(all(unlist(lapply(x, function(element) is(element, "data.frame")), use.names=FALSE)))

   tbl.reg <- x[[sources[[1]]]]
   checkTrue(all(colnames(tbl.reg) == getRegulatoryTableColumnNames(prep)))
   checkTrue(nrow(tbl.reg) > 20)
   checkEquals(length(grep("AQP4.dhs", tbl.reg$id)), nrow(tbl.reg))

   closeAllPostgresConnections()

} # test_getRegulatoryRegions_encodeDHS
#------------------------------------------------------------------------------------------------------------------------
test_getRegulatoryRegions_twoFootprintSources_oneDHS <- function()
{
   printf("--- test_getRegulatoryRegions_twoFootprintSources_oneDHS")
   aqp4.tss <- 26865884

   sources <- list("postgres://whovian/brain_hint_20",
                   "encodeHumanDHS",
                   "postgres://whovian/brain_wellington_16")

   prep <- TrenaPrep("AQP4", aqp4.tss, "chr18", aqp4.tss-100, aqp4.tss+100, regulatoryRegionSources=sources)

      # get, but don't combined, footprints from both sources
   x <- getRegulatoryRegions(prep, combine=FALSE)  # the default
   checkTrue(is(x, "list"))
   checkEquals(length(x), 3)
   checkEquals(sort(names(x)), sort(as.character(sources)))
   checkTrue(all(unlist(lapply(x, function(element) is(element, "data.frame")), use.names=FALSE)))

   tbl.reg <- x[[sources[[1]]]]
   checkTrue(all(colnames(tbl.reg) == getRegulatoryTableColumnNames(prep)))

   summary <- lapply(x, nrow)
   checkEquals(summary[["postgres://whovian/brain_hint_20"]], 13)
   checkEquals(summary[["postgres://whovian/brain_wellington_16"]], 11)
   checkEquals(summary[["encodeHumanDHS"]], 23)


      # get AND combine footprints and DHS regions  all three  sources
   x <- getRegulatoryRegions(prep, combine=TRUE)  # the default
   checkTrue(is(x, "list"))
   checkEquals(length(x), 4)
   checkEquals(sort(names(x)), sort(c(as.character(sources), "all")))
   checkTrue(all(unlist(lapply(x, function(element) is(element, "data.frame")), use.names=FALSE)))

       # the "all" table should have rowcount the sum of the other two
   summary <- lapply(x, nrow)
   checkEquals(summary[["postgres://whovian/brain_hint_20"]], 13)
   checkEquals(summary[["postgres://whovian/brain_wellington_16"]], 11)
   checkEquals(summary[["encodeHumanDHS"]], 23)
   checkEquals(summary[["all"]], 47)

   checkTrue(all(colnames(x[[1]]) == getRegulatoryTableColumnNames(prep)))
   checkTrue(all(colnames(x[[2]]) == getRegulatoryTableColumnNames(prep)))
   checkTrue(all(colnames(x[[3]]) == getRegulatoryTableColumnNames(prep)))
   checkTrue(all(colnames(x[[4]]) == getRegulatoryTableColumnNames(prep)))

   closeAllPostgresConnections()

} # test_getRegulatoryRegions_twoFootprintSources_oneDHS
#------------------------------------------------------------------------------------------------------------------------
test_createGeneModel <- function()
{
   printf("--- test_createGeneModel")
   aqp4.tss <- 26865884
   fp.source <- "postgres://whovian/brain_hint_20"
   sources <- list(fp.source)

   prep <- TrenaPrep("AQP4", aqp4.tss, "chr18", aqp4.tss-100, aqp4.tss+100, regulatoryRegionSources=sources)
   x <- getRegulatoryRegions(prep)
   closeAllPostgresConnections()
   tbl.regulatoryRegions <- x[[fp.source]]

   tbl.geneModel <- createGeneModel(prep, "randomForest", tbl.regulatoryRegions, mtx)
   checkTrue(all(getGeneModelTableColumnNames(prep) == colnames(tbl.geneModel)))
   tbl.strong <- subset(tbl.geneModel, randomForest > 3)
   checkTrue(all(c("TEAD1", "SP3", "KLF3", "NEUROD2") %in% tbl.strong$tf))

} # test_createGeneModel
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
