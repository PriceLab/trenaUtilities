library(trenaUtilities)
library(RPostgreSQL)   # TODO: part of the nasty hack.  remove
library(RUnit)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)   # load here so that it is ready when needed
library(MotifDb)
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
   test_buildMultiModelGraph_oneModel()
   test_buildMultiModelGraph_fiveModels()
   test_buildMultiModelGraph_twoModels_15k_span()

   test_geneModelLayout();

   test_assessSnp()
   test_assessSnp_mutOnly_wtOnly()
   test_assessSnpMotifDbMatrices()

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

   tbl.fp <- trenaUtilities:::.callFootprintFilter(prep, source, chromosome, chromStart, chromEnd,
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

   tbl.fp <- trenaUtilities:::.callHumanDHSFilter(prep, source, chromosome, chromStart, chromEnd,
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

   invisible(list(tbl.regulatoryRegions=tbl.regulatoryRegions,
                  tbl.geneModel=tbl.geneModel))

} # test_createGeneModel
#------------------------------------------------------------------------------------------------------------------------
test_buildMultiModelGraph_oneModel <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_oneModel")
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   fp.source <- "postgres://whovian/brain_hint_20"
   sources <- list(fp.source)

   prep <- TrenaPrep(targetGene, aqp4.tss, "chr18", aqp4.tss-100, aqp4.tss+100, regulatoryRegionSources=sources)
   x <- getRegulatoryRegions(prep)
   closeAllPostgresConnections()
   tbl.regulatoryRegions <- expandRegulatoryRegionsTableByTF(prep, x[[fp.source]])

   tbl.geneModel <- createGeneModel(prep, "randomForest", tbl.regulatoryRegions, mtx)
   tbl.geneModel.strong <- subset(tbl.geneModel, randomForest > 3)
   tbl.regulatoryRegions.strong <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.strong$tf)

   models <- list(rf3=list(tbl.regulatoryRegions=tbl.regulatoryRegions.strong, tbl.geneModel=tbl.geneModel.strong))

   #save(models, file="testModel.RData")
   #load("testModel.RData")

   g <- buildMultiModelGraph(prep, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(models[[1]]$tbl.regulatoryRegions$id)
   tfNodes <- unique(models[[1]]$tbl.regulatoryRegions$tf)
   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))
   tbl.reg <- models[[1]]$tbl.regulatoryRegions
   checkEquals(length(edgeNames(g)), nrow(tbl.reg) + length(unique(tbl.reg$id)))

   g.lo <- addGeneModelLayout(prep, g, xPos.span=1500)
   min.xPos <- min(as.numeric(nodeData(g.lo, attr="xPos")))
   max.xPos <- max(as.numeric(nodeData(g.lo, attr="xPos")))
   checkEquals(abs(max.xPos - min.xPos), 1500)

   if(display){
     tViz <<- TrenaViz()
     httpAddGraph(tViz, g.lo, names(models))
     loadStyle(tViz, system.file(package="trenaUtilities", "extdata", "style.js"))
     Sys.sleep(3); fit(tViz)
     browser()
     }

} # test_buildMultiModelGraph_oneModel
#------------------------------------------------------------------------------------------------------------------------
test_buildMultiModelGraph_fiveModels <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_fiveModels")
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   fp.source <- "postgres://whovian/brain_hint_20"
   sources <- list(fp.source)

   prep <- TrenaPrep(targetGene, aqp4.tss, "chr18", aqp4.tss-1000, aqp4.tss+1000, regulatoryRegionSources=sources)
   x <- getRegulatoryRegions(prep)
   closeAllPostgresConnections()
   tbl.regulatoryRegions <- expandRegulatoryRegionsTableByTF(prep, x[[fp.source]])

   tbl.geneModel <- createGeneModel(prep, "randomForest", tbl.regulatoryRegions, mtx)

      # two get multiple models, filter on randomForest score
   tbl.geneModel.rf10 <- subset(tbl.geneModel, randomForest > 10)
   tbl.regulatoryRegions.rf10 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf10$tf)

   tbl.geneModel.rf5 <- subset(tbl.geneModel, randomForest > 5)
   tbl.regulatoryRegions.rf5 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf5$tf)


   tbl.geneModel.rf3 <- subset(tbl.geneModel, randomForest > 3)
   tbl.regulatoryRegions.rf3 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf3$tf)

   tbl.geneModel.rf2 <- subset(tbl.geneModel, randomForest > 2)
   tbl.regulatoryRegions.rf2 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf2$tf)

   tbl.geneModel.rf1 <- subset(tbl.geneModel, randomForest > 1)
   tbl.regulatoryRegions.rf1 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf1$tf)

   models <- list(rf01=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf1,  tbl.geneModel=tbl.geneModel.rf1),
                  rf02=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf2,  tbl.geneModel=tbl.geneModel.rf2),
                  rf03=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf3,  tbl.geneModel=tbl.geneModel.rf3),
                  rf05=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf5,  tbl.geneModel=tbl.geneModel.rf5),
                  rf10=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf10, tbl.geneModel=tbl.geneModel.rf10)
                  )

   #save(models, file="testModel.RData")
   #load("testModel.RData")

   g <- buildMultiModelGraph(prep, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(models[[1]]$tbl.regulatoryRegions$id)
   tfNodes <- unique(models[[1]]$tbl.regulatoryRegions$tf)
   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))
   tbl.reg <- models[[1]]$tbl.regulatoryRegions
   #checkEquals(length(edgeNames(g)), nrow(tbl.reg) + length(unique(tbl.reg$id)))

   g.lo <- addGeneModelLayout(prep, g, xPos.span=1500)
   min.xPos <- min(as.numeric(nodeData(g.lo, attr="xPos")))
   max.xPos <- max(as.numeric(nodeData(g.lo, attr="xPos")))
   checkEquals(abs(max.xPos - min.xPos), 1500)

   if(display){
     tViz <<- TrenaViz()
     httpAddGraph(tViz, g.lo, names(models))
     loadStyle(tViz, system.file(package="trenaUtilities", "extdata", "style.js"))
     Sys.sleep(3); fit(tViz)
     browser()
     }

} # test_buildMultiModelGraph_fiveModels
#------------------------------------------------------------------------------------------------------------------------
test_buildMultiModelGraph_twoModels_15k_span <- function(display=FALSE)
{
   printf("--- test_buildMultiModelGraph_twoModels_15k_span")
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   fp.source <- "postgres://whovian/brain_hint_20"
   sources <- list(fp.source)

   prep <- TrenaPrep(targetGene, aqp4.tss, "chr18", aqp4.tss-5000, aqp4.tss+10000, regulatoryRegionSources=sources)
   x <- getRegulatoryRegions(prep)
   closeAllPostgresConnections()
   tbl.regulatoryRegions <- expandRegulatoryRegionsTableByTF(prep, x[[fp.source]])

   tbl.geneModel <- createGeneModel(prep, "randomForest", tbl.regulatoryRegions, mtx)

      # two get multiple models, filter on randomForest score
   tbl.geneModel.rf10 <- subset(tbl.geneModel, randomForest > 10)
   tbl.regulatoryRegions.rf10 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf10$tf)

   tbl.geneModel.rf1 <- subset(tbl.geneModel, randomForest > 1)
   tbl.regulatoryRegions.rf1 <- subset(tbl.regulatoryRegions, tf %in% tbl.geneModel.rf1$tf)


   models <- list(rf01=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf1,  tbl.geneModel=tbl.geneModel.rf1),
                  rf10=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rf10, tbl.geneModel=tbl.geneModel.rf10)
                  )

   #save(models, file="models.2.big.RData")
   #load("models.2.big.RData")


   g <- buildMultiModelGraph(prep, models)
   nodesInGraph <- nodes(g)
   regionNodes <- unique(models[[1]]$tbl.regulatoryRegions$id)
   tfNodes <- unique(models[[1]]$tbl.regulatoryRegions$tf)
   checkEquals(length(nodesInGraph), length(regionNodes) + length(tfNodes) + length(targetGene))
   tbl.reg <- models[[1]]$tbl.regulatoryRegions
   #checkEquals(length(edgeNames(g)), nrow(tbl.reg) + length(unique(tbl.reg$id)))

   g.lo <- addGeneModelLayout(prep, g, xPos.span=1500)
   min.xPos <- min(as.numeric(nodeData(g.lo, attr="xPos")))
   max.xPos <- max(as.numeric(nodeData(g.lo, attr="xPos")))
   checkEquals(abs(max.xPos - min.xPos), 1500)

   if(display){
     tViz <<- TrenaViz()
     httpAddGraph(tViz, g.lo, names(models))
     loadStyle(tViz, system.file(package="trenaUtilities", "extdata", "style.js"))
     Sys.sleep(3); fit(tViz)
     browser()
     }

} # test_buildMultiModelGraph_fiveModels
#------------------------------------------------------------------------------------------------------------------------
test_geneModelLayout <- function()
{
   printf("--- test_geneModelLayout")
   filename <- "testModel.RData";
   if(!file.exists(filename)){
      printf("   skipping test_geneModelLayout, no file named '%s'", filename)
      return()
      }

   load(filename)
   print(load("ohsu.aqp4.graphLayoutNaNBug.input.RData"))
   # g <- buildMultiModelGraph(prep, models)
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   fp.source <- "postgres://whovian/brain_hint_20"
   sources <- list(fp.source)

   prep <- TrenaPrep(targetGene, aqp4.tss, "chr18", aqp4.tss-5000, aqp4.tss+10000, regulatoryRegionSources=sources)

   g.lo <- addGeneModelLayout(prep, g, xPos.span=1500)

} # test_geneModelLayout
#------------------------------------------------------------------------------------------------------------------------
test_assessSnp <- function()
{
   printf("--- test_assessSnp")
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   fp.source <- "postgres://whovian/brain_hint_20"
   sources <- list(fp.source)
   prep <- TrenaPrep(targetGene, aqp4.tss, "chr18", aqp4.tss-5000, aqp4.tss+10000, regulatoryRegionSources=sources)

   tbl.assay <- assessSnp(prep, pfms=list(), "rs3875089", 10, pwmMatchMinimumAsPercentage=80)
   checkEquals(ncol(tbl.assay), 13)
      # check a few of the columns including the last-added, "variant"
   checkTrue(all(c("status", "assessed", "delta", "tf", "variant") %in% colnames(tbl.assay)))
   tbl.changed <- subset(tbl.assay, assessed != "in.both")

      # only 3 motifs, all in the same location (+/- 1 base) are broken,
   checkEquals(nrow(tbl.changed), 3)
   checkEquals(length(grep("TEAD1", tbl.changed$tf)), 3)
   checkEquals(sort(tbl.changed$signature), c("MA0090.2;26865465;+", "MA0808.1;26865466;+", "MA0809.1;26865465;+"))

      # now lets look deeper and see how bad the mut sequence scored
   tbl.assay.deeper <- assessSnp(prep, pfms=list(), "rs3875089", 10, pwmMatchMinimumAsPercentage=50)
   tbl.toCompare <- subset(tbl.assay.deeper, signature %in% tbl.changed$signature)
   wt.mean.match.score  <- mean(subset(tbl.toCompare, status=="wt")$motifRelativeScore)  # [1] 0.8201358
   mut.mean.match.score <- mean(subset(tbl.toCompare, status=="mut")$motifRelativeScore) # [1] 0.6746936
      # > 20% drop off in match
   dropOff <- abs((mut.mean.match.score - wt.mean.match.score)/mut.mean.match.score)
   checkTrue(dropOff > 0.20)

   suppressWarnings(tbl.assay.short <- assessSnp(prep, pfms=list(), "rs3875089", 3, pwmMatchMinimumAsPercentage=95))
   checkEquals(nrow(tbl.assay.short), 0)

} # test_assessSnp
#------------------------------------------------------------------------------------------------------------------------
# in preparation for adding, and ongoing testing of, a delta column for all entries, here we use a snp which at the 80%
# match level # returns all three kinds of match: in.both. wt.only, mut.only
test_assessSnp_mutOnly_wtOnly <- function()
{
   printf("--- test_assessSnp_mutOnly_wtOnly")

   snp <- "rs3763043"
   targetGene <- "AQP4"
   aqp4.tss <- 26865884
   fp.source <- "postgres://whovian/brain_hint_20"
   sources <- list(fp.source)
   prep <- TrenaPrep(targetGene, aqp4.tss, "chr18", aqp4.tss-5000, aqp4.tss+10000, regulatoryRegionSources=sources)

   tbl.assay <- assessSnp(prep, pfms=list(), snp, shoulder=10, pwmMatchMinimumAsPercentage=80)

   checkEquals(ncol(tbl.assay), 13)
   checkTrue("delta" %in% colnames(tbl.assay))
   checkEqualsNumeric(min(tbl.assay$delta), -0.0487, tol=1e-3)
   checkEqualsNumeric(max(tbl.assay$delta),  0.165, tol=1e-2)

} # test_assessSnp_mutOnly_wtOnly
#------------------------------------------------------------------------------------------------------------------------
test_assessSnpMotifDbMatrices <- function()
{
   snp <- "rs9357347"
   snp.loc <-  list(chrom="chr6", start=41182853L, end=41182853L)  #  A>C
   targetGene <- "TREM2"
   tss <- 41163190   # - strand
   chrom <- "chr6"
   fp.source <- "postgres://whovian/brain_hint_20"
   sources <- list(fp.source)
   prep <- TrenaPrep(targetGene, tss, chrom, tss-5000, tss+10000, regulatoryRegionSources=sources)

   pfms <- as.list(query(MotifDb, "MESP1"))
   tbl.assay <- assessSnp(prep, pfms=pfms, snp, shoulder=8, pwmMatchMinimumAsPercentage=70)

   hocomoco.human <- query(query(MotifDb, "sapiens"), "hocomoco")
   tbl.assay.hohu <- assessSnp(prep, pfms=as.list(hocomoco.human), snp, shoulder=8, pwmMatchMinimumAsPercentage=80)

   mesp1.all <- query(MotifDb, "MESP1")
   tbl.assay.mesp1 <- assessSnp(prep, pfms=as.list(mesp1.all), snp, shoulder=8, pwmMatchMinimumAsPercentage=70)


} # test_assessSnpMotifDbMatrices
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
