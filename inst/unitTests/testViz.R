library(TReNA)
library(TrenaHelpers)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function(display=FALSE)
{
   if(display){
      test_multiple.geneRegulatoryModelsToGraph(display=display);
      browser();
      xyz <- 99;
      }
   else{
      test_multiple.geneRegulatoryModelsToGraph(display=FALSE);
      }

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_multiple.geneRegulatoryModelsToGraph <- function(display=FALSE)
{
   printf("--- test_multiple.geneRegulatoryModelToGraph")
   load(system.file(package="TReNA", "extdata", "twoAQP4modelsForTesting.RData"))
   tViz <- TrenaViz()

   models <- list(wt=list(tbl.model=x.wt$tbl.model, tbl.regulatoryRegions=x.wt$tbl.regulatoryRegions),
                  rs3875089=list(tbl.model=x.mut$tbl.model, tbl.regulatoryRegions=x.mut$tbl.regulatoryRegions))
   g <- buildMultiModelGraph(tViz, "AQP4", models)
   checkEquals(length(nodes(g)), 34)
   checkEquals(length(edgeNames(g)), 46)

   noa.names <- sort(names(nodeDataDefaults(g)))
   checkEquals(length(noa.names), 33)
   checkEquals(length(grep("rs3875089.", noa.names, fixed=TRUE)), 11)
   checkEquals(length(grep("wt.", noa.names, fixed=TRUE)), 11)
   g.lo <- addGeneModelLayout(tViz, g)

   if(display){
      addGraph(tViz, g.lo, names(models))
      loadStyle(tViz, system.file(package="TrenaHelpers", "extdata", "style.js"))
      Sys.sleep(3); fit(tViz)
      browser();
      xyz <- 99
      }

   return(TRUE)


} # test_multiple.geneRegulatoryModelToGraph
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests(display=FALSE)
