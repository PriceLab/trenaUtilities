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
# 25 jun 2017:  updating httpAddGraph to use fetch rather than $.getScript
test_httpAddGraph_tiny <- function()
{
   printf("--- test_httpAddGraph_tiny")
   node.name <- "x"
   g <- graphNEL(nodes=node.name, edgemode="directed")
   g <- graph::addEdge(node.name, node.name, g)
   tViz <- TrenaViz()
   httpAddGraph(tViz, g);
   selectNodes(tViz, node.name)
   checkEquals(getSelectedNodes(tViz), node.name)

} # test_httpAddGraph_tiny
#------------------------------------------------------------------------------------------------------------------------
# 25 jun 2017:  updating httpAddGraph to use fetch rather than $.getScript
test_httpAddGraph_1000nodes <- function()
{
   printf("--- test_httpAddGraph_1000nodes")
   max <- 1000
   node.names <- as.character(1:max)
   tbl.edges <- data.frame(a=sample(node.names, max, replace=TRUE),
                           b=sample(node.names, max, replace=TRUE),
                           stringsAsFactors=FALSE)
   g <- graphNEL(nodes=node.names, edgemode="directed")
   g <- graph::addEdge(tbl.edges$a, tbl.edges$b, g)
   tViz <- TrenaViz()
   httpAddGraph(tViz, g);
   layout(tViz, "grid")

   clearSelection(tViz)
   nodes.to.select <- sample(node.names, 100)
   selectNodes(tViz, nodes.to.select)
   checkEquals(sort(getSelectedNodes(tViz)), sort(nodes.to.select))

} # test_httpAddGraph_1000nodes
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
