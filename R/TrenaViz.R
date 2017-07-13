#----------------------------------------------------------------------------------------------------
cyjsBrowserFile <- system.file(package="trenaUtilities", "scripts", "trenaViz.html")
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
.TrenaViz <- setClass ("TrenaViz",
                        representation = representation(model="data.frame"),
                        contains = "BrowserVizClass",
                        prototype = prototype (uri="http://localhost", 9000)
                        )

#----------------------------------------------------------------------------------------------------
setGeneric('addGraph',         signature='obj', function(obj, graph, modelNames=list()) standardGeneric ('addGraph'))
setGeneric('httpAddGraph',     signature='obj', function(obj, graph, modelNames=list()) standardGeneric ('httpAddGraph'))
setGeneric('httpAddStructureGraph',     signature='obj', function(obj, graph, modelNames=list()) standardGeneric ('httpAddStructureGraph'))
setGeneric('loadStyle',        signature='obj', function(obj, filename) standardGeneric ('loadStyle'))
setGeneric('loadStructureStyle',        signature='obj', function(obj, filename) standardGeneric ('loadStructureStyle'))
setGeneric('raiseTab',         signature='obj', function(obj, tabName) standardGeneric('raiseTab'))
setGeneric('fit',              signature='obj', function(obj, padding=30) standardGeneric('fit'))
setGeneric('fitSelected',      signature='obj', function(obj, padding=30) standardGeneric('fitSelectedContent'))
setGeneric('selectNodes',      signature='obj', function(obj, nodeIDs) standardGeneric('selectNodes'))
setGeneric('setNodeAttributes',   signature='obj', function(obj, attribute, nodes, values, cyInstance="default") standardGeneric('setNodeAttributes'))
setGeneric('selectStructureNodes',      signature='obj', function(obj, nodeIDs) standardGeneric('selectStructureNodes'))
setGeneric('getSelectedNodes', signature='obj', function(obj) standardGeneric('getSelectedNodes'))
setGeneric('clearSelection',   signature='obj', function(obj) standardGeneric('clearSelection'))
setGeneric('sfn',              signature='obj', function(obj) standardGeneric('sfn'))
setGeneric('addBedTrackFromDataFrame',   signature='obj', function(obj, trackName, tbl.bed, displayMode="COLLAPSED", color)
                                  standardGeneric('addBedTrackFromDataFrame'))
setGeneric('addBedTrackFromHostedFile',   signature='obj',
                        function(obj, trackName, uri, index.uri=NA, displayMode="COLLAPSED", color)
                                  standardGeneric('addBedTrackFromHostedFile'))
setGeneric('showGenomicRegion',   signature='obj', function(obj, regionString) standardGeneric('showGenomicRegion'))
setGeneric('getGenomicRegion',    signature='obj', function(obj, regionString) standardGeneric('getGenomicRegion'))
setGeneric('layout',              signature='obj', function(obj, strategy) standardGeneric('layout'))
setGeneric('layoutStrategies',    signature='obj', function(obj) standardGeneric('layoutStrategies'))
setGeneric('geneRegulatoryModelToGraph',    signature='obj', function(obj, target.gene, tbl.gm, tbl.reg)
                                  standardGeneric('geneRegulatoryModelToGraph'))
#setGeneric('buildMultiModelGraph', signature='obj', function(obj, target.gene, models) standardGeneric('buildMultiModelGraph'))
#setGeneric('addGeneModelLayout', signature='obj', function(obj, g, xPos.span) standardGeneric('addGeneModelLayout'))
#----------------------------------------------------------------------------------------------------
# constructor
trenaViz = function(portRange=11000:11025, host="localhost", title="TReNA-Viz", quiet=TRUE)
{

   model <- data.frame()
   obj <- .TrenaViz(BrowserViz(portRange=portRange, host=host, title=title,
                                quiet=quiet, browserFile=cyjsBrowserFile,
                                httpQueryProcessingFunction=myQP),
                     model=model)

  while (!browserResponseReady(obj)){
      Sys.sleep(.1)
      }
   if(!quiet) {
      message(sprintf("BrowserViz ctor called from TReNA-Viz ctor got browser response"))
      print(getBrowserResponse(obj))
      }

   obj

} # trenaViz constructor
#----------------------------------------------------------------------------------------------------
setMethod('addGraph', 'TrenaViz',

  function (obj, graph, modelNames=list()) {
     printf("TrenaViz::addGraph");
     print(graph)
     printf("--- calling .graphToJSON");
     g.json <- .graphToJSON(graph)
     payload <- list(graph=g.json, modelNames=modelNames)
     printf("about to send g.json: %d chars", nchar(g.json));
     send(obj, list(cmd="addGraph", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('httpAddGraph', 'TrenaViz',

  function (obj, graph, modelNames=list()) {
     printf("TrenaViz::httpAddGraph");
     print(graph)
     printf("--- converting graph to JSON");
     g.json <- .graphToJSON(graph)
     printf("--- conversion complete");
     #g.json <- paste("network = ", .graphToJSON(graph))
     #g.json <- paste("network = ", as.character(biocGraphToCytoscapeJSON(graph)))
     filename <- "g.json"
     payload <- list(graph=filename, modelNames=modelNames)
     printf("--- about to write file 'g.json' with %d characters", nchar(g.json))
     printf("--- first few characters: %s", substr(g.json, 1, 20))
     write(g.json, file=filename)
     printf("--- file writing complete")
     send(obj, list(cmd="httpAddGraph", callback="handleResponse", status="request",
                    payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('httpAddStructureGraph', 'TrenaViz',

  function (obj, graph, modelNames=list()) {
     printf("TrenaViz::httpAddStructureGraph");
     print(graph)
     printf("--- converting graph to JSON");
     g.json <- .graphToJSON(graph)
     printf("--- conversion complete");
     #g.json <- paste("network = ", .graphToJSON(graph))
     #g.json <- paste("network = ", as.character(biocGraphToCytoscapeJSON(graph)))
     filename <- "g.json"
     payload <- list(graph=filename, modelNames=modelNames)
     printf("--- about to write file 'g.json' with %d characters", nchar(g.json))
     printf("--- first few characters: %s", substr(g.json, 1, 20))
     write(g.json, file=filename)
     printf("--- file writing complete")
     send(obj, list(cmd="httpAddStructureGraph", callback="handleResponse", status="request",
                    payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('loadStyle', 'TrenaViz',

  function (obj, filename) {
     printf("TrenaViz::loadStyle");
     send(obj, list(cmd="httpSetStyle", callback="handleResponse", status="request", payload=filename))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('loadStructureStyle', 'TrenaViz',

  function (obj, filename) {
     printf("TrenaViz::loadStructureStyle");
     send(obj, list(cmd="httpSetStructureStyle", callback="handleResponse", status="request", payload=filename))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('raiseTab', 'TrenaViz',

  function (obj, tabName) {
     send(obj, list(cmd="raiseTab", callback="handleResponse", status="request", payload=tabName))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('addBedTrackFromDataFrame', 'TrenaViz',

  function (obj, trackName, tbl.bed, displayMode="COLLAPSED", color) {
     printf("TrenaViz::addBedTrackFromDataFrame");
     temp.filename <- "tmp.bed"
     write.table(tbl.bed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
     payload <- list(name=trackName, bedFileName=temp.filename, displayMode=displayMode, color=color)
     send(obj, list(cmd="addBedTrackFromDataFrame", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('addBedTrackFromHostedFile', 'TrenaViz',

  function (obj, trackName, uri, index.uri, displayMode="COLLAPSED", color) {
     printf("TrenaViz::addBedTrackFromHostedFile");
     payload <- list(name=trackName, uri=uri, indexUri=index.uri, displayMode=displayMode, color=color)
     send(obj, list(cmd="addBedTrackFromHostedFile", callback="handleResponse",
                    status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('showGenomicRegion', 'TrenaViz',

   function (obj, regionString) {
     payload <- list(regionString=regionString)
     send(obj, list(cmd="showGenomicRegion", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('getGenomicRegion', 'TrenaViz',

   function (obj, regionString) {
     payload <- ""
     send(obj, list(cmd="getGenomicRegion", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('fit', 'TrenaViz',

  function (obj, padding=30) {
     send(obj, list(cmd="fit", callback="handleResponse", status="request", payload=padding))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('fitSelected', 'TrenaViz',

  function (obj, padding=30) {
     send(obj, list(cmd="fitSelected", callback="handleResponse", status="request", payload=padding))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('selectNodes', 'TrenaViz',

  function (obj, nodeIDs) {
     payload <- list(nodeIDs=nodeIDs)
     send(obj, list(cmd="selectNodes", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('selectStructureNodes', 'TrenaViz',

  function (obj, nodeIDs) {
     payload <- list(nodeIDs=nodeIDs)
     send(obj, list(cmd="selectStructureNodes", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('getSelectedNodes', 'TrenaViz',

  function (obj) {
     payload <- ""
     send(obj, list(cmd="getSelectedNodes", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     fromJSON(getBrowserResponse(obj))$id;
     })

#----------------------------------------------------------------------------------------------------
setMethod('sfn', 'TrenaViz',

  function (obj) {
     payload <- ""
     send(obj, list(cmd="sfn", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('clearSelection', 'TrenaViz',

  function (obj) {
     payload <- ""
     send(obj, list(cmd="clearSelection", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('layoutStrategies', 'TrenaViz',

  function (obj) {
     send(obj, list(cmd="layoutStrategies", callback="handleResponse", status="request",
                                  payload=""))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj)
     })

#----------------------------------------------------------------------------------------------------
setMethod('layout', 'TrenaViz',

  function (obj, strategy="random") {
     send(obj, list(cmd="doLayout", callback="handleResponse", status="request",
                                  payload=strategy))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj)
     })

#----------------------------------------------------------------------------------------------------
setMethod('setNodeAttributes', 'TrenaViz',

   function(obj, attribute, nodes, values, cyInstance="default"){

     if (length (nodes) == 0)
       return ()

     if(length(values) == 1)
        values <- rep(values, length(nodes))

     payload <- list(attribute=attribute, nodes=nodes, values=values, cyInstance=cyInstance)
     send(obj, list(cmd="setNodeAttributes", callback="handleResponse", status="request",
                    payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     result <- getBrowserResponse(obj)
     if(nchar(result) > 0)
       return(fromJSON(getBrowserResponse(obj)))
     else
       invisible("")
     }) # setNodeAttributes

#------------------------------------------------------------------------------------------------------------------------
# {elements: [
#    {data: {id: 'a', score:5}, position: {x: 100, y: 200}},
#    {data: {id: 'b', score:100}, position: {x: 200, y: 200}},
#    {data: {id: 'e1', source: 'a', target: 'b'}}
#    ],  // elements array
# layout: { name: 'preset'},
# style: [{selector: 'node', style: {'content': 'data(id)'}}]
# }
.graphToJSON <- function(g)
{
    x <- '{"elements": [';
    nodes <- nodes(g)
    edgeNames <- edgeNames(g)
    edges <- strsplit(edgeNames, "~")  # a list of pairs
    edgeNames <- sub("~", "->", edgeNames)
    names(edges) <- edgeNames

    noa.names <- names(nodeDataDefaults(g))
    eda.names <- names(edgeDataDefaults(g))
    nodeCount <- length(nodes)
    edgeCount <- length(edgeNames)

    for(n in 1:nodeCount){
       node <- nodes[n]
       x <- sprintf('%s {"data": {"id": "%s"', x, node);
       nodeAttributeCount <- length(noa.names)
       for(i in seq_len(nodeAttributeCount)){
          noa.name <- noa.names[i];
          value <-  nodeData(g, node, noa.name)[[1]]
          if(is.numeric(value))
             x <- sprintf('%s, "%s": %s', x, noa.name, value)
          else
             x <- sprintf('%s, "%s": "%s"', x, noa.name, value)
          } # for i
       x <- sprintf('%s}', x)     # close off this node data element
       if(all(c("xPos", "yPos") %in% noa.names)){
           xPos <- as.integer(nodeData(g, node, "xPos"))
           yPos <- as.integer(nodeData(g, node, "yPos"))
           x <- sprintf('%s, "position": {"x": %d, "y": %d}', x, xPos, yPos)
           } # add position element
       x <- sprintf('%s}', x)     # close off this node data element
       if(n != nodeCount)
           x <- sprintf("%s,", x)  # another node coming, add a comma
       } # for n

    for(e in seq_len(edgeCount)) {
       edgeName <- edgeNames[e]
       edge <- edges[[e]]
       sourceNode <- edge[[1]]
       targetNode <- edge[[2]]
       x <- sprintf('%s, {"data": {"id": "%s", "source": "%s", "target": "%s"', x, edgeName, sourceNode, targetNode);
       edgeAttributeCount <- length(eda.names)
       for(i in seq_len(edgeAttributeCount)){
          eda.name <- eda.names[i];
          value <-  edgeData(g, sourceNode, targetNode, eda.name)[[1]]
          if(is.numeric(value))
             x <- sprintf('%s, "%s": %s', x, eda.name, value)
          else
             x <- sprintf('%s, "%s": "%s"', x, eda.name, value)
          } # for each edgeAttribute
       x <- sprintf('%s}}', x)     # close off this edge data element
       } # for e

    x <- sprintf("%s]}", x)

    x

} # .graphToJSON
#------------------------------------------------------------------------------------------------------------------------
myQP <- function(queryString)
{
   printf("=== TReNA-Viz::myQP");
   #print(queryString)
     # for reasons not quite clear, the query string comes in with extra characters
     # following the expected filename:
     #
     #  "?sampleStyle.js&_=1443650062946"
     #
     # check for that, cleanup the string, then see if the file can be found

   ampersand.loc <- as.integer(regexpr("&", queryString, fixed=TRUE))
   #printf("ampersand.loc: %d", ampersand.loc)

   if(ampersand.loc > 0){
      queryString <- substring(queryString, 1, ampersand.loc - 1);
      }

   questionMark.loc <- as.integer(regexpr("?", queryString, fixed=TRUE));
   #printf("questionMark.loc: %d", questionMark.loc)

   if(questionMark.loc == 1)
      queryString <- substring(queryString, 2, nchar(queryString))

   filename <- queryString;
   #printf("myQP filename: '%s'", filename)
   #printf("       exists?  %s", file.exists(filename));

   stopifnot(file.exists(filename))

   printf("--- about to scan %s", filename);
      # reconstitute linefeeds though collapsing file into one string, so json
      # structure is intact, and any "//" comment tokens only affect one line
   text <- paste(scan(filename, what=character(0), sep="\n", quiet=TRUE), collapse="\n")
   printf("%d chars read from %s", nchar(text), filename);

   return(text);

} # myQP
#----------------------------------------------------------------------------------------------------
setMethod('geneRegulatoryModelToGraph', 'TrenaViz',

  function (obj, target.gene, tbl.gm, tbl.reg) {

     required.geneModelColumnNames <- c("tf", "pearson", "spearman", "betaLasso", "randomForest", "pcaMax", "concordance")
     required.regulatoryRegionsColumnNames <- c("motifName", "chrom", "motifStart", "motifEnd", "strand",
                                                "motifScore", "motifRelativeScore", "match",
                                                "distance.from.tss",
                                                "chromStart", "chromEnd", "seq", "status", "tf")
     stopifnot(all(required.geneModelColumnNames %in% colnames(tbl.gm)))
     stopifnot(all(required.regulatoryRegionsColumnNames %in% colnames(tbl.reg)))

     printf("genes: %d, %d occurences of %d motifs", length(tbl.gm$tf), length(tbl.reg$motifName),
            length(unique(tbl.reg$motifName)))

     g <- graphNEL(edgemode = "directed")

     nodeDataDefaults(g, attr = "type") <- "undefined"             # targetGene, tf, footprint
     nodeDataDefaults(g, attr = "label") <- "default node label"
     nodeDataDefaults(g, attr = "distance") <- 0
     nodeDataDefaults(g, attr = "pearson") <- 0
     nodeDataDefaults(g, attr = "randomForest") <- 0
     nodeDataDefaults(g, attr = "pcaMax") <- 0
     nodeDataDefaults(g, attr = "concordance") <- 0
     nodeDataDefaults(g, attr = "betaLasso") <- 0
     nodeDataDefaults(g, attr = "motif") <- ""
     nodeDataDefaults(g, attr = "xPos") <- 0
     nodeDataDefaults(g, attr = "yPos") <- 0

     edgeDataDefaults(g, attr = "edgeType") <- "undefined"

     tfs <- tbl.gm$tf

     regRegions.names <- unlist(lapply(1:nrow(tbl.reg), function(i){
         distance.from.tss <- tbl.reg$distance.from.tss[i]
         region.size <- nchar(tbl.reg$match[i])
         motif.name <- tbl.reg$motifName[i]
         if(distance.from.tss < 0)
            sprintf("%s.fp.downstream.%05d.L%d.%s", target.gene, abs(distance.from.tss), region.size, motif.name)
          else
            sprintf("%s.fp.upstream.%05d.L%d.%s", target.gene, abs(distance.from.tss), region.size, motif.name)
          }))

   tbl.reg$regionName <- regRegions.names
   all.nodes <- unique(c(target.gene, tfs, regRegions.names))
   g <- addNode(all.nodes, g)

   nodeData(g, target.gene, "type") <- "targetGene"
   nodeData(g, tfs, "type")         <- "TF"
   nodeData(g, regRegions.names, "type")  <- "regulatoryRegion"
   nodeData(g, all.nodes, "label")  <- all.nodes
   nodeData(g, regRegions.names, "label") <- tbl.reg$motifName
   nodeData(g, regRegions.names, "distance") <- tbl.reg$distance
   nodeData(g, regRegions.names, "motif") <- tbl.reg$motifName

   nodeData(g, tfs, "pearson") <- tbl.gm$pearson
   nodeData(g, tfs, "betaLasso") <- tbl.gm$betaLasso
   nodeData(g, tfs, "randomForest") <- tbl.gm$randomForest
   nodeData(g, tfs, "pcaMax") <- tbl.gm$pcaMax
   nodeData(g, tfs, "concordance") <- tbl.gm$concordance

   #browser()
   g <- addEdge(tbl.reg$tf, tbl.reg$regionName, g)
   edgeData(g,  tbl.reg$tf, tbl.reg$regionName, "edgeType") <- "bindsTo"

   g <- graph::addEdge(tbl.reg$regionName, target.gene, g)
   edgeData(g, tbl.reg$regionName, target.gene, "edgeType") <- "regulatorySiteFor"

   g

   }) # geneRegulatoryModelToGraph

#------------------------------------------------------------------------------------------------------------------------
# {elements: [
#    {data: {id: 'a', score:5}, position: {x: 100, y: 200}},
#    {data: {id: 'b', score:100}, position: {x: 200, y: 200}},
#    {data: {id: 'e1', source: 'a', target: 'b'}}
#    ],  // elements array
# layout: { name: 'preset'},
# style: [{selector: 'node', style: {'content': 'data(id)'}}]
# }
graphToJSON <- function(g)
{
    x <- '{"elements": [';
    nodes <- nodes(g)
    edgeNames <- edgeNames(g)
    edges <- strsplit(edgeNames, "~")  # a list of pairs
    edgeNames <- sub("~", "->", edgeNames)
    names(edges) <- edgeNames

    noa.names <- names(nodeDataDefaults(g))
    eda.names <- names(edgeDataDefaults(g))
    nodeCount <- length(nodes)
    edgeCount <- length(edgeNames)

    for(n in 1:nodeCount){
       #printf("---- node %d", n)
       node <- nodes[n]
       x <- sprintf('%s {"data": {"id": "%s"', x, node);
       nodeAttributeCount <- length(noa.names)
       for(i in seq_len(nodeAttributeCount)){
          noa.name <- noa.names[i];
          #printf("node %s, noa.name: %s", node, noa.name)
          value <-  nodeData(g, node, noa.name)[[1]]
          if(is.numeric(value))
             x <- sprintf('%s, "%s": %s', x, noa.name, value)
          else
             x <- sprintf('%s, "%s": "%s"', x, noa.name, value)
          #browser();
          #xyz <- 99
          } # for i
       x <- sprintf('%s}', x)     # close off this node data element
       #printf("-- x partway: %s", x)
       if(all(c("xPos", "yPos") %in% noa.names)){
           xPos <- as.integer(nodeData(g, node, "xPos"))
           yPos <- as.integer(nodeData(g, node, "yPos"))
           x <- sprintf('%s, "position": {"x": %d, "y": %d}', x, xPos, yPos)
           #xyz <- 99
           } # add position element
       x <- sprintf('%s}', x)     # close off this node data element
       if(n != nodeCount)
           x <- sprintf("%s,", x)  # another node coming, add a comma
       } # for n

    #browser()
    #xyz <- 99

    for(e in seq_len(edgeCount)) {
       edgeName <- edgeNames[e]
       edge <- edges[[e]]
       sourceNode <- edge[[1]]
       targetNode <- edge[[2]]
       x <- sprintf('%s, {"data": {"id": "%s", "source": "%s", "target": "%s"', x, edgeName, sourceNode, targetNode);
       edgeAttributeCount <- length(eda.names)
       for(i in seq_len(edgeAttributeCount)){
          eda.name <- eda.names[i];
          value <-  edgeData(g, sourceNode, targetNode, eda.name)[[1]]
          if(is.numeric(value))
             x <- sprintf('%s, "%s": %s', x, eda.name, value)
          else
             x <- sprintf('%s, "%s": "%s"', x, eda.name, value)
          } # for each edgeAttribute
       x <- sprintf('%s}}', x)     # close off this edge data element
       } # for e

    x <- sprintf("%s]}", x)

    x

} # graphToJSON
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
.validModelData <- function(models)
{
   required.geneModelColumnNames <- c("tf", "pearson", "spearman", "betaLasso", "randomForest", "pcaMax", "concordance")
   required.regulatoryRegionsColumnNames <- c("motifName", "chrom", "motifStart", "motifEnd", "strand",
                                              "score", "length", "label")

   valid <- TRUE;  # be optimistic
   model.names <- names(models)
   checkTrue(length(model.names) > 0)
      # for simplicity, and especially for use in graph node attributes, we want model names with NO spaces
   noSpacesInNames <- all(grepl(" ", model.names) == FALSE)
   if(!noSpacesInNames){
      valid <- FALSE
      warning("model names should contain no spaces")
      }

   for(model.name in names(models)){
      stopifnot(all(names(models[[model.name]]) %in% c("tbl.model", "tbl.regulatoryRegions")))
      tbl.model <- models[[model.name]]$tbl.model

      missing.in.model <- setdiff(required.geneModelColumnNames, colnames(tbl.model))
      valid.0 <- length(missing.in.model) == 0;
      if(!valid.0)
         warning(sprintf("missing columns in tbl.model for %s: %s", model.name, paste(missing.in.model, collapse=", ")))
      valid <- valid & valid.0
      missing.in.regions <- setdiff(required.geneModelColumnNames, colnames(tbl.model))
      valid.1 <- length(missing.in.regions) == 0;
      if(!valid.1)
         warning(sprintf("missing columns in tbl.model for %s: %s", model.name, paste(missing.in.regions, collapse=", ")))
      valid <- valid & valid.1
      } # for model.name

    valid

} # .validModelData
#------------------------------------------------------------------------------------------------------------------------
test_validModelData <- function()
{
   printf("--- test_validModelData")

   tbl.model.test <- data.frame(tf=c("a"),
                                pearson=c(0),
                                spearman=c(0),
                                betaLasso=c(0),
                                randomForest=c(0),
                                pcaMax=c(0),
                                concordance=c(0),
                                stringsAsFactors=FALSE)

   tbl.regions.test <- data.frame(motifName=c("a"),
                                  chrom=c("a"),
                                  motifStart=c(0),
                                  motifEnd=c(0),
                                  strand=c("+"),
                                  score=c(0),
                                  length=c(10),
                                  label=c("xyz"),
                                  stringsAsFactors=FALSE
                                  )

   models <- list(wt=list(tbl.model=tbl.model.test,  tbl.regulatoryRegions=tbl.regions.test),
                  mut=list(tbl.model=tbl.model.test, tbl.regulatoryRegions=tbl.regions.test),
                  m3=list(tbl.model=tbl.model.test,  tbl.regulatoryRegions=tbl.regions.test))

   checkTrue(validModelData(models))

   names(models)[1] <- "wild type"
   suppressWarnings(checkTrue(!validModelData(models)))

} # test_validModelData
#------------------------------------------------------------------------------------------------------------------------
# setMethod('buildMultiModelGraph', 'TrenaViz',
#
#   function (obj, target.gene, models){
#
#     stopifnot(.validModelData(models))
#     g <- graphNEL(edgemode = "directed")
#     model.names <- names(models)
#
#     node.attribute.specs <- list(type="undefined",
#                                  label="default node label",
#                                  distance=0,
#                                  pearson=0,
#                                  randomForest=0,
#                                  pcaMax=0,
#                                  concordance=0,
#                                  betaLasso=0,
#                                  motif="",
#                                  xPos=0,
#                                  yPos=0)
#     edge.attribute.spec <- list(edgeType="undefined")
#     attribute.classes <- c("", model.names)  # "" (no prefix) is the currently displayed set of attibutes
#
#       # create current version of these attributes, and then
#       # per-model versions, which get mapped to current
#       # in response to user's interactive choice on the cyjs user interface
#       # the "current version" is, e.g., "distance".
#       # per-model ("wt" and "mut" versions) become "wt.distance" and "mut.distance"
#       # and are used by copying e.g. all wt.xxx attributes into the current (non-prefixed)
#       # attribute, upon which the cyjs style is defined
#
#     for(class.name in attribute.classes){
#        class.name.prefix <- class.name  # with possible "." appended, permits standard and model-specific attributes
#        if(nchar(class.name) > 0)
#           class.name.prefix <- sprintf("%s.", class.name)
#        noa.names.without.prefix <- names(node.attribute.specs)
#        noa.names <- sprintf("%s%s", class.name.prefix, noa.names.without.prefix)
#        noa.count <- length(node.attribute.specs)
#        for(i in 1:noa.count){
#           nodeDataDefaults(g, attr=noa.names[i]) <- node.attribute.specs[[noa.names.without.prefix[i]]]
#           }
#        } # for class
#
#     edgeDataDefaults(g, attr = "edgeType") <- "undefined"
#
#     tfs <- c()
#     regulatoryRegions <- c()
#
#     for(model in models){  # collect all the tf and regulatory region nodes
#        tbl.model <- model$tbl.model
#        tfs <- unique(c(tfs, tbl.model$tf))
#        tbl.reg <- model$tbl.regulatoryRegions
#        regulatoryRegions <- unique(c(regulatoryRegions, tbl.reg$id))
#        } # for model
#
#     all.nodes <- unique(c(target.gene, tfs, regulatoryRegions))
#     g <- addNode(all.nodes, g)
#
#     browser()
#     nodeData(g, target.gene, "type") <- "targetGene"
#     nodeData(g, tfs, "type")         <- "TF"
#     nodeData(g, regulatoryRegions, "type")  <- "regulatoryRegion"
#     nodeData(g, all.nodes, "label")  <- all.nodes
#
#       # add edges, edge attribute, and the constant attributes for all of the regulatoryRegion nodes
#
#     for(model in models){
#        tfs <- model$tbl.regulatoryRegions$tf
#        regRegions <- model$tbl.regulatoryRegions$id
#        suppressWarnings(g <- addEdge(tfs, regRegions, g))
#        edgeData(g,  tfs, regRegions, "edgeType") <- "bindsTo"
#        suppressWarnings(g <- addEdge(regRegions, target.gene, g))
#        edgeData(g, regRegions, target.gene, "edgeType") <- "regulatorySiteFor"
#        nodeData(g, tbl.reg$id, "label") <- tbl.reg$motifName
#        nodeData(g, tbl.reg$id, "distance") <- tbl.reg$distance.from.tss
#        nodeData(g, tbl.reg$id, "motif") <- tbl.reg$motifName
#        } # for model
#
#       # now copy in the first model's tf node data
#
#     tbl.model <- models[[1]]$tbl.model
#     nodeData(g, tbl.model$tf, attr="randomForest") <- tbl.model$randomForest
#     nodeData(g, tbl.model$tf, attr="pearson") <- tbl.model$pearson
#
#      # now copy in each of the model's tf node data in turn
#     model.names <- names(models)
#     for(model.name in model.names){
#        tbl.model <- models[[model.name]]$tbl.model
#        noa.name <- sprintf("%s.%s", model.name, "randomForest")
#        nodeData(g,  tbl.model$tf, attr=noa.name) <- tbl.model$randomForest
#        noa.name <- sprintf("%s.%s", model.name, "pearson")
#        nodeData(g,  tbl.model$tf, attr=noa.name) <- tbl.model$pearson
#       } # for model.name
#
#     g
#
#     }) # buildMultiModelGraph
#
# #------------------------------------------------------------------------------------------------------------------------
