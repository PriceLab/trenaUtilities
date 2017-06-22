.TrenaPrep <- setClass ("TrenaPrep",
                        representation = representation(
                           targetGene="character",
                           chromosome="character",
                           chromStart="numeric",
                           chromEnd="numeric",
                           regulatoryRegionSources="list",
                           targetGeneTSS="numeric")
                        )
#------------------------------------------------------------------------------------------------------------------------
setGeneric('getRegulatoryRegions',  signature='obj', function(obj, combine=FALSE) standardGeneric ('getRegulatoryRegions'))
setGeneric('getRegulatoryTableColumnNames',  signature='obj', function(obj) standardGeneric ('getRegulatoryTableColumnNames'))
setGeneric('getGeneModelTableColumnNames',  signature='obj', function(obj) standardGeneric ('getGeneModelTableColumnNames'))
setGeneric('expandRegulatoryRegionsTableByTF', signature='obj', function(obj, tbl.reg) standardGeneric('expandRegulatoryRegionsTableByTF'))
setGeneric('createGeneModel', signature='obj', function(obj,  solvers, tbl.regulatoryRegions, mtx)
              standardGeneric('createGeneModel'))
setGeneric('buildMultiModelGraph', signature='obj', function(obj, models) standardGeneric('buildMultiModelGraph'))
setGeneric('addGeneModelLayout', signature='obj', function(obj, g, xPos.span=1500) standardGeneric('addGeneModelLayout'))
#------------------------------------------------------------------------------------------------------------------------
# a temporary hack: some constants
genome.db.uri <- "postgres://bddsrds.globusgenomics.org/hg38"   # has gtf and motifsgenes tables
#------------------------------------------------------------------------------------------------------------------------
TrenaPrep = function(targetGene, targetGeneTSS, chromosome, chromStart, chromEnd, regulatoryRegionSources=list(), quiet=TRUE)
{

   obj <- .TrenaPrep(targetGene=targetGene,
                     targetGeneTSS=targetGeneTSS,
                     chromosome=chromosome,
                     chromStart=chromStart,
                     chromEnd=chromEnd,
                     regulatoryRegionSources=regulatoryRegionSources)

   obj

} # constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod('getRegulatoryTableColumnNames', 'TrenaPrep',

      function(obj){
         c("chrom", "motifStart", "motifEnd", "motifName", "strand", "score", "length", "distance.from.tss", "id", "tf")
         })

#------------------------------------------------------------------------------------------------------------------------
setMethod('getGeneModelTableColumnNames', 'TrenaPrep',

      function(obj){
         c("tf", "randomForest", "pearson", "spearman", "betaLasso", "pcaMax", "concordance")
         })

#------------------------------------------------------------------------------------------------------------------------
.callFootprintFilter <- function(obj, source, chromosome, chromStart, chromEnd, targetGene, targetGeneTSS)
{
    chromLocString <- sprintf("%s:%d-%d", chromosome, chromStart, chromEnd)
    fpFilter <- FootprintFilter(genome.db.uri, source,  geneCenteredSpec=list(), regionsSpec=chromLocString)
    x.fp <- getCandidates(fpFilter)
    tbl.fp <- x.fp$tbl
    colnames(tbl.fp) <- c("chrom", "motifStart", "motifEnd", "motifName", "length", "strand", "score1", "score", "score3", "tf")
    distance <- tbl.fp$motifStart - targetGeneTSS
    direction <- rep("upstream", length(distance))
    direction[which(distance < 0)] <- "downstream"
    tbl.fp$distance.from.tss <- distance
    tbl.fp$id <- sprintf("%s.fp.%s.%06d.%s", targetGene, direction, abs(distance), tbl.fp$motifName)
    tbl.fp <- tbl.fp[, getRegulatoryTableColumnNames(obj)]

    tbl.fp

} # .callFootprintFilter
#------------------------------------------------------------------------------------------------------------------------
.callHumanDHSFilter <- function(obj, source, chromosome, chromStart, chromEnd, targetGene, targetGeneTSS)
{
    chromLocString <- sprintf("%s:%d-%d", chromosome, chromStart, chromEnd)
    dhsFilter <- HumanDHSFilter(genome="hg38",
                                encodeTableName="wgEncodeRegDnaseClustered",
                                pwmMatchPercentageThreshold=85L,
                                geneInfoDatabase.uri=genome.db.uri,
                                geneCenteredSpec=list(),
                                regionsSpec=chromLocString)

    x.dhs <- getCandidates(dhsFilter)
    tbl.dhs <- x.dhs$tbl
    tbl.dhs$length <- nchar(tbl.dhs$match)
    distance <- tbl.dhs$motifStart - targetGeneTSS
    direction <- rep("upstream", length(distance))
    direction[which(distance < 0)] <- "downstream"

    colnames(tbl.dhs)[grep("motifRelativeScore", colnames(tbl.dhs))] <- "score"
    colnames(tbl.dhs)[grep("tfs", colnames(tbl.dhs))] <- "tf"
    tbl.dhs$distance.from.tss <- distance
    tbl.dhs$id <- sprintf("%s.dhs.%s.%06d.%s", targetGene, direction, abs(distance), tbl.dhs$motifName)

    tbl.dhs <- tbl.dhs[, getRegulatoryTableColumnNames(obj)]

    tbl.dhs

} # .callHumanDHSFilter
#------------------------------------------------------------------------------------------------------------------------
setMethod('getRegulatoryRegions', 'TrenaPrep',

      function(obj, combine=FALSE){
         tbl.combined <- data.frame()
         result <- list()

         sources <- obj@regulatoryRegionSources
         encodeDHS.source.index <- grep("encodeHumanDHS", sources)

         if(length(encodeDHS.source.index)){
            sources <- sources[-encodeDHS.source.index]
            tbl.dhs <- .callHumanDHSFilter(obj, source, obj@chromosome, obj@chromStart, obj@chromEnd,
                                           obj@targetGene, obj@targetGeneTSS)
            result[["encodeHumanDHS"]] <- tbl.dhs
            if(combine)
               tbl.combined <- rbind(tbl.combined, tbl.dhs)
            } # if encode DSH source requested

         for(source in sources){  # don't use the object slot.  use locally scoped, possibly modified variable
            tbl.fp <- .callFootprintFilter(obj, source, obj@chromosome, obj@chromStart, obj@chromEnd,
                                           obj@targetGene, obj@targetGeneTSS)
            if(combine)
               tbl.combined <- rbind(tbl.combined, tbl.fp)
            result[[source]] <- tbl.fp
            } # for source
         if(combine)
            result[["all"]] <- tbl.combined
         result
         }) # getRegulatoryRegions

#------------------------------------------------------------------------------------------------------------------------
setMethod('expandRegulatoryRegionsTableByTF', 'TrenaPrep',

     function(obj, tbl.reg){
        tbl.trimmed <- subset(tbl.reg, nchar(tf) != 0)
        tfs.split <- strsplit(tbl.trimmed$tf, ";")
        #length(tfs.split) # [1] 36929
        counts <- unlist(lapply(tfs.split, length))
        tfs.split.vec <- unlist(tfs.split)
        tbl.expanded <- expandRows(tbl.trimmed, counts, count.is.col=FALSE, drop=FALSE)
        checkEquals(length(tfs.split.vec), nrow(tbl.expanded))
        tbl.expanded$tf <- tfs.split.vec
        tbl.expanded
        }) # expandRegulatoryRegionsTableByTF

#------------------------------------------------------------------------------------------------------------------------
setMethod('createGeneModel', 'TrenaPrep',

      function(obj, solvers, tbl.regulatoryRegions, mtx){

         stopifnot(solvers=="randomForest")  # more solvers to come
         tbl.small <- subset(tbl.regulatoryRegions,
                             chrom==obj@chromosome & motifStart >= obj@chromStart & motifEnd <= obj@chromEnd)

         tfs <- sort(unique(unlist(strsplit(tbl.small$tf, ";"))))
         tfs <- intersect(tfs, rownames(mtx))
         printf("tf candidate count: %d", length(tfs))
         solver.wt <- RandomForestSolver(mtx, targetGene=obj@targetGene, candidateRegulators=tfs)
         model.wt  <-run(solver.wt)
         count <- nrow(model.wt$edges)
         tbl.model <- data.frame(tf=rownames(model.wt$edges),
                                 randomForest=model.wt$edges$IncNodePurity,
                                 pearson=model.wt$edges$gene.cor,
                                 spearman=rep(0, count),
                                 betaLasso=rep(0, count),
                                 pcaMax=rep(0, count),
                                 concordance=rep(0, count),
                                 stringsAsFactors=FALSE)
         tbl.model
      }) # createGeneModel

#------------------------------------------------------------------------------------------------------------------------
setMethod('buildMultiModelGraph', 'TrenaPrep',

  function (obj, models){

    g <- graphNEL(edgemode = "directed")
    model.names <- names(models)

    node.attribute.specs <- list(type="undefined",
                                 label="default node label",
                                 distance=0,
                                 pearson=0,
                                 randomForest=0,
                                 pcaMax=0,
                                 concordance=0,
                                 betaLasso=0,
                                 motif="",
                                 xPos=0,
                                 yPos=0)
    edge.attribute.spec <- list(edgeType="undefined")
    attribute.classes <- c("", model.names)  # "" (no prefix) is the currently displayed set of attibutes

      # create current version of these attributes, and then
      # per-model versions, which get mapped to current
      # in response to user's interactive choice on the cyjs user interface
      # the "current version" is, e.g., "distance".
      # per-model ("wt" and "mut" versions) become "wt.distance" and "mut.distance"
      # and are used by copying e.g. all wt.xxx attributes into the current (non-prefixed)
      # attribute, upon which the cyjs style is defined

    for(class.name in attribute.classes){
       class.name.prefix <- class.name  # with possible "." appended, permits standard and model-specific attributes
       if(nchar(class.name) > 0)
          class.name.prefix <- sprintf("%s.", class.name)
       noa.names.without.prefix <- names(node.attribute.specs)
       noa.names <- sprintf("%s%s", class.name.prefix, noa.names.without.prefix)
       noa.count <- length(node.attribute.specs)
       for(i in 1:noa.count){
          nodeDataDefaults(g, attr=noa.names[i]) <- node.attribute.specs[[noa.names.without.prefix[i]]]
          }
       } # for class

    edgeDataDefaults(g, attr = "edgeType") <- "undefined"

    tfs <- c()
    regulatoryRegions <- c()

    for(model in models){  # collect all the tf and regulatory region nodes
       tbl.model <- model$tbl.geneModel
       tfs <- unique(c(tfs, tbl.model$tf))
       tbl.reg <- model$tbl.regulatoryRegions
       regulatoryRegions <- unique(c(regulatoryRegions, tbl.reg$id))
       } # for model

    all.nodes <- unique(c(obj@targetGene, tfs, regulatoryRegions))
    g <- addNode(all.nodes, g)

    nodeData(g, obj@targetGene, "type") <- "targetGene"
    nodeData(g, tfs, "type")         <- "TF"
    nodeData(g, regulatoryRegions, "type")  <- "regulatoryRegion"
    nodeData(g, all.nodes, "label")  <- all.nodes

      # add edges, edge attribute, and the constant attributes for all of the regulatoryRegion nodes

    for(model in models){
       tfs <- model$tbl.regulatoryRegions$tf
       regRegions <- model$tbl.regulatoryRegions$id
       suppressWarnings(g <- addEdge(tfs, regRegions, g))
       edgeData(g,  tfs, regRegions, "edgeType") <- "bindsTo"
       suppressWarnings(g <- addEdge(regRegions, obj@targetGene, g))
       edgeData(g, regRegions, obj@targetGene, "edgeType") <- "regulatorySiteFor"
       nodeData(g, tbl.reg$id, "label") <- tbl.reg$motifName
       nodeData(g, tbl.reg$id, "distance") <- tbl.reg$distance.from.tss
       nodeData(g, tbl.reg$id, "motif") <- tbl.reg$motifName
       } # for model

      # now copy in the first model's tf node data

    tbl.model <- models[[1]]$tbl.geneModel
    nodeData(g, tbl.model$tf, attr="randomForest") <- tbl.model$randomForest
    nodeData(g, tbl.model$tf, attr="pearson") <- tbl.model$pearson

     # now copy in each of the model's tf node data in turn
    model.names <- names(models)
    for(model.name in model.names){
       tbl.model <- models[[model.name]]$tbl.geneModel
       noa.name <- sprintf("%s.%s", model.name, "randomForest")
       nodeData(g,  tbl.model$tf, attr=noa.name) <- tbl.model$randomForest
       noa.name <- sprintf("%s.%s", model.name, "pearson")
       nodeData(g,  tbl.model$tf, attr=noa.name) <- tbl.model$pearson
      } # for model.name

    g

    }) # buildMultiModelGraph

#------------------------------------------------------------------------------------------------------------------------
setMethod('addGeneModelLayout', 'TrenaPrep',

  function (obj, g, xPos.span=1500){
    printf("--- addGeneModelLayout")
    all.distances <- sort(unique(unlist(nodeData(g, attr='distance'), use.names=FALSE)))
    print(all.distances)

    fp.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "regulatoryRegion")]
    tf.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "TF")]
    targetGene.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "targetGene")]

     # add in a zero in case all of the footprints are up or downstream of the 0 coordinate, the TSS
    span.endpoints <- range(c(0, as.numeric(nodeData(g, fp.nodes, attr="distance"))))
    span <- max(span.endpoints) - min(span.endpoints)
    footprintLayoutFactor <- 1
    printf("initial:  span: %d  footprintLayoutFactor: %f", span, footprintLayoutFactor)

    footprintLayoutFactor <- xPos.span/span

    #if(span < 600)  #
    #   footprintLayoutFactor <- 600/span
    #if(span > 1000)
    #   footprintLayoutFactor <- span/1000

    printf("corrected:  span: %d  footprintLayoutFactor: %f", span, footprintLayoutFactor)

    xPos <- as.numeric(nodeData(g, fp.nodes, attr="distance")) * footprintLayoutFactor
    yPos <- 0
    nodeData(g, fp.nodes, "xPos") <- xPos
    nodeData(g, fp.nodes, "yPos") <- yPos

    adjusted.span.endpoints <- range(c(0, as.numeric(nodeData(g, fp.nodes, attr="xPos"))))
    printf("raw span of footprints: %d   footprintLayoutFactor: %f  new span: %8.0f",
           span, footprintLayoutFactor, abs(max(adjusted.span.endpoints) - min(adjusted.span.endpoints)))

    tfs <- names(which(nodeData(g, attr="type") == "TF"))

    for(tf in tfs){
       footprint.neighbors <- edges(g)[[tf]]
       footprint.positions <- as.integer(nodeData(g, footprint.neighbors, attr="xPos"))
       new.xPos <- mean(footprint.positions)
       #printf("%8s: %5d", tf, new.xPos)
       nodeData(g, tf, "xPos") <- new.xPos
       nodeData(g, tf, "yPos") <- sample(300:1200, 1)
       } # for tf

    nodeData(g, targetGene.nodes, "xPos") <- 0
    nodeData(g, targetGene.nodes, "yPos") <- -200

    g

    }) # addGeneModelLayout

#------------------------------------------------------------------------------------------------------------------------
