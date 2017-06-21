library(TReNA)
library(TrenaHelpers)
library(RPostgreSQL)
library(RUnit)
library(splitstackshape)
#------------------------------------------------------------------------------------------------------------------------
tbl.snp <- data.frame(target.gene=rep("AQP4", 5),
                      chromosome=rep("chr18", 5),
                      loc=c(26864410, 26865469, 26855623, 26855854, 26850565),
                      snp=c("rs3763040", "rs3875089", "rs335929", "rs3763043", "rs9951307"),
                      genome=rep("hg38", 5),
                      stringsAsFactors=FALSE)

tbl.reg.colnames <- c("chrom", "motifStart", "motifEnd", "motifName", "strand", "score", "length", "distance.from.tss", "id", "tf")
tbl.model.colnames <- c("tf", "randomForest", "pearson", "spearman", "betaLasso", "pcaMax", "concordance")

#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx")){
  load("~/github/projects/examples/microservices/trenaGeneModel/datasets/coryAD/rosmap_counts_matrix_normalized_geneSymbols_25031x638.RData")
  #load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
  mtx <- asinh(mtx)
  mtx.var <- apply(mtx, 1, var)
  deleters <- which(mtx.var < 0.01)
  if(length(deleters) > 0)   # 15838 x 638
     mtx <- mtx[-deleters,]
  mtx <<- mtx
  }
#------------------------------------------------------------------------------------------------------------------------
init <- function()
{
   tv <- TrenaViz()
   tv

} # init
#------------------------------------------------------------------------------------------------------------------------
getAndDisplayFootprintsFromDatabases <- function(tv, target.gene, tss, chrom, start, end)
{
      # make sure all 4 brain footprint databases are available
   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", host="whovian")
   dboi <- c("brain_wellington_16", "brain_hint_16", "brain_hint_20", "brain_wellington_20")
   checkTrue(all(dboi %in% dbGetQuery(db, "select datname from pg_database")$datname))
   dbDisconnect(db)

     # genome information
   genome.db.uri <- "postgres://bddsrds.globusgenomics.org/hg38"   # has gtf and motifsgenes tables
   tbl.combined <- data.frame()

   chromLocString <- sprintf("%s:%d-%d", chrom, start, end)

   for(dbName in dboi[1]){
     footprint.db.uri <- sprintf("postgres://whovian/%s", dbName)
     fpFilter <- FootprintFilter(genome.db.uri, footprint.db.uri, list(), chromLocString)
     x.fp <- getCandidates(fpFilter)
     tbl.fp <- x.fp$tbl
     colnames(tbl.fp) <- c("chrom", "motifStart", "motifEnd", "motifName", "length", "strand", "score1", "score", "score3", "tf")
     distance <- tbl.fp$motifStart - tss
     direction <- rep("upstream", length(distance))
     direction[which(distance < 0)] <- "downstream"
     tbl.fp$distance.from.tss <- distance
     tbl.fp$id <- sprintf("%s.fp.%s.%06d.%s", target.gene, direction, abs(distance), tbl.fp$motifName)
     tbl.fp <- tbl.fp[, tbl.reg.colnames]
     tbl.combined <- rbind(tbl.combined, tbl.fp)
     addBedTrackFromDataFrame(tv, dbName, tbl.fp[, c("chrom", "motifStart", "motifEnd", "motifName", "score")], color="blue")
     }

    dhsFilter <- HumanDHSFilter(genome="hg38",
                                encodeTableName="wgEncodeRegDnaseClustered",
                                pwmMatchPercentageThreshold=85L,
                                geneInfoDatabase.uri=genome.db.uri,
                                geneCenteredSpec=list(),
                                regionsSpec=chromLocString)
    x.dhs <- getCandidates(dhsFilter)
    tbl.dhs <- x.dhs$tbl
    tbl.dhs$length <- nchar(tbl.dhs$match)
    distance <- tbl.dhs$motifStart - tss
    direction <- rep("upstream", length(distance))
    direction[which(distance < 0)] <- "downstream"

    colnames(tbl.dhs)[grep("motifRelativeScore", colnames(tbl.dhs))] <- "score"
    colnames(tbl.dhs)[grep("tfs", colnames(tbl.dhs))] <- "tf"
    tbl.dhs$distance.from.tss <- distance
    tbl.dhs$id <- sprintf("%s.fp.%s.%06d.%s", target.gene, direction, abs(distance), tbl.dhs$motifName)

    tbl.dhs <- tbl.dhs[, tbl.reg.colnames]
    tbl.combined <- rbind(tbl.combined, tbl.dhs)

   printf("%d dhs regions in %d bases from %s", nrow(tbl.dhs), 1  + end - start, "DHS")
   addBedTrackFromDataFrame(tv, "DHS", tbl.dhs[, c("chrom", "motifStart", "motifEnd", "motifName", "score")], color="blue")

   tbl.combined.bed <- unique(tbl.combined[,c(1,2,3,4)])
   tbl.combined.bed$score <- 1
   addBedTrackFromDataFrame(tv, "all", tbl.combined.bed, color="magenta")

   browser()
   printf("----- expanding tfs, starting row count: %d", nrow(tbl.combined))
   tbl.trimmed <- subset(tbl.combined, nchar(tf) != 0)
   tfs.split <- strsplit(tbl.trimmed$tf, ";")
   length(tfs.split) # [1] 36929
   counts <- unlist(lapply(tfs.split, length))
   tfs.split.vec <- unlist(tfs.split)
   tbl.expanded <- expandRows(tbl.trimmed, counts, count.is.col=FALSE, drop=FALSE)
   checkEquals(length(tfs.split.vec), nrow(tbl.expanded))
   tbl.expanded$tf <- tfs.split.vec
   printf("----- expanding tfs, done, concluding row count: %d", nrow(tbl.expanded))

   invisible(tbl.expanded)

} # getAndDisplayFootprintsFromDatabases
#------------------------------------------------------------------------------------------------------------------------
getDHSregions <- function(chrom, start, end)
{
   hdf <- HumanDHSFilter(genome="hg38",
                         encodeTableName="wgEncodeRegDnaseClustered",
                         pwmMatchPercentageThreshold=85L,
                         geneInfoDatabase.uri="postgres://whovian/gtf",
                         geneCenteredSpec=list(),
                         regionsSpec=sprintf("%s:%d-%d", chrom, start, end))

   tbl.dhs <- getRegulatoryRegions(hdf, "wgEncodeRegDnaseClustered", chrom, start, end)

   colnames(tbl.dhs) <- c("chrom", "start", "end", "name", "score")  # count -> name, good for display in igv
   tbl.dhs

} # displayDHSregions
#------------------------------------------------------------------------------------------------------------------------
# brain_hint_16 has one fp that intersects with snp2(26865469): chr18:26865467-26865487
# brain_hint_20: chr18:26865467-26865487
findAndDisplayRegulatoryRegions <- function()
{
   tv <- init()

   target.gene <- "AQP4"
   tss <- 26865884
   chrom <- unique(tbl.snp$chromosome)
   region.start <- min(tbl.snp$loc[2] - 100)
   region.end   <- max(tbl.snp$loc[2] + 100)
   snp.start <- min(tbl.snp$loc) - 1000;
   snp.end   <- max(tbl.snp$loc) + 10000;

   start <- snp.start
   end   <- snp.end

   tbl.combined <- getAndDisplayFootprintsFromDatabases(tv, target.gene, tss, chrom, start, end)
   browser()
   region <- sprintf("%s:%d-%d", tbl.snp$chromosome[1], start-1000, end+1000)
   showGenomicRegion(tv, region)
   addBedTrackFromDataFrame(tv, "snp",
                            tbl.snp[, c("chromosome", "loc", "loc", "snp")],
                            color="purple")



   return(list(tv=tv, tbl=tbl.combined))

} # findAndDisplayRegulatoryRegions
#------------------------------------------------------------------------------------------------------------------------
create.wildtype.model <- function(tbl.regulatoryRegions, target.gene, chrom.name, chrom.start, chrom.end)
{
   tbl.mg <- read.table(system.file(package="TReNA", "extdata", "motifGenes.tsv"), sep="\t", as.is=TRUE, header=TRUE)
   tbl.small <- subset(tbl.regulatoryRegions, chrom==chrom.name & motifStart >= chrom.start & motifEnd <= chrom.end)
   tfs <- sort(unique(unlist(strsplit(tbl.small$tf, ";"))))
   tfs <- intersect(tfs, rownames(mtx))
   printf("tf candidate count: %d", length(tfs))
   solver.wt <- RandomForestSolver(mtx, targetGene=target.gene, candidateRegulators=tfs)
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
   #tbl.model <- subset(tbl.model, randomForest >= 2)

   #aqp4.tss <- 26865884
   #required.regulatoryRegionsColumnNames <- c("chrom", "motifStart", "motifEnd", "motifName",
   #                                           "strand", "score", "length", "id")

   #tfs <- lapply(tbl.regulatoryRegions$motifName, function(m) subset(tbl.mg, motif==m)$tf.gene)
   #tbl.motifs <- with(tbl.regulatoryRegions, data.frame(motifName=name, chrom=chrom, motifStart=start, motifEnd=end,
   #                                                     strand=stran, score=score, tf=tfs))
   #print(head(model.wt$edges, n=20))
   #browser()
   #xyz <- 99;

   tbl.model

} # create.wildtype.model
#------------------------------------------------------------------------------------------------------------------------
test.findRegulatoryRegions <- function()
{
   printf("--- test.findRegulatoryRegions")

   x <- findAndDisplayRegulatoryRegions()
   tbl.reg <- x$tbl
   checkEquals(colnames(tbl.reg), tbl.reg.colnames)

      # should be no semi-colon-separated tf names
   checkEquals(length(grep(";", tbl.reg$tf)), 0)
      # no empty tf names either
   length(which(nchar(tbl.reg$tf)==0))

      # with current profligate assignment of tfs to motifs, the
      # expanded tbl.reg - with one tf per row - is big
   checkTrue(nrow(tbl.reg) > 10000)

   save(tbl.reg, file="tbl.reg.4k.sample.RData")

   #showGenomicRegion(x$tv, "18:26,862,000-26,866,000")
   #chromLoc <<- TReNA:::.parseChromLocString(getGenomicRegion(x$tv))
   #model <<- create.wildtype.model(x$tbl, target.gene="AQP4", chromLoc$chrom, chromLoc$start, chromLoc$end)
   #browser()
   #xyz <- 99;

} # test.findRegulatoryRegions
#------------------------------------------------------------------------------------------------------------------------
test.create.wildtype.model <- function()
{
   printf("--- test.create.wildtype.model")

   load("tbl.reg.4k.sample.RData")
   tbl.model <- create.wildtype.model(tbl.reg, target.gene="AQP4", "chr18", 26865188, 26866330)
   tbl.model <- subset(tbl.model, randomForest > 10)  # just a few rows
   checkTrue(nrow(tbl.model) >= 4)
   checkEquals(colnames(tbl.model), tbl.model.colnames)
   save(tbl.model, file="tbl.model.4k.sample.RData")

} # test.create.wildtype.model
#------------------------------------------------------------------------------------------------------------------------
test.createMultiModelGraph <- function()
{
   printf("--- test.create.wildtype.model")

   load("tbl.reg.4k.sample.RData")    # 42545 x 10
   load("tbl.model.4k.sample.RData")  # 4 x 7
   tbl.reg.reduced <- subset(tbl.reg, tf %in% tbl.model$tf)

   tViz <- TrenaViz()

   models <- list(wt=list(tbl.model=tbl.model, tbl.regulatoryRegions=tbl.reg.reduced))

   g <- buildMultiModelGraph(tViz, "AQP4", models)
   checkTrue(all(c("AQP4", tbl.reg.reduced$id, tbl.model$tf) %in% nodes(g)))
   save(g, file="g.400.RData")

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

} # test.createMultiModelGraph
#------------------------------------------------------------------------------------------------------------------------
test.layout <- function()
{
   printf("--- test.layout")
   load("g.4k.sample.RData")
   tViz <- TrenaViz()
   g.lo <- addGeneModelLayout(tViz, g, xPos.span=1500) # scale all genome distances so that they fit into 1500 pixels
   distances <- as.numeric(nodeData(g.lo, attr="distance"))   # -13608 .. 9430
   xPos <- as.numeric(nodeData(g.lo, attr="xPos"))
   actual.xPos.span <- abs(max(xPos) - min(xPos))
   checkEqualsNumeric(actual.xPos.span, 1500)

   addGraph(tViz, g.lo, "AQP4")
   loadStyle(tViz, system.file(package="TrenaHelpers", "extdata", "style.js"))
   Sys.sleep(3); fit(tViz)

} # test.layout
#------------------------------------------------------------------------------------------------------------------------
