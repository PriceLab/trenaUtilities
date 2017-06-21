library(TReNA)
library(TrenaHelpers)
library(RPostgreSQL)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
tbl.snp <- data.frame(target.gene=rep("AQP4", 5),
                      chromosome=rep("chr18", 5),
                      loc=c(26864410, 26865469, 26855623, 26855854, 26850565),
                      snp=c("rs3763040", "rs3875089", "rs335929", "rs3763043", "rs9951307"),
                      genome=rep("hg38", 5),
                      stringsAsFactors=FALSE)
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

   coi <- c("chrom", "motifStart", "motifEnd", "motifName", "strand", "score", "length", "id")
   chromLocString <- sprintf("%s:%d-%d", chrom, start, end)

   for(dbName in dboi[1]){
     footprint.db.uri <- sprintf("postgres://whovian/%s", dbName)
     fpFilter <- FootprintFilter(genome.db.uri, footprint.db.uri, list(), chromLocString)
     x.fp <- getCandidates(fpFilter)
     tbl.fp <- x.fp$tbl
     colnames(tbl.fp) <- c("chrom", "motifStart", "motifEnd", "motifName", "length", "strand", "score1", "score", "score3")
     distance <- tbl.fp$motifStart - tss
     direction <- rep("upstream", length(distance))
     direction[which(distance < 0)] <- "downstream"
     tbl.fp$id <- sprintf("%s.fp.%s.%06d.%s", target.gene, direction, abs(distance), tbl.fp$motifName)
     tbl.fp <- tbl.fp[, coi]
     tbl.combined <- rbind(tbl.combined, tbl.fp)
     addBedTrackFromDataFrame(tv, dbName, tbl.fp[, c("chrom", "motifStart", "motifEnd", "motifName", "score")], color="blue")
     #closeDatabaseConnections(fpFilter)
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

    #browser()
    colnames(tbl.dhs)[grep("motifRelativeScore", colnames(tbl.dhs))] <- "score"
    tbl.dhs$id <- sprintf("%s.fp.%s.%06d.%s", target.gene, direction, abs(distance), tbl.dhs$motifName)

    tbl.dhs <- tbl.dhs[, coi]
    tbl.combined <- rbind(tbl.combined, tbl.dhs)

   printf("%d dhs regions in %d bases from %s", nrow(tbl.dhs), 1  + end - start, "DHS")
   addBedTrackFromDataFrame(tv, "DHS", tbl.dhs[, c("chrom", "motifStart", "motifEnd", "motifName", "score")], color="blue")

   tbl.combined.bed <- unique(tbl.combined[,c(1,2,3,4)])
   tbl.combined.bed$score <- 1
   addBedTrackFromDataFrame(tv, "all", tbl.combined.bed, color="magenta")

   invisible(tbl.combined)

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
   browser();
   tbl.small <- subset(tbl.regulatoryRegions, chrom==chrom.name & start >= chrom.start & end <= chrom.end)
   tfs <- unique(unlist(lapply(tbl.small$name, function(m) subset(tbl.mg, motif==m)$tf.gene)))
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
   tbl.model <- subset(tbl.model, randomForest >= 2)

   aqp4.tss <- 26865884
   required.regulatoryRegionsColumnNames <- c("chrom", "motifStart", "motifEnd", "motifName",
                                              "strand", "score", "length", "id")

   tfs <- lapply(tbl.regulatoryRegions$motifName, function(m) subset(tbl.mg, motif==m)$tf.gene)
   tbl.motifs <- with(tbl.regulatoryRegions, data.frame(motifName=name, chrom=chrom, motifStart=start, motifEnd=end,
                                                        strand=stran, score=score, tf=tfs))

   print(head(model.wt$edges, n=20))
   browser()
   xyz <- 99;

} # create.wildtype.model
#------------------------------------------------------------------------------------------------------------------------
run.all <- function()
{
   x <<- findAndDisplayRegulatoryRegions()
   showGenomicRegion(x$tv, "18:26,862,000-26,866,000")
   chromLoc <<- TReNA:::.parseChromLocString(getGenomicRegion(x$tv))
   #model <<- create.wildtype.model(x$tbl, target.gene="AQP4", chromLoc$chrom, chromLoc$start, chromLoc$end)
   browser()
   xyz <- 99;

} # run.all
#------------------------------------------------------------------------------------------------------------------------
