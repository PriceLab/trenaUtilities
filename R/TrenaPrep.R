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
   browser()

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
    tbl.dhs$id <- sprintf("%s.fp.%s.%06d.%s", targetGene, direction, abs(distance), tbl.dhs$motifName)

    tbl.dhs <- tbl.dhs[, getRegulatoryTableColumnNames(obj)]

    tbl.dhs

} # .callHumanDHSFilter
#------------------------------------------------------------------------------------------------------------------------
setMethod('getRegulatoryRegions', 'TrenaPrep',

      function(obj, combine=FALSE){
         tbl.combined <- data.frame()
         result <- list()
         encodeDHS.source.index <- grep("encodeHumanDHS", sources)

         if(length(encodeDHS.source.index)){
            sources <- sources[-encodeDHS.source.index]
            tbl.dhs <- .callHumanDHSFilter(obj, source, obj@chromosome, obj@chromStart, obj@chromEnd,
                                           obj@targetGene, obj@targetGeneTSS)
            if(combine)
               tbl.combined <- rbind(tbl.combined, tbl.dhs)
            } # if encode DSH source requested

         for(source in obj@regulatoryRegionSources){
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
