library(TReNA)
library(TrenaHelpers)


# slide 1
tv <- TrenaViz()
showGenomicRegion(tv, "18:26,845,093-26,876,905")

# slide 2
tbl.snp <- data.frame(chromosome=rep("chr18", 5),
                      loc=c(26864410, 26865469, 26855623, 26855854, 26850565),
                      rsid=c("rs3763040", "rs3875089", "rs335929", "rs3763043", "rs9951307"),
                      genome=rep("hg38", 5),
                      stringsAsFactors=FALSE)

addBedTrackFromDataFrame(tv, "snp", tbl.snp[, c("chromosome", "loc", "loc", "rsid")], color="red")


# slide 3: load 5 tracks of regulatory regions
chrom <- "chr18"
start <- 26848475
end   <- 26869165
aqp4.tss <- 26865884
db.names <- c("brain_hint_16", "brain_hint_20")
sources <- list("postgres://whovian/brain_hint_16",
                "postgres://whovian/brain_hint_20",
                "postgres://whovian/brain_wellington_16",
                "postgres://whovian/brain_wellington_20")

trackNames <- unlist(lapply(sources, function(uri) sub("postgres://whovian/", "", uri, fixed=TRUE)))
utils <- TrenaPrep("AQP4", aqp4.tss, chrom, start, end, regulatoryRegionSources=sources)
x <- getRegulatoryRegions(utils)
for(i in 1:length(x)){
   tbl.fp <- x[[i]]
   addBedTrackFromDataFrame(tv, trackNames[i],
                            tbl.fp[, c("chrom", "motifStart", "motifEnd", "motifName", "score")],
                            color="blue")
   }

# slide 4: add the dhs track (not yet part of TrenaUtils)
genome.db.uri <- "postgres://bddsrds.globusgenomics.org/hg38"   # has gtf and motifsgenes tables
dhsFilter <- HumanDHSFilter(genome="hg38",
                            encodeTableName="wgEncodeRegDnaseClustered",
                            pwmMatchPercentageThreshold=85L,
                            geneInfoDatabase.uri=genome.db.uri,
                            geneCenteredSpec=list(),
                            regionsSpec=sprintf("%s:%d-%d", chrom, start, end))
x.dhs <- getCandidates(dhsFilter)
addBedTrackFromDataFrame(tv, "DHS",
                         x.dhs$tbl[, c("chrom", "motifStart", "motifEnd", "motifName", "motifScore")],
                         color="darkgreen")


# slide 5: load an expression matrix

load("~/github/projects/examples/microservices/trenaGeneModel/datasets/coryAD/rosmap_counts_matrix_normalized_geneSymbols_25031x638.RData")
mtx <- asinh(mtx)
mtx.var <- apply(mtx, 1, var)
deleters <- which(mtx.var < 0.01)
if(length(deleters) > 0)   # 15838 x 638
   mtx <- mtx[-deleters,]

load("mayo.rnaSeq.cer.and.tcx.matrices.RData")

dim(mtx); fivenum(mtx)

# slide 6: build a model assuming wild type sequence, using igv display for coordinates
region <- TReNA:::.parseChromLocString(getGenomicRegion(tv))
targetGene <- "AQP4"
aqp4.tss <- 26865884
fp.source <- "postgres://whovian/brain_hint_20"
sources <- list(fp.source)

prep <- TrenaPrep(targetGene, aqp4.tss, "chr18", region$start, region$end, regulatoryRegionSources=sources)
x <- getRegulatoryRegions(prep)
tbl.regulatoryRegions <- expandRegulatoryRegionsTableByTF(prep, x[[fp.source]])

dim(tbl.regulatoryRegions)  # 11014 10
tbl.regulatoryRegions[sample(1:nrow(tbl.regulatoryRegions), 5),]   # a few rows at random

tbl.geneModel <- createGeneModel(prep, "randomForest", tbl.regulatoryRegions, mtx)
fivenum(tbl.geneModel$randomForest)
subset(tbl.geneModel, randomForest > 1)[, c(1:3)]

# slide 7: create graph, display
tbl.strong <- subset(tbl.geneModel, randomForest > 1)
tbl.regulatoryRegions.strong <- subset(tbl.regulatoryRegions, tf %in% tbl.strong$tf)
dim(tbl.regulatoryRegions.strong)  # 423 10


model.for.viewing <- list(wt=list(tbl.regulatoryRegions=tbl.regulatoryRegions.strong,
                          tbl.geneModel=tbl.strong))

g <- buildMultiModelGraph(prep, model.for.viewing)   # 254 nodes, 602 edges
g.lo <- addGeneModelLayout(prep, g, xPos.span=1500)
httpAddGraph(tv, g.lo, "wt")
loadStyle(tv, system.file(package="TrenaHelpers", "extdata", "style.js"))
fit(tv)

# slide 8: explore rs3875089
tbl.assay <- assessSnp(prep, "rs3875089", 10, pwmMatchMinimumAsPercentage=80)
subset(tbl.assay, assessed != "in.both")[, c(1,2,3,4,10,11)]



# slide 9: how low are the TEAD-related motifs in the mutant?
tbl.assay <- assessSnp(prep, "rs3875089", 10, pwmMatchMinimumAsPercentage=60)
subset(tbl.assay[grep("TEAD", tbl.assay$tf), c(1,2,3,4,10,11)], assessed == "in.both")

broken.motifs <- c("MA0809.1", "MA0808.1", "MA0090.2")

# slide 10: filter these motifs out of tbl.regulatoryRegions, create a new gene model
tbl.regulatoryRegions.rs3875089 <- subset(tbl.regulatoryRegions, !motifName %in% broken.motifs)
nrow(tbl.regulatoryRegions) - nrow(tbl.regulatoryRegions.rs3875089)   # just 12 lost
tbl.geneModel.rs3875089 <- createGeneModel(prep, "randomForest", tbl.regulatoryRegions.rs3875089, mtx)

fivenum(tbl.geneModel.rs3875089$randomForest)
subset(tbl.geneModel.rs3875089, randomForest > 1)[, c(1:3)]

tbl.strong.rs3875089 <- subset(tbl.geneModel.rs3875089, randomForest > 1)
tbl.regulatoryRegions.rs.strong <- subset(tbl.regulatoryRegions.rs3875089, tf %in% tbl.strong.rs3875089$tf)
dim(tbl.regulatoryRegions.rs.strong)  # 423 10


dual.model <- list(wt=list(tbl.regulatoryRegions=tbl.regulatoryRegions.strong,tbl.geneModel=tbl.strong),
                   rs3875089=list(tbl.regulatoryRegions=tbl.regulatoryRegions.rs.strong,
                                  tbl.geneModel=tbl.strong.rs3875089))


g2 <- buildMultiModelGraph(prep, dual.model)   # 254 nodes, 602 edges
g2.lo <- addGeneModelLayout(prep, g2, xPos.span=1500)
addGraph(tv, g2.lo, c("wt", "rs3875089"))
loadStyle(tv, system.file(package="TrenaHelpers", "extdata", "style.js"))


# slide 11: create two contrasting models of mef2c using different expression matrices:
#   mtx.tcx, mtx.cer

mef2c.chrom <- "chr5"
mef2c.tss <- 88826234

showGenomicRegion(tv, sprintf("%s:%d-%d", mef2c.chrom, mef2c.tss-1000, mef2c.tss+1000))
chrom <- mef2c.chrom
start <- mef2c.tss - 1000
end   <- mef2c.tss + 2000
tss <- mef2c.tss
db.names <- c("brain_hint_16", "brain_hint_20")
sources <- list("postgres://whovian/brain_hint_16",
                "postgres://whovian/brain_hint_20",
                "postgres://whovian/brain_wellington_16",
                "postgres://whovian/brain_wellington_20")

trackNames <- unlist(lapply(sources, function(uri) sub("postgres://whovian/", "", uri, fixed=TRUE)))
utils <- TrenaPrep("MEF2C", tss, chrom, start, end, regulatoryRegionSources=sources)
x <- getRegulatoryRegions(utils)
for(i in 1:length(x)){
   tbl.fp <- x[[i]]
   addBedTrackFromDataFrame(tv, trackNames[i],
                            tbl.fp[, c("chrom", "motifStart", "motifEnd", "motifName", "score")],
                            color="blue")
   }

# slide 4: add the dhs track (not yet part of TrenaUtils)
genome.db.uri <- "postgres://bddsrds.globusgenomics.org/hg38"   # has gtf and motifsgenes tables
dhsFilter <- HumanDHSFilter(genome="hg38",
                            encodeTableName="wgEncodeRegDnaseClustered",
                            pwmMatchPercentageThreshold=85L,
                            geneInfoDatabase.uri=genome.db.uri,
                            geneCenteredSpec=list(),
                            regionsSpec=sprintf("%s:%d-%d", chrom, start, end))
x.dhs <- getCandidates(dhsFilter)
addBedTrackFromDataFrame(tv, "DHS",
                         x.dhs$tbl[, c("chrom", "motifStart", "motifEnd", "motifName", "motifScore")],
                         color="darkgreen")

print(load("mayo.rnaSeq.cer.and.tcx.matrices.RData"))

region <- TReNA:::.parseChromLocString(getGenomicRegion(tv))
targetGene <- "MEF2C"
fp.source <- "postgres://whovian/brain_hint_20"
sources <- list(fp.source)

prep <- TrenaPrep(targetGene, tss, "chr5", region$start, region$end, regulatoryRegionSources=sources)
x <- getRegulatoryRegions(prep)
tbl.regulatoryRegions <- expandRegulatoryRegionsTableByTF(prep, x[[fp.source]])

dim(tbl.regulatoryRegions)  # 11014 10
tbl.regulatoryRegions[sample(1:nrow(tbl.regulatoryRegions), 5),]   # a few rows at random

mtx <- mtx.tcx
mtx <- mtx.cer
tbl.geneModel.cer <- createGeneModel(prep, "randomForest", tbl.regulatoryRegions, mtx.cer)
fivenum(tbl.geneModel.cer$randomForest)
fivenum(tbl.geneModel.tcx$randomForest)
subset(tbl.geneModel.cer, randomForest > 5)[, c(1:3)]
subset(tbl.geneModel.tcx, randomForest > 5)[, c(1:3)]

# slide 7: create graph, display
tbl.strong.cer <- subset(tbl.geneModel.cer, randomForest > 5)
tbl.strong.tcx <- subset(tbl.geneModel.tcx, randomForest > 5)

tbl.regulatoryRegions.cer.strong <- subset(tbl.regulatoryRegions, tf %in% tbl.strong.cer$tf)
tbl.regulatoryRegions.tcx.strong <- subset(tbl.regulatoryRegions, tf %in% tbl.strong.tcx$tf)
dim(tbl.regulatoryRegions.cer.strong)  # 423 10
dim(tbl.regulatoryRegions.tcx.strong)  # 423 10


model.for.viewing <- list(cer=list(tbl.regulatoryRegions=tbl.regulatoryRegions.cer.strong,
                                   tbl.geneModel=tbl.strong.cer),
                          tcx=list(tbl.regulatoryRegions=tbl.regulatoryRegions.tcx.strong,
                                   tbl.geneModel=tbl.strong.tcx))

g <- buildMultiModelGraph(prep, model.for.viewing)   # 254 nodes, 602 edges
g.lo <- addGeneModelLayout(prep, g, xPos.span=1500)
addGraph(tv, g.lo, c("cer", "tcx"))
loadStyle(tv, system.file(package="TrenaHelpers", "extdata", "style.js"))
fit(tv)
