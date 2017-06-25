library(TReNA)
library(TrenaHelpers)
library(RUnit)

# use 15kb region centered on the 5 snps: 18:
aqp4.tss <- 26865884
chrom <- "chr18"
start <- 26849640
end   <- 26869007

tv <- TrenaViz()
showGenomicRegion(tv, sprintf("%s:%d-%d", chrom, start, end))

# define the snps
tbl.snp <- data.frame(chromosome=rep("chr18", 5),
                      loc=c(26864410, 26865469, 26855623, 26855854, 26850565),
                      rsid=c("rs3763040", "rs3875089", "rs335929", "rs3763043", "rs9951307"),
                      genome=rep("hg38", 5),
                      stringsAsFactors=FALSE)

addBedTrackFromDataFrame(tv, "snp", tbl.snp[, c("chromosome", "loc", "loc", "rsid")], color="red")

fp.source <- "postgres://whovian/brain_hint_20"
sources <- list(fp.source)

tUtils <- TrenaPrep("AQP4", aqp4.tss, chrom, start, end, regulatoryRegionSources=sources)

# assess the snps: 2/5 snps destory TEAD1 binding
snp.assays <- list()
for(snp in tbl.snp$rsid){
   tbl.assay <- assessSnp(tUtils, snp, 10, pwmMatchMinimumAsPercentage=80)
   snp.assays[[snp]] <- tbl.assay
   }

lapply(snp.assays, function(x) length(grep("TEAD", x$tf)))
# $rs3763040:  0
# $rs3875089:  3
# $rs335929:   0
# $rs3763043:  3
# $rs9951307:  0

# ---- create 3 gene models for aqp4: rosman, tcx, cer

x <- getRegulatoryRegions(tUtils)

addBedTrackFromDataFrame(tv, names(x)[1],
                            x[[1]][, c("chrom", "motifStart", "motifEnd", "motifName", "score")],
                            color="blue")

tbl.regulatoryRegions <- expandRegulatoryRegionsTableByTF(tUtils, x[[1]])

load("~/github/projects/examples/microservices/trenaGeneModel/datasets/coryAD/rosmap_counts_matrix_normalized_geneSymbols_25031x638.RData")
mtx <- asinh(mtx)
mtx.var <- apply(mtx, 1, var)
deleters <- which(mtx.var < 0.01)
if(length(deleters) > 0)   # 15838 x 638
   mtx <- mtx[-deleters,]
mtx.rosmap <- mtx

print(load("mayo.rnaSeq.cer.and.tcx.matrices.RData"))  # mtx.tcx, mtx.cer

models <- list()
i <- 1
for(mtx in list(mtx.rosmap, mtx.cer, mtx.tcx)){
  models[[i]] <- createGeneModel(tUtils, "randomForest", tbl.regulatoryRegions, mtx)
  i <- i + 1
  }


tfs.strong <- c()
rf.threshold <- 1

strongGeneModels <- list()
i <- 1

for(model in models){
   tbl.geneModel.strong <- subset(model, randomForest >= rf.threshold)
   new.tfs <- tbl.geneModel.strong$tf
   tfs.strong <- c(tfs.strong, new.tfs)
   strongGeneModels[[i]] <- tbl.geneModel.strong
   i <- i + 1
   }

names(strongGeneModels) <- c("rosmap", "cer", "tcx")

tfs.strong <- unique(tfs.strong)
printf("tfs.strong: %d", length(tfs.strong))

tbl.regulatoryRegions.strong <- subset(tbl.regulatoryRegions, tf %in% tfs.strong)
dim(tbl.regulatoryRegions.strong)  # 315 10

model.for.viewing <- list(rosmap=list(tbl.regulatoryRegions=tbl.regulatoryRegions.strong, tbl.geneModel=strongGeneModels[["rosmap"]]),
                             cer=list(tbl.regulatoryRegions=tbl.regulatoryRegions.strong, tbl.geneModel=strongGeneModels[["cer"]]),
                             tcx=list(tbl.regulatoryRegions=tbl.regulatoryRegions.strong, tbl.geneModel=strongGeneModels[["tcx"]]))

g <- buildMultiModelGraph(tUtils, model.for.viewing)   # 254 nodes, 602 edges
g.lo <- addGeneModelLayout(tUtils, g, xPos.span=1500)
httpAddGraph(tv, g.lo, names(strongGeneModels))
#loadStyle(tv, system.file(package="TrenaHelpers", "extdata", "style.js"))
loadStyle(tv, "style.js")
fit(tv)

lapply(snp.assays, function(x) length(grep("TEAD", x$tf)))

rosmap.tfs <- strongGeneModels[["rosmap"]]$tf   # 27
cer.tfs    <- strongGeneModels[["cer"]]$tf      # 32
tcx.tfs    <- strongGeneModels[["tcx"]]$tf      # 17

all.tfs <- list(rosmap=rosmap.tfs, cer=cer.tfs, tcx=tcx.tfs)

for(mtx.name in names(all.tfs)){
   tfs <- all.tfs[[mtx.name]]
   printf("--- mtx: %s", mtx.name)
   for(snp in names(snp.assays)){
      printf("    snp: %s", snp)
      snp.tf.losses <- sapply(tfs, function(tf) length(grep(tf, subset(snp.assays[[snp]], assessed=="wt.only")$tf)))
      snp.tf.losses <- snp.tf.losses[snp.tf.losses != 0]
      snp.tf.gains <- sapply(tfs, function(tf) length(grep(tf, subset(snp.assays[[snp]], assessed=="mut.only")$tf)))
      snp.tf.gains <- snp.tf.gains[snp.tf.gains != 0]
      printf(" tfs with lost motifs:")
      print(snp.tf.losses)
      printf(" tfs with gained motifs:")
      print(snp.tf.gains)
      browser();
      xyz <- 99;
      }
   }


#                           tfs w/motifs lost                    tfs w/motifs gained
# rosmap rs3763040              SP3 KLF3
#        rs3875089              TEAD1
#        rs335929               NR2E1                              NR2E1  OTX1
#        rs3763043              TEAD1 NEUROD2 TCF12 OTX1
#        rs9951307              NEUROD2   STAT1                    ATF7 NFE2L2 ID4 ID2 STAT1
#
# cer    rs335929               MSX2                               PRRX1  MSX2  HOPX
#        rs3763043              PRRX1  HOPX
#        rs9951307              TRPS1 MEIS2 ZGLP1 STAT3            ID4 ID2 ID1 STAT3
#
# tcx    rs3875089              TEAD1
#        rs335929               EMX2                               EMX2
#        rs3763043              TEAD1  ELK3
#        rs9951307              TRPS1 STAT3                        STAT3
#
#
#
# cer rs335929
#     tfs with motif loss: MSX2
#     tfs with motif gain: PRRX1  MSX2  HOPX
# cer model  rs3763043
#     tfs with motif loss: PRRX1  HOPX
# cer model rs9951307
#     tfs with motif loss: TRPS1 MEIS2 ZGLP1 STAT3
#     tfs with motif gain:   ID4   ID2   ID1 STAT3
#
# tcx model rs3875089
#    tfs with motif loss: TEAD1
# tcx model rs335929
#    tfs with motif loss: EMX2
#    tfs with motif gain: EMX2
# tcx model  rs3763043
#    tfs with motif loss: TEAD1  ELK3
# tcx model  rs9951307
#    tfs with motif loss:  TRPS1 STAT3
#    tfs with motif gain: ID4    STAT3
