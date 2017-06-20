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
init <- function()
{
   tv <- TrenaViz()
   tv

} # init
#------------------------------------------------------------------------------------------------------------------------
displayFootprintsFromDatabases <- function(tv, chrom, start, end)
{
   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", host="whovian")
   dboi <- c("brain_wellington_16", "brain_hint_16", "brain_hint_20", "brain_wellington_20")
   checkTrue(all(dboi %in% dbGetQuery(db, "select datname from pg_database")$datname))
   dbDisconnect(db)

   genome.db.uri    <- "postgres://bddsrds.globusgenomics.org/hg38"   # has gtf and motifsgenes tables
   for(dbName in dboi){
     footprint.db.uri <- sprintf("postgres://whovian/%s", dbName)
     fpf <- FootprintFinder(genome.db.uri, footprint.db.uri, quiet=FALSE)
     tbl.fp <- getFootprintsInRegion(fpf, chrom, start, end)
     printf("%d footprints in %d bases from %s", nrow(tbl.fp), 1  + end - start, dbName)
     addBedTrackFromDataFrame(tv, dbName, tbl.fp[, c("chrom", "start", "endpos", "name", "score1")], color="blue")
     closeDatabaseConnections(fpf)
     }

} # displayFootprintsFromDatabases
#------------------------------------------------------------------------------------------------------------------------
# brain_hint_16 has one fp that intersects with snp2(26865469): chr18:26865467-26865487
# brain_hint_20: chr18:26865467-26865487
runDemo <- function()
{
   tv <- init()

   tbl.snp.display <- data.frame(chromosome=tbl.snp$chromosome, start=tbl.snp$loc, end=tbl.snp$loc, name=tbl.snp$snp,
                                 stringsAsFactors=FALSE)
   addBedTrackFromDataFrame(tv, "snp", tbl.snp.display[, c("chromosome", "start", "end", "name")], color="purple")

   chrom <- unique(tbl.snp$chromosome)
   region.start <- min(tbl.snp$loc[2] - 100)
   region.end   <- max(tbl.snp$loc[2] + 100)
   snp.start <- min(tbl.snp$loc) - 1000;
   snp.end   <- max(tbl.snp$loc) + 1000;

   #tbl.snp.display <- data.frame(chromosome=tbl.snp$chromosome, start=tbl.snp$loc, end=tbl.snp$loc+1, name=tbl.snp$snp)
   #addBedTrackFromDataFrame(tv, "snp", tbl.snp.display, color="purple")

   start <- snp.start
   end   <- snp.end

   tbl.fp <- displayFootprintsFromDatabases(tv, chrom, start, end)
   region <- sprintf("%s:%d-%d", tbl.snp$chromosome[1], min(tbl.snp$loc) - 1000, max(tbl.snp$loc) + 1000)
   showGenomicRegion(tv, region)

   tv

} # runDemo
#------------------------------------------------------------------------------------------------------------------------
