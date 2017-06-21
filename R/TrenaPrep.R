.TrenaPrep <- setClass ("TrenaPrep",
                        representation = representation(
                           targetGene="character",
                           chromosome="character",
                           chromStart="numeric",
                           chromEnd="numeric",
                           regulatoryRegionSources="list")
                        )
#------------------------------------------------------------------------------------------------------------------------
setGeneric('addGraph',         signature='obj', function(obj, graph, modelNames=list()) standardGeneric ('addGraph'))
setGeneric('httpAddGraph',     signature='obj', function(obj, graph) standardGeneric ('httpAddGraph'))
setGeneric('loadStyle',        signature='obj', function(obj, filename) standardGeneric ('loadStyle'))
#------------------------------------------------------------------------------------------------------------------------
TrenaPrep = function(targetGene, chromosome, chromStart, chromEnd, regulatoryRegionSources=list(), quiet=TRUE)
{

   obj <- .TrenaPrep(targetGene=targetGene,
                     chromosome=chromosome,
                     chromStart=chromStart,
                     chromEnd=chromEnd,
                     regulatoryRegionSources=regulatoryRegionSources)

   obj

} # constructor
#------------------------------------------------------------------------------------------------------------------------

