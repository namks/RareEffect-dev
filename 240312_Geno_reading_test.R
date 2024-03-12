# R code to read genotype from PLINK file
# Should be run in SAIGE_1.3.3 env

library(SAIGE, lib.loc = "~/utils/SAIGE")
library(data.table)

setwd("~/utils/SAIGE/extdata/input")
bedFile <- "plinkforGRM_1000samples_10kMarkers.bed"
bimFile <- "plinkforGRM_1000samples_10kMarkers.bim"
famFile <- "plinkforGRM_1000samples_10kMarkers.fam"
AlleleOrder <- "alt-first"
fam <- fread(famFile)
sampleID <- fam$V2

objGeno = setGenoInput(bgenFile = "",
    bgenFileIndex = "",
    vcfFile = "",   #not activate yet
    vcfFileIndex = "",
    vcfField = "",
    savFile = "",
    savFileIndex = "",
    sampleFile = "",
    bedFile=bedFile,
    bimFile=bimFile,
    famFile=famFile,
    idstoIncludeFile = "",
    rangestoIncludeFile = "",
    chrom = "",
    AlleleOrder = AlleleOrder,
    sampleInModel = sampleID)

n <- length(sampleID)
t_GVec <- rep(0, n)

SAIGE::Unified_getOneMarker(t_genoType = objGeno$genoType,
        t_gIndex_prev = objGeno$markerInfo$genoIndex_prev[1],
        t_gIndex = objGeno$markerInfo$genoIndex[1],
        t_ref = "2",
        t_alt = "1",
        t_marker = objGeno$markerInfo$ID[1],
        t_pd = objGeno$markerInfo$POS[1],
        t_chr = toString(objGeno$markerInfo$CHROM[1]),
        t_altFreq = 0,
        t_altCounts = 0,
        t_missingRate = 0,
        t_imputeInfo = 0,
        t_isOutputIndexForMissing = TRUE,
        t_indexForMissing = 0,
        t_isOnlyOutputNonZero = FALSE,
        t_indexForNonZero = 0,
        t_GVec = t_GVec,
        t_isImputation = FALSE)
