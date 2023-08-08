#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument('input', type="character", help='bedfile')
parser$add_argument('-s', dest='species', type="character", default=NULL, help='hg38 or mm10')
parser$add_argument('-g', dest='gtf', type="character", default=NULL, help='gtf file')
parser$add_argument('-t', dest='tss', type="integer", nargs="+", default=c(-3000, 3000), help='tssRegion')
parser$add_argument('-o', dest='outdir', type="character", required=TRUE, help='output dir')
parser$add_argument('--strand', dest='strand', action='store_true', help='same strand')
parser$add_argument('-f', dest='fig', default=NULL, help='pdf or png')

args <- parser$parse_args()
bedfile <- args$input
species <- args$species
gtf <- args$gtf
tss <- args$tss
outdir <- args$outdir
samestrand <- args$strand
fig <- args$fig
ylab <- args$ylab

stopifnot(!(is.null(species) & is.null(gtf)))

library(ChIPseeker)

if(!is.null(gtf)) {
    library(GenomicFeatures)
    txdb <- makeTxDbFromGFF(file=gtf)
    annoDb <- NULL
}else {
    if(species == "hg38") {
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    annoDb <- "org.Hs.eg.db"
    }else if(species == "mm10") {
    library(org.Mm.eg.db)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    annoDb <- "org.Mm.eg.db"
    }else stop("species not support")
}

peak <- readPeakFile(bedfile)
peakAnno <- annotatePeak(
    peak,
    tssRegion = tss,
    TxDb = txdb,
    level = "transcript",
    assignGenomicAnnotation = TRUE,
    genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
    annoDb = annoDb,
    addFlankGeneInfo = FALSE,
    flankDistance = 5000,
    sameStrand = samestrand,
    ignoreOverlap = FALSE,
    ignoreUpstream = FALSE,
    ignoreDownstream = FALSE,
    overlap = "TSS",
    verbose = TRUE
)

if(!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

setwd(outdir)

anno <- as.data.frame(peakAnno)
write.table(anno,file="anno.xls",sep='\t',quote=F)

if(!is.null(fig)) {

    if(fig=="pdf") {
        pdf("AnnoPie.pdf")
    }else {
        png("AnnoPie.png")
    }
    plotAnnoPie(peakAnno)
    dev.off()

    if(fig=="pdf") {
        pdf("VennPie.pdf")
    }else {
        png("VennPie.png")
    }
    vennpie(peakAnno)
    dev.off()

    library(ggplot2)

    p <- plotAnnoBar(peakAnno, xlab="", ylab='Percentage(%)', title="Feature Distribution")
    ggsave(p, file=paste0("AnnoBar.", fig))

    p <- plotDistToTSS(peakAnno, distanceColumn="distanceToTSS", xlab="", ylab="Binding sites (%) (5'->3')",
        title="Distribution of transcription factor-binding loci relative to TSS")
    ggsave(p, file=paste0("dis2TSS.", fig))

    p <- upsetplot(peakAnno, vennpie = TRUE)
    ggsave(p, file=paste0("upset.", fig), width = 8, height = 10, bg="white")

    p <- covplot(peak)
    ggsave(p, file=paste0("cov.", fig), width = 20, height = 20)
}
