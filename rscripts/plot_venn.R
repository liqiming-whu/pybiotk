#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument('input', type="character", nargs="+", help='input file.')
parser$add_argument('-s', dest='samples', nargs="+", default=NULL, type="character", help='sample names.')
parser$add_argument('-t', dest='title', type="character", default='', help='title.')
parser$add_argument('--noheader', dest='noheader', action='store_true', help='if no header.')
parser$add_argument('-c', dest='column', type="integer", default=1, help='use column to calculate.')
parser$add_argument('-o', dest='outdir', type="character", required=TRUE, help='output dir.')
parser$add_argument('-p', dest='pdf', action='store_true', help='pdf or png')
parser$add_argument('--name_size', dest='name_size', type="integer", default=5, help='name size.')
parser$add_argument('--text_size', dest='text_size', type="integer", default=4, help='text size.')
parser$add_argument('--show_percentage', dest='show_percentage', action='store_true', help='show percentage')

args <- parser$parse_args()

library(ggvenn)
library(ggplot2)

outdir <- args$outdir
pdf <- args$pdf

filesnum <- length(args$input)

samples <- args$samples

if(is.null(samples)) {
    samples <- 1:filesnum
}

stopifnot(filesnum==length(samples))

if(args$noheader) {
    header <- FALSE
}else {
    header <- TRUE
}

data <- list()

for(idx in 1:filesnum) {
    table <- read.table(args$input[idx], sep="\t", header=header)
    table <- table[, args$column]
    data[[samples[idx]]] <- table
}

if(!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

if(pdf) {
    pdf(file.path(outdir, "venn.pdf"), width=540, height=480)
}else {
    png(file.path(outdir, "venn.png"), width=540, height=480)
}

ggvenn(data, samples, show_percentage = args$show_percentage,
    stroke_color = "black",
    fill_color = c("#E41A1C","#1E90FF","#FF8C00", "#4DAF4A","#984EA3"),
    set_name_size = args$name_size,
    text_size = args$text_size
    ) + labs(title=args$title) + theme(plot.title = element_text(face="bold", hjust=0.5, color="black", size=14))

dev.off()
