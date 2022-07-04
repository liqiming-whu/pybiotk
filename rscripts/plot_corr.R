#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument('input', type="character", nargs='+', help='input files.')
parser$add_argument('-s', dest='samples', type="character", nargs='+', default=NULL, help='sample names.')
parser$add_argument('--noheader', dest='noheader', action='store_true', help='if no header.')
parser$add_argument('-b', dest='by_columns', type="integer", nargs='+', required=TRUE, help='use columns to merge.')
parser$add_argument('-c', dest='column', type="integer", required=TRUE, help='use column to calculate.')
parser$add_argument('-o', dest='outdir', type="character", required=TRUE, help='output dir')
parser$add_argument('-p', dest='pdf', action='store_true', help='pdf or png')
parser$add_argument('-m', dest='min_lable', type="double", default=NULL, help='min lable')
parser$add_argument('--break', dest='breaklen', type="double", default=0.05, help='break length')

args <- parser$parse_args()

library(pheatmap)          
library(corrplot)                                                                                                                                                               
library(RColorBrewer)

input_files <- args$input
samples_name <- args$samples
by_columns <- args$by_columns
column <- args$column
columns <- append(by_columns, column)
outdir <- args$outdir
pdf <- args$pdf
min_lable <- args$min_lable

if(!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

if(args$noheader) {
    header <- FALSE
}else {
    header <- TRUE
}

filename_preix <- function(file_path) {
    name <- basename(file_path)
    name_vec <- strsplit(name, split=".", fixed=T)[[1]]
    return(name_vec[1])
}

if(!is.null(samples_name)) {
    stopifnot(length(input_files) == length(samples_name))
}else {
    samples_name <- c()
    for(filename in input_files) {
        samples_name <- append(samples_name, filename_preix(filename))   
    }
}

by_columns_len <- length(by_columns)
bycolumns <- 1:by_columns_len

read_data <- function(filename, sample) {
    table <- read.table(filename, sep="\t", header=header)
    table <- table[, columns]
    rawcolnames <- colnames(table)[1:by_columns_len]
    names(table) <- c(rawcolnames, sample)
    table <- table[complete.cases(table),]
    return(table)
}

all_data <- read_data(input_files[1], samples_name[1])
for (idx in seq_along(input_files)[-1]) {
    filename <- input_files[idx]
    sample_name <- samples_name[idx]
    table <- read_data(filename, sample_name)
    all_data <- merge(all_data, table, by=bycolumns, all=T)
}

us_count <- subset(all_data,select=-bycolumns)
if(by_columns_len == 1){
    rownames(us_count) <- all_data[,1]
}
us_count[is.na(us_count)] <- 0
us_count <- us_count[rowSums(us_count) > 0,]
us_count <- cor(us_count, use='everything', method='pearson')
write.table(us_count, file=file.path(outdir, "corr.tsv"), sep="\t")

if(pdf) {
    pdf(file.path(outdir, "heatmap.pdf"))
}else {
    png(file.path(outdir, "heatmap.png"))
}

if (is.null(min_lable)) {
pheatmap(
  us_count,
  clustering_distance_rows = "correlation", 
  clustering_distance_cols = "correlation", 
  )
}else {
    breaksList <- seq(min_lable, 1, by=args$breaklen)
pheatmap(
  us_count,
  clustering_distance_rows = "correlation", 
  clustering_distance_cols = "correlation", 
  color = colorRampPalette(rev(brewer.pal(n=7,name = "RdYlBu")))(length(breaksList)),
  breaks = breaksList
  )
}
dev.off()

if(pdf) {
    pdf(file.path(outdir, "corr.pdf"))
}else {
    png(file.path(outdir, "corr.png"))
}
col <- colorRampPalette(c("green", "white", "red"))
corrplot(us_count, tl.col = "black", order="origin", col=col(100),)
dev.off()
