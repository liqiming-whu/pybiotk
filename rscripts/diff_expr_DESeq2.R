#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument('input', type="character", nargs='+', help='control and case')
parser$add_argument('-s', dest='samples', type="character", nargs='+', default=NULL, help='sample names.')
parser$add_argument('-g', dest='group', type="character", nargs='+', default=NULL, help='sample groups.')
parser$add_argument('-b', dest='by', type="integer", default=1, help='use columns to merge.')
parser$add_argument('-c', dest='column', type="integer", default=2, help='use column to calculate.')
parser$add_argument('-o', dest='outdir', type="character", required=TRUE, help='output dir')
parser$add_argument('-l', dest='log2fc', type="double", default=1.0, help='log2fc [default=1].')
parser$add_argument("-p", dest="padj_value", type="double", default=0.05, help="cutoff of padj_value [default=0.05].")
parser$add_argument('--spike-in', dest="spike_in", type="character", nargs='+', default=NULL, help='control and case spike_in.')

args <- parser$parse_args()
input_files <- args$input
samples_name <- args$samples
columns <- c(args$by, args$column)
outdir <- args$outdir
log2fc <- args$log2fc
padj_value <- args$padj_value
spike_in <- args$spike_in

input_files_len <- length(input_files)

filename_preix <- function(file_path) {
    name <- basename(file_path)
    name_vec <- strsplit(name, split=".", fixed=T)[[1]]
    return(name_vec[1])
}

if(!is.null(samples_name)) {
    stopifnot(input_files_len == length(samples_name))
}else {
    samples_name <- c()
    for(filename in input_files) {
        samples_name <- append(samples_name, filename_preix(filename))   
    }
}

library(DESeq2)

read_data <- function(filename, sample) {
    table <- read.table(filename, sep="\t", header=T)
    table <- table[, columns]
    names(table) <- c("id", sample)
    table <- table[complete.cases(table),]
    return(table)
}

table <- read_data(input_files[1], samples_name[1])
for (idx in seq_along(input_files)[-1]) {
    filename <- input_files[idx]
    sample_name <- samples_name[idx]
    data <- read_data(filename, sample_name)
    table <- merge(table, data, by="id", all=T)
}

table[is.na(table)] <- 0 

if (is.null(args$group)) {
    group_length <- input_files_len / 2
    group <- rep(c('control', 'treat'), each=group_length)
}else {
    group <- args$group
}

counts_table <- table[,-1]
rownames(counts_table) <- table[,1]

counts_table <- round(counts_table, digits=0) 
counts_table <- as.matrix(counts_table)
condition <- factor(group, levels=c("control", "treat"))
coldata <- data.frame(row.names=colnames(counts_table), condition)

dds <- DESeqDataSetFromMatrix(counts_table, coldata, design=~condition)

if(!is.null(spike_in)) {

    spike_table <- read_data(spike_in[1], samples_name[1])
    for (idx in seq_along(spike_in)[-1]) {
        filename <- spike_in[idx]
        sample_name <- samples_name[idx]
        data <- read_data(filename, sample_name)
        spike_table <- merge(spike_table, data, by="id", all=T)
    }

    spike_table[is.na(spike_table)] <- 0 
    spike_counts_table <- spike_table[,-1]
    rownames(spike_counts_table) <- spike_table[,1]

    spike_counts_table <- round(spike_counts_table, digits=0) 
    spike_counts_table <- as.matrix(spike_counts_table)
    spike_coldata <- data.frame(row.names=colnames(spike_counts_table), condition)
    spike_dds <- DESeqDataSetFromMatrix(spike_counts_table, spike_coldata, design=~condition)
    spike_dds <- estimateSizeFactors(spike_dds)

    sizeFactors(dds) <- sizeFactors(spike_dds)
}

dds <- DESeq(dds) # This function performs a default analysis through the steps:
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
# dds <- nbinomWaldTest(dds)

res <- results(dds)

if(!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

setwd(outdir)

res <- res[order(res$padj),]
res <- merge(as.data.frame(res), as.data.frame(counts(dds, normalize=T)),
             by="row.names", sort=F)
deseq_res <- data.frame(res)

up_diff <- subset(deseq_res, (padj < padj_value) & (log2FoldChange > log2fc))
down_diff <- subset(deseq_res, (padj < padj_value) & (log2FoldChange < -log2fc))
sig_result <- subset(deseq_res, (padj < padj_value) & (abs(log2FoldChange) > log2fc))
all_result <- subset(deseq_res, baseMean != 0)

fpm_table <- fpm(dds)

write.table(up_diff, "up.xls", quote=F, row.names=F, col.names=T, sep="\t")
write.table(down_diff, "down.xls", quote=F, row.names=F, col.names=T, sep="\t")
write.table(sig_result, "sig.xls", quote=F, row.names=F, col.names=T, sep="\t")
write.table(all_result, "all.xls", quote=F, row.names=T, col.names=T, sep="\t")
write.table(fpm_table, "cpm.xls", quote=F, row.names=T, col.names=T, sep="\t")