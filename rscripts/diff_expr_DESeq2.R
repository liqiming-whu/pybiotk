#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument('input', type="character", nargs='+', help='control and case')
parser$add_argument('--info', dest='info', type="character", default=NULL, help='gene info table with a header, first column must be same as -b option.')
parser$add_argument('-s', dest='samples', type="character", nargs='+', default=NULL, help='sample names.')
parser$add_argument('-g', dest='group', type="character", nargs='+', default=NULL, help='sample groups.')
parser$add_argument('--level', dest='level', type="character", nargs='+', default=NULL, help='sample groups.')
parser$add_argument('-b', dest='by', type="integer", default=1, help='use columns to merge.')
parser$add_argument('-c', dest='column', type="integer", default=2, help='use column to calculate.')
parser$add_argument('-o', dest='outdir', type="character", required=TRUE, help='output dir')
parser$add_argument('-l', dest='log2fc', type="double", default=1.0, help='log2fc [default=1].')
parser$add_argument("-p", dest="padj_value", type="double", default=0.05, help="cutoff of padj_value [default=0.05].")
parser$add_argument('--spike-in', dest="spike_in", type="character", nargs='+', default=NULL, help='control and case spike_in.')
parser$add_argument('--coef', dest="coef", type="character", nargs='+', default=NULL, help='control and case group pair, like control_vs_treat')

args <- parser$parse_args()
input_files <- args$input
info <- args$info
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

info.table <- NULL
if(!is.null(info)) {
    raw.table <- read.table(info, header=T)
    info.table <- raw.table[,-1]
    rownames(info.table) <- raw.table[,1]
}


add.info <- function(data) {
    if(!is.null(info.table)) {
        data <- merge(data, info.table, by="row.names", sort=F)
        stopifnot(colnames(data)[1] == "Row.names")
        rownames(data) <- data[, 1]
        data <- data[,-1]
    }
    return(data)
}


read_data <- function(filename, sample) {
    table <- read.table(filename, sep="\t", header=T)
    table <- table[, columns]
    names(table) <- c("id", sample)
    table <- na.omit(table)
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
if(is.null(args$level)) {
    condition <- factor(group)
}else {
    condition <- factor(group, levels=args$level)
}

coldata <- data.frame(row.names=colnames(counts_table), condition)
dds <- DESeqDataSetFromMatrix(counts_table, coldata, design=~condition)


min.group <- min(as.data.frame(table(group))$Freq)

keep <- rowSums(counts(dds) >= 10) >= min.group
dds <- dds[keep,]

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

raw.counts <- as.data.frame(counts(dds, normalize=F))
dds <- DESeq(dds) # This function performs a default analysis through the steps:
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
# dds <- nbinomWaldTest(dds)

norm.counts <- as.data.frame(counts(dds, normalize=T))

if(!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

setwd(outdir)
fpm_table <- fpm(dds)

raw.counts <- add.info(raw.counts)
norm.counts <- add.info(norm.counts)
fpm_table <- add.info(fpm_table)

info.columns <- c()
if(ncol(norm.counts)-input_files_len > 0) info.columns <- colnames(norm.counts)[(input_files_len+1):ncol(norm.counts)]


write.table(raw.counts, "raw.counts.xls", quote=F, row.names=T, col.names=T, sep="\t")
write.table(norm.counts, "norm.counts.xls", quote=F, row.names=T, col.names=T, sep="\t")
write.table(fpm_table, "cpm.xls", quote=F, row.names=T, col.names=T, sep="\t")

group.diff <- function(control, treat) {
res <- results(dds, contrast=c("condition", treat, control))
res <- res[order(res$padj),]

res.dir = file.path(paste(control, "vs", treat, sep="-"))
if(!dir.exists(res.dir)) dir.create(res.dir, recursive=TRUE)

res <- merge(as.data.frame(res), norm.counts[,append(rownames(coldata)[coldata$condition %in% c(control, treat)], info.columns)], by="row.names", sort=F)
deseq_res <- data.frame(res)
up_diff <- subset(deseq_res, (padj < padj_value) & (log2FoldChange > log2fc))
down_diff <- subset(deseq_res, (padj < padj_value) & (log2FoldChange < -log2fc))
sig_result <- subset(deseq_res, (padj < padj_value) & (abs(log2FoldChange) > log2fc))
all_result <- subset(deseq_res, baseMean != 0)

if(nrow(sig_result) > 0){
    sig_result$states <- "Up"
    sig_result[sig_result$log2FoldChange > 0,]$states <- "Up"
    sig_result[sig_result$log2FoldChange < 0,]$states <- "Down"
}
write.table(up_diff, file.path(res.dir, "up.xls"), quote=F, row.names=F, col.names=T, sep="\t")
write.table(down_diff, file.path(res.dir, "down.xls"), quote=F, row.names=F, col.names=T, sep="\t")
write.table(sig_result, file.path(res.dir, "sig.xls"), quote=F, row.names=F, col.names=T, sep="\t")
write.table(all_result, file.path(res.dir, "all.xls"), quote=F, row.names=F, col.names=T, sep="\t")
}

if(is.null(args$coef)) {
    combo <- t(combn(levels(condition), 2))
    for(i in 1:nrow(combo)) {
        group.diff(combo[i, 1], combo[i, 2])
        }
}else {
for(i in seq_along(args$coef)) {
    pair <- unlist(strsplit(args$coef[i], "_vs_", fixed=T))
    group.diff(pair[1], pair[2])
}
}
