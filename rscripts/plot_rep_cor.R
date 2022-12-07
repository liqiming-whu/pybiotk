#!/usr/bin/env Rscript
library(argparse)
library(ggplot2)
library(viridis)

parser <- ArgumentParser()
parser$add_argument(dest="input", type="character", nargs='+', help='input files.')
parser$add_argument('-s', dest='samples', type="character", nargs='+', default=NULL, help='sample names.')
parser$add_argument('--noheader', dest='noheader', action='store_true', help='if no header.')
parser$add_argument('-b', dest='by_columns', type="integer", nargs='+', required=TRUE, help='use columns to merge.')
parser$add_argument('-c', dest='column', type="integer", required=TRUE, help='use column to calculate.')
parser$add_argument('-o', dest='outfig', type="character", required=TRUE, help='output fig')
parser$add_argument('--drop_low', dest='drop', action="store_true", help='remove < 1 rows.')
parser$add_argument('--density', dest='density', action="store_true", help='plot density')

args <- parser$parse_args()

input_files <- args$input
samples_name <- args$samples
by_columns <- args$by_columns
column <- args$column
columns <- append(by_columns, column)


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

stopifnot(length(input_files) == 2)

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

data1 <- read_data(input_files[1], "x")
data2 <- read_data(input_files[2], "y")

all_data <- merge(data1, data2, by=bycolumns, all=T)

data <- subset(all_data,select=-bycolumns)

data[is.na(data)] <- 0
data <- data[rowSums(data) > 0,]
if(args$drop) data <- data[apply(data, 1, min) >= 1,]
correlation <- round(cor(data$x, data$y, use="everything", method="pearson"), 2)

label <- paste0("R = ", as.character(correlation))

log10_x <- log10(data$x)[!is.infinite(log10(data$x))]

p <- ggplot(data) +
    geom_point(aes(x=log10(x), y=log10(y)), color="blue") +
    geom_text(aes(x=min(log10_x) + 0.1*(max(log10_x) - min(log10_x)), y=max(log10(y))), label=label, size = 2, check_overlap=T) +
    labs(x=paste0(samples_name[1]," log10(FPKM)"), y=paste0(samples_name[2], " log10(FPKM)")) +
    theme_bw() +
    theme(panel.grid=element_blank(), legend.position="none", text = element_text(size = 8))

if(args$density) {
    p <- p + stat_density2d(geom = 'polygon', aes(x=log10(x), y=log10(y), fill=..level..)) +
        scale_fill_viridis(begin = 0.1, end = 0.9, option = "turbo", alpha = 0.8)
}

ggsave(args$outfig, p, width=5, height=5, dpi=300, units="cm")
