#!/usr/bin/env Rscript
library(argparse)
library(ggplot2)

parser <- ArgumentParser()
parser$add_argument('input', type="character", help='input file.')
parser$add_argument('-x', dest='x_lable', type="character", default='', help='x lab.')
parser$add_argument('-y', dest='y_lable', type="character", default='', help='y lab.')
parser$add_argument('-t', dest='title', type="character", default='', help='title.')
parser$add_argument('--noheader', dest='noheader', action='store_true', help='if no header.')
parser$add_argument('-o', dest='output', type="character", required=TRUE, help='output.')
parser$add_argument('-w', dest='width', type="integer", default=12, help='outfig width.')
parser$add_argument('-l', dest='length', type="integer", default=9, help='outfig length.')
parser$add_argument('-b', dest='bins', type="integer", default=30, help='bin nums.')
parser$add_argument('-f', dest='font_size', type="integer", default=12, help='font size.')

args <- parser$parse_args()

options(scipen=3)

if(args$noheader) {
    header <- FALSE
}else {
    header <- TRUE
}

df <- read.table(args$input, sep="\t", header=header)
names(df) <- c("Length")
stats <- boxplot.stats(df$Length)
df <- subset(df, Length < stats$stats[5] & Length > stats$stats[1])

p <- ggplot(df, aes(x=Length)) +
    geom_histogram(aes(y=..density..), position = "identity", fill="#FF000080", color="black", bins=args$bins, alpha=0.5) +
    geom_density() +
    labs(x=args$x_lable, y=args$y_lable, title=args$title) +
    theme_bw(base_size = 12) +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust=0.5, face="bold", color="black", size=14)) +
    theme(axis.text.x = element_text(color='black',size=args$font_size)) +
    theme(axis.text.y = element_text(color='black',size=args$font_size)) +
    theme(panel.grid = element_blank())

ggsave(p, filename = args$output, width = args$width, height = args$length)
