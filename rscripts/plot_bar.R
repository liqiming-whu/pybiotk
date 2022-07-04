#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument('input', type="character", help='input file.')
parser$add_argument('-x', dest='x_lable', type="character", default='', help='x lab.')
parser$add_argument('-y', dest='y_lable', type="character", default='', help='y lab.')
parser$add_argument('-t', dest='title', type="character", default='', help='title.')
parser$add_argument('--noheader', dest='noheader', action='store_true', help='if no header.')
parser$add_argument('-c', dest='column', type="integer", nargs="+", default=c(1,2), help='use column to calculate.')
parser$add_argument('-o', dest='output', type="character", required=TRUE, help='output.')
parser$add_argument('-w', dest='width', type="integer", default=12, help='outfig width.')
parser$add_argument('-l', dest='length', type="integer", default=9, help='outfig length.')
parser$add_argument('-b', dest='bin_width', type="double", default=0.5, help='bin width.')
parser$add_argument('-f', dest='font_size', type="integer", default=12, help='font size.')
parser$add_argument('--line', dest='line', action='store_true', help='add a line.')
parser$add_argument('--text', dest='text', action='store_true', help='add text.')
parser$add_argument('--ylog10', dest='ylog10', action='store_true', help='log10(y).')

args <- parser$parse_args()

library(ggplot2)
options(scipen=3)

if(args$noheader) {
    header <- FALSE
}else {
    header <- TRUE
}

stopifnot(length(args$column)==2)

df <- read.table(args$input, sep="\t", header=header)
df <- df[, args$column]
df <- df[complete.cases(df),]
df[,1] <- factor(df[,1], levels=df[, 1], ordered=TRUE)

if(args$ylog10) plot.y <- log10(df[,2]) else plot.y <- df[,2]

p <- ggplot(df) +
    geom_bar(aes(x=df[,1], y=plot.y), stat="identity", width=args$bin_width, fill="#FF000080") +
    labs(x=args$x_lable, y=args$y_lable, title=args$title) +
    theme_bw(base_size = 12) +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust=0.5, face="bold", color="black", size=14)) +
    theme(axis.text.x = element_text(color='black',size=args$font_size)) +
    theme(axis.text.y = element_text(color='black',size=args$font_size)) +
    theme(panel.grid = element_blank())

if(args$line) {
    p <- p + geom_line(aes(x=df[,1], y=plot.y, group=1), stat="identity", color="red")
}

if(args$text) {
    p <- p + geom_text(aes(label=df[,2], x=df[,1], y=1.02*plot.y), color="black")
}

ggsave(p, filename = args$output, width = args$width, height = args$length)
