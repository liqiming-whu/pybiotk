#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument('input', type="character", help='inputfile')
parser$add_argument('-c', dest='columns', type="integer", nargs="+", help='use columns to plot, order:gene, log2fc, other')
parser$add_argument('-o', dest='output', type="character", required=TRUE, help='output')
parser$add_argument('-t', dest='title', type="character", default="", help='title')
parser$add_argument('-l', dest='log2fc', type="double", default=1.0, help='cutoff of log2fc')
parser$add_argument('-p', dest='pvalue', type="double", default=0.05, help='cutoff of pvalue')
parser$add_argument('-f', dest='font_size', type="integer", default=8, help='font size')
parser$add_argument('--nudge_x', dest='nudge_x', type="double", default=0, help='nudge_x')
parser$add_argument('--nudge_y', dest='nudge_y', type="double", default=0, help='nudge_y')
parser$add_argument('--nopvalue', dest='nopvalue', action='store_true', help='plot log2fc and reads')
parser$add_argument('--rankplot', dest='rankplot', action='store_true', help='plot log2fc and rank')
parser$add_argument('--labels', dest='labels', type="character", nargs="+", default=NULL, help="labels list")

args <- parser$parse_args()

data_table <- args$input
columns <- args$columns
output <- args$output
title <- args$title
log2fc <- args$log2fc
p_value <- args$pvalue
font_size <- args$font_size
nudge_x <- args$nudge_x
nudge_y <- args$nudge_y
labels <- args$labels

# library(ggthemes)
library(ggplot2)
library(ggrepel)

data <- read.table(file=data_table, header=T)
if(is.null(columns)) data <- data[c("Row.names", "log2FoldChange", "pvalue")] else data <- data[, columns]
data <- data[complete.cases(data),]

if(args$nopvalue) {
names(data) <- c("gene_id", "log2FoldChange", "reads")
data$threshold <- as.factor(ifelse(data$reads > 3 & abs(data$log2FoldChange) > log2fc,
                            ifelse((data$log2FoldChange) > log2fc ,'up','down'),'not'))

if(!is.null(labels)){
need_label <- data[data$gene_id%in%labels,]
need_label$color = "black"
}else{
label_data <- data[!data$threshold=='not',]
label_data$pos <- 'middle'
label_data <- label_data[order(label_data$log2FoldChange), ]
label_data$rank <- 1:nrow(label_data)
dnrow <- nrow(label_data[label_data$threshold=='down', ])
if(dnrow > 5) down_cut <- 5 else down_cut <- dnrow
label_data[1:dnrow, 'pos'] <- 'down'
label_data[1:dnrow, 'color'] <- '#546de5'

unrow <- nrow(label_data[label_data$threshold=='up', ])
if(unrow > 5) up_cut <- 5 else up_cut <- unrow

label_data[(nrow(label_data)-up_cut+1):nrow(label_data), 'pos'] <- 'up'
label_data[(nrow(label_data)-up_cut+1):nrow(label_data), 'color'] <- '#ff4757'

need_label <- label_data[c(1:down_cut, (nrow(label_data)-up_cut+1):nrow(label_data)),]
}
if(args$rankplot){
rank_data <- data[data$reads > 3]
rank_data$pos <- 'not'
rank_data <- rank_data[order(rank_data$log2FoldChange), ]
rank_data$rank <- 1:nrow(rank_data)
dnrow <- nrow(rank_data[rank_data$threshold=='down', ])
if(dnrow > 5) down_cut <- 5 else down_cut <- dnrow
rank_data[1:dnrow, 'pos'] <- 'down'
rank_data[1:dnrow, 'color'] <- '#546de5'

unrow <- nrow(rank_data[rank_data$threshold=='up', ])
if(unrow > 5) up_cut <- 5 else up_cut <- unrow

rank_data[(nrow(rank_data)-unrow+1):nrow(rank_data), 'pos'] <- 'up'
rank_data[(nrow(rank_data)-unrow+1):nrow(rank_data), 'color'] <- '#ff4757'

if(!is.null(labels)){
need_rank <- rank_data[rank_data$gene_id%in%labels,]
need_rank$color = "black"
}else{
need_rank <- rank_data[c(1:down_cut, (nrow(rank_data)-up_cut+1):nrow(rank_data)),]
}

up.count <- nrow(rank_data[rank_data$threshold=="up",])
down.count <- nrow(rank_data[rank_data$threshold=="down",])
max.y <- max(rank_data$log2FoldChange)
min.y <- min(rank_data$log2FoldChange)
max.x <- max(rank_data$rank)
plot <- ggplot(data=rank_data, aes(x=rank, y=log2FoldChange,color=pos, fill=pos)) +
    scale_color_manual(values=c("down"="#546de5", "not"="#d2dae2", "up"="#ff4757")) +
    geom_point(size=0.1) +
    geom_hline(yintercept=c(-log2fc, log2fc), linetype=1, col="grey", linewidth=0.6) +
    geom_hline(yintercept=0, linetype=1, col="grey", linewidth=0.6) +
    geom_text(x=max.x*0.9, y=min.y, label= paste0("down: ", as.character(down.count)), size = 2, color="#546de5", check_overlap=T) +
    geom_text(x=max.x*0.1, y=max.y, label= paste0("up: ", as.character(up.count)), size = 2, color="#ff4757", check_overlap=T) +
    geom_text_repel(inherit.aes=F, data=need_rank,
        aes(x=rank, y=log2FoldChange, label=gene_id), color=need_rank$color, na.rm=T, segment.size = 0.3, segment.colour="black", 
        box.padding = unit(0.35, "lines"), size=2, direction='both', nudge_x=nudge_x, nudge_y=nudge_y, min.segment.length = 0, max.overlaps=Inf) +
    # scale_y_continuous(breaks=c(-6, -3, -1, 0, 1, 3, 6)) +
    theme_bw(base_size=8) +
    theme(legend.position="right",
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text= element_text(color="black", size=8),
          plot.title = element_text(hjust=0.5, color="black", size=font_size),
          axis.text.x = element_text(color="black", size=font_size),
          axis.text.y = element_text( color="black", size=font_size),
          axis.title.x = element_text(color="black", size=font_size),
          axis.title.y = element_text(color="black", size=font_size))+
    labs(x="rank",y="log2(fold change)",title=title)

}else{

up.count <- nrow(data[data$threshold=="up",])
down.count <- nrow(data[data$threshold=="down",])
max.y <- max(log10(data$reads))
max.x <- max(data$log2FoldChange)
min.x <- min(data$log2FoldChange)

plot <- ggplot(data=data, 
               aes(x=log2FoldChange, y=log10(reads), 
                   color=threshold, fill=threshold)) +
    scale_color_manual(values=c("down"="#546de5", "not"="#d2dae2", "up"="#ff4757")) +
    # geom_point(alpha=0.4, size=0.2) +
    geom_jitter(size=0.1, width=0.05, height=0.05) +
    geom_vline(xintercept=c(-log2fc, log2fc), linetype=1, col="grey", linewidth=0.6) +
    geom_hline(yintercept=log10(3.162), linetype=2, col="grey", linewidth=0.6) +
    geom_text(x=min.x/2, y=max.y, label= paste0("down: ", as.character(down.count)), size = 2, color="#546de5", check_overlap=T) +
    geom_text(x=max.x/2, y=max.y, label= paste0("up: ", as.character(up.count)), size = 2, color="#ff4757", check_overlap=T) +
    geom_text_repel(inherit.aes=F, data=need_label,
        aes(x=log2FoldChange, y=log10(reads), label=gene_id), color=need_label$color, na.rm=T, segment.size = 0.3, segment.colour="black",
        box.padding = unit(0.35, "lines"), size=2, direction='both', nudge_x=nudge_x, nudge_y=nudge_y, min.segment.length = 0, max.overlaps=Inf) +
    # scale_x_continuous(breaks=c(-6, -3, -1, 0, 3, 1, 6)) +
    theme_bw(base_size=8) +
    theme(legend.position="right",
          panel.grid=element_blank(),
          legend.title = element_blank(),
          legend.text= element_text(color="black", size=8),
          plot.title = element_text(hjust=0.5, color="black", size=font_size),
          axis.text.x = element_text(color="black", size=font_size),
          axis.text.y = element_text( color="black", size=font_size),
          axis.title.x = element_text(color="black", size=font_size),
          axis.title.y = element_text(color="black", size=font_size))+
    labs(x="log2(fold change)",y="log10(Reads)",title=title)
}

}else {
names(data) <- c("gene_id", "log2FoldChange", "pvalue")
data$threshold <- as.factor(ifelse(data$pvalue < p_value & abs(data$log2FoldChange) > log2fc,
                            ifelse((data$log2FoldChange) > log2fc ,'up','down'),'not'))

if(!is.null(labels)){
need_label <- data[data$gene_id%in%labels,]
need_label$color = "black"
}else{
label_data <- data[!data$threshold=='not',]
label_data$pos <- 'middle'
label_data <- label_data[order(label_data$log2FoldChange), ]
label_data$rank <- 1:nrow(label_data)
dnrow <- nrow(label_data[label_data$threshold=='down', ])
if(dnrow > 5) down_cut <- 5 else down_cut <- dnrow
label_data[1:dnrow, 'pos'] <- 'down'
label_data[1:dnrow, 'color'] <- '#546de5'

unrow <- nrow(label_data[label_data$threshold=='up', ])
if(unrow > 5) up_cut <- 5 else up_cut <- unrow

label_data[(nrow(label_data)-up_cut+1):nrow(label_data), 'pos'] <- 'up'
label_data[(nrow(label_data)-up_cut+1):nrow(label_data), 'color'] <- '#ff4757'

need_label <- label_data[c(1:down_cut, (nrow(label_data)-up_cut+1):nrow(label_data)),]
}
if(args$rankplot){

rank_data <- data[abs(data$log2FoldChange) <= log2fc | data$pvalue < p_value,]
rank_data$pos <- 'not'
rank_data <- rank_data[order(rank_data$log2FoldChange), ]
rank_data$rank <- 1:nrow(rank_data)
dnrow <- nrow(rank_data[rank_data$threshold=='down', ])
if(dnrow > 5) down_cut <- 5 else down_cut <- dnrow
rank_data[1:dnrow, 'pos'] <- 'down'
rank_data[1:dnrow, 'color'] <- '#546de5'

unrow <- nrow(rank_data[rank_data$threshold=='up', ])
if(unrow > 5) up_cut <- 5 else up_cut <- unrow

rank_data[(nrow(rank_data)-unrow+1):nrow(rank_data), 'pos'] <- 'up'
rank_data[(nrow(rank_data)-unrow+1):nrow(rank_data), 'color'] <- '#ff4757'

if(!is.null(labels)){
need_rank <- rank_data[rank_data$gene_id%in%labels,]
need_rank$color = "black"
}else{
need_rank <- rank_data[c(1:down_cut, (nrow(rank_data)-up_cut+1):nrow(rank_data)),]
}

up.count <- nrow(rank_data[rank_data$threshold=="up",])
down.count <- nrow(rank_data[rank_data$threshold=="down",])
max.y <- max(rank_data$log2FoldChange)
min.y <- min(rank_data$log2FoldChange)
max.x <- max(rank_data$rank)

plot <- ggplot(data=rank_data, aes(x=rank, y=log2FoldChange,color=pos, fill=pos)) +
    scale_color_manual(values=c("down"="#546de5", "not"="#d2dae2", "up"="#ff4757")) +
    geom_point(size=0.1) +
    geom_hline(yintercept=c(-log2fc, log2fc), linetype=1, color="grey", linewidth=0.6) +
    geom_hline(yintercept=0, linetype=1, color="grey", linewidth=0.6) +
    geom_text(x=max.x*0.9, y=min.y, label= paste0("down: ", as.character(down.count)), size = 2, color="#546de5", check_overlap=T) +
    geom_text(x=max.x*0.1, y=max.y, label= paste0("up: ", as.character(up.count)), size = 2, color="#ff4757", check_overlap=T) +
    geom_text_repel(inherit.aes=F, data=need_rank,
        aes(x=rank, y=log2FoldChange, label=gene_id), color=need_rank$color, na.rm=T, segment.size = 0.3, segment.colour="black",
        box.padding = unit(0.35, "lines"), size=2, direction='both', nudge_x=nudge_x, nudge_y=nudge_y, min.segment.length = 0, max.overlaps=Inf) +
    # scale_y_continuous(breaks=c(-6, -3, -1, 0, 1, 3, 6)) +
    theme_bw(base_size=8) +
    theme(legend.position="right",
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(color="black", size=8),
          plot.title = element_text(hjust=0.5, color="black", size=font_size),
          axis.text.x = element_text(color="black", size=font_size),
          axis.text.y = element_text( color="black", size=font_size),
          axis.title.x = element_text(color="black", size=font_size),
          axis.title.y = element_text(color="black", size=font_size))+
    labs(x="rank",y="log2(fold change)",title=title)
}else{

up.count <- nrow(data[data$threshold=="up",])
down.count <- nrow(data[data$threshold=="down",])
max.y <- max(-log10(data$pvalue))
max.x <- max(data$log2FoldChange)
min.x <- min(data$log2FoldChange)

plot <- ggplot(data=data, 
               aes(x=log2FoldChange, y=-log10(pvalue), 
                   color=threshold, fill=threshold)) +
    scale_color_manual(values=c("down"="#546de5", "not"="#d2dae2", "up"="#ff4757")) +
    # geom_point(alpha=0.4, size=0.2) +
    geom_jitter(size=0.1, width=0.05, height=0.05) +
    geom_vline(xintercept=c(-log2fc, log2fc), linetype=2, color="grey", linewidth=0.6) +
    geom_hline(yintercept=-log10(p_value), linetype=2, color="grey", linewidth=0.6) +
    geom_text(x=min.x/2, y=max.y, label= paste0("down: ", as.character(down.count)), size = 2, color="#546de5", check_overlap=T) +
    geom_text(x=max.x/2, y=max.y, label= paste0("up: ", as.character(up.count)), size = 2, color="#ff4757", check_overlap=T) +
    geom_text_repel(inherit.aes=F, data=need_label,
        aes(x=log2FoldChange, y=-log10(pvalue), label=gene_id), color=need_label$color, na.rm=T, segment.size = 0.3, segment.colour="black",
        box.padding = unit(0.35, "lines"), size=2, direction='both', nudge_x=nudge_x, nudge_y=nudge_y, min.segment.length = 0, max.overlaps=Inf) +
    theme_bw(base_size=8) +
    theme(legend.position="right",
         panel.grid=element_blank(),
         legend.title = element_blank(),
         legend.text= element_text(color="black", size=8),
         plot.title = element_text(hjust=0.5, color="black", size=font_size),
         axis.text.x = element_text(color="black", size=font_size),
         axis.text.y = element_text(color="black", size=font_size),
         axis.title.x = element_text(color="black", size=font_size),
         axis.title.y = element_text(color="black", size=font_size)) +
    labs(x="log2(fold change)",y="-log10(pvalue)",title=title)
}
}

ggsave(output, width=10, height=7, unit="cm", dpi=300)
