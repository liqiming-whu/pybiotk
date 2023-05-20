#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument('input', type="character", help='inputfile')
parser$add_argument('-c', dest='columns', type="integer", nargs="+", help='use columns to plot, order:gene, basemean, log2fc, padj')
parser$add_argument('-o', dest='output', type="character", required=TRUE, help='output')
parser$add_argument('-t', dest='title', type="character", default="", help='title')
parser$add_argument('-l', dest='log2fc', type="double", default=1.0, help='cutoff of log2fc')
parser$add_argument('-p', dest='padj_value', type="double", default=0.05, help='cutoff of padj_value, or pvalue')
parser$add_argument('-f', dest='font_size', type="integer", default=8, help='font size')
parser$add_argument('--nudge_x', dest='nudge_x', type="double", default=0, help='nudge_x')
parser$add_argument('--nudge_y', dest='nudge_y', type="double", default=0, help='nudge_y')
parser$add_argument('--box_padding', dest='box_padding', type="double", default=0.35, help='box_padding')
parser$add_argument('--labels', dest='labels', type="character", nargs="+", default=NULL, help="labels list")
parser$add_argument('--plot-type', dest='plot_type', type="character", default="volcano", help='plot type, volcano, ma, rank')


args <- parser$parse_args()

data_table <- args$input
columns <- args$columns
output <- args$output
title <- args$title
log2fc <- args$log2fc
padj_value <- args$padj_value
font_size <- args$font_size
nudge_x <- args$nudge_x
nudge_y <- args$nudge_y
box_padding <- args$box_padding
labels <- args$labels
plot_type <- args$plot_type


library(ggplot2)
library(ggrepel)
library(ggpubr)

.lwd <- ggplot2::.pt / ggplot2::.stroke
gg_lwd_convert <- function(value, unit = "pt") {

  # convert input type to mm
  value_out <- grid::convertUnit(grid::unit(value, unit), "mm", valueOnly = TRUE)

  # return with conversion factor
  return(
    value_out / .lwd
  )
}

volcano.plot <- function(data, path, title="", labels=NULL, up_label=5, down_label=5) {
  names(data) <- c("gene_id", "baseMean", "log2FoldChange", "padj")
  data <- na.omit(data)
  data$threshold <- as.factor(ifelse(data$padj <= padj_value & abs(data$log2FoldChange) >= log2fc,
                                     ifelse((data$log2FoldChange) >= log2fc ,'up','down'),'not'))

  up.count <- nrow(data[data$threshold=="up",])
  down.count <- nrow(data[data$threshold=="down",])

  data$label <- "NS"
  if(up.count > 0) data[data$threshold=="up",]$label <- paste0("Up: ", up.count)
  if(down.count > 0) data[data$threshold=="down",]$label <- paste0("Down: ", down.count)

  data$label <- factor(data$label, levels = c(paste0("Up: ", up.count), paste0("Down: ", down.count), "NS"))

  if(nrow(data[data$padj==0,]) > 0) data[data$padj==0,]$padj <- .Machine$double.xmin
  max.y <- max(-log10(data$padj))
  max.x <- max(data$log2FoldChange)
  min.x <- min(data$log2FoldChange)

  if(!is.null(labels)){
    need_label <- data[data$gene_id%in%labels,]
    need_label$color = "black"
  }else{
    label_data <- data[!data$threshold=='not',]
    if(nrow(label_data) > 0) label_data$pos <- 'middle'
    label_data <- label_data[order(label_data$log2FoldChange), ]
    dnrow <- nrow(label_data[label_data$threshold=='down', ])
    if(dnrow > down_label) down_cut <- down_label else down_cut <- dnrow
    label_data[1:dnrow, 'pos'] <- 'down'
    label_data[1:dnrow, 'color'] <- '#546de5'

    unrow <- nrow(label_data[label_data$threshold=='up', ])
    if(unrow > up_label) up_cut <- up_label else up_cut <- unrow

    label_data[(nrow(label_data)-up_cut+1):nrow(label_data), 'pos'] <- 'up'
    label_data[(nrow(label_data)-up_cut+1):nrow(label_data), 'color'] <- '#ff4757'

    need_label <- label_data[c(1:down_cut, (nrow(label_data)-up_cut+1):nrow(label_data)),]
    need_label <- need_label[nchar(need_label$gene_id) < 10,]
  }

  max.overlaps = getOption("ggrepel.max.overlaps", default = Inf)

  cols = setNames(c("#B31B21", "#1465AC", "darkgray"), c(paste0("Up: ", up.count), paste0("Down: ", down.count), "NS"))
  p <- ggplot(data=data, aes(x=log2FoldChange, y=-log10(padj), color=label, fill=label)) +
    scale_color_manual(values=cols) +
    geom_point(size=0.4, alpha=1) +
    geom_vline(xintercept=c(-log2fc, log2fc), linetype=2, color="black", linewidth=gg_lwd_convert(1)) +
    geom_hline(yintercept=-log10(padj_value), linetype=2, color="black", linewidth=gg_lwd_convert(1)) +
    # geom_text(x=0.5*min.x, y=1.02*max.y, label= paste0("down: ", as.character(down.count)), size = font_size/3, color="black", check_overlap=T) +
    # geom_text(x=0.5*max.x, y=1.02*max.y, label= paste0("up: ", as.character(up.count)), size = font_size/3, color="black", check_overlap=T) +
    geom_text_repel(data=need_label, aes(x=log2FoldChange, y=-log10(padj), label=gene_id), color="black", na.rm=T, force = 1, seed=42, segment.colour="black",
                    box.padding = unit(box_padding, "lines"), point.padding = unit(0.3, "lines"), size=font_size/3, direction="both", nudge_x=nudge_x, nudge_y=nudge_y,
                    min.segment.length = 0, max.overlaps=max.overlaps) +
    theme_bw(base_size=font_size) +
    theme(
      panel.grid=element_blank(),
      panel.border = element_rect(color="black", linewidth=gg_lwd_convert(1)),
      axis.line = element_blank(),
      legend.position="right",
      legend.title = element_blank(),
      legend.text= element_text(color="black", size=font_size),
      plot.title = element_text(hjust=0.5, color="black", size=font_size),
      axis.text = element_text(color="black", size=font_size),
      axis.title = element_text(color="black", size=font_size),
    ) +
    labs(x="Log2 fold change",y="-Log10 p-adjust",title=title)

  ggsave(path, p, width=12, height=9, unit="cm", dpi=300)

}


.levels <- function(x){
  if(!is.factor(x)) x <- as.factor(x)
  levels(x)
}

.parse_font <- function(font){
  if(is.null(font)) res <- NULL
  else if(inherits(font, "list")) res <- font
  else{
    # matching size and face
    size <- grep("^[0-9]+$", font, perl = TRUE)
    face <- grep("plain|bold|italic|bold.italic", font, perl = TRUE)
    if(length(size) == 0) size <- NULL else size <- as.numeric(font[size])
    if(length(face) == 0) face <- NULL else face <- font[face]
    color <- setdiff(font, c(size, face))
    if(length(color) == 0) color <- NULL
    res <- list(size=size, face = face, color = color)
  }
  res
}

ggmaplot <- function (data, fdr = 0.05, log2fc = 1.5, genenames = NULL,
                      detection_call = NULL, size = NULL, alpha = 1,
                      seed = 42,
                      font.label = c(12, "plain", "black"), label.rectangle = FALSE,
                      palette = c("#B31B21", "#1465AC", "darkgray"),
                      top = 15, select.top.method = c("padj", "fc"),
                      label.select = NULL,
                      main = NULL, xlab = "Log2 mean expression",  ylab = "Log2 fold change",
                      ggtheme = theme_classic(),...)
{

  if(!base::inherits(data, c("matrix", "data.frame", "DataFrame", "DE_Results", "DESeqResults")))
    stop("data must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  if(!is.null(detection_call)){
    if(nrow(data)!=length(detection_call))
      stop("detection_call must be a numeric vector of length = nrow(data)")
  }
  else if("detection_call" %in% colnames(data)){
    detection_call <- as.vector(data$detection_call)
  }
  else detection_call = rep(1, nrow(data))

  # Legend position
  if(is.null(list(...)$legend)) legend <- c(0.12, 0.9)
  # If basemean logged, we'll leave it as is, otherwise log2 transform
  is.basemean.logged <- "baseMeanLog2" %in% colnames(data)
  if(is.basemean.logged){
    data$baseMean <- data$baseMeanLog2
  }
  else if("baseMean" %in% colnames(data)){
    data$baseMean <- log2(data$baseMean +1)
  }

  # Check data format
  ss <- base::setdiff(c("baseMean", "log2FoldChange", "padj"), colnames(data))
  if(length(ss)>0) stop("The colnames of data must contain: ",
                        paste(ss, collapse = ", "))

  if(is.null(genenames)) genenames <- rownames(data)
  else if(length(genenames)!=nrow(data))
    stop("genenames should be of length nrow(data).")

  sig <- rep(3, nrow(data))
  sig[which(data$padj <= fdr & data$log2FoldChange < 0 & abs(data$log2FoldChange) >= log2fc & detection_call ==1)] = 2
  sig[which(data$padj <= fdr & data$log2FoldChange > 0 & abs(data$log2FoldChange) >= log2fc & detection_call ==1)] = 1
  data <- data.frame(name = genenames, mean = data$baseMean, lfc = data$log2FoldChange,
                     padj = data$padj, sig = sig)

  # Change level labels
  . <- NULL
  data$sig <- as.factor(data$sig)
  .lev <- .levels(data$sig) %>% as.numeric()
  palette <- palette[.lev]
  new.levels <- c(
    paste0("Up: ", sum(sig == 1)),
    paste0("Down: ", sum(sig == 2)),
    "NS"
  ) %>% .[.lev]

  data$sig <- factor(data$sig, labels = new.levels)

  # Ordering for selecting top gene
  select.top.method <- match.arg(select.top.method)
  if(select.top.method == "padj") data <- data[order(data$padj), ]
  else if(select.top.method == "fc") data <- data[order(abs(data$lfc), decreasing = TRUE), ]
  # select data for top genes
  complete_data <- stats::na.omit(data)
  labs_data <- subset(complete_data, padj <= fdr & name!="" & abs(lfc) >= log2fc)
  labs_data <- utils::head(labs_data, top)
  # Select some specific labels to show
  if(!is.null(label.select)){
    selected_labels  <- complete_data %>%
      subset(complete_data$name  %in% label.select, drop = FALSE)
    labs_data <- dplyr::bind_rows(labs_data, selected_labels) %>%
      dplyr::distinct(.data$name, .keep_all = TRUE)
  }

  font.label <- .parse_font(font.label)
  font.label$size <- ifelse(is.null(font.label$size), 12, font.label$size)
  font.label$color <- ifelse(is.null(font.label$color), "black", font.label$color)
  font.label$face <- ifelse(is.null(font.label$face), "plain", font.label$face)

  # Plot
  mean <- lfc <- sig <- name <- padj <-  NULL
  p <- ggplot(data, aes(x = mean, y = lfc)) +
    geom_point(aes(color = sig), size = size, alpha = alpha)

  max.overlaps = getOption("ggrepel.max.overlaps", default = Inf)

  if(label.rectangle){
    p <- p + ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = name),
                                       box.padding = unit(box_padding, "lines"),
                                       point.padding = unit(0.3, "lines"),
                                       force = 1, seed = seed, fontface = font.label$face,
                                       size = font.label$size/3, color = font.label$color,
                                       direction="both", nudge_x=nudge_x, nudge_y=nudge_y,
                                       max.overlaps = max.overlaps)
  }
  else{
    p <- p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = name),
                                      box.padding = unit(box_padding, "lines"),
                                      point.padding = unit(0.3, "lines"),
                                      force = 1, seed = seed, fontface = font.label$face,
                                      size = font.label$size/3, color = font.label$color,
                                      direction="both", nudge_x=nudge_x, nudge_y=nudge_y,
                                      max.overlaps = max.overlaps)
  }

  p <- p + scale_x_continuous(breaks=seq(0, max(data$mean), 2))+
    labs(x = xlab, y = ylab, title = main, color = "")+ # to remove legend title use color = ""
    geom_hline(yintercept = c(0, -log2fc, log2fc), linetype = c(1, 2, 2),
               color = c("black", "black", "black"), linewidth=gg_lwd_convert(1))

  p <- ggpar(p, palette = palette, ggtheme = ggtheme, ...)
  p
}


ma.plot <-function(data, path, title="", labels=NULL) {
  names(data) <- c("gene_id", "baseMean", "log2FoldChange", "padj")
  data <- na.omit(data)
  p <- ggmaplot(
    data,
    fdr = padj_value,
    log2fc = log2fc,
    genenames = as.vector(data$gene_id),
    size = 0.4,
    alpha = 1,
    seed = 42,
    font.label = c(font_size, "plain", "black"),
    label.rectangle = FALSE,
    palette = c("#B31B21", "#1465AC", "darkgray"),
    select.top.method = "padj",
    top = 10,
    label.select = labels,
    main = title,
    xlab = "Log2 mean expression",
    ylab = "Log2 fold change",
    ggtheme = theme_bw(base_size = font_size)+
      theme(
        panel.grid=element_blank(),
        panel.border = element_rect(color="black", linewidth=gg_lwd_convert(1)),
        axis.line = element_blank(),
        legend.position="right",
        legend.title = element_blank(),
        legend.text= element_text(color="black", size=font_size),
        plot.title = element_text(hjust=0.5, color="black", size=font_size),
        axis.text.x = element_text(color="black", size=font_size),
        axis.text.y = element_text(color="black", size=font_size),
        axis.title.x = element_text(color="black", size=font_size),
        axis.title.y = element_text(color="black", size=font_size),
      )
  )
  ggsave(path, p, width=12, height=9, unit="cm", dpi=300)
}


rank.plot <- function(data, path, title="", labels=NULL) {
  names(data) <- c("gene_id", "baseMean", "log2FoldChange", "padj")
  data$threshold <- as.factor(ifelse(data$padj <= padj_value & abs(data$log2FoldChange) >= log2fc,
                                     ifelse((data$log2FoldChange) >= log2fc ,'up','down'),'not'))
  up.count <- nrow(data[data$threshold=="up",])
  down.count <- nrow(data[data$threshold=="down",])
  rank_data <- data[data$padj < padj_value,]
  rank_data$pos <- 'NS'
  rank_data <- rank_data[order(rank_data$log2FoldChange), ]
  rank_data$rank <- 1:nrow(rank_data)
  dnrow <- nrow(rank_data[rank_data$threshold=='down', ])
  marks <- c()
  if(dnrow > 5) down_cut <- 5 else down_cut <- dnrow
  if(dnrow > 0) {
    rank_data[1:dnrow, 'pos'] <- paste0("Down: ", dnrow)
    marks <- append(marks, rank_data[1:down_cut,]$gene_id)
  }

  unrow <- nrow(rank_data[rank_data$threshold=='up', ])
  if(unrow > 5) up_cut <- 5 else up_cut <- unrow
  if(unrow > 0) {
    rank_data[(nrow(rank_data)-unrow+1):nrow(rank_data), 'pos'] <- paste0("UP: ", unrow)
    marks <- append(marks, rank_data[(nrow(rank_data)-up_cut+1):nrow(rank_data),]$gene_id)
  }

  rank_data$pos <- factor(rank_data$pos, levels = c(paste0("UP: ", unrow), paste0("Down: ", dnrow), "NS"))

  if(!is.null(labels)){
    need_rank <- rank_data[rank_data$gene_id%in%labels,]
   }else{
    need_rank <- rank_data[rank_data$gene_id%in%marks,]
   }

  max.overlaps = getOption("ggrepel.max.overlaps", default = Inf)

  cols = setNames(c("#B31B21", "#1465AC", "darkgray"), c(paste0("Up: ", up.count), paste0("Down: ", down.count), "NS"))
  plot <- ggplot(data=rank_data, aes(x=rank, y=log2FoldChange, fill=pos)) +
        scale_color_manual(values=cols) +
        geom_point(size=0.4, alpha=1) +
        geom_hline(yintercept=c(-log2fc, log2fc), linetype=1, color="grey", linewidth=0.6) +
        geom_hline(yintercept=0, linetype=1, color="grey", linewidth=0.6) +
        geom_text_repel(inherit.aes=F, data=need_rank,
            aes(x=rank, y=log2FoldChange, label=gene_id), color="black", na.rm=T, force = 1, seed=42, segment.colour="black",
                    box.padding = unit(box_padding, "lines"), point.padding = unit(0.3, "lines"), size=font_size/3, direction="both", nudge_x=nudge_x, nudge_y=nudge_y,
                    min.segment.length = 0, max.overlaps=max.overlaps) +
        theme_bw(base_size=font_size) +
        theme(panel.grid=element_blank(),
              panel.border = element_rect(color="black", linewidth=gg_lwd_convert(1)),
              axis.line = element_blank(),
              legend.position="right",
              legend.title = element_blank(),
              legend.text= element_text(color="black", size=font_size),
              plot.title = element_text(hjust=0.5, color="black", size=font_size),
              axis.text = element_text(color="black", size=font_size),
              axis.title = element_text(color="black", size=font_size),
            ) +
        labs(x="rank",y="Log2 fold change", title=title)

}


data <- read.table(file=data_table, header=T)
if(is.null(columns)) data <- data[,c("Row.names", "baseMean", "log2FoldChange", "padj")] else data <- data[, columns]
data <- na.omit(data)
columns <- colnames(data)
if((ncol(data) < 4) & (!"padj" %in% colnames(data))) data$padj <- 0
data <- data[,append(columns, "padj")]

if(plot_type == "volcano") {
    volcano.plot(data, output, title, labels)
}else if(plot_type == "ma") {
    ma.plot(data, output, title, labels)
}else if(plot_type == "rank") {
    rank.plot(data, output, title, labels)
}else volcano.plot(data, output, title, labels)
