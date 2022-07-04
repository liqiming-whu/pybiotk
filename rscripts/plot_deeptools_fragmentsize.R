#!/usr/bin/env Rscript
library(argparse)
library(ggplot2)
library(viridis)
library(ggpubr)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument('input', type="character", help='input file.')
parser$add_argument('-o', dest='outfig', type="character", required=TRUE, help='output fig')
parser$add_argument('-s', dest='samples', type="character", nargs='+', required=TRUE, help='sample names.')
parser$add_argument('-g', dest='group', type="character", nargs='+', required=TRUE, help='sample groups.')
parser$add_argument('-r', dest='replicate', type="character", nargs='+', required=TRUE, help='sample replicate.')
args <- parser$parse_args()

data <- read.table(args$input, header=TRUE)

samples <- args$samples
group <- args$group
replicate <- args$replicate

data_info <- data.frame("Sample"=samples, "group"=group, "replicate"=replicate)

data <- merge(data, data_info, by="Sample")
data$weight <- data$Occurrences/sum(data$Occurrences)
data$Sample <- factor(data$Sample, levels=samples)
data$group <- factor(data$group, levels=unique(group))

max_col <- data[which(data$Occurrences==max(data$Occurrences),arr.ind=TRUE),]

data2 <- data %>% group_by(Size) %>% summarise("Occ"=sum(Occurrences))
max_Occ <- data2[which(data2$Occ==max(data2$Occ), arr.ind=TRUE),]

fig1A <- ggplot(data, aes(x = Sample, y = Size, weight = weight, fill = group)) +
  geom_violin(bw=5) +
  scale_y_continuous(breaks = seq(0, 700, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("") +
  ggtitle("1A. Fragment Length Distribution(violin plot)")

fig1B <- ggplot(data, aes(x = Size, y = Occurrences, group=Sample, color = Sample)) +
  geom_line(size = 1) +
  geom_text(label=max_col$Size, x=max_col$Size+30, y=max_col$Occurrences, color="black", size=2, check_overlap=T) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "turbo") +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(50, 700))+
  ggtitle("1B. Fragment Length Distribution(line plot)")

fig1C <- ggplot(data, aes(x = Sample, y = Size, weight = weight, fill = group)) +
  geom_boxplot(outlier.shape=NA) +
  scale_y_continuous(breaks = seq(0, 700, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 11) +
  theme(panel.grid=element_blank()) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("") +
  ggtitle("1C. Fragment Length Distribution(box plot)")

fig1D <- ggplot(data, aes(x = Size, y = Occurrences, group=Sample, fill = Sample)) +
  geom_col(position = "stack") +
  geom_text(label=max_Occ$Size, x=max_Occ$Size+30, y=max_Occ$Occ, color="black", size=2, check_overlap=T) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "turbo") +
  theme_bw(base_size = 11) +
  theme(panel.grid=element_blank()) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(50, 700))+
  ggtitle("1D. Fragment Length Distribution(bar plot)")

fig <- ggarrange(fig1A, fig1B, fig1C, fig1D, ncol = 2, nrow=2)

ggsave(args$outfig, fig, width=12, height=9, bg="white")
extension <- tail(strsplit(args$outfig, split=".", fixed=T)[[1]], n=1)
prefix <- sub(extension, "", args$outfig)
ggsave(paste0(prefix, "violin.", extension), fig1A, width=10, height=7, bg="white")
ggsave(paste0(prefix, "line.", extension), fig1B, width=10, height=7, bg="white")
ggsave(paste0(prefix, "box.", extension), fig1C, width=10, height=7, bg="white")
ggsave(paste0(prefix, "bar.", extension), fig1D, width=10, height=7, bg="white")
