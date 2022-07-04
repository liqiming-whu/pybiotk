#!/usr/bin/env Rscript
# Plot Fragment proportion in Peaks regions (FRiPs).
library(argparse)
library(ggplot2)
library(viridis)
library(ggpubr)
library(ggrepel)
parser <- ArgumentParser()
parser$add_argument('input', type="character", help='inputfile, summary.xls')
parser$add_argument('-o', dest='outfig', type="character", required=TRUE, help='output fig')
args <- parser$parse_args()


data <- read.table(args$input, sep="\t", header=TRUE)

data$FRiPs <- data$reads_in_peaks*100/data$filtered_reads
data$samples <- factor(data$samples, levels=data$samples)
fig1A <- ggplot(data, aes(x=group, y=peaks_count, fill=group))+
  geom_boxplot() +
  geom_jitter(aes(color = replicates), position = position_jitter(0.15))+
  geom_text_repel(aes(label=peaks_count, x=group, y=peaks_count), color="black", segment.color="grey", size=2, direction='both') +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Number of Peaks") +
  xlab("") +
  ggtitle("1A. Number of Peaks")

fig1B <- ggplot(data, aes(x=group, y=reads_in_peaks/2000000, fill=group))+
  geom_boxplot() +
  geom_jitter(aes(color = replicates), position = position_jitter(0.15))+
  geom_text_repel(aes(label=round(reads_in_peaks/2), x=group, y=reads_in_peaks/2000000), color="black", segment.color="grey", size=2, direction='both') +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragments in Peaks regions per Million") +
  xlab("") +
  ggtitle("1B. Fragments in Peaks regions")

fig1C <- ggplot(data, aes(x=group, y=FRiPs, fill=group)) + 
  geom_boxplot() +
  geom_jitter(aes(color = replicates), position = position_jitter(0.15))+
  geom_text_repel(aes(label=paste0(as.character(round(FRiPs, 1)), "%"), x=group, y=FRiPs), color="black", segment.color="grey", size=2, direction='both') +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("% of Fragments in Peaks") +
  xlab("") +
  ggtitle("1C. Fragment proportion in Peaks regions (FRiPs)")

fig <- ggarrange(fig1A, fig1B, fig1C, ncol = 3, common.legend = TRUE, legend="right", widths=c(1,1,1.1))

ggsave(args$outfig, fig, width=12, height=9, bg="white")
extension <- tail(strsplit(args$outfig, split=".", fixed=T)[[1]], n=1)
prefix <- sub(extension, "", args$outfig)
ggsave(paste0(prefix, "Peaks.", extension), fig1A, width=10, height=7, bg="white")
ggsave(paste0(prefix, "FiPs.", extension), fig1B, width=10, height=7, bg="white")
ggsave(paste0(prefix, "FRiPs.", extension), fig1C, width=10, height=7, bg="white")
