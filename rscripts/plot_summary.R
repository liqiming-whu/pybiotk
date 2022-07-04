#!/usr/bin/env Rscript
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
data$samples <- factor(data$samples, levels=data$samples)

trim <- data.frame("samples"=rep(data$samples,2),
                   "withadapter"=c(data$read1_with_adapter, data$read2_with_adapter),
                   "withadapter_percent"=c(data$read1_with_adapter_percent, data$read2_with_adapter_percent),
                   "termini"=c(rep("read1", nrow(data)), rep("read2", nrow(data))))
                   
fig1A <- ggplot(data) + 
  geom_col(aes(x=samples, y=total_read_pairs/1000000, fill=samples)) +
  geom_jitter(data=trim, aes(x=samples, y=withadapter/1000000, color = termini), position = position_jitter(0.15))+
  # guides(color = guide_legend(override.aes = list(size=2))) + 
  geom_text_repel(aes(label=total_read_pairs, x=samples, y=total_read_pairs/1000000), color="black", segment.color="grey", size=2, direction='y') +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank(), legend.position = "bottom") +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Read pairs per Million") +
  xlab("") +
  ggtitle("1A. Rawdata Read Pairs Count")

fig1B <- ggplot(trim, aes(x=samples, y=withadapter/1000000, fill=samples)) + 
  geom_boxplot() +
  geom_jitter(aes(color = termini), position = position_jitter(0.15))+
  geom_text_repel(aes(label=withadapter, x=samples, y=withadapter/1000000), color="black", segment.color="grey", size=2, direction='both') +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Reads with Adapter per Million") +
  xlab("") +
  ggtitle("1B. Reads with Adapter")

fig1C <- ggplot(trim, aes(x=samples, y=as.numeric(sub("%", "", withadapter_percent)), fill=samples)) + 
  geom_boxplot() +
  geom_jitter(aes(color = termini), position = position_jitter(0.15))+
  geom_text_repel(aes(label=withadapter_percent, x=samples, y=as.numeric(sub("%", "", withadapter_percent))), color="black", segment.color="grey", size=2, direction='both') +
  # geom_text(aes(label=withadapter_percent, x=samples, y=1.02*as.numeric(sub("%", "", withadapter_percent))), color="black")+
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("% of Reads with Adapter") +
  xlab("") +
  ggtitle("1C. Reads with Adapter(%)")

fig1 <- ggarrange(fig1A, fig1B, fig1C, ncol = 3, common.legend = TRUE, legend="bottom")

alignres <- data.frame("group"=data$group, "replicate"=data$replicates,
                       "sequencedepth"=data$input_read_pairs,
                       "mappedfragments"=data$mapped_read_pairs,
                       "alignmentrate"=data$alignment_rate)

fig2A <- ggplot(alignres, aes(x=group, y=sequencedepth/1000000, fill=group))+
  geom_boxplot() +
  geom_jitter(aes(color = replicate), position = position_jitter(0.15))+
  geom_text_repel(aes(label=sequencedepth, x=group, y=sequencedepth/1000000), color="black", segment.color="grey", size=2, direction='both') +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Sequencing Depth per Million") +
  xlab("") +
  ggtitle("2A. Sequencing Depth")


fig2B <- ggplot(alignres, aes(x=group, y=mappedfragments/1000000, fill=group))+
  geom_boxplot() +
  geom_jitter(aes(color = replicate), position = position_jitter(0.15))+
  geom_text_repel(aes(label=mappedfragments, x=group, y=mappedfragments/1000000), color="black", segment.color="grey", size=2, direction='both') +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Mapped Fragments per Million") +
  xlab("") +
  ggtitle("2B. Alignable Fragments")


fig2C <- ggplot(alignres, aes(x=group, y=as.numeric(sub("%", "", alignmentrate)), fill=group))+
  geom_boxplot() +
  geom_jitter(aes(color = replicate), position = position_jitter(0.15))+
  geom_text_repel(aes(label=alignmentrate, x=group, y=as.numeric(sub("%", "", alignmentrate))), color="black", segment.color="grey", size=2, direction='both') +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("2C. Alignment Rate")

fig2 <- ggarrange(fig2A, fig2B, fig2C, ncol = 3, common.legend = TRUE, legend="bottom")


dupdata <- data.frame("group"=data$group, "replicate"=data$replicates,
                      "duprate"=data$duplication_rate,
                      "uniquefragments"=data$mapped_read_pairs * (100 - as.numeric(sub("%", "", data$duplication_rate)))/100,
                      "libarysize"=data$filtered_reads)


fig3A <- ggplot(dupdata, aes(x=group, y=as.numeric(sub("%", "", duprate)), fill=group))+
  geom_boxplot() +
  geom_jitter(aes(color = replicate), position = position_jitter(0.15))+
  geom_text_repel(aes(label=duprate, x=group, y=as.numeric(sub("%", "", duprate))), color="black", segment.color="grey", size=2, direction='both') +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Duplication Rate (*100%)") +
  xlab("") +
  ggtitle("3A. Duplication Rate")


fig3B <- ggplot(dupdata, aes(x=group, y=uniquefragments/1000000, fill=group))+
  geom_boxplot() +
  geom_jitter(aes(color = replicate), position = position_jitter(0.15))+
  geom_text_repel(aes(label=round(uniquefragments), x=group, y=uniquefragments/1000000), color="black", segment.color="grey", size=2, direction='both') +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Unique Fragments per Million") +
  xlab("") +
  ggtitle("3B. Unique Fragments")

fig3C <- ggplot(dupdata, aes(x=group, y=libarysize/1000000, fill=group))+
  geom_boxplot() +
  geom_jitter(aes(color = replicate), position = position_jitter(0.15))+
  geom_text_repel(aes(label=libarysize, x=group, y=libarysize/1000000), color="black", segment.color="grey", size=2, direction='both') +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.85, option = "turbo", alpha = 0.8) +
  # scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Reads per Million") +
  xlab("") +
  ggtitle("3C. Final Library Size")


fig3 <- ggarrange(fig3A, fig3B, fig3C, ncol = 3, common.legend = TRUE, legend="bottom")

fig <- ggarrange(fig1, fig2, fig3, ncol=1)

ggsave(args$outfig, fig, width=12, height=12, bg="white")
extension <- tail(strsplit(args$outfig, split=".", fixed=T)[[1]], n=1)
prefix <- sub(extension, "", args$outfig)
ggsave(paste0(prefix, "fig1.", extension), fig1, width=10, height=7, bg="white")
ggsave(paste0(prefix, "fig2.", extension), fig2, width=10, height=7, bg="white")
ggsave(paste0(prefix, "fig3.", extension), fig3, width=10, height=7, bg="white")
