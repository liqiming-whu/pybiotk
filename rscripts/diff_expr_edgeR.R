#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument('input', type="character", nargs='+', help='control and case count tables')
parser$add_argument('-s', dest='samples', type="character", nargs='+', default=NULL, help='sample names.')
parser$add_argument('-g', dest='group', type="character", nargs='+', default=NULL, help='sample groups.')
parser$add_argument('-b', dest='by', type="integer", default=1, help='use columns to merge.')
parser$add_argument('-c', dest='column', type="integer", default=2, help='use column to calculate.')
parser$add_argument('-o', dest='outdir', type="character", required=TRUE, help='output dir')
parser$add_argument('-l', dest='log2fc', type="double", default=1.0, help='log2fc [default=1].')
parser$add_argument("-p", dest="padj_value", type="double", default=0.05, help="cutoff of padj_value [default=0.05].")
parser$add_argument('--bcv', dest='bcv', type="double", default=0.1, help='bcv value. [default=0.1].')
parser$add_argument('--housekeeping', dest='housekeeping', type="character", default=NULL, help='housekeeping genes for estimate dispersions.')
parser$add_argument('--spkie-in', dest="spike_in", type="character", nargs='+', default=NULL, help='control and case spike_in count tables.')

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

library(edgeR)

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

# 构建DGEList对象
dgelist <- DGEList(counts = table[,-1], genes = table[,1], group=group)

# 过滤低表达基因，两种方式

# 1.直接依据选定的count值过滤
# keep <- rowSums(dgelist$counts) >= 50

# 2.CPM标准化
# keep <- rowSums(cpm(dgelist)>1) >= 1

# 3.自动过滤
keep<- filterByExpr(deglist, group=group)

dgelist <- dgelist[keep, , keep.lib.sizes=FALSE]

if(!is.null(spike_in)) {

    spike_table <- read_data(spike_in[1], samples_name[1])
    for (idx in seq_along(spike_in)[-1]) {
        filename <- spike_in[idx]
        sample_name <- samples_name[idx]
        data <- read_data(filename, sample_name)
        spike_table <- merge(spike_table, data, by="id", all=T)
    }

    spike_table[is.na(spike_table)] <- 0
    spike_dgelist <- DGEList(counts = spike_table[,-1], genes = spike_table[,1], group=group)
    spike_dgelist <- calcNormFactors(dgelist, method='TMM')

    deglist$samples$norm.factors <- spike_dgelist$samples$norm.factors
}else{
# 计算标准化因子，TMM
dgelist <- calcNormFactors(dgelist, method='TMM')
}

if(is.null(args$housekeeping)) {
    if (args$bcv){
        # 1.自定义离散值
        bcv <- args$bcv
        et <- exactTest(dgelist, dispersion = bcv ^ 2)
        genes <- decideTestsDGE(et, p.value=0.05, lfc=0)
        summary(genes)
        result <- cbind(dgelist$genes,dgelist$counts,et$table)
    }else{
        # 2.估算离散值(有重复)
        design <- model.matrix(~group)
        deglist <- estimateDisp(dgelist, design, robust=TRUE)
        fit <- glmFit(dgelist, design, robust=TRUE)
        lrt <- glmLRT(fit)
        genes <- decideTestsDGE(lrt, p.value=0.05, lfc = 0)
        summary(genes)
        result <- cbind(dgelist$genes,dgelist$counts,lrt$table)
    }
}else {
# 3.根据已知一些不会发生改变的基因推测离散值
    housekeeping_list <- read.table(args$housekeeping, header=F)
    housekeeping_list <- housekeeping_list[,1]
    genes <- dgelist$genes[, 1]
    housekeeping <- which(genes %in% housekeeping_list)
    dge <- dgelist
    dge$samples$group <- 1
    rownames(dge$counts) <- dgelist$genes[, 1]
    dge <- estimateDisp(dge[housekeeping,], trend="none", tagwise=FALSE, robust=TRUE)
    dgelist$common.dispersion <- dge$common.dispersion
    design <- model.matrix(~group)
    fit <- glmFit(dgelist, design, robust=TRUE)
    lrt <- glmLRT(fit)
    genes <- decideTestsDGE(lrt, p.value=0.05,lfc = 0)
    summary(genes)
    result <- cbind(dgelist$genes,dgelist$counts,lrt$table)
}

if(!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

setwd(outdir)

result$Count <- rowSums(result[, 2:(2+input_files_len-1)])
up_diff <- subset(result, logFC > log2fc)
down_diff <- subset(result, logFC < -log2fc)
sig_result <- subset(result, abs(logFC) > log2fc)

write.table(up_diff, "up.xls", quote=F, row.names=F, col.names=T, sep="\t")
write.table(down_diff, "down.xls", quote=F, row.names=F, col.names=T, sep="\t")
write.table(sig_result, "sig.xls", quote=F, row.names=F, col.names=T, sep="\t")
write.table(result, "all.xls", quote=F, row.names=F, col.names=T, sep="\t")

