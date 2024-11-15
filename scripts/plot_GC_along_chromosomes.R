#!/usr/bin/env Rscript

require('argparse')

parser <- ArgumentParser()
parser$add_argument("-i", "--coverage_filename", type = 'character',
                    help="input bed file;")
parser$add_argument("-o", "--output_file",
                    help="filename of the output plot")
parser$add_argument("-t", "--title",
                    help="name printed at the top of the plot [default is the filename]")
parser$add_argument("-c", "--chr_size_filter", type = "double",
                    help="Plot chromosomes only", default = 0)

args <- parser$parse_args()

# coverage_filename <- 'data/GCA_910594885.2.GC.1k.bedGraph'
# chromosome_size_threshold <- 0.5
# title <- 'GCA_910594885'
coverage_filename <- args$coverage_filename
output_fileame <- args$output_file
chromosome_size_threshold <- args$chr_size_filter
title <- args$title


GC_tab <- read.table(coverage_filename)

chomoromes <- unique(GC_tab[, 1])
chomorome_sizes <- sapply(chomoromes, function(chr){ max(GC_tab[GC_tab[, 1] == chr, 3]) })

retained_chromosomes <- chomorome_sizes > (chromosome_size_threshold * 1e6)

chomoromes <- chomoromes[retained_chromosomes]
chomorome_sizes <- chomorome_sizes[retained_chromosomes]
GC_tab <- GC_tab[GC_tab[, 1] %in% chomoromes, ]

adjustment <- c(0, cumsum(chomorome_sizes[1:(length(chomorome_sizes) - 1)]))
names(adjustment) <- chomoromes

GC_tab[, 2] <- GC_tab[, 2] + adjustment[GC_tab[, 1]]
GC_tab[, 3] <- GC_tab[, 3] + adjustment[GC_tab[, 1]]

cols <- rep(c('black', 'gray'), ceiling(length(chomoromes) / 2))
names(cols) <- chomoromes
GC_tab[, 'col'] <- cols[GC_tab[, 1]]

output_suffix <- gsub(".*\\.","", output_fileame)

if (output_suffix == 'pdf'){
    pdf(output_fileame)
} else {
    png(output_fileame, width = 1040, height = 800, units = "px")
}

    plot(GC_tab[, 4] ~ GC_tab[, 2], col = GC_tab[, 'col'], cex = 0.2, xlab = 'Position in a genome', ylab = 'GC content', main = title)
dev.off()
