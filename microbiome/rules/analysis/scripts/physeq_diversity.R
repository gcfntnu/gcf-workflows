

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep("--file=", initial.options)])
script.basename <- dirname(script.name)
source(file=file.path(script.basename, "physeq_common.R"))

suppressWarnings(library(ggplot2))
library(phyloseq)
library(stringr)

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

## prune empty
PF <- prune_taxa(taxa_sums(physeq) > 0, physeq)
R <- estimate_richness(PF)
colnames(R) <- str_replace_all(colnames(R), "\\.", "_")
R <- cbind(Sample_ID=sample_names(PF), R)
R <- cbind(R, sample_data(PF))

## diversity fig
if (!is.null(args$condition)){
    x <- make.names(args$condition)
} else{
    x <- NULL
}

if (!is.null(args$block)){
    color <- make.names(args$block)
} else{
    color <- NULL
}

p <- plot_richness(PF, x=x, color=color, measures=c("Chao1", "Shannon"))



if (is.null(args$subset)){
    prefix <- paste0(args$output_dir, "/")
} else{
    prefix <- paste0(file.path(args$output_dir, args$subset_keep), "_")
}

rich.fn <- paste0(prefix, "diversity.txt")
write.table(R, file=rich.fn, sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("Richness table: %s", rich.fn), "\n")
p.fn <- rich.fn <- paste0(prefix, "diversity.pdf")
ggplot2::ggsave(p.fn, plot=p)
cat(sprintf("Richness fig: %s", p.fn), "\n")




