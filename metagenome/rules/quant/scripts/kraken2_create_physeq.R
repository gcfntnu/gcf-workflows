#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(library(phyloseq))

biom.fn <- args[1]
tree.fn <- args[2]
out.fn <- args[3]

OTU <- import_biom(biom.fn, tree.fn, parseFunction=parse_taxonomy_default)
TAX <- tax_table(OTU)
colnames(TAX) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_table(OTU) <- TAX

print(OTU)
print(head(TAX))
saveRDS(OTU, out.fn)
