library(magrittr)
library(phyloseq)
library(breakaway)
library(DivNet)

source("physeq_common.R")


tax = tax_glom(physeq, taxrank="Phylum")

res <- divnet(tax, X=arg.condition, ncores = 4)

test <- testDiversity(res)
