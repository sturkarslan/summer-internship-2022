library(tidyverse)
library(here)
library(DT)
library(gridExtra)
library(jsonlite)
library(plyr)
library(purrr)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ReactomePA)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(data.table)
source("code/regulon_enrichment.R")

############### function  for pathway enrichment ###############
drug_pathway_enrichment <- function(input.df, regulons) {
  # drug pathway enrichment
  drugs_all_regulon_pathways <- data.frame()

  total.targets <- length(unique(input.df$`Target Gene(s)`))
  count = 0

  for (i in unique(input.df$`Target Gene(s)`)) {
    count = count + 1
    cat("\nProcessing ",
        i,
        " ",
        count,
        "/",
        total.targets,
        "\n")

    # select the drug row
    my.row = input.df %>%
      filter(`Target Gene(s)` == i)

    # get all the overactive regulons
    my.regulons1 = input.df %>%
      filter(`Target Gene(s)` == i) %>%
      pull(`Overactive Regulon(s)`) %>%
      unique()


    if (is.na(my.regulons1)) {
      cat("No Overactive regulons for :", i, "\n")
      drugs_all_regulon_pathways <-
        rbind(drugs_all_regulon_pathways,
              cbind(
                my.row,
                Reactome = NA,
                GO = NA,
                Hallmarks =  NA
              ))
    } else {
      # split regulons if multiple
      my.regulons2 <- strsplit(my.regulons1, split = ",")

      enrichment.list <- list()

      ## For each regulon get the regulon genes
      all.genes <- c()
      ii = 0
      total = length(my.regulons2[[1]])
      cat("Total ", total, "regulons are overactive\n")

      for (regulon in my.regulons2[[1]]) {
        ii = ii + 1

        regulon.genes <- regulons[[regulon]]
        cat("Total ",
            length(regulon.genes),
            "genes in regulon",
            regulon,
            "\n")

        # combine all regulon genes
        all.genes <- append(all.genes, regulon.genes)
      }

      # remve duplicate genes
      all.genes = unique(all.genes)

      ## Convert gene symbols
      gene.entrez <-
        AnnotationDbi::select(
          org.Hs.eg.db,
          keys = all.genes,
          column = "ENTREZID",
          keytype = "ENSEMBL",
          multiVals = "first"
        )$ENTREZID
      gene.symbols <-
        AnnotationDbi::select(
          org.Hs.eg.db,
          keys = all.genes,
          column = "SYMBOL",
          keytype = "ENSEMBL",
          multiVals = "first"
        )$SYMBOL

      enrichment.list <-
        try(regulon_enrichment(my.genes = gene.symbols))
      #reactome.enrich <- enrichPathway(gene = gene.entrez, pvalueCutoff = 0.05, pAdjustMethod = "BH")

      drugs_all_regulon_pathways <-
        rbind(
          drugs_all_regulon_pathways,
          cbind(
            my.row,
            Reactome = paste0(enrichment.list$table8$Description[1:5], collapse =
                                ":"),
            GO = paste0(enrichment.list$table6$Description[1:5], collapse =
                          ":"),
            Hallmarks =  paste0(enrichment.list$table1$Description, collapse =
                                  ":"),
            KEGG =  paste0(enrichment.list$table7$Description, collapse =
                                  ":")
          )
        )
    }



  }

  return(drugs_all_regulon_pathways)

}
