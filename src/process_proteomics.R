proteomics <- read.csv(
    file.path(PROJECT_PATH, 'data/1602426_TableS4.csv'), 
    na.strings = ""
)
ensembl <- useDataset("celegans_gene_ensembl", mart = useMart("ensembl"))
ids <- getBM(
    attributes = c("uniprotswissprot", "uniprotsptrembl", "wormbase_gene"),
    mart = ensembl
)
ids$protID <- paste0(ids$uniprotswissprot, ids$uniprotsptrembl)
proteomics2 <- proteomics[!is.na(match(proteomics$WB, ids$wormbase_gene)),]
proteomics2$id <- ids$protID[match(proteomics2$WB, ids$wormbase_gene)]
proteomics2 <- proteomics2[!duplicated(proteomics2$id),]
proteomics <- proteomics2[, c(1, 3:(ncol(proteomics2)-1))]
proteomics$UniProtKB <- proteomics2$id
write.csv(proteomics, 'data/1602426_TableS4_processed.csv', row.names = FALSE, quote = F)
