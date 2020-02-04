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
#
#
#
# Process ATAC files
ATACage <- read.table('~/_20180427_devATAC-JJ-paper/_DESeq2-files/age_DESeq2-analysis.txt', header = T, row.names = 1)[,1:5]
all$is.new <- read_tsv('~/shared/data/reg_elements_dev_ageing_tissues_20190824.tsv')$atac_source == 'atac_tissues'
clusters <- read.table('~/20200115_thesis_snippets//20180530.supp-table.txt', header = TRUE, sep = '\t')
all_old <- all[!all$is.new,]
all_old$dev_cluster <- factor(clusters$clustersATAC.dev, levels = c(".", "G1", "G2", "G3", "G4", "I1", "I2", "H", "N+M", "Mix1", "Mix2", "Mix3", "Mix4", "Mix5", "Mix6", "Mix7", "Mix8", "stable"))
g <- data.frame(
    all_old[c(1:3, 6, 9, 10, 11)], 
    log2(ATACage+1),
    aging_cluster = clusters$clustersATAC.age
) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
REs <- g
#
res <- read.table('~/_devATAC-JJ-paper_2/all-regulatory-elements.cv-domains-ATACcov.txt', header = TRUE, sep = '\t')
df <- read.table('~/20200115_thesis_snippets/elife-37344-fig4-data1-v2.txt', header = TRUE, sep = '\t')
ATACage <- read.table('~/_20180427_devATAC-JJ-paper/_DESeq2-files/age_DESeq2-analysis.txt', header = T, row.names = 1)[,1:5]
ATACage_meancentered <- apply(ATACage, 1, function(ROW) {log2(ROW/mean(ROW))}) %>% t() %>% data.frame()
tab <- cbind(
    locus = res$locus, 
    WormBaseID = res$WormBaseID, 
    regulatory_class = res$regulatory_class, 
    ATACage,
    ageing_prom_cluster_label = df$ageing_prom_cluster_label
)
g <- granges(REs)
mcols(g) <- tab
export(g, 'data/ATAC-seq_aging.gff3')
