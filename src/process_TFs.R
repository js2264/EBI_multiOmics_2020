TFs.bed.folder <- '~/../nh337/data/modencode/ce11/idr_narrowpeak/'
TFs <- list.files('~/../nh337/data/modencode/ce11/idr_narrowpeak/') %>% gsub('-GFP.*', '', .) %>% gsub('-1_.*', '-1', .) %>% unique() %>% grep('_LAST_UPDATED', ., value = TRUE, invert = TRUE)
TFs.peaks <- parallel::mclapply(mc.cores = 50, TFs, function(TF) {
    TF.peaks.files <- list.files(TFs.bed.folder, pattern = TF, full.names = T)
    lapply(TF.peaks.files, function(FILE) {read.table(FILE) %>% makeGRangesFromDataFrame(seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', keep.extra.columns = T)}) %>% GRangesList() %>% unlist() %>% GenomicRanges::reduce()
}) %>% setNames(TFs)
lapply(names(TFs.peaks), function(NAME) {export(TFs.peaks[[NAME]], paste0('../data/bed/', NAME, '.bed'))})
