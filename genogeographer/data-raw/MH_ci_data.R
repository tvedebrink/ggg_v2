## Builds inst/extdata/MH_ci_data.Rds: the reference haplotype counts and graph margins needed for
## the first-order (delta method) confidence interval of the MH profile log-likelihood.
## Run from the package root. Source: the full (non-lean) out_metas*.Rds of the 'microhaps' analysis.

src <- "../../pedro/rds_files_meta"
mhs <- readRDS(file.path(src, "out_metas_MHs.Rds"))
full <- readRDS(file.path(src, "out_metas.Rds"))
snp_map <- readRDS("inst/extdata/MH_snp_map.Rds")

ci_data <- lapply(names(mhs), function(meta){
  A_full <- full[[meta]]$A$A1
  out <- lapply(names(mhs[[meta]]), function(mh){
    m <- mhs[[meta]][[mh]]
    A <- m$A$A1
    snps <- snp_map$snp[snp_map$MH == mh][order(snp_map$position[snp_map$MH == mh])]
    stopifnot(identical(colnames(A), snps), all(nchar(A) == 1))
    ## rows must be the same haplotypes (same order) in every marker: the joint variance relies on it
    stopifnot(identical(unname(A), unname(A_full[, snps, drop = FALSE])))
    motif <- do.call(paste0, lapply(seq_len(ncol(A)), function(j) A[, j]))
    motifs <- unique(motif)
    cl <- lapply(m$cms$cms1, attr, "vars")
    sp <- Filter(Negate(is.null), lapply(m$sms$sms1, function(s) if(is.null(s)) NULL else attr(s, "vars")))
    list(
      motifs = motifs, row_code = match(motif, motifs),
      pos = lapply(c(cl, sp), match, table = snps),
      sign = c(rep(1L, length(cl)), rep(-1L, length(sp)))
    )
  })
  setNames(out, names(mhs[[meta]]))
})
names(ci_data) <- names(mhs)

saveRDS(ci_data, "inst/extdata/MH_ci_data.Rds", compress = "xz")
