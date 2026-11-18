# Simulated example files for the "Add reference populations" tab (not real data).
# Writes inst/deployable_app/EUROFORGEN_global_db.csv and EUROFORGEN_global_info.csv.
set.seed(2026)
load("data/ggg_allele_list.rda")

info <- data.frame(
  pop = c("DK", "NL", "CN", "JP", "NG", "KE"),
  population = c("Danish", "Dutch", "Chinese", "Japanese", "Nigerian", "Kenyan"),
  meta = c("EUR", "EUR", "EAS", "EAS", "AFR", "AFR"),
  metapopulation = c("Europe", "Europe", "East Asia", "East Asia", "Africa", "Africa"),
  lat = c(56.0, 52.2, 35.0, 36.0, 9.1, 0.0),
  lon = c(10.0, 5.3, 105.0, 138.0, 8.7, 37.9)
)

n_per_pop <- 80
loci <- ggg_allele_list$locus
meta_freq <- sapply(unique(info$meta), function(m) rbeta(length(loci), 0.6, 0.6))
rownames(meta_freq) <- loci

samples <- do.call(rbind, lapply(seq_len(nrow(info)), function(i) {
  p <- pmin(pmax(meta_freq[, info$meta[i]] + rnorm(length(loci), 0, 0.03), 0.02), 0.98)
  geno <- sapply(seq_along(loci), function(l) {
    a <- ifelse(runif(n_per_pop * 2) < p[l], ggg_allele_list$x1[l], ggg_allele_list$x2[l])
    paste0(a[seq_len(n_per_pop)], a[n_per_pop + seq_len(n_per_pop)])
  })
  colnames(geno) <- loci
  data.frame(sample = sprintf("%s%02d", info$pop[i], seq_len(n_per_pop)),
             pop = info$pop[i], meta = info$meta[i], geno, check.names = FALSE)
}))

write.csv(samples, "inst/deployable_app/EUROFORGEN_global_db.csv", row.names = FALSE)
write.csv(info, "inst/deployable_app/EUROFORGEN_global_info.csv", row.names = FALSE)
