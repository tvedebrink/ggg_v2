library(devtools)

db_ggg <- readRDS("../EUROFORGEN_120_SNP.rds")
## db_mh <- readRDS("../../pedro/rds_files_meta/MH_profiles_probs.Rds")

mh_list <- list(
  probs     = readRDS("inst/extdata/MH_profiles_probs.Rds"),
  snp_map   = readRDS("inst/extdata/"),
  out_metas = readRDS("inst/extdata/MH_out_metas.Rds")
)

genogeographer.devel::app_genogeo(db_list = db_list, reporting_panel = TRUE, mh_data = mh_data)

load_all()
app_genogeo(db_list = db_ggg, mh_data = mh_list)

library(tidyr)
library(readr)

hw <- readr::read_csv2("inst/extdata/mh_sample_wide.csv", col_types = cols(.default = "c"))

hw2 <- hw %>%
  separate_wider_regex(cols = sample_name, patterns = c("sample_name" = ".*", g = ".{1}$")) %>%
  arrange(g, .by = sample_name) %>%
  summarise(across(-g, ~ paste0(.x, collapse = "")), .by = sample_name)

hw2 %>% write_csv2("inst/extdata/mh_sample_wide2.csv")

hl2 <- hw2 %>% pivot_longer(cols = -sample_name, names_to = "locus", values_to = "allele")

hl2 %>% write_csv2("inst/extdata/mh_sample_long2.csv")

hl <- hl2 %>% separate_wider_regex(cols = allele, patterns = c(A = ".{1}", B = ".{1}"))

hl %>% write_csv2("inst/extdata/mh_sample_long.csv")
