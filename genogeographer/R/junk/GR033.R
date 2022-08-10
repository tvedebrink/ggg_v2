
## GR033

library(tidyverse)

gr33 <- read_csv("~/work/artikler/AIMs/data/valid_v2/data/Pakstis/GR033.csv") %>%
  rename(sample = Sample, locus = `Target ID`, genotype = Genotype) %>%
  ## left_join(
  genotype_x0(.) #, by = c("locus")) %>%

admix1_v2 <- read_rds("meta_DB_1admix2020-09-16Rds")
meta_v2 <- read_rds("meta_DB2020-09-16Rds")

##
DB_v1 <- read_rds("~/work/artikler/AIMs/reference_DB/GGG_v1_DB_14122020.Rds")
source("R/helper_db.R")
ggg_allele_list <- read_rds("ggg_allele_list.Rds")

meta_v1 <- make_x1(DB_v1, groups = "meta", exclude = c(1, 2)) %>%
  mutate(meta = case_when(
    meta == "EURO" ~ "EUROPE",
    meta == "SAHARA" ~ "AFRICA",
    meta == "SOMALI" ~ "HORN_AFRICA",
    meta == "N AFRC" ~ "N AFRICA",
    ## meta == "" ~ "",
    TRUE ~ meta
  ))

##

meta_v2 %>% select(1) %>% mutate(v2 = TRUE) %>%
  full_join(meta_v1 %>% select(1) %>% mutate(v1 = TRUE), by = "meta")

meta_compare <- meta_v1 %>% unnest(col = data) %>%
  inner_join(meta_v2 %>% unnest(col = data),
             by = c("meta", "locus", "x0"), suffix = c("_v1", "_v2"))

meta_compare_gr33 <- meta_compare %>%
  left_join(gr33 %>% mutate(GR033 = TRUE), by = c("locus", "x0"))

meta_compare_gr33 %>% filter(!is.na(GR033)) %>%
  group_by(meta) %>%
  summarise(
    z_v1 = sum(z_raw_v1)/sqrt(sum(z_var_v1)),
    z_v2 = sum(z_raw_v2)/sqrt(sum(z_var_v2))
  )

## meta_compare %>%
meta_compare_gr33 %>% filter(meta == "EUROPE") %>%
  mutate(GR033 = !is.na(GR033)) %>%
  ggplot(aes(x = z_raw_v1, z_raw_v2, colour = GR033)) +
  geom_abline() +
  geom_point() +
  facet_wrap(~meta)

meta_compare %>%
  ggplot(aes(x = logP_v1, logP_v2)) +
  geom_point() +
  facet_wrap(~meta)
