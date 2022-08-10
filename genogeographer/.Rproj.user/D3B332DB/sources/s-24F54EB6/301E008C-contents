## Z-score

if(FALSE){

  dane1_raw <- read_csv(paste0("~/work/artikler/AIMs/data/helle_populationer/RGA profiles/Danes_113432.csv"), col_types = cols())
  dane1 <- dane1_raw %>% genotype_x0()

  dane_res <- ggg_score(dane1, pop_DB)

  dane_res_loo <- ggg_score(dane1, pop_DB, LOO = "DK")
  list("In" = anti_join(dane_res, dane_res_loo), "LOO" = anti_join(dane_res_loo, dane_res)) %>%
    bind_rows(.id = "type") %>%
    arrange(pop) # %>% mutate(across(where(is.numeric), ~sprintf("%0.8f",.x)))

  dane_res_meta <- ggg_score(dane1, meta_DB)
    dane_res_meta_loo <- ggg_score(dane1, meta_DB, LOO = "EUROPE")

  list("In" = anti_join(dane_res_meta, dane_res_meta_loo), "LOO" = anti_join(dane_res_meta_loo, dane_res_meta)) %>%
    bind_rows(.id = "type") %>%
    arrange(meta) # %>% mutate(across(where(is.numeric), ~sprintf("%0.8f",.x)))

  dane_res_meta_loo <- ggg_score(dane1, meta_DB, LOO = "CDP")

  dane_res_1admix <- ggg_score(dane1, pop_DB_1admix)
  dane_res_meta_1admix <- ggg_score(dane1, meta_DB_1admix)

  dane_res_meta_ <- bind_rows(dane_res_meta, dane_res_meta_1admix)
  attr(dane_res_meta_, "info") <- list(dane_res_meta, dane_res_meta_1admix) %>% map_dfr(~attr(.x,"info"))
  dane_res_ <- bind_rows(dane_res, dane_res_1admix)

### LR table

  LR_table(dane_res)
  LR_table(dane_res_meta)

  dane_lr <- LR_table(bind_rows(dane_res_meta, dane_res_meta_1admix))
  saveRDS(dane_lr, file = "../dane_lr.Rds")

  meta_info <- attr(db_list$`Precision ID`$meta$db, "info") %>% select(starts_with("meta"), n) %>% mutate(n = paste(n))
  meta_1admix_info <- attr(db_list$`Precision ID`$meta$admix, "info")
  meta_info_ <- bind_rows(meta_info, meta_1admix_info)

  pop_info <- attr(db_list$`Precision ID`$pop$db, "info") %>% select(starts_with("pop"), n) %>% mutate(n = paste(n))
  pop_1admix_info <- attr(db_list$`Precision ID`$pop$admix, "info")
  pop_info_ <- bind_rows(pop_info, pop_1admix_info)

  LR_plot(dane_res_meta, db_info = meta_info)

## DB API

  data_file <- read_rds("~/work/artikler/AIMs/reference_DB/DB_17062020.Rds")
  data_file <- data_file %>% mutate(meta = ifelse(meta == "AFRICA", "SUB-SAHARAN AFRICA", meta))
  data_file <- data_file %>% filter(pop != "POR")
  ## data_file_v1 <- read_rds(paste0("~/", work, "artikler/AIMs/reference_DB/GGG_v1_DB_14122020.Rds"))

  data_file %>% write_csv(file = "../data_csv.csv")
  data_file %>% write_csv2(file = "../data_csv2.csv")
  data_file %>% rio::export(file = "../data_excel.xlsx")

  pop_latlon <- readRDS("~/work/R/ggg_v2/genogeographer/data/population_latlon.Rds")

  pop_latlon %>% write_csv(file = "../info_csv.csv")
  pop_latlon %>% write_csv2(file = "../info_csv2.csv")
  pop_latlon %>% rio::export(file = "../info_excel.xlsx")

  pop_count <- count_alleles(data_file, groups = "pop", exclude = c("sample", "meta"))

  ### HIRO DATA

  hiro_file_all <- "~/work/artikler/AIMs/hiro2022/data/Data_Hiro_1000G_revised_allsets_FINAL.xlsx"
  hiro_sheets_all <- readxl::excel_sheets(hiro_file_all) %>% set_names(.)
  hiro_data_all <- hiro_sheets_all %>%
    imap(~ readxl::read_excel(path = hiro_file_all, sheet = .y, range = readxl::cell_limits(c(3,2), c(601, NA))))

  hiro_data_all$`75 SNPs`$rs10512572 %>% table()
  hiro_DB <- hiro_data_all$`75 SNPs` %>% mutate(across(starts_with("rs"), as.integer)) %>%
    make_x1(df = ., latlon = pop_latlon, groups = "Origin", exclude = c("Sample"))

  ## GGG v2
  pop_DB <- make_x1(df = data_file, latlon = pop_latlon, groups = "pop", exclude = c("sample", "meta"))
  meta_DB <- make_x1(df = data_file, latlon = pop_latlon, groups = "meta", exclude = c("sample", "pop"),
                     out_of_place = c("CEU", "GIH", "ITU", "STU"))

  pop_DB_1admix <- admix_dbs(pop_DB, no_cores = 8)
  meta_DB_1admix <- admix_dbs(meta_DB, no_cores = 8)

  db_list <- list("Precision ID" =
                    list(pop = list(db = pop_DB, admix = pop_DB_1admix),
                         meta = list(db = meta_DB, admix = meta_DB_1admix)))
  db_list <- c(db_list,
               list("Seldin loci" = db_set(db_list, seldin_loci),
                    "Kidd loci" = db_set(db_list, kidd_loci))
               )

  dane_res_pop <- ggg_score(dane1, db_list$`Precision ID`$pop$db)
  dane_res_pop_admix <- ggg_score(dane1, db_list$`Precision ID`$pop$admix)
  dane_res_meta <- ggg_score(dane1, db_list$`Precision ID`$meta$db)
  dane_res_meta_admix <- ggg_score(dane1, db_list$`Precision ID`$meta$admix)

  dane_res_meta_ <- bind_rows(dane_res_meta, dane_res_meta_admix)
  attr(dane_res_meta_, "info") <- list(dane_res_meta, dane_res_meta_admix) %>% map_dfr(~attr(.x,"info"))
  dane_res_pop_ <- bind_rows(dane_res_pop, dane_res_pop_admix)
  attr(dane_res_pop_, "info") <- list(dane_res_pop, dane_res_pop_admix) %>% map_dfr(~attr(.x,"info"))


  db_list %>% write_rds(file = "data/db_list.Rds")

  db_list <- read_rds(file = "data/db_list.Rds")
  pop_DB <- db_list$`Precision ID`$pop$db
  pop_DB_1admix <- db_list$`Precision ID`$pop$admix
  meta_DB <- db_list$`Precision ID`$meta$db
  meta_DB_1admix <- db_list$`Precision ID`$meta$admix

  app_genogeo(db_list = db_list)

  ## GGG v1
  # pop_DB_v1 <- make_x1(data_file_v1, groups = "pop", exclude = c(1, 3))
  # meta_DB_v1 <- make_x1(data_file_v1, groups = "meta", exclude = c(1, 2))
  #
  # meta_DB_1admix_v1 <- admix_dbs(meta_DB_v1)

}
