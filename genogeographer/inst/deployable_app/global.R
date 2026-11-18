## REQUIRED: named list of reference databases (the names appear in the app's database selector).
## Each entry is created by make_x1()/admix_dbs(), e.g. downloaded from the
## "Add reference populations" tab.
db_list <- readRDS("db_list.Rds")

## REQUIRED: TRUE shows the "Download HTML report" panel after an analysis.
reporting_panel <- TRUE

## OPTIONAL: override the microhaplotype data bundled with the package.
# mh_data <- list(
#   probs     = readRDS("MH_profiles_probs.Rds"),
#   snp_map   = readRDS("MH_snp_map.Rds"),
#   out_metas = readRDS("MH_out_metas.Rds")
# )
