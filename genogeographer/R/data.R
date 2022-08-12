#' List of alleles in the Applied Biosystems Precision ID Ancestry Panel
#'
#' The data consists of the 164 AIMs included in genogeographer's reference data.
#' The 164 AISNPs are from the Applied Biosystems Precision ID Ancestry Panel excluding rs10954737
#' (due to high degree of missingness in the genotyped reference populations, see Mogensen et al, 2020)
#'
#' Two columns indicate whether the SNP is included in the Kidd and Seldin panels, respectively.
#'
#' @format A dataframe of 164 rows and 5 columns.
#' \describe{
#'   \item{locus}{The rs identifier for the SNP}
#'   \item{x1}{The first allelic variant at the locus. In the genogeographer setup, this is
#'   the first allele in lexicographic order. When a profile is converted to a integer valued
#'   representation the number of `x1` alleles is counted.}
#'   \item{x1}{The second allelic variant}
#'   \item{Kidd}{Boolian value indicating whether the locus is included in the Kidd panel}
#'   \item{Seldin}{Boolian value indicating whether the locus is included in the Seldin panel}
#' }
"ggg_allele_list"
