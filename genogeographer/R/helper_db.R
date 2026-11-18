## Helper functions

dist_x0_cond_x <- function(n,x) {
  x1 <- n+2-x
  cbind(p0 = x1*(x1-1), p1 = 2*x*x1, p2 = x*(x-1))/((n+2)*(n+1))
}

dist_cond_x <- function(n,x) {
  x1 <- n+2-x
  cbind(p0 = x1*(x1-1), p1 = 2*x*x1, p2 = x*(x-1))/((n+2)*(n+1))
}


xlx <- function(x){
  x[x<=0] <- 1
  log(x)*x
}

ss_moments <- function(x0, x1, N){
  x <- x1+x0
  p <- dist_x0_cond_x(N,x)
  x <- matrix(rep(x, each = 3), ncol = 3, byrow = TRUE)
  N <- matrix(rep(N, each = 3), ncol = 3, byrow = TRUE)
  m02 <- matrix(0:2, ncol = 3, nrow = nrow(x), byrow = TRUE)
  xx <- x - m02*(p>0)
  stat <- xlx(xx)+xlx(N-xx)-(m02*(p>0) == 1)*2*log(2)
  Estat <- rowSums((p*stat)*(p>0))
  stat <- stat-Estat
  Vstat <- rowSums((p*stat^2)*(p>0))
  list(z_exp = Estat, z_var = Vstat)
}

log_P <- function(x0, p){
  case_when(
    x0 == 0 ~ 2*log(1-p),
    x0 == 1 ~ log(2) + log(p) + log(1-p),
    x0 == 2 ~ 2*log(p),
    TRUE ~ NA_real_
    )
}

varlog_P <- function(x0, n, p){
  q <- 1-p
  case_when(
    x0 == 0 ~ 4*p/(q*n) + 2*p*(3-5*q)/(q*n)^2 + p*(1-6*q*p)/(q*n)^3,
    x0 == 1 ~ ( (1-4*p*q)/n + (10*p*q-2)/n^2 - (6*p*q-1)/n^3 )/(p*q),
    x0 == 2 ~ 4*q/(p*n) + 2*q*(3-5*p)/(p*n)^2 + q*(1-6*p*q)/(p*n)^3,
    TRUE ~ NA_real_
  )
}

##

nullset <- function(x, ...){
  x <- subset(x, ...)
  if(length(x) == 0) return(NULL)
  x
}

locus_guess <- function(profile, loci = "^rs.*", tol = 0.9){
  nullset(names(profile), map_lgl(profile, ~ mean(grepl(pattern = loci, .x)) > tol))
}

genotype_guess <- function(profile, DNA_chars = c("A", "C", "G", "T", "N", "-"), tol = 0.9){
  nullset(names(profile), map_lgl(profile, ~ all(nchar(.x)<=2) & mean(unlist(strsplit(paste(.x),"")) %in% DNA_chars) > tol))
}

profile_guess <- function(profile, tol = 0.9){
  list(locus = locus_guess(profile, tol = tol), genotype = genotype_guess(profile, tol = tol))
}
