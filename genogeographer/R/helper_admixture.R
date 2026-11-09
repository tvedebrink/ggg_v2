# xlx <- function(x) ifelse(x<=0, 0, x*log(x)) ## x*log(x), per definition 0 for x = 0.
## it also catches the case of x negative in loops below.

EM <- function(x,n, tol = 1e-8){ ## x = (xj, xk), n = (nj, nk)
  dev <- Inf
  ## dev_trace <- dev
  p <- x/n                  ## Init estimate of pj and pk: xj/(2n)
  ## p_trace <- list(p)
  ## pY_trace <- NULL
  p[p == 0] <- 1e-6
  p[p == 1] <- 1 - 1e-6
  if(sum(p)==0) p <- p+1e-6 ## if both 0, set p to (eps, eps)
  if(sum(p)==2) p <- p-1e-6 ## if both 1, set p to (1-eps, 1-eps)
  while(dev > tol){
      p0 <- p                         ##
      p. <- prod(p)                   ## p. = pj*pk
      pY_x0_1 <- (p-p.)/(sum(p)-2*p.) ## P(Y | x0 = 1)
      p <- (x+pY_x0_1)/(n+1)          ## p_j^{(t+1)} = {x_j + E(y|x0 = 1, xj,xk)} / {2nj + 1}
      dev <- mean(abs(p-p0))          ## difference between p^{(t)} and p^{(t+1)}
      ## dev_trace <- c(dev_trace, dev)
      ## p_trace <- c(p_trace, list(p))
      ## pY_trace <- c(pY_trace, list(pY_x0_1))
  }
  # list(pY_x0_1, # er fordelingen af Y givet at x0 = x0^1 + x0^2 = 1
  #     dev = dev_trace, p_trace = p_trace, pY_trace = pY_trace)
  pY_x0_1
}

## EM(x = c(1238, 150), n = c(1244, 150))

logQ <- function(y, x, n){ ## -2 log-likelihood ratio
    -2*(xlx(n)-xlx(n+1) + xlx(x+y)-xlx(x) + xlx(n-x+1-y)-xlx(n-x) - xlx(y)-xlx(1-y))
}

het_marker <- function(x,n, tol){
  pY_x0_1 <- EM(x, n, tol = tol)
  dev <- 0
  for(y in 1:0){
    dev <- dev + pY_x0_1[2-y]*(logQ(y,x[1],n[1])+logQ(1-y,x[2],n[2]))
  }
  list(dev = dev, pY_x0_1 = pY_x0_1)
}

dev_marker <- function(x0, x, n, tol){
  if(x0==1){ ## if heterozygote - do this;
    return(het_marker(x, n, tol = tol))
  }
  else{ ## else add contributions from x0=0 and x0=2 (in cases: a/2 equals 0 or 1)
    return(list(dev = logQ(x0/2,x[1],n[1]) + logQ(x0/2,x[2],n[2])))
  }
}

p_x0 <- function(x, n){
  p <- x/n
  p. <- prod(p)
  c(prod(1-p), sum(p)-2*p., p.)
}

log_Padmix <- function(x0, f1, f2){
  case_when(
    x0 == 0 ~ log(1-f1) + log(1-f2),
    x0 == 1 ~ log(f1*(1-f2) + (1-f1)*f2),
    x0 == 2 ~ log(f1) + log(f2),
    TRUE ~ NA_real_
  )
}

varlog_Padmix <- function(x0, f1, f2, n1, n2){
  case_when(
    x0 == 0 ~ 1/(2*(1-f1)*(1-f2))*( f1*f2/(2*n1*n2) + (1-f1)*f2/n1 + (1-f2)*f1/n2 ),
    x0 == 1 ~ (n2*f1*(1-f1)*(1-2*f2) + n1*f2*(1-f2)*(1-2*f1) + 2*f1*f2*( (1-f1)*(1-f2) + 2*n2*f2*(1-f1) + 2*n1*f1*(1-f2) )) /
      (2*n2*n2*(f1*(1-f2) + f2*(1-f1))^2),
    x0 == 2 ~ 1/(2*f1*f2)*( (1-f1)*(1-f2)/(2*n1*n2) + (1-f1)*f2/n1 + (1-f2)*f1/n2 ),
    TRUE ~ NA_real_
  )
}

marker_z <- function(x0, xj, nj, xk, nk, return_p = TRUE, tol = 1e-8){  ## we assume that n is 2*n
  # x0: count main allele (0:2)
  # x: count af main allele in pop_j and pop_k
  # n: n1 and n2 in pop_j and pop_k ## ASSUMES n is 2*n_individuals
    #pe <- pest(a,p,n)
  marker_stats <- lapply(0:2, dev_marker, x=c(xj,xk), n=c(nj,nk), tol = tol)
  pY_x0_1 <- marker_stats[[2]]$pY_x0_1
  devs <- unlist(lapply(marker_stats, "[[", "dev"))
  obs <- devs[x0+1]
  p_x0 <- p_x0(c(xj,xk), c(nj,nk))
  E_z <- sum(p_x0*devs)
  V_z <- sum(p_x0*(devs-E_z)^2)
  # list(tibble(z_center = obs-E_z, V_z = V_z, pY_x0_1 = list(pY_x0_1)))
  if(!return_p) return(list(tibble(z_raw = obs-E_z, z_var = V_z)))
  list(tibble(z_raw = obs-E_z, z_var = V_z, p1 = pY_x0_1[1], p2 = pY_x0_1[2]))
}

markers_z <- Vectorize(marker_z, vectorize.args = list("x0", "xj", "nj", "xk", "nk"))

