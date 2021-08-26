#' ---
#' title: Clusters, Correlations and Inference
#' author: Rafael F. Bressan
#' date: "`r Sys.Date()`"
#' ---
#' 
#+ warnings=FALSE, message=FALSE
library(data.table)
library(fixest)
library(kableExtra)

setDTthreads(100)
setFixest_nthreads(4)

nsim <- 1000
n <- c(10, 50, 500) # individuals per cluster
J <- c(5, 50, 500) # number of clusters
#' 
#' # Individual RCT inference
#' 
#' $y_{ij}=\beta W_{ij}+v_{j}+e_{ij}$
#' 
cl_rct <- function(n, J, alpha = 0.05) {
  while (TRUE) {
    W <- rbinom(n*J, 1, 0.5)
    if (!all(W == 1) & !all(W == 0)) break
  }
  e <- rnorm(n*J, sd = 0.2)
  cl <- rep(seq_len(J), each = n)
  dt <- data.table(W, e, cl = factor(cl))
  dt[, by = cl, `:=`(v = rnorm(1))]
  dt[, `:=`(y = 0.1*W + v + e,
            y0 = v + e)]
  #'
  reg <- feols(y~W, data = dt)
  reg0 <- feols(y0~W, data = dt)
  pvals <- c(y_sd = pvalue(reg)["W"],
             y_cl = pvalue(reg, cluster = "cl")["W"],
             y0_sd = pvalue(reg0)["W"],
             y0_cl = pvalue(reg0, cluster = "cl")["W"])
  lapply(pvals, function(x) x < alpha)
}
#'
#+ cache=TRUE
sim_dt <- CJ(n, J, sim = seq_len(nsim))
set.seed(123)
sim_dt[, rej := purrr::map2(n, J, cl_rct)]
sim_dt[, by = .(n, J, sim), 
       c("Pwr Std", "Pwr CRVE", "Size Std", "Size CRVE") := rej[[1]]]
sim_dt[, rej := NULL]

size <- sim_dt[, by = .(n, J), lapply(.SD, mean),
                .SDcols = c("Size Std", "Size CRVE")] |> 
  dcast(n~J, value.var = c("Size Std", "Size CRVE"))
power <- sim_dt[, by = .(n, J), lapply(.SD, mean),
                .SDcols = c("Pwr Std", "Pwr CRVE")] |> 
  dcast(n~J, value.var = c("Pwr Std", "Pwr CRVE"))
#'
#'
kbl(size, col.names = c("n", rep(c("5", "50", "500"), 2))) |> 
  kable_classic(full_width = FALSE) |> 
  add_header_above(c(" " = 1, "Standard" = 3, "CRVE" = 3))
#' The test size is adequate using regular inference even in the presence of
#' unobserved cluster-specific characteristics. 
#'
kbl(power, col.names = c("n", rep(c("5", "50", "500"), 2))) |> 
  kable_classic(full_width = FALSE) |> 
  add_header_above(c(" " = 1, "Standard" = 3, "CRVE" = 3))
#' 
#' We have more power using the regular inference method.
#'
#' # Treatment designation by cluster
#' 
#' $y_{ij}=\beta W_{j}+v_{j}+e_{ij}$
#' 
cl_grp <- function(n, J, alpha = 0.05) {
  e <- rnorm(n*J, sd = 0.2)
  cl <- rep(seq_len(J), each = n)
  while (TRUE) {
    W <- rep(rbinom(J, 1, 0.5), each = n)
    if (!all(W == 1) & !all(W == 0)) break
  }
  
  dt <- data.table(W, e, cl = factor(cl))
  dt[, by = cl, `:=`(v = rnorm(1)) ]
  dt[, `:=`(y = 0.4*W + v + e,
            y0 = v + e)]
  #'
  reg <- feols(y~W, data = dt)
  reg0 <- feols(y0~W, data = dt)
  pvals <- c(y_sd = pvalue(reg)["W"],
             y_cl = pvalue(reg, cluster = "cl")["W"],
             y0_sd = pvalue(reg0)["W"],
             y0_cl = pvalue(reg0, cluster = "cl")["W"])
  lapply(pvals, function(x) x < alpha)
}
#'
#+ cache=TRUE
sim1_dt <- CJ(n, J, sim = seq_len(nsim))
set.seed(123)
sim1_dt[, rej := purrr::map2(n, J, cl_grp)]
sim1_dt[, by = .(n, J, sim), 
       c("Pwr Std", "Pwr CRVE", "Size Std", "Size CRVE") := rej[[1]]]
sim1_dt[, rej := NULL]

size1 <- sim1_dt[, by = .(n, J), lapply(.SD, mean),
               .SDcols = c("Size Std", "Size CRVE")] |> 
  dcast(n~J, value.var = c("Size Std", "Size CRVE"))
power1 <- sim1_dt[, by = .(n, J), lapply(.SD, mean),
                .SDcols = c("Pwr Std", "Pwr CRVE")] |> 
  dcast(n~J, value.var = c("Pwr Std", "Pwr CRVE"))
#'
#'
kbl(size1, col.names = c("n", rep(c("5", "50", "500"), 2))) |> 
  kable_classic(full_width = FALSE) |> 
  add_header_above(c(" " = 1, "Standard" = 3, "CRVE" = 3))
#' 
#' The test size is completely **inadequate** (oversized) when using regular 
#' inference with cluster randomization. The researcher is better off using CRVE
#' inference in this case. 
#'
kbl(power1, col.names = c("n", rep(c("5", "50", "500"), 2))) |> 
  kable_classic(full_width = FALSE) |> 
  add_header_above(c(" " = 1, "Standard" = 3, "CRVE" = 3))
#' 


