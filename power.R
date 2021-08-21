#' ---
#' title: Power Calculation Example
#' subtitle: using `paramtest` package
#' author: Rafael F. Bressan
#' date: "`r Sys.Date()`"
#' ---
#'
#' This is an old script I used to benchmark some power test algorithms, 
#' including the usage of the `paramtest` package and its parallelism 
#' capabilities. You will notice I don't use `data.table` and `fixest` here, 
#' up to you to adapt the following code, :-)
#' 
#+ message = FALSE, warning = FALSE
library(dplyr)
library(sandwich)
library(paramtest)
library(microbenchmark)
library(ggplot2)
#' Load data
data <- read.csv("Data/escolas.csv")
# Clean up the data. Never forget it!!
data <- data %>%
  filter(!is.na(logico)) %>%
  mutate(escola = if_else(escola == "S/escola", "99", escola)) %>%
  mutate(escola = as.integer(escola))
  
#' We'll consider a fine grid of treatment effects
grid = seq(0, 0.5, 0.05)

#' 
pt_func <- function(simNum, eff) {
  #Drawing a sample, with replacement, of schools from our data. The idea here is that our sample well approximates the population distribution.
  #and that the sampling scheme is close to random sampling from the population of interest. 
  sample.draw = sample(list.schools, replace = TRUE)
  
  #Construct the artificial data from a draw. Note that repeated schools should be treated as distinct so as to mimic random sampling.
  list.new.data = lapply(1:length(sample.draw), function(x){
    extract = data[data$escola == sample.draw[x],]
    return(extract)
  })
  
  #Concatenating data on lists
  # data.artificial ends up with lower number of observations that original 
  # data
  data.artificial = do.call(rbind, list.new.data)
  
  #Next, we select that half of the schools will be treated
  treat.draw = sample(unique(data.artificial$escola),size = length(unique(data.artificial$escola))/2, replace = F)
  data.artificial$treat = 1*(data.artificial$escola %in% treat.draw)
  
  #Create outcome
  data.artificial$y = data.artificial$logico + data.artificial$treat*eff
  
  #Running models and storing whether we reject the null that effect is 0
  model1 = lm(y~treat, data = data.artificial)
  se = sqrt(diag(sandwich::vcovCL(model1, cluster = data.artificial$escola)))
  rej.col.1 = 1*(abs(model1$coefficients["treat"]/se["treat"]) > qnorm(0.975) )
  
  model2 = lm(y~treat+mulher, data = data.artificial)
  se = sqrt(diag(sandwich::vcovCL(model2, cluster = data.artificial$escola)))
  rej.col.2 = 1*(abs(model2$coefficients["treat"]/se["treat"]) > qnorm(0.975) )
  
  model3 = lm(y~treat+mulher+idade, data = data.artificial)
  se = sqrt(diag(sandwich::vcovCL(model3, cluster = data.artificial$escola)))
  rej.col.3 = 1*(abs(model3$coefficients["treat"]/se["treat"]) > qnorm(0.975) )
  
  return(c(M1 = rej.col.1, M2 = rej.col.2, M3 = rej.col.3))
}

#' Passing data frame along
pt_func2 <- function(simNum, df, eff) {
  list.schools <- unique(df$escola)
  #Drawing a sample, with replacement, of schools from our data. The idea here is that our sample well approximates the population distribution.
  #and that the sampling scheme is close to random sampling from the population of interest. 
  sample.draw = sample(list.schools, replace = T)
  
  #Construct the artificial data from a draw. Note that repeated schools should be treated as distinct so as to mimic random sampling.
  list.new.data = lapply(1:length(sample.draw), function(x){
    extract = df[df$escola == sample.draw[x],]
    #extract$escola = x
    return(extract)
  })
  
  #Concatenating data on lists
  # data.artificial ends up with lower number of observations that original 
  # data
  data.artificial = do.call(rbind, list.new.data)
  
  #Next, we select that half of the schools will be treated
  treat.draw = sample(unique(data.artificial$escola),
                      size = length(unique(data.artificial$escola))/2, 
                      replace = F)
  data.artificial$treat = 1*(data.artificial$escola %in% treat.draw)
  
  #Create outcome
  data.artificial$y = data.artificial$logico + data.artificial$treat*eff
  
  #Running models and storing whether we reject the null that effect is 0
  model1 = lm(y~treat, data = data.artificial)
  se = sqrt(diag(sandwich::vcovCL(model1, cluster = data.artificial$escola)))
  rej.col.1 = 1*(abs(model1$coefficients["treat"]/se["treat"]) > qnorm(0.975) )
  
  model2 = lm(y~treat+mulher, data = data.artificial)
  se = sqrt(diag(sandwich::vcovCL(model2, cluster = data.artificial$escola)))
  rej.col.2 = 1*(abs(model2$coefficients["treat"]/se["treat"]) > qnorm(0.975) )
  
  model3 = lm(y~treat+mulher+idade, data = data.artificial)
  se = sqrt(diag(sandwich::vcovCL(model3, cluster = data.artificial$escola)))
  rej.col.3 = 1*(abs(model3$coefficients["treat"]/se["treat"]) > qnorm(0.975) )
  
  return(c(M1 = rej.col.1, M2 = rej.col.2, M3 = rej.col.3))
}


#' _Testing paramtest_
#' 
#' First let as global variables: data and list.schools
#' 
list.schools <- unique(data$escola)
pt_power <- run_test(pt_func, n.iter = 1000, output = "data.frame",
                     eff = 0.5)
#' Test power results for one effect level only
colMeans(results(pt_power))[2:4]

#' Now test passing data to pt_func2
pt_power2 <- run_test(pt_func2, n.iter = 1000, output = "data.frame",
                                 df = data, eff = 0.5)
#' Test power results
colMeans(results(pt_power2))[2:4]
#' 
#' WARNING: Benchmarking may take a long time to run depending on your number of
#' simulations and number of replication in the benchmark (parameter "times")
#' 
#' Benchmarking. Looks like there is no harm in passing the data.frame to the
#' function. 
#+ echo = TRUE
pt_bench <- microbenchmark(
  pt1 = run_test(pt_func, n.iter = 1000, output = "data.frame", eff = 0.5),
  pt2 = run_test(pt_func2, n.iter = 1000, output = "data.frame", df = data, 
                 eff = 0.5),
  times = 5
)
#'
#+ warning = FALSE
autoplot(pt_bench)

#' Now let's vary the eff argument, to capture the Minimum Detectable Effect
#'
#+ cache = TRUE 
pt_pwr2 <- grid_search(pt_func2, n.iter = 1000, output = "data.frame",
                       df = data, params = list(eff = grid))
pwr2_df <- results(pt_pwr2) %>% 
  group_by(eff.test) %>% 
  summarise(across(contains("treat"), mean))
pwr2_df

#'
#+ cache = TRUE
pt_par <- grid_search(pt_func2, n.iter = 1000, output = "data.frame",
                      df = data, params = list(eff = grid),
                      parallel = "multicore", ncpus = 4)
par_df <- results(pt_par) %>% 
  group_by(eff.test) %>% 
  summarise(across(contains("treat"), mean))
par_df
#' Is it worth to use parallelization? 
#' 
#' Benchmark it! Parallelized simulation is almost twice as fast.
#' 
#+ cache = TRUE
par_bench <- microbenchmark(
  no_par = grid_search(pt_func2, n.iter = 1000, output = "data.frame",
                       df = data, params = list(eff = grid)),
  par = grid_search(pt_func2, n.iter = 1000, output = "data.frame",
                    df = data, params = list(eff = grid),
                    parallel = "multicore", ncpus = 4),
  times = 5
)
par_bench

