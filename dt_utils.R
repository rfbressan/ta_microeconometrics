#' Data.Table utils
#' 
library(magrittr) # Use the pipe operator
# unnest_dt function --------------------------------------------
unnest_dt <- function(dt, col, id) {
  stopifnot(data.table::is.data.table(dt))
  by <- substitute(id)
  col <- substitute(unlist(col, recursive = FALSE))
  dt[, eval(col), by = eval(by)]
}

#' Summaries by binary treatment 0/1
#' 
dt_summary_bin <- function(dt, treat, cols) {
  stopifnot(data.table::is.data.table(dt))
  if (dt[, uniqueN(get(treat))] != 2L)
    stop("Treatment column must be binary")
  else if (!all(dt[, unique(get(treat))] %in% c(0,1)))
    stop("Treatment levels must be either 0 or 1")
  
  d <- dt[, by = .(get(treat)),
             list(
               variable = rep(cols, uniqueN(get(treat))),
               nobs = .N,
               mean = sapply(.SD, function(x) mean(x, na.rm = TRUE)),
               sd = sapply(.SD, function(x) sd(x, na.rm = TRUE))
             ),
             .SDcols = cols] 
  desc <- dcast(d, variable~get, 
                  value.var = c("nobs", "mean", "sd"))
  setcolorder(desc,
              c("variable", 
                grep("_0", colnames(desc), value = TRUE),
                grep("_1", colnames(desc), value = TRUE)))
  return(desc)
}

# Three part table function -----------------------------------------------
threeparttable_notes <- function(file_in, file_out = NULL) {
  if (is.null(file_out)) file_out <- file_in
  
  char_vec <- readLines(file_in) 
  char_vec <- gsub("(\\\\begin\\{tablenotes\\}).*", "\\1[flushleft]", char_vec)
  char_vec <- gsub("(\\\\begin\\{TableNotes\\}).*", "\\1[flushleft]", char_vec)
  char_vec <- gsub("(?<=& )(X)", "$\\\\checkmark$", char_vec, perl = TRUE)
  
  writeLines(char_vec, con = file_out)
}

# Plot time series - control and treatment --------------------------------
dt_plot_treat_ts <- function(dt, xaxis, yaxis, treat_group, 
                             xlabel = "", ylabel = "") {
  stopifnot(data.table::is.data.table(dt))
  
  ggplot2::ggplot(dt, aes(get(xaxis))) +
    ggplot2::geom_line(aes(y = get(yaxis), color = factor(get(treat_group)))) +
    ggplot2::scale_x_date(breaks = "years", date_labels = "%Y") +
    ggplot2::labs(x = xlabel, y = ylabel) +
    ggplot2::scale_color_discrete(
      labels = c("Controle", "Tratamento"),
      guide = guide_legend(title = "Grupo")) +
    ggplot2::theme_classic()
}

# Cumulative paste of strings -------------------------------------------- 
# Useful to create incremental controls in a regression
cumpaste <- function(x, .sep = "+") {
  Reduce(function(x1, x2) paste(x1, x2, sep = .sep), 
         x, 
         accumulate = TRUE)
  }

# Balance of covariates table -------------------------------------------

balance_tbl <- function(dt, treat, alpha = 0.05) {
  stopifnot(data.table::is.data.table(dt))
  # Compute the normalized difference
  norm_diff <- function(xt, xc, sdt, sdc) {
    (xt - xc)/(sqrt((sdt^2 + sdc^2)/2))
  }
  # Coverages
  coverage <- function(dt, treat_var, variables, values, alpha = 0.05) {
    # stopifnot(is.data.table(dt))
    
    quantiles <- dt[, by = .(get(treat_var), get(variables)),
                         lapply(.SD, function(x){
                           c(quantile(x, 1 - alpha/2),
                             quantile(x, alpha/2)
                           )
                         }),
                         .SDcols = values
    ][, q := rep(c("high", "low"), .N/2)]
    
    quantiles <- data.table::dcast(quantiles, get+get.1~q, value.var = "value")
    data.table::setnames(quantiles, c("get", "get.1"), c(treat_var, variables))
    quantiles[, not_treat_var := 1*(!get(treat_var))]
    
    pi_dt <- merge(dt_long, quantiles[, !c(..treat_var)], 
                   by.x = c(treat_var, variables),
                   by.y = c("not_treat_var", variables),
                   all.x = TRUE)
    pi_dt <- pi_dt[, by = .(get(treat_var), get(variables)),
                   .(pi = (1 - mean(value <= high)) + mean(value <= low))]
    pi_dt <- data.table::dcast(pi_dt, get.1~get, value.var = "pi")
    data.table::setnames(pi_dt, c("get.1", "0", "1"), c(variables, "pi_cont", "pi_treat"))
    return(pi_dt)
  }
  # Select the numerical columns and compute the normalized difference and log
  # ratio
  num_cols <- colnames(dt)[sapply(dt, is.numeric)]
  
  dt_long <- data.table::melt(dt[, ..num_cols], id.vars = treat,
                  variable.factor = FALSE)
  
  dt_summa <- dt_long[, by = .(variable, get(treat)),
                      .(N = .N,
                        mean = mean(value, na.rm = TRUE),
                        sd = sd(value, na.rm = TRUE)
                      )]
  data.table::setnames(dt_summa, "get", treat)
  dt_wide <- data.table::dcast(dt_summa, as.formula(sprintf("variable~%s", treat)), 
                               value.var = c("mean", "sd", "N"))
  dt_wide[,  `:=`(
    diff_mean = mean_1 - mean_0,
    std_err = sqrt(sd_1^2/N_1 + sd_0^2/N_0),
    t_stat = (mean_1 - mean_0)/sqrt(sd_1^2/N_1 + sd_0^2/N_0),
    norm_diff = norm_diff(mean_1, mean_0, sd_1, sd_0),
    log_ratio = log(sd_1) - log(sd_0)
  ) ]
  
  pi_dt <- coverage(dt_long, treat_var = "union", 
                    variables = "variable", 
                    values = "value",
                    alpha = alpha)
  # Join coverages to dt_wide
  bal_dt <- dt_wide[pi_dt, on = "variable"]
  # Reorder the columns by their index
  vari <- which(colnames(dt_wide) == "variable")
  treats <- grep("_1", colnames(dt_wide))
  cont <- grep("_0", colnames(dt_wide))
  data.table::setcolorder(bal_dt, c(vari, treats, cont))
  
  return(bal_dt)
}

# Print the balance table using kableExtra ------------------------------

balance_kbl <- function(bal_dt, digits = 2, alpha = 0.05) {
  stopifnot(data.table::is.data.table(bal_dt))
  Nt <- bal_dt[1, N_1]
  Nc <- bal_dt[1, N_0]
  t_title <- sprintf("Treatment ($N_t = %s$)", Nt)
  c_title <- sprintf("Control ($N_c = %s$)", Nc)
  kbl_head <- c(1, 2, 2, 5)
  names(kbl_head) <- c(" ", t_title, c_title, "Overlap Measures")
  
  kableExtra::kbl(bal_dt[, -c("N_1", "N_0", "diff_mean", "std_err")], 
                  digits = digits,
                  col.names = c("Variable", "Mean", "Std.Dev",
                                "Mean", "Std.Dev",
                                "t-stat", "Norm. Diff.", "Log Ratio", 
                                sprintf("$\\pi_c^{%s}$", alpha), 
                                sprintf("$\\pi_t^{%s}$", alpha)),
                  escape = FALSE) %>% 
    kableExtra::kable_classic(full_width = FALSE) %>% 
    kableExtra::add_header_above(header = kbl_head)
}
