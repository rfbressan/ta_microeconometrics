#' ---
#' title: Assignment I
#' subtitle: THE IMPACT OF HIGH SCHOOL FINANCIAL EDUCATION EVIDENCE FROM A LARGE-SCALE EVALUATION IN BRAZIL
#' author: Rafael F. Bressan
#' date: "`r Sys.Date()`"
#' output:
#'   html_document:
#'     code_folding: show
#' ---
#' 
#+ message = FALSE, warning = FALSE
library(data.table)
library(haven) # Load Stata's dta files
library(fixest)
library(kableExtra)

setDTthreads(75) # data.table uses 3/4 of available cores
setFixest_nthreads(4) # fixest uses 4 cores

#' Reading the data
#+ cache = TRUE
df <- read_stata("Data/school_intervention_panel_final.dta", 
                 encoding = "latin1")
setDT(df) # set `df` as data.table by reference

#' Auxiliar function to extract regression statistics
sum_stats <- function(reg) {
  stopifnot(class(reg) == "fixest")
  
  n <- reg$nobs
  coef_tbl <- coeftable(reg)
  
  Control <- coef_tbl["(Intercept)", "Estimate"] 
  sdC <- coef_tbl["(Intercept)", "Std. Error"]
  Treatment <- coef_tbl["treatment", "Estimate"] + Control
  sdT <- coef_tbl["treatment", "Std. Error"] 
  pv <- coef_tbl["treatment", "Pr(>|t|))"] 
  return(data.table(nobs = n, avg_c = Control, sd_c = sdC, 
                    avg_t = Treatment, sd_t = sdT, pv = pv))
}

#' # Replicating the paper
#'
#' ## TABLE 1: Sum Stats and Balance
base_vars <- c("cd_escola", "id_geral", "treatment")
test <- c("female", "dumm_rp_08_bl", "dumm_rp_09_bl", "dumm_rp_24_bl", 
          "dumm_rp_14_bl", "dumm_rp_23_bl", "vl_proficiencia_bl")
aluno <- c("dumm_rp_49_bl", "dumm_rp_50_bl", "dumm_rp_65A_bl", "poupar_final2_bl",
           "dumm_rp_64A_bl", "dumm_negotiates_bl", "autonomia_final2_bl")
school <- c("matriculas", "docentes", "abandonona1sriemdio", "aprovaona1sriemdio")
dt1 <- df[round == 1 & !is.na(treatment)][order(cd_escola, id_geral)]
#' Count number of schools. column (1)
ns <- dt1[, sapply(c(..school, ..test, ..aluno), function(x){
  .SD[!is.na(get(x)), uniqueN(cd_escola)]
})]
#' School-level
dt1sc <- dt1[, c(..base_vars, ..school)]
dt1sc <- dt1sc[, by = cd_escola, lapply(.SD, first)]
sc_dt <- data.table(desc = sapply(dt1sc, attr, which = 'label'),
                     dep_var = colnames(dt1sc))[-(1:3)]

sc_regs <- lapply(sc_dt$dep_var, function(x){
  feols(as.formula(sprintf("%s~treatment", x)),
        cluster = ~cd_escola,
        lean = TRUE,
        data = dt1sc)
})
#' Extract summary statistics from regressions, creates a DT and column bind it
#' to sc_dt
sc_dt <- cbind(sc_dt, rbindlist(lapply(sc_regs, sum_stats)))
#' Standard Deviations
sc_std <- dt1sc[, by = treatment, lapply(.SD, sd, na.rm = TRUE),
                .SDcols = school] |> 
  melt(id.vars = "treatment") |> 
  dcast(variable~treatment, value.var = "value")

#' Student background and financial characteristics
dt1sum <- dt1[, c(..base_vars, ..test, ..aluno)]
st_dt <- data.table(desc = sapply(dt1sum, attr, which = 'label'),
                     dep_var = colnames(dt1sum))[-(1:3)]
                                            
st_regs <- lapply(st_dt$dep_var, function(x){
  feols(as.formula(sprintf("%s~treatment", x)),
        cluster = ~cd_escola,
        lean = TRUE,
        data = dt1sum)
})

#' Extract summary statistics from regressions, creates a DT and column bind it
#' to bal_dt
st_dt <- cbind(st_dt, rbindlist(lapply(st_regs, sum_stats)))
#' Standard Deviations
st_std <- dt1sum[, by = treatment, lapply(.SD, sd, na.rm = TRUE),
                .SDcols = c(test, aluno)] |> 
  melt(id.vars = "treatment") |> 
  dcast(variable~treatment, value.var = "value")
std <- rbind(sc_std, st_std)
#' Row bind data.tables and create the column with number of schools
sum_dt <- rbind(sc_dt, st_dt)
sum_dt[, nschools := ns[dep_var]]
#' Reordering columns
setcolorder(sum_dt, c("desc", "dep_var", "nschools"))
#' Show tables
kbl(sum_dt, digits = 2) |> 
  kable_classic(full_width = TRUE)

kbl(std, digits = 2) |> 
  kable_classic(full_width = FALSE)