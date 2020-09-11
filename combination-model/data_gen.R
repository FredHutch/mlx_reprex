# ---------- Functions for simulating bolus PK starting with standard pk parameterization -------

# C(t) = dose * (Aexp(-alpha t) + Bexp(-beta t)
# note that 1/Vc scaling already in A and B
convert_2cmp_pkparms = function (CL, Vc, Q, Vp)
{
  k_sum = Q / Vc + Q / Vp + CL / Vc
  beta_parm = 0.5 * (k_sum - sqrt(k_sum ^ 2 - 4 * Q * CL / (Vp * Vc)))
  alpha_parm = Q * CL / (beta_parm * Vp * Vc)
  A = 1 / Vc * (alpha_parm - Q / Vp) / (alpha_parm - beta_parm)
  B = 1 / Vc * (beta_parm - Q / Vp) / (beta_parm - alpha_parm)
  data.frame(A, B, alpha_parm, beta_parm)
}

# parms = data.frame(A, B, alpha_parm, beta_parm)
bolus_model = function(mtime, parms, dose)
{
  with(parms, {
    dose * (A * exp(-alpha_parm * mtime) + (B * exp(-beta_parm * mtime)))
  })
}

simulate_multivariate = function(sim_id, times, parms1, parms2, dose1, dose2, constant_error){
  
  data.frame(
    ID = sim_id,
    TIME = times,
    DV1 = bolus_model(times, parms = parms1, dose = dose1) + rnorm(length(times), 0, constant_error),
    DV2 = bolus_model(times, parms = parms2, dose = dose2) + rnorm(length(times), 0, constant_error)
    ,stringsAsFactors = F)
  
}

# -------- simulation ----------

times = c(0, 0.05, 1:10/10, seq(2, 7) * 7)
pk_parms1 = convert_2cmp_pkparms(CL = 0.1, Vc = 1, Q = 1, Vp = 1)
pk_parms2 = convert_2cmp_pkparms(CL = 0.2, Vc = 2, Q = 1, Vp = 1)
dose1 = 1500
dose2 = 1000

library(purrr)

set.seed(100)
sim_data = map_df(1:25, simulate_multivariate, times = times, 
                         parms1 = pk_parms1, parms2 = pk_parms2,
                         dose1 = dose1, dose2 = dose2,
                         constant_error = 10
                         ) 


library(dplyr)
library(tidyr)
library(readr)

# reformat for mlx input
sim_data_long = sim_data %>%
  gather(DVID, DV, DV1, DV2) %>%
  mutate(
    DVID = parse_number(DVID),
    ID_UNIQUE = paste(DVID, ID, sep = "."),
    AMT = case_when(TIME == 0 & DVID == 1 ~ dose1,
                    TIME == 0 & DVID == 2 ~ dose2,
                    TRUE ~ 0)
  )

write_csv(sim_data_long, "sim_test_data.csv")


