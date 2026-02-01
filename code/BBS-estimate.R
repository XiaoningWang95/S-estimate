library(beta4s)
library(dplyr)
library(tibble)
library(reshape2)
library(data.table)
library(furrr)

routes=read.csv("data/routes.csv")
bird_data=read.csv("data/BBS-2023.csv")[,-1]
source('code/run_beta4s_bootstrap.R')

routes=routes[routes$RouteTypeID==1,]
routes=routes[routes$RouteTypeDetailID==1,]

routes$StateNum=as.numeric(routes$StateNum)
bird_data <- bird_data[bird_data$RPID==101,]

bird_total=apply(bird_data[,8:57],1,sum)
bird_abundance=data.frame(RouteDataID=bird_data$RouteDataID,state=bird_data$CountryNum,Route=bird_data$Route,
                          Year=bird_data$Year,AOU=bird_data$AOU,abundance=bird_total)

bird_abundance=bird_abundance[bird_abundance$abundance>0,]
bird_abundance=bird_abundance[bird_abundance$state>0,]

data2017 <- bird_abundance[bird_abundance$Year=='2017',]
data2018 <- bird_abundance[bird_abundance$Year=='2018',]
data2019 <- bird_abundance[bird_abundance$Year=='2019',]
rm(bird_abundance,bird_data,routes,bird_total)

sampe_size <- c(20,50,100,500,1000,1500,2000,3000)

plan(multisession,workers=10)

S_2017 <- future_map_dfr(sampe_size,~run_beta4s_bootstrap(data2017,sample_size = .x,n_boot = 30),.id = 'sampe_size')
save(S_2017,file = 'result/S_2017.Rdata')

S_2018 <- future_map_dfr(sampe_size,~run_beta4s_bootstrap(data2018,sample_size = .x,n_boot = 30),.id = 'sampe_size')
save(S_2018,file = 'result/S_2018.Rdata')

S_2019 <- future_map_dfr(sampe_size,~run_beta4s_bootstrap(data2019,sample_size = .x,n_boot = 30),.id = 'sampe_size')
save(S_2019,file = 'result/S_2019.Rdata')


## =========================================================
##  BBS abundance matrix → beta4s scale-up of regional S
## =========================================================

## ---------------------------------------------------------
## 0. BBS abundance matrix
##    rows   = sampling units (routes or route × year)
##    cols   = species
##    values = abundances
## ---------------------------------------------------------

# bbs_mat <- your abundance matrix (already prepared)

## ---------------------------------------------------------
## 1. Sampling units
## ---------------------------------------------------------

mat <- as.matrix(
  xtabs(
    abundance ~ RouteDataID + AOU,
    data = data2017
  )
)

bbs_mat <- sample_routes_mat(mat, x = 50)

bbs_mat <- bbs_mat[rowSums(bbs_mat) > 0, ]

M_smpl <- nrow(bbs_mat)

## ---------------------------------------------------------
## 2. Construct ab_occ object (required by beta4s)
## ---------------------------------------------------------

ab_occ_bbs <- data.frame(
  species = colnames(bbs_mat),
  n.ind   = colSums(bbs_mat),
  n.plots = colSums(bbs_mat > 0)
)

# Remove species never observed
ab_occ_bbs <- ab_occ_bbs[ab_occ_bbs$n.ind > 0, ]

## ---------------------------------------------------------
## 3. Regional constraints
## ---------------------------------------------------------

# Observed total individuals
J_obs <- sum(bbs_mat)

# Sampling effort (user-defined or externally estimated)
# Example: BBS covers ~5% of the region

# Regional total individuals
J_tot <-  7.5e8
effort <- J_obs/J_tot

# Regional total sampling units
M_tot <- M_smpl / effort

## ---------------------------------------------------------
## 4. Observed anchors
## ---------------------------------------------------------

# Beta diversity (primary anchor)
beta_obs <- betadiv(comp = t(bbs_mat),
  ab_occ = ab_occ_bbs,
  M      = M_smpl
)

# Spatial singletons (Q1) – recommended for NB model
Q1_obs <- sum(ab_occ_bbs$n.plots == 1)

# Simpson dominance (alternative anchor)
D_obs <- simpdom(ab_occ_bbs$n.ind)

## ---------------------------------------------------------
## 5. k–μ scaling parameter estimation
## ---------------------------------------------------------

kmu_fit <- kmuscaling(
  comp = t(bbs_mat),
  ab_occ = ab_occ_bbs,
  method = "abundance"
)

# General model (non-neutral)
c_gen <- kmu_fit$coef["c_General", "estimate"]
d_gen <- kmu_fit$coef["d_General", "estimate"]

# Neutral model (d fixed to 1)
c_neu <- kmu_fit$coef["c_Neutral", "estimate"]

## ---------------------------------------------------------
## 6. Scale-up using Logseries SAD
## ---------------------------------------------------------

# LS – Neutral
res_ls_neu <- upscaleS(
  beta  = beta_obs,
  J_tot = J_tot,
  M_tot = M_tot,
  c     = c_neu,
  d     = 1,
  SAD   = "ls"
)

# LS – General
res_ls_gen <- upscaleS(
  beta  = beta_obs,
  J_tot = J_tot,
  M_tot = M_tot,
  c     = c_gen,
  d     = d_gen,
  SAD   = "ls"
)

## ---------------------------------------------------------
## 7. Scale-up using Negative Binomial SAD
## ---------------------------------------------------------

# NB with Q1 anchor (recommended for BBS)
res_nb_q1 <- upscaleS(
  beta     = beta_obs,
  J_tot    = J_tot,
  M_tot    = M_tot,
  c        = c_gen,
  d        = d_gen,
  SAD      = "nb",
  Q1       = Q1_obs,
  M_sample = M_smpl
)

# NB with Simpson dominance (alternative)
res_nb_D <- upscaleS(
  beta  = beta_obs,
  J_tot = J_tot,
  M_tot = M_tot,
  c     = c_gen,
  d     = d_gen,
  SAD   = "nb",
  D     = D_obs
)

## ---------------------------------------------------------
## 8. Results
## ---------------------------------------------------------

cat("LS (Neutral) Estimated S:", round(res_ls_neu$S), "\n")
cat("LS (General) Estimated S:", round(res_ls_gen$S), "\n")
cat("NB (Q1 anchor) Estimated S:", round(res_nb_q1$S), "\n")
cat("NB (Q1) Estimated r:", round(res_nb_q1$params$r, 4), "\n")
cat("NB (D anchor) Estimated S:", round(res_nb_D$S), "\n")


## ---------------------------------------------------------
## 9. 汇总结果到一个表
## ---------------------------------------------------------




# 创建空数据框
results_df <- data.frame(
  Method = c("LS_Neutral", "LS_General", "NB_Q1", "NB_D"),
  Estimated_S = c(
    round(res_ls_neu$S),
    round(res_ls_gen$S),
    round(res_nb_q1$S),
    round(res_nb_D$S)
  ),
  Estimated_r = c(
    NA,                   # LS 方法没有 r
    NA,
    round(res_nb_q1$params$r, 4),
    NA                    # NB D 方法没有 r
  ),
  c_gen = c(
    NA,                   # LS 方法没有 r
    c_gen,
    NA,
    NA                    # NB D 方法没有 r
  ),
  d_gen = c(
    NA,                   # LS 方法没有 r
    d_gen,
    NA,
    NA                    # NB D 方法没有 r
  ),
  c_neu = c(
    c_neu,                   # LS 方法没有 r
    NA,
    NA,
    NA                    # NB D 方法没有 r
  )
)

# 查看结果

# 可选：写出 CSV
# write.csv(results_df, "beta4s_results_summary.csv", row.names = FALSE)

