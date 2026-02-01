run_beta4s_bootstrap <- function(
    bbs_mat,
    n_boot = 20,
    sample_size = nrow(bbs_mat),
    J_tot = 7.5e8,
    seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # ------------------------------
  # helper: sample routes
  # ------------------------------
  sample_routes_mat <- function(mat, x) {
    if (x > nrow(mat)) stop("sample_size > number of routes")
    mat[sample(seq_len(nrow(mat)), x, replace = FALSE), , drop = FALSE]
  }
  
  # ------------------------------
  # results container
  # ------------------------------
  res_list <- vector("list", n_boot)
  
  # ------------------------------
  # bootstrap loop
  # ------------------------------
  
  bbs_mat <- as.matrix(
    xtabs(
      abundance ~ RouteDataID + AOU,
      data = bbs_mat
    )
  )
  
  
  for (i in seq_len(n_boot)) {
    
    # ---- 1. resample routes ----
    mat_i <- sample_routes_mat(bbs_mat, sample_size)
    mat_i <- mat_i[rowSums(mat_i) > 0, ]
    
    M_smpl <- nrow(mat_i)
    if (M_smpl < 2) next
    
    # ---- 2. ab_occ ----
    ab_occ <- data.frame(
      species = colnames(mat_i),
      n.ind   = colSums(mat_i),
      n.plots = colSums(mat_i > 0)
    )
    ab_occ <- ab_occ[ab_occ$n.ind > 0, ]
    
    # ---- 3. regional constraints ----
    J_obs <- sum(mat_i)
    effort <- J_obs / J_tot
    M_tot <- M_smpl / effort
    
    # ---- 4. anchors ----
    beta_obs <- betadiv(
      comp = t(mat_i),
      ab_occ = ab_occ,
      M = M_smpl
    )
    
    Q1_obs <- sum(ab_occ$n.plots == 1)
    D_obs  <- simpdom(ab_occ$n.ind)
    
    # ---- 5. kmu scaling ----
    kmu_fit <- tryCatch(
      kmuscaling(comp = t(mat_i), ab_occ = ab_occ, method = "abundance"),
      error = function(e) NULL
    )
    if (is.null(kmu_fit)) next
    
    c_gen <- kmu_fit$coef["c_General", "estimate"]
    d_gen <- kmu_fit$coef["d_General", "estimate"]
    c_neu <- kmu_fit$coef["c_Neutral", "estimate"]
    
    # ---- 6. scale-up ----
    res_ls_neu <- upscaleS(
      beta = beta_obs, J_tot = J_tot, M_tot = M_tot,
      c = c_neu, d = 1, SAD = "ls"
    )
    
    res_ls_gen <- upscaleS(
      beta = beta_obs, J_tot = J_tot, M_tot = M_tot,
      c = c_gen, d = d_gen, SAD = "ls"
    )
    
    res_nb_q1 <- upscaleS(
      beta = beta_obs, J_tot = J_tot, M_tot = M_tot,
      c = c_gen, d = d_gen,
      SAD = "nb", Q1 = Q1_obs, M_sample = M_smpl
    )
    
    res_nb_D <- upscaleS(
      beta = beta_obs, J_tot = J_tot, M_tot = M_tot,
      c = c_gen, d = d_gen,
      SAD = "nb", D = D_obs
    )

  
    # ---- 7. store ----
    res_list[[i]] <- data.frame(
      iter = i,
      sample_size = sample_size,
      
      LS_Neutral_S = round(res_ls_neu$S),
      LS_General_S = round(res_ls_gen$S),
      NB_Q1_S      = round(res_nb_q1$S),
      NB_D_S       = round(res_nb_D$S),
      NB_Q1_r      = res_nb_q1$params$r,
      c_gen = c_gen,
      d_gen = d_gen,
      c_neu = c_neu
    )
  }
  do.call(rbind, res_list)
}
