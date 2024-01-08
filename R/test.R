# LCGA - MSM
 obsdata_msm <- genobsdatamsm(n=5000,t=10, seed = 345)
 V <- "V"
 L <- "L"
 K = 10
 time <- 1:10
 A <- "A"
 Y <- "Y"
 Adat_long <- widetolong(obsdata = obsdata_msm[, paste0(A,time)], varying = 1:10)
 res.traj = buildtraj(Rdat = Adat_long, J=3,formula = cbind(A,1-A) ~ time, id= "id")
 table(res.traj$dpost$class)
 obsdata_msm2 = merge(obsdata_msm,res.traj$dpost, by = "id")
 obsdata_msm2[,"class"] <- relevel(as.factor(obsdata_msm2[,"class"]), ref = 3)

 #IPW
 trajMSM_IPW(formula1 = as.formula("Y ~ class"), numerator = "stabilized",
     id = id, V = V,L = L,A = A,Y=Y,J=3,
           C = FALSE,K = K, time = time, family = "binomial",
          obsdata = obsdata_msm2)

 #g-computation


trajMSM_gform(formula = paste0("Y~", paste0("A", 1:K,collapse = "+"), "+",
                               paste0("L", 1:K,collapse = "+"),"+",  V, collapse = "+"),rep = 50,
id = id, V = V,L = L,A = A,Y=Y,K=K,time = time,timevar= "time",J=3,
 trajmodel = res.traj$model, ref = "2", obsdata = obsdata)


#Pooled-ltmle

trajMSM_pltmle(formula = paste0("Y~", paste0("A", 1:K,collapse = "+"), "+",
                              paste0("L", 1:K,collapse = "+"),"+",
                            V, collapse = "+"),
             id = id, V = V,L = L,A = A,Y=Y,K = K,time = time,timevar = "time",J=3,
             trajmodel = res.traj$model, obsdata = obsdata, ref = "2")
# LCGA - HRMSM

obsdata =  gendatahrmsm(5000)
res_IPW = trajHRMSM_IPW(formulaY = as.formula("Y ~ factor(traj) + factor(Interv)"),
                       numerator = "stabilized", degree_traj = "linear",
                       idvar = "id", V = V,L = L,A = A,Y=Y,
                       C = FALSE,s = 3, K = 5, J=3, family = "poisson",
                       obsdata = obsdata, v.names = v.names)

res_gform = trajHRMSM_gform(obsdata = obsdata,Rep=50,degree_traj = "linear",
                           A = A,L = L,V=V,Y = Y, s = 3,K = 5, timevar = timevar, idvar = "id",
                           J = 3, family = "poisson")

res_pltmle = trajHRMSM_pltmle(obsdata = obsdata,degree_traj = "linear",
                              A = A,L = L,V=V,Y = Y, s = 3,K = 5,
                              timevar = timevar, idvar = idvar,
                              J = 3, family = "poisson")
res_IPW[[1]]
res_gform[[1]]
res_pltmle[[1]]
