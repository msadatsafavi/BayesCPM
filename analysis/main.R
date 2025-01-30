evidence2 <- list(prev=list(type="beta", var=0.000872595702409054, mean=0.427966984132821),
                 cstat=list(type="beta", mean=0.760730571861002, var=3.82247885993882e-05),
                 cal_mean=list(type="norm", mean=-0.00934717199436785, var=0.0155046339663481),
                 cal_slp=list(type="norm", mean=0.995017759715243, var=0.000563011700688932))

evidence3 <- list(prev=list(type="beta", m=0.47, cih=0.65),
                  cstat=list(type="beta", m=0.75, sd=0.1),
                  cal_slp=list(type="norm", m=1, sd=0.1),
                  cal_mean=list(type="norm", m=0, sd=0.1)
                  )

N<- c(50, 100, 200, 500, 1000, 2000)
#N <- 350
n_sim=1000

target_ciws = list(cstat=0.1, cal_oe=0.22, cal_slp=0.3)
rules=list(fciw=T, eciw=T, qciw=0.9, nb_voi=T, nb_assurance=T)
threshold=0.2
dist_type="logitnorm"
impute_cor <- T
method="sample"

library(bayescpm)
set.seed(1)
res <- BayesCPM(N, evidence2, dist_type=dist_type, method=method, target_ciws=target_ciws, rules=rules, n_sim=n_sim, impute_cor=impute_cor, threshold=threshold)


