evidence_isaric <- list(prev=list(type="beta", var=0.000872595702409054, mean=0.427966984132821),
                 cstat=list(type="beta", mean=0.760730571861002, var=3.82247885993882e-05),
                 cal_mean=list(type="norm", mean=-0.00934717199436785, var=0.0155046339663481),
                 cal_slp=list(type="norm", mean=0.995017759715243, var=0.000563011700688932))

evidence3 <- list(prev=list(type="beta", m=0.47, cih=0.65),
                  cstat=list(type="beta", m=0.75, sd=0.1),
                  cal_slp=list(type="norm", m=1, sd=0.1),
                  cal_mean=list(type="norm", m=0, sd=0.1)
                  )

targets_ss <- list(eciw.cstat=0.1, eciw.cal_oe=0.1, qciw.cal_oe=c(0.15, 0.9), assurance.nb=0.9)
targets_pow <- list(eciw.cstat=T, eciw.cal_oe=T, qciw.cal_oe=c(0.9), assurance.nb=T, voi.nb=T)

targets_riley <- list(eciw.cstat=0.1, eciw.cal_slp=0.15, qciw.cstat=c(0.1, 0.9), qciw.cal_slp=c(0.15,0.9), assurance.nb=0.9)


N<- c(50, 100, 200, 500, 1000, 2000)

n_sim=1000

threshold=0.2
dist_type="logitnorm"
impute_cor <- T
method="sample"

library(bayescpm)
set.seed(1)
res <- bpm_valsamp(evidence=evidence_isaric, dist_type=dist_type, method=method, targets=targets_ss, n_sim=n_sim, impute_cor=impute_cor, threshold=threshold)

res_riley <- bpm_valsamp(evidence=evidence_isaric, dist_type=dist_type, method=method, targets=targets_riley, n_sim=n_sim, impute_cor=impute_cor, threshold=threshold)


