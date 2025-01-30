#' @import fastLogisticRegressionWrap
#' @import mcmapper
#' @import pROC
# #' @import evsiexval
#' @import mc2d
#' @import logitnorm
#' @import cobs
#' @import quantreg

#' @importFrom stats binomial coef cor dbeta glm integrate optim pbeta pnorm predict qnorm quantile rbeta rbinom rnorm t.test uniroot



logit <- function(x) {log(x/(1-x))}

expit <- function(x) {1/(1+exp(-x))}




#'@export
moments <- function(type, parms) {
  # Check if the input values are valid
  if(type=="norm")
  {
    return(c(m=parms[1], v=sqrt(parms[2])))
  }
  if(type=="beta")
  {
    a <- parms[1]; b <- parms[2]
    return(c(m=a/(a+b), v=a*b/(((a+b)^2)*(a+b+1))))
  }
  if(type=="logitnorm")
  {
    return(momentsLogitnorm(parms[[1]], parms[[2]]))
  }
  if(type=="probitnorm")
  {
    #TODO
  }
}





#'@export
inv_moments <- function(type, moments) {
  # Check if the input values are valid
  m <- moments[[1]]
  v <- moments[[2]]
  
  if(type=="norm")
  {
    sd <- sqrt(v)
    return(c(mean=m, sd=sd))
  }
  
  if(type=="logitnorm")
  {
    if (m <= 0 || m >= 1) {
      stop("Mean must be between 0 and 1 (exclusive).")
    }
    if (v <= 0) {
      stop("Variance must be positive.")
    }
    
    objective <- function(params) {
      mu <- params[1]
      sigma <- params[2]
      
      tmp <- logitnorm::momentsLogitnorm(mu, sigma)
      mean_computed <- tmp[1]
      var_computed <- tmp[2]
      
      # Sum of squared differences between target and computed values
      (mean_computed - m)^2 + (var_computed - v)^2
    }
    
    # Initial guesses for mu and sigma
    initial_guess <- c(0, 1)
    
    # Use optim() to minimize the objective function
    result <- optim(initial_guess, objective, method = "L-BFGS-B", lower = c(-Inf, 1e-6), upper = c(Inf, Inf))
    
    # Extract mu and sigma
    mu <- result$par[1]
    sigma <- result$par[2]
    
    return(c(mu = mu, sigma = sigma))
  }
  
  if(type=="beta")
  {
    if (m <= 0 || m >= 1) {
      stop("Mean must be between 0 and 1 (exclusive).")
    }
    if (v <= 0) {
      stop("Variance must be positive.")
    }
    
    # Calculate alpha and beta
    alpha <- m * ((m * (1 - m) / v) - 1)
    beta <- (1 - m) * ((m * (1 - m) / v) - 1)
    
    # Return the alpha and beta as a named list
    return(c(shape1=alpha, shape2=beta))
  }
}

#'@export
inv_mean_quantile <- function(type, m, q, p)
{
  if(type=="logitnorm")
  {
    res <- tryCatch(
      {logitnorm::twCoefLogitnormE(mean=m, quant=q, perc=p)},
      error=(function(cond) {logitnorm::twCoefLogitnormE(mean=m, quant=q, perc=p,theta0=logitnorm::twCoefLogitnorm(m,q,p))}))
    
    out <- c(mu=res[1], sigma=res[2])
  }
  if(type=="beta") 
  {
    res <- uniroot(function(x) {pbeta(q,x,x*(1-m)/m)-p}, interval=c(0.0001,10000))
    out <- c(shape1=res$root, shape2=res$root*(1-m)/m)
  }
  if(type=="probitnorm") #TODO
  {
    
  }
  if(type=="norm")
  {
    res <- uniroot(function(x) {pnorm(q,m,x)-p}, interval=c(0.0001,10))
    out <- c(mean=m, sd=res$root)
  }
  
  out
}




# #@export
# ENBS <- function(parms)
# {
#   if(parms$evidence_type=="ind_beta")
#   {
#     VoI <- EVSI_ag(parms$evidence, parms$z, n_sim=parms$n_sim_inner, future_sample_sizes=parms$n_stars)
#   }
#   if(parms$evidence_type=="prev_cs_A_B")
#   {
#     VoI <- EVSI_gf(parms$sample[,c('prev','se','sp')], parms$z, n_sim=parms$n_sim_inner, future_sample_sizes=parms$n_stars)
#   }
#   if(params$evidence_type=="sample")
#   {
#     VoI <- EVSI_gf(global_vars$params$evidence[,c('prev','se','sp')], global_vars$params$z, n_sim=global_vars$params$n_sim_inner, future_sample_sizes=global_vars$params$n_stars)
#   }
#
#   EVPI <- VoI$EVPI
#   EVSIs <- c(VoI$EVSI)
#   ENBS <- EVSIs*parms$N-parms$n_stars/parms$lambda
#
#   data.frame("n_star"=as.integer(parms$n_stars),
#              "EVPI"=EVPI,
#              "EVSI"=EVSIs,
#              "Population EVSI"=as.double(EVSIs*parms$N),
#              "ENBS"=as.double(ENBS)
#   )
# }




# export
# find_n_star <- function(parms, n_iter=200)
# {
#   n_min <- min(parms$n_stars)
#   n_max <- max(parms$n_stars)
#
#   require(OOR)
#   res <- OOR::StoSOO(fn=function(x){
#       parms$n_stars<-round(n_min+x*(n_max-n_min),0);
#       ENBS(parms)$ENBS
#     },
#     par=c(0.5),
#     nb_iter=n_iter,
#     control=list(verbose=1, max=T))
#   round(n_min+res$par*(n_max-n_min),0);
# }



#'@export
infer_cal_int_from_mean <- function(dist_type, dist_parms, cal_mean, cal_slp, prev=NULL)
{
  if(dist_type=="logitnorm")
  {
    f <- function(intercept)
    {
      integrate(function(x, intercept, slope){expit((logit(x)-intercept)/cal_slp)*mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])}, 0, 1,intercept=intercept, slope=cal_slp)$value
    }
    res <- optim(0, fn=function(x) {((prev-f(x))-cal_mean)^2}, lower=-2, upper=2, method="Brent")
  }

  if(dist_type=="beta")
  {
    f <- function(intercept)
    {
      integrate(function(x, intercept, slope){expit((logit(x)-intercept)/cal_slp)*dbeta(x,dist_parms[1],dist_parms[2])}, 0, 1,intercept=intercept, slope=cal_slp)$value
    }
    res <- optim(0, fn=function(x) {((prev-f(x))-cal_mean)^2}, lower=-2, upper=2, method="Brent")
  }

  if(dist_type=="probitnorm")
  {
    f <- function(intercept)
    {
      integrate(function(x, intercept, slope){expit((logit(x)-intercept)/cal_slp)*mcmapper::dprobitnorm(x,dist_parms[1],dist_parms[2])}, 0, 1,intercept=intercept, slope=cal_slp)$value
    }
    res <- optim(0, fn=function(x) {((prev-f(x))-cal_mean)^2}, lower=-4, upper=4, method="Brent")
  }

  unname(res$par)
}





#'@export
calc_se_sp <- function(dist_type, dist_parms, cal_int, cal_slp, threshold, prev=NULL)
{

  hz <- expit(logit(threshold)*cal_slp+cal_int)

  if(dist_type=="logitnorm")
  {
    tp <- integrate(f=function(x) {x*mcmapper::dlogitnorm(x, dist_parms[1], dist_parms[2])}, hz, 1)$value
    se <- unname(tp/prev)
    sp <- unname((mcmapper::plogitnorm(hz, dist_parms[1], dist_parms[2])-(prev-tp))/(1-prev))
  }

  if(dist_type=="beta")
  {
    tp <- integrate(f=function(x) {x*dbeta(x, dist_parms[1], dist_parms[2])}, hz, 1)$value
    se <- unname(tp/prev)
    sp <- unname((pbeta(hz, dist_parms[1], dist_parms[2])-(prev-tp))/(1-prev))
  }

  if(dist_type=="probitnorm")
  {
    tp <- integrate(f=function(x) {x*mcmapper::dprobitnorm(x, dist_parms[1], dist_parms[2])}, hz, 1)$value
    se <- unname(tp/prev)
    sp <- unname((mcmapper::pprobitnorm(hz, dist_parms[1], dist_parms[2])-(prev-tp))/(1-prev))
  }

  return(c(se=max(0,min(se,1)), sp=max(0,min(sp,1))))
}










#'@export
infer_correlation <- function(dist_type, dist_parms, cal_int, cal_slp, n, n_sim)
{
  require(pROC)
  require(mcmapper)

  out <- matrix(NA, nrow=n_sim, ncol=5)
  colnames(out) <- c("prev", "cstat", "cal_mean", "cal_int", "cal_slp")

  for(i in 1:n_sim)
  {
    p <- do.call(paste0("r",dist_type), args=as.list(c(n,unname(dist_parms))))
    Y <- rbinom(n,1,p)
    pi <- expit((logit(p)-cal_int)/cal_slp)
    df <- cbind(pi=pi,Y=Y)
    out[i,]<-c(mean(df[2]),
               roc(df[,2]~df[,1], quiet=TRUE)$auc,
               mean(df[2]-df[1]),
               coef(glm(df[,2]~logit(df[,1]), family=binomial(link="logit")))[1:2])
  }

  cor(out)
}




















#N can be vector!
calc_vars <- function(N, parms)
{
  prev <- parms$prev
  C <- parms$cstat
  dist_type <- parms$dist_type
  dist_parms <- c(parms$dist_parm1,parms$dist_parm2)
  cal_int <- parms$cal_int
  cal_slp <- parms$cal_slp

  if(dist_type=="logitnorm")
  {
    #tryCatch({
      E_pi <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(expit((logit(x)-cal_int)/cal_slp))},0,1)$value
      I_a <- integrate(function(x){ mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(x*(1-x))}, 0, 1)$value
      I_b <- integrate(function(x){ mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)^2*(x*(1-x)))}, 0, 1)$value
      I_ab <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)*(x*(1-x)))}, 0, 1)$value
    #}, error=function(e){bad_place$parms <<- parms})
  }
  if(dist_type=="beta")
  {
    E_pi <- integrate(function(x){x*dbeta(expit(logit(x)*cal_slp+cal_int),dist_parms[1],dist_parms[2])}, 0, 1,intercept=cal_int, slope=cal_slp)$value
    I_a <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*(exp(cal_int+x*cal_slp)/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
    I_b <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*((x^2*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
    I_ab <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*((x*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
  }
  if(dist_type=="probitnorm")
  {
    E_pi <- integrate(function(x){x*mcmapper::dprobitnorm(expit(logit(x)*cal_slp+cal_int),dist_parms[1],dist_parms[2])}, 0, 1,intercept=cal_int, slope=cal_slp)$value
    I_a <- N*integrate(function(x){mcmapper::dprobitnorm(expit(logit(x)), dist_parms[1], dist_parms[2])*(exp(cal_int+x*cal_slp)/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
    I_b <- N*integrate(function(x){mcmapper::dprobitnorm(expit(logit(x)), dist_parms[1], dist_parms[2])*((x^2*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
    I_ab <- N*integrate(function(x){mcmapper::dprobitnorm(expit(logit(x)), dist_parms[1], dist_parms[2])*((x*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
  }

  v_prev <- prev*(1-prev)/N
  v_cstat <- (C*(1-C)*(1+(N/2-1)*((1-C)/(2-C))+(N/2-1)*C/(1+C)))/(N*N*prev*(1-prev))
  v_cal_mean <- v_prev #TODO
  v_cal_oe <- (1-prev)/(N*E_pi)
  det <- N*(I_a*I_b-I_ab^2)
  v_cal_int <- I_b/det
  v_cal_slp <- I_a/det
  cov_int_slp <- -I_ab/det

  list(prev=v_prev,cstat=v_cstat, cal_mean=v_cal_mean, cal_oe=v_cal_oe,  cal_int=v_cal_int, cal_slp=v_cal_slp, cov_int_slp=cov_int_slp)
}




#N can be vector
#DEPRECATED
# calc_ciw_1s <- function(N, parms)
# {
#   vars <- calc_vars(N, parms)
#   list(oe=sqrt(vars$oe), cstat=sqrt(vars$cstat), cal_slp=sqrt(vars$cal_slp))
# }


#N can be vector
#'@export
calc_ciw_2s <- function(N, parms)
{
  ciw_cstat <- ciw_cal_oe <- ciw_cal_mean <- ciw_cal_int <- ciw_cal_slp <- rep(NA, length(N))

  l <- length(N)

  s1 <- calc_vars(N, parms=parms)
  
  hat_cals <- rbnorm(l,
                rep(parms$cal_int,l),
                rep(parms$cal_slp,l),
                s1$cal_int,
                s1$cal_slp,
                s1$cov_int_slp)

  k <- 2*qnorm(0.975)
  
  for(i in 1:l)
  {
    hat_prev <- inv_moments("beta",c(parms$prev, s1$prev[i]))
    hat_cstat <- inv_moments("beta",c(parms$cstat, s1$cstat[i]))

    f <- function()
    {
      repeat{
        x<-rbeta(1, hat_cstat[1], hat_cstat[2])
        if(x>=0.51 & x<0.99) break;
        }
      x
    }

    new_parms <- list(prev=rbeta(1, hat_prev[1], hat_prev[2]),
                    cstat=f(),
                    cal_int=hat_cals[i,1],
                    cal_slp=hat_cals[i,2])


    dist_parms <- mcmapper::mcmap_logitnorm(c(new_parms$prev, new_parms$cstat)) #TODO


    s2_vars <- calc_vars(N[i], parms=list(prev=new_parms$prev,
                                        cstat=new_parms$cstat,
                                        cal_int=new_parms$cal_int,
                                        cal_slp=new_parms$cal_slp,
                                        dist_type="logitnorm",
                                        dist_parm1=dist_parms[1],
                                        dist_parm2=dist_parms[2]
                                        )
                        )
    ciw_cstat[i] <- k*sqrt(s2_vars$cstat)
    ciw_cal_oe[i] <- k*sqrt(s2_vars$cal_oe)
    ciw_cal_mean[i] <- k*sqrt(s2_vars$cal_mean)
    ciw_cal_int[i] <- k*sqrt(s2_vars$cal_int)
    ciw_cal_slp[i] <- k*sqrt(s2_vars$cal_slp)
  }

  if(is.infinite(sum(ciw_cal_oe))) {browser()}

  list(cstat=ciw_cstat,
       cal_oe=ciw_cal_oe,
       cal_mean=ciw_cal_mean,
       cal_int=ciw_cal_int,
       cal_slp=ciw_cal_slp)
}



#N can be vector
#'@export
calc_ciw_sample <- function(N, parms)
{
  require(fastLogisticRegressionWrap)
  se_cstat <- se_cal_oe <- se_cal_mean <- se_cal_int <- se_cal_slp <- rep(NA, length(N))

  prev <- parms$prev
  C <- parms$cstat
  dist_type <- parms$dist_type
  dist_parms <- c(parms$dist_parm1, parms$dist_parm2)
  cal_int <- parms$cal_int
  cal_slp <- parms$cal_slp

  N_max <- max(N)
  p <- do.call(paste0("r",dist_type), as.list(c(N_max, dist_parms)))
  pi <- expit((logit(p)-cal_int)/cal_slp)
  #logit_pi <- logit(pi)
  logit_pi <- cbind(1,logit(pi))
  Y <- rbinom(N_max, 1, p)

  for(i in 1:length(N))
  {
    n <- N[i]
    O <- sum(Y[1:n])
    se_cal_oe[i] <- sqrt((1-O/n)/O)
    se_cal_mean[i] <- t.test(Y[1:n]-pi[1:n])$stderr
    # reg <- glm(Y[1:n]~logit_pi[1:n], family = binomial())
    # se_cal_int[i] <- sqrt(vcov(reg)[1,1])
    # se_cal_slp[i] <- sqrt(vcov(reg)[2,2])
    reg <- fast_logistic_regression(logit_pi[1:n,],Y[1:n],do_inference_on_var="all")
    se_cal_int[i] <- reg$se[1]
    se_cal_slp[i] <- reg$se[2]

    #tmp <- ci.auc(Y[1:n]~logit_pi[1:n], quiet=TRUE)
    tmp <- ci.auc(Y[1:n]~logit_pi[1:n,2], quiet=TRUE)
    se_cstat[i] <- (tmp[3]-tmp[2])/qnorm(0.975)
  }

  k <- 2*qnorm(0.975)
  list(cstat=k*se_cstat, cal_oe=k*se_cal_oe, cal_mean=k*se_cal_mean, cal_int=k*se_cal_int, cal_slp=k*se_cal_slp)
}







#'@export
calc_ciw_mc <- function(N, parms_sample, method)
{
  n <- nrow(parms_sample)

  ciw_cstat <- ciw_cal_oe <- ciw_cal_mean <- ciw_cal_int <- ciw_cal_slp <- matrix(NA, nrow=n, ncol=length(N))

  if(method=="sample")
  {
    for(i in 1:n)
    {
      tmp <- calc_ciw_sample(N, parms=parms_sample[i,])
      ciw_cstat[i, ] <- tmp$cstat
      ciw_cal_oe[i, ] <- tmp$cal_oe
      ciw_cal_mean[i, ] <- tmp$cal_mean
      ciw_cal_int[i, ] <- tmp$cal_int
      ciw_cal_slp[i, ] <- tmp$cal_slp
    }
  } else if(method=="2s")
  {
    for(i in 1:n)
    {
      tmp <- calc_ciw_2s(N, parms=parms_sample[i,])
      ciw_cstat[i, ] <- tmp$cstat
      ciw_cal_oe[i, ] <- tmp$cal_oe
      ciw_cal_mean[i, ] <- tmp$cal_mean
      ciw_cal_int[i, ] <- tmp$cal_int
      ciw_cal_slp[i, ] <- tmp$cal_slp
    }
  }

  list(cstat=ciw_cstat, cal_oe=ciw_cal_oe, cal_mean=ciw_cal_mean, cal_int=ciw_cal_int, cal_slp=ciw_cal_slp)
}


#Convext monotonically decreasing root finder
find_n_mean <- function(target, N, ciws, decreasing=T, convex=T)
{
  require(cobs)
  
  X <- rep(N, each=nrow(ciws))
  Y <- as.vector(ciws)
  
  constraint <- c(ifelse(decreasing,"decreasing",""),ifelse(convex,"convex",""))
  
  S <- cobs(X,Y, constraint=c("decrease","convex"))
  #plot(predict(S, seq(min(N),max(N),length.out=100)), type='l')
  
  f <- function(x) {predict(S, z=x)[,2] - target}
  round(uniroot(f, c(min(N),max(N)))$root,0)
}



#Convext monotonically decreasing root finder
find_n_quantile <- function(target, N, q, ciws)
{
  require(quantreg)
  
  X <- rep(N, each=nrow(ciws))
  Y <- as.vector(ciws)
  
  
  fit <- rqss(Y ~ qss(X, lambda = 1), tau=q)
  
  # plot(X, Y, main = "Non-parametric Quantile Regression", xlab = "X", ylab = "Y")
  # for (i in 1:length(tau)) {
  #   lines(sort(X), predict(fit, newdata = data.frame(X = sort(X)), tau = tau[i]), col = i + 1)
  # }
  # legend("topleft", legend = paste("tau =", tau), col = 2:(length(tau) + 1), lty = 1)
  
  f <- function(x) {predict(fit, newdata=data.frame(X=x)) - target}
  round(uniroot(f, c(min(N),max(N)))$root,0)
}






rbnorm <- function(n, mu1, mu2, var1, var2, cov) {

  # Calculate standard deviations and correlation from variances and covariance
  sigma1 <- sqrt(var1)
  sigma2 <- sqrt(var2)
  rho <- cov / (sigma1 * sigma2)

  x1 <- rnorm(n, mu1, sigma1)

  conditional_mean <- mu2 + rho * (sigma2 / sigma1) * (x1 - mu1)
  conditional_sd <- sigma2 * sqrt(1 - rho^2)

  x2 <- rnorm(n, conditional_mean, conditional_sd)

  return(cbind(x1, x2))
}



















# Creates a clean list including type, parms, moments
process_evidence_element <- function(element)
{
  e <- list()
  e$type <- element$type
  
  if(any(nchar(names(element))==0)) {stop("Unnamed objects found in ...")}
  ##Should have any of the following members: (mu, var), (mu, sd), (alpha, beta)
  nms <- names(element)
  possible_args <- list(c("mean","var"), c("m","v"), c("mean","sd"), c("m","sd"), c("alpha","beta"), c("mean","cih"), c("m","cih"))
  renamed_args <-  list(c("m","v"),   c("m","v"), c("m","s"), c("m","s"), c("shape1","shape2"), c("m","cih"), c("m","cih"))
  parms <- c()
  for(i in 1:length(possible_args))
  {
    res <- match(possible_args[[i]], nms)
    if(!any(is.na(res)))
    {
      if(length(parms)>0) stop("Multiple arguments matched")
      parms <- element[res]
      names(parms) <- renamed_args[[i]]
    }
  }
  if(length(parms)==0) stop("No valid parameter specification")
  ##IF (m,s), apply the method of moments. Note that for normal it is a stupid tail chasing but OK!
  if(names(parms)[1]=='m')
  {
    m <- parms[[1]]
    if(names(parms)[2]=="cih")
    {
      e$parms <- inv_mean_quantile(element$type, m, q=parms[[2]], p=0.975)
      e$moments <- c(m=m, v=0)
      e$moments[2] <- moments(element$type, e$parms)[[2]]
    }
    else
    {
      v <- parms[[2]]
      if(names(parms)[2]=="s") v <- v^2
      e$moments <- c(m=m,v=v)
      e$parms <- inv_moments(element$type, list(m,v))
    }
  }
  else
  {
    e$parms <- parms
    e$moments <- moments(element$type, parms)
  }

  e
}


#' @export
process_evidence <- function(evidence)
{
  if(is.null(evidence$prev)) {stop("evidence object must have a prev (prevalence) member")}
  if(is.null(evidence$prev$type)) {evidence$prev$type<-"beta"; message("Assuming prev has a beta distribution") }
  evidence$prev <- process_evidence_element(evidence$prev)
  
  if(is.null(evidence$cstat)) {stop("evidence object must have a cstat (c-statistic) member")}
  if(is.null(evidence$cstat$type)) {evidence$cstat$type<-"beta"; message("Assuming cstat has a beta distribution") }
  evidence$cstat <- process_evidence_element(evidence$cstat)
  
  nms <- names(evidence)
  possible_args <- list(c("cal_mean","cal_slp"), c("cal_oe","cal_slp"), c("cal_int","cal_slp"))
  renamed_args <-  list(c("cal_mean","cal_slp"), c("cal_oe","cal_slp"), c("cal_int","cal_slp"))
  cal_parms <- c()
  for(i in 1:length(possible_args))
  {
    res <- match(possible_args[[i]], nms)
    if(!any(is.na(res)))
    {
      if(length(cal_parms)>0) stop("Multiple arguments matched")
      cal_parms <- evidence[res]
      names(cal_parms) <- renamed_args[[i]]
    }
  }
  if(length(cal_parms)==0) stop("No valid parameter specification for calibration (calibration slope AND at lease one of intercept, mean calibration, or O/E ratio are required)")
  
  cal_mean <- match(c("cal_mean"),names(cal_parms))
  cal_int <- match(c("cal_int"),names(cal_parms))
  cal_slp <- match(c("cal_slp"),names(cal_parms))
  cal_oe <- match(c("cal_oe"),names(cal_parms))
  
  if(!is.na(cal_mean))
  {
    if(is.na(cal_parms[[cal_mean]]$type)) { cal_parms[[cal_mean]]$type<-"normal"; message("Assuming normal distirbution for calibration mean")}
    evidence$cal_mean <- process_evidence_element(cal_parms[[cal_mean]])
  }
  if(!is.na(cal_int))
  {
    if(is.na(cal_parms[[cal_int]]$type)) { cal_parms[[cal_int]]$type<-"normal"; message("Assuming normal distirbution for calibration mean")}
    evidence$cal_int <- process_evidence_element(cal_parms[[cal_int]])
  }
  if(!is.na(cal_slp))
  {
    if(is.na(cal_parms[[cal_slp]]$type)) { cal_parms[[cal_slp]]$type<-"normal"; message("Assuming normal distirbution for calibration mean")}
    evidence$cal_slp <- process_evidence_element(cal_parms[[cal_slp]])
  }
  if(!is.na(cal_oe))
  {
    if(is.na(cal_parms[[cal_oe]]$type)) { cal_parms[[cal_oe]]$type<-"normal"; message("Assuming normal distirbution for calibration mean")}
    evidence$cal_oe <- process_evidence_element(cal_parms[[cal_oe]])
  }
  
  evidence
}



#' @export
BayesCPM <- function(N, evidence, dist_type="logitnorm", method="sample", target_ciws, rules=list(fciw=FALSE, eciw=TRUE, qciw=0.9, nb_voi=TRUE, nb_assurance=0.9), n_sim, impute_cor=FALSE, threshold=NULL, ex_args=NULL)
{
  out <- list(N=N)
  
  if(is.function(ex_args$f_progress))
  {
    f_progress <- ex_args$f_progress
  }else
  {
    f_progress <- base::message
  }
  
  ##Step 1: Process evidence
  f_progress("Processing evidence...")
  evidence <- process_evidence(evidence)
  out$evidence <- evidence
  
  ##Generate base parms
  base <- list()
  base$dist_type <- dist_type
  base$prev <- evidence$prev$moments[[1]]
  base$cstat <- evidence$cstat$moments[[1]]
  tmp <- mcmapper::mcmap(c(base$prev, base$cstat), type=dist_type)$valu
  base$dist_parm1 <- tmp[1]; base$dist_parm2 <- tmp[2]
  base$cal_slp <- evidence$cal_slp$moments[[1]]
    
  #Cal intercept is not provided. Need to derive it for base params for correlation induction
  if(is.na(match("cal_int", names(evidence))))
  {
    base$cal_int <- infer_cal_int_from_mean(base$dist_type, c(base$dist_parm1, base$dist_parm2), cal_mean=evidence$cal_mean$moments[[1]], cal_slp=base$cal_slp, prev=base$prev)
  }
  else
  {
    base$cal_int <- evidence$cal_int$moments[[1]]
  }
  
    

  #Step 2: generate sample of marginals
  f_progress("Generating Monte Carlo sample...")
  sample <- NULL
  for(element in evidence)
  {
    sample <- cbind(sample, do.call(paste0("r",element$type), args=as.list(c(n=n_sim,element$parms))))
  }
  colnames(sample) <- names(evidence)
  
  #TODO: replace bad c-statistic values
  n_bads <- 0
  repeat{
    bads <- which(sample[,'cstat']<0.51 |  sample[,'cstat']>0.99)
    if(length(bads)==0) break;
    n_bads <- n_bads+length(bads)
    subsample <- NULL
    for(element in evidence)
    {
      subsample <- cbind(subsample, do.call(paste0("r",element$type), args=as.list(c(n=length(bads),element$parms))))
    }
    sample[bads,] <- subsample
  }
  if(n_bads>0) warning(paste("in step 'Generating MOnte Carlo sample' - ", n_bads, "observations were replaced due to bad value of c-statistic."))
  
  #Step 3: induce correlation (if asked)
  if(impute_cor)
  {
    f_progress("Imputing correlation...")
    eff_n <- round(evidence$prev$moments[[1]]*(1-evidence$prev$moments[[1]])/evidence$prev$moments[[2]],0)
    f_progress(paste("Based on effective sample size:", eff_n))
    base_cor <- infer_correlation(base$dist_type, c(base$dist_parm1,base$dist_parm2), base$cal_int, base$cal_slp, eff_n, 1000)
    good_rows <- match(colnames(sample), colnames(base_cor))
    base_cor <- base_cor[good_rows,good_rows]
    sample <- mc2d::cornode(sample, target=base_cor)
  }


  #Step 4: if intercept is missing, impute it for the whole sample
  f_progress("Infering calibration intercept...")
  sample <- as.data.frame(sample)
  sample$dist_type <- dist_type
  sample$dist_parm1 <- 0
  sample$dist_parm2 <- 0

  if(is.na(match("cal_int",names(evidence))))
  {
    sample$cal_int <- NA
    for(i in 1:nrow(sample))
    {
      prev <- unname(sample[i,'prev'])
      cstat<- unname(sample[i,'cstat'])
      cal_mean <- unname(sample[i,'cal_mean']) #TODO
      cal_slp <- unname(sample[i,'cal_slp'])

      parms <- mcmap(c(prev, cstat), dist_type)$value
      sample$dist_parm1[i] <- parms[1]
      sample$dist_parm2[i] <- parms[2]

      cal_int <- infer_cal_int_from_mean(dist_type=dist_type, dist_parms=parms, cal_mean=cal_mean, cal_slp=cal_slp, prev=prev)

      sample[i,'cal_int'] <- cal_int
    }
  }


  # Step 5: Bayesian Riley
  f_progress("Computing CI widths...")
  nms <- names(target_ciws)
  
  if(!is.na(match("fciw",names(rules)))) #Frequentist CIWs
  {
    fv <- calc_vars(N, parms=base)
  }
  
  ciws <- calc_ciw_mc(N, sample, method=method)
  out$ciws <- ciws
  
  
  for(i in 1:length(target_ciws))
  {
    ciw <- target_ciws[[i]]
    nm <- nms[i]

    if(!is.na(match("fciw",names(rules)))) 
    {
      out$fciw[[nm]] <- sqrt(fv[[nm]])*2*qnorm(0.975)
    }
    if(!is.na(match("eciw",names(rules))))
    {
      out$eciw[[nm]] <- colMeans(ciws[[nm]])
    }
    if(!is.na(match("qciw",names(rules))))
    {
      out$qciw[[nm]] <- apply(ciws[[nm]],2,FUN=quantile, rules$qciw)
    }
    if(!is.na(match("assurance",names(rules))))
    {
      out$assurance[[nm]] <- colMeans(ciws[[nm]]<ciw)
    }
  }





  # Step 6: Calculate se and sp
  f_progress("Computing se/sp...")
  sample$sp <- sample$se <- NA
  for(i in 1:nrow(sample))
  {
    se_sp <- calc_se_sp(sample$dist_type[i],
                               c(sample$dist_parm1[i], sample$dist_parm2[i]),
                               sample$cal_int[i],
                               sample$cal_slp[i],
                               threshold,
                               sample$prev[i]
    )
    sample[i,c('se','sp')] <- se_sp
  }


  #Step 7: VoI
  b_voi <- !is.na(match("nb_voi",names(rules)))
  b_assurance <- !is.na(match("nb_assurance",names(rules)))
  if(b_voi | b_assurance)
  {
    f_progress("VoI and NB assuraance...")
    
    #require(evsiexval)
    #res <- evsiexval::EVSI_gf(sample[,c('prev','se','sp')], future_sample_sizes=N,  ignore_prior=TRUE, z=threshold)
    
    tnb1 <- sample$prev*sample$se - (1-sample$prev)*(1-sample$sp)*threshold/(1-threshold)
    tnb2 <- sample$prev - (1-sample$prev)*threshold/(1-threshold)
    tnbs <- cbind(0,tnb1,tnb2)
    maxnbs <- apply(tnbs, 1, max)
    maxenb <- max(colMeans(tnbs))
    evsi <- assurance <- rep(0, length(N))
    evpi <- mean(maxnbs)-maxenb
    assurance0 <- mean(apply(tnbs, 1, which.max)==which.max(colMeans(tnbs)))
    time <- Sys.time()
    for(I in 1:n_sim)
    {
      for(i in 1:length(N))
      {
        n <- N[i]
        nd <- rbinom(rep(n,n_sim), n, sample$prev)
        ntp <- rbinom(rep(n,n_sim), nd, sample$se)
        nfp <- rbinom(rep(n,n_sim), n-nd, 1-sample$sp)
        
        nb1 <- ntp/n - nfp/n*threshold/(1-threshold)
        nb2 <- nd/n - (1-nd/n)*threshold/(1-threshold)
        
        winners <- apply(cbind(0,nb1,nb2), 1, which.max)
        winnertnbs <- tnbs[cbind(1:n_sim,winners)]
        evsi[i] <- evsi[i]+mean(winnertnbs)-maxenb
        assurance[i] <- assurance[i]+mean(winnertnbs==maxnbs)
      }
    }
    message(Sys.time()-time)
    if(b_voi)
    {
      out$nb$evpi <- evpi #res$EVPI
      out$nb$evsi <- evsi/n_sim #res$EVSI
    }
    if(b_assurance)
    {
      out$nb$assurance0 <- assurance0
      out$nb$assurance <- assurance/n_sim #res$EVSIp
    }
  }
  
  out$sample <- sample
  out
}

