make.data.generator <- function(
  estimates = cbind("Intercept" = 0, "TS.vs.MT" = 0, "MT.vs.ST" = 0, "ST.vs.CNTL" = 0), 
  n.subj.perCondition = 80,
  n.subj = n.subj.perCondition * ncond, # Total number of subjects (for between subject design)
  n.obs = 16,                        # Number of items per subject
  subject.ranef.covar = matrix(
    rep(0, length(estimates)**2), 
    nrow = length(estimates),
    dimnames = list(colnames(estimates), colnames(estimates))
  ),
  item.ranef.covar = matrix(
    rep(0, length(estimates)**2), 
    nrow = length(estimates),
    dimnames = list(colnames(estimates), colnames(estimates))
  )
) {
  require(mvtnorm)
  require(magrittr)
  require(dplyr)
  
  ncond = length(estimates)
  n = n.obs * n.subj 
  d = expand.grid(Subject = 1:n.subj,
                  Item = 1:n.obs,
                  Intercept = 1)
  d = d[order(d$Subject),]
  d$Cond = factor(
    sort(rep(
      c("Control", "Single talker", "Multi-talker", "Talker-specific"),
      n / ncond)),
    levels = c("Talker-specific", "Multi-talker", "Single talker", "Control"))
  contrasts(d$Cond) = cbind("TS.vs.MT" = c(3/4, -1/4, -1/4, -1/4),
                            "MT.vs.ST" = c(1/2, 1/2, -1/2, -1/2),
                            "ST.vs.CNTL" = c(1/4, 1/4, 1/4, -3/4))
  
  d %<>%
    mutate(
      TS.vs.MT = case_when(
        Cond == "Talker-specific" ~ 3/4,
        Cond == "Multi-talker" ~ -1/4,
        Cond == "Single talker" ~ -1/4,
        T ~ -1/4
      ),
      MT.vs.ST = case_when(
        Cond == "Talker-specific" ~ 1/2,
        Cond == "Multi-talker" ~ 1/2,
        Cond == "Single talker" ~ -1/2,
        T ~ -1/2
      ),
      ST.vs.CNTL = case_when(
        Cond == "Talker-specific" ~ 1/4,
        Cond == "Multi-talker" ~ 1/4,
        Cond == "Single talker" ~ 1/4,
        T ~ -3/4
      )
    )
  
  generate.data <- function() {
    if (all(dim(item.ranef.covar) == 1)) {
      item.adjustment <- data.frame(sapply(1:n.obs, function(x) { rnorm(n = 1, mean = 0, eval(item.ranef.covar**.5)) } ))
    } else {
      item.adjustment <- as.data.frame(t(sapply(1:n.obs, function(x) t(rmvnorm(n = 1, sigma = eval(item.ranef.covar))))))
    }
    names(item.adjustment) <- colnames(eval(item.ranef.covar))
    item.adjustment$Item <- 1:n.obs
    
    if (all(dim(subject.ranef.covar) == 1)) {
      subject.adjustment <- data.frame(sapply(1:n.subj, function(x) { rnorm(n = 1, mean = 0, subject.ranef.covar**.5) } ))
    } else{
      subject.adjustment <- as.data.frame(t(sapply(1:n.subj, function(x) t(rmvnorm(n = 1, sigma = subject.ranef.covar)))))
    }
    names(subject.adjustment) <- colnames(subject.ranef.covar)
    subject.adjustment$Subject <- 1:n.subj

    d$y <- rowSums(d[,colnames(estimates)] * (
      matrix(
        rep(as.numeric(estimates), n),
        byrow = T, 
        nrow = n
      ) +
        item.adjustment[d$Item, colnames(estimates)] + 
        subject.adjustment[d$Subject, colnames(estimates)]
    ))
    d$Item = factor(d$Item)
    d$Subject = factor(d$Subject)
    
    d$Outcome = sapply(d$y, FUN = function(x) { rbinom(1,1,plogis(x))})
    
    return(d)
  }
  
  return(generate.data)
}


fit.models <- function(d.sim) {
  require(lme4)
  
  m.tmp <- glmer(Outcome ~ 1 + Cond + 
                   (1 | Subject) + (1 + Cond | Item), # only by-subject intercepts; no item random effects 
                 d.sim, 
                 family = "binomial")
  m.tmp.coef <- coef(summary(m.tmp))
  
  simulation <- data.frame(
    Int.estimate = m.tmp.coef[1,1],
    Int.z = m.tmp.coef[1,3],
    TS.vs.MT.estimate = m.tmp.coef[2,1],
    TS.vs.MT.z = m.tmp.coef[2,3],
    MT.vs.ST.estimate = m.tmp.coef[3,1],
    MT.vs.ST.z = m.tmp.coef[3,3],
    ST.vs.CNTL.estimate = m.tmp.coef[4,1],
    ST.vs.CNTL.z = m.tmp.coef[4,3],
    Subj.sd = VarCorr(m.tmp)[[1]][1]^.5,
    Item.sd = VarCorr(m.tmp)[[2]][1]^.5,
    ConvergenceFailure = any(grepl("failed to converge", m.tmp@optinfo$conv$lme4$messages)),
    IsSingular = isSingular(m.tmp),
    n.subj = length(unique(d.sim$Subject)),
    n.item = length(unique(d.sim$Item))
  )
  
  return(simulation)
}


data_grid = function (data, ..., .model = NULL) 
{
  require(modelr)
  expanded <- tidyr::expand(data, ...)
  if (is.null(.model)) {
    return(expanded)
  }
  requested <- unlist(strsplit(
    gsub("[ ]+", " ", 
         gsub(" 1 ", " ", 
              gsub("[()|+]", "", .model$formula$formula[3]))), " ", fixed = F))
  
  # print(requested)
  
  needed <- setdiff(requested, names(expanded))
  typical_vals <- lapply(data[needed], typical)
  typical_df <- tidyr::crossing(!!!typical_vals)
  tidyr::crossing(expanded, typical_df)
}

xhypothesis = function(hypothesis, model = NULL, label = NULL, format = "markdown", ...) {
  require(knitr)
  
  if (is.null(hypothesis)) 
    stop("ERROR: Must specify one hypothesis.\n")
  else if (!is.character(hypothesis) & class(hypothesis) != "brmshypothesis")
    stop("ERROR: Argument hypothesis must be a brmshypothesis or character.\n")
  else if (is.character(hypothesis) & is.null(model)) 
    stop("ERROR: Must specify a model, if hypothesis is a character.\n")
  
  if (class(hypothesis) == "brmshypothesis")
    k = hypothesis$hypothesis
  else 
    k = hypothesis(x = model, hypothesis = hypothesis, ...)$hypothesis
    
  if (length(label) != 1)
      stop("ERROR: Number of row labels does not match number of hypotheses (=1).\n")
    else 
      hypotheses[, 1] = labels
  
  return(kable(k,
               digits = c(0,2,3,2,2,1,3,0),
               col.names = c("Hypothesis","Est.","SE","CI~L~","CI~U~","BF","p~postior~",""),
               caption = "Summary of non-linear Bayesian hypothesis testing",
               format = format)
  )
}

xhypotheses = function(hypotheses = list(), model = NULL, labels = NULL, format = "markdown", ...) {
  require(knitr)
  require(dplyr)
  require(brms)
  
  if (is.null(hypotheses)) 
    stop("ERROR: Must specify at least one hypothesis.\n")
  else if (!is.list(hypotheses))
    xhypothesis(hypotheses, model, labels, ...)
  else if (!is.character(hypotheses[[1]]) & class(hypotheses[[1]]) != "brmshypothesis")
    stop("ERROR: Argument hypotheses is a list, but must contain elements that are characters or brmshypothesis objects.\n")
  else if (is.character(hypotheses[[1]]) & is.null(model)) 
    stop("ERROR: Must specify a model, if hypothesis is a list of characters.\n")
  else if (class(hypotheses[[1]]) == "brmshypothesis" & !is.null(model)) 
    cat("WARNING: Ignoring argument model because argument hypotheses contains brmsfithypothesis objects.\n")
  else if (is.character(hypotheses[[1]]))
    hypotheses = lapply(
      hypotheses,
      FUN = function(h) brms::hypothesis(
        x = model,
        hypothesis = h, 
        ...
      )
    )

  # extract the data frame
  hypotheses = 
    do.call(
      rbind.data.frame,
      lapply(
        hypotheses,
        function(h) h$hypothesis
      )
    )
  
  if (!is.null(labels)) {
    if (length(labels) != nrow(hypotheses))
      stop(paste0("ERROR: Number of row labels does not match number of hypotheses (=", nrow(hypotheses), ").\n"))
    else 
      hypotheses[, 1] = labels
  }
  
  return(kable(hypotheses, 
               digits = c(0,2,3,2,2,1,3,0), 
               col.names = c("Hypothesis","Est.","SE","CI~L~","CI~U~","BF","p~post~(h)",""),
               caption = "Summary of non-linear Bayesian hypothesis testing",
               format = format)
  )
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

bracketsGrob <- function(...){
  l <- list(...)
  e <- new.env()
  e$l <- l
  grid:::recordGrob(  {
    do.call(grid.brackets, l)
  }, e)
}

myGplot.defaults = function(
  type = c("paper","poster","slides")[1],
  base_size = if (type == "paper") { 10 } else if (type == "slides") { 32 } else if (type == "poster") { 36 } else { 10 }, 
  margin=c(0.6,0.5,0.5,0.3)
)
{
  require(ggplot2)
  
  theme_set(theme_bw(base_size=base_size))
  theme_update(
    axis.text.x = element_text(size=base_size, vjust=1),
    axis.text.y = element_text(size=base_size, hjust=1, vjust=.5),
    axis.title.x = element_text(size=base_size+1, vjust=0, hjust=0.5, face = "bold"), 
    axis.title.y = element_text(angle=90, size=base_size+1, hjust= 0.4, vjust=0.5, face = "bold"), 
    legend.title = element_text(size=base_size+1, face = "bold", hjust= 0), 
    legend.text = element_text(size=base_size)
    #  1) top, 2) right, 3) bottom, 4) left
    # plot.margin = unit(margin, "lines")
  )
}

whichDFMethod <- function(model) {
  l = length(fitted(model))
  
  # This threshold was set more or less arbitrarily: we noticed that in some cases the smaller [F13] data set(s) 
  # yielded much inflated DF estimates under the Satterthwaite approximation, presumably because of variance 
  # estimates for the random slope associated with the relevant predictor that were indistinguishable form zero.
  # Using Kenward-Roger's approximation instead seems to solve the problem and yield reasonable DF estimates for
  # all tests.
  if (l < 3000) {
    cat(paste0("Using Kenward-Roger's approximation of DFs for t-test over coefficient because the data set is small: ", l, "\n"))
    return("Kenward-Roger") 
  } else { return("Satterthwaite") }
}

get_stat = function(model, coef) {
  s = summary(model)
  
  if (is.brmsfit(model)) {
    return(s$fixed[coef, "Estimate"] / s$fixed[coef, "Est.Error"])
  } else if (is.GLMM(model)) {
    return(coef(s)[coef,"z value"])
  } else {
    return(coef(s)[coef,"t value"])
  }
}

# using Kenward-Rogers approximation since it's more accurate and there was one instance (for 
# the smaller data set [F13] for which Satterthwaite approximation resulted in largely inflated
# DFs.)
get_df = function(model, coef, printModel) {
  if (is.brmsfit(model)) {
    s = summary(model)
    df = s$fixed[coef, "Eff.Sample"]
  } else if (is.GLMM(model)) {
    s = summary(model)
    df = coef(s)[coef,"df"] + 1 # Since the ReplicationBayesfactorNCT function subtracts 1 for the DFs (for sample == 1)
  } else {
    require(lmerTest)
    s = summary(model, ddf = whichDF(model))
    df = coef(s)[coef,"df"] + 1 # Since the ReplicationBayesfactorNCT function subtracts 1 for the DFs (for sample == 1)
  }
  
  if (printModel) {
    cat("\n===================================================================\n")
    print(s)
    cat("\n===================================================================\n")
  } 
  
  return(df)
}

replicationBFs = function(model1, model2, df1 = NULL, df2 = NULL, coefs) {
  p.rep = list()
  
  for (c in coefs) {
    TEMP = replicationBF(model1, model2, df1, df2, coef = c)
    p <- last_plot()
    
    # Store legend in list element '0' (which shouldn't be used for other purposes)
    # if (c == coefs[1]) p.rep[[0]] = get_legend(p)
    p.rep[[c]] = p + theme(legend.position = "none")
  }
  
  # returns an array of plots
  return(p.rep)
}

# wrapper function that calls ReplicationBayesfactorNCT with t-values and n1, n2 observations
# based on the output of mixed model analyses. The test is conducted for the k=th t-test within
# the mixed model, specified by the "coef" argument. 
#
# *** This function has been adapted for the use with glmer logistic models, which use the z-
#     instead of t-statistics. (and therefore also do not require approximation of the DFs for
#     the t-statistics)
replicationBF <- function(model1, model2, df1 = NULL, df2 = NULL, yhigh = 0, limits.x = NULL, coef = 1, M = 500000, printModels = F) {
  require(pbkrtest)
  
  # cat("\n\nCalculting replication Bayes Factor:\n")
  out <- ReplicationBayesfactorNCT(
    get_stat(model1, coef), 
    get_stat(model2, coef),
    if (is.null(df1)) get_df(model1, coef, printModels) else df1, 
    if (is.null(df2)) get_df(model1, coef, printModels) else df2, 
    sample = 1,
    plot = 1, 
    post = 1,
    yhigh = yhigh,
    limits.x = limits.x,
    M = M
  )
  
  return (out)
}

# Replication Bayes Factor taken from Josine Verhagen's webpage
# (http://josineverhagen.com/?page_id=76).
# See Verhagen and Wagenmakers (2014) for details.
ReplicationBayesfactorNCT <- function(
  tobs,                  # t value in first experiment
  trep,                  # t value in replicated experiment
  n1,                    # first experiment: n in group 1 or total n  
  n2,                    # second experiment: n in group 1 or total n 
  m1       = 1,          # first experiment: n in group 2 or total n
  m2       = 1,          # second experiment: n in group 2 or total n
  sample   = 1,        # 1 = one sample t-test (or within), 2 = two sample t-test	
  wod      = dir,	      # working directory 
  plot = 0,  # 0 = no plot 1 = replication   
  post = 0,  # 0 = no posterior, 1 = estimate posterior
  M = 500000, 
  yhigh = 0,
  limits.x = NULL
)	
{
  require(MCMCpack)
  require(ggplot2) # added by TFJ
  
  ##################################################################
  #STEP 1: compute the prior for delta based on the first experiment
  ################################### ###############################
  
  D    <- tobs 
  
  if (sample==1) {
    sqrt.n.orig  <-  sqrt(n1)
    df.orig <- n1 -1  
  }
  
  if (sample==2) {
    sqrt.n.orig  <- sqrt( 1/(1/n1+1/m1)) # two sample alternative for sqrt(n)
    df.orig <- n1 + m1 -2   #degrees of freedom
  }
  
  # To find out quickly the lower and upper bound areas, the .025 and .975 quantiles of tobs at a range of values for D are computed. 
  # To make the algorithm faster, only values in a reasonable area around D are computed.    
  # Determination of area, larger with large D and small N
  range.D <- 4 + abs(D/(2*sqrt.n.orig)) 
  sequence.D <- seq(D,D + range.D,.01) #make sequence with range 
  # determine which D gives a quantile closest to tobs with an accuracy of .01
  options(warn=-1)
  # TFJ:  qt() is the quantile function for the Student t-distribution.
  #       For each value of the sequence of non-centrality parameters (ncp) given in sequence.D, we are looking for the the .025 quantile 
  #       for a t-distribution with the DFs from the original data. We select the value in sequence.D that gives us the lowest .025 quantile.
  #
  # TFJ:  The non-central t-distribution is used when testing hypotheses that t is different from 0 (i.e., the hypothesis that t = ncp).
  #       Here the non-centrality parameter is referred to as D. In the R documentation for the t-distribution, it is referred to as ncp (and as "Del" for Delta).
  #       On wikipedia, it's simply mu.
  approximatelow.D <- sequence.D[which( abs(qt(.025,df.orig,sequence.D)-tobs)==min(abs(qt(.025,df.orig,sequence.D)-tobs)) )]  
  options(warn=0)
  # TFJ:  now we will use increasingly more fine-grained smaller step sizes to generate new sequences of D. Repeating the same step as above.
  #       This gives us low.D and sdlow.D
  #
  # Then a more accurate interval is computed within this area
  # Make sequence within .01 from value found before 
  sequenceappr.D <- seq((approximatelow.D-.01),(approximatelow.D+.01),.00001) 
  # determine which D gives a quantile closest to tobs with an accuracy of .00001
  low.D <- sequenceappr.D[which( abs(qt(.025,df.orig,sequenceappr.D)-tobs)==min(abs(qt(.025,df.orig,sequenceappr.D)-tobs)) )]  
  
  # Compute standard deviation for the corresponding normal distribution.
  sdlow.D <- (D-low.D)/qnorm(.025) 
  
  # compute prior mean and as for delta 
  prior.mudelta <- D/sqrt.n.orig
  prior.sdelta <-  sdlow.D/sqrt.n.orig
  
  ##################################################################
  #STEP 2: Compute Replication Bayes Factor
  ##################################################################
  # For one sample t-test: within (group2==1) or between (group2==vector)
  
  if (sample==1)
  {  
    df.rep <- n2-1 
    sqrt.n.rep <- sqrt(n2)
    # TFJ:  Get the density of t = 0 for a central t-distribution with the t-value from the replication data and the degrees of 
    #       freedom from the replication data. If the same density under the original data is divided by this density (Likelihood.Y.H0) 
    #       under the replication data, this gives the replication Bayes factor (BF_r0).
    Likelihood.Y.H0 <- dt(trep,df.rep)
    
    # TFJ:  Sample M samples from the normal prior with mean prior.mudelta and standard deviation prior.sdelta. This mean is the 
    #       t-value observed in the original data (D or tobs) divided by the square root of the original data points (for sample == 1).
    #       Below that mean that is multiplied by the square root of the replication data points.
    # TFJ:  Each of these samples corresponds to a possible effect size (delta(i) in equation 6, drawn from the posterior distribution
    #       of effect sizes delta, given the original data [given the original t-value and DFs]). THIS IS THE ONLY NON-DETERMINISTIC
    #       (SAMPLING) STEP UP TO THIS POINT.
    sample.prior <- rnorm(M,prior.mudelta,prior.sdelta) 
    options(warn=-1)
    # TFJ:  What is the mean density of the t-value observed in the replication, given the replication's degrees of freedom and
    #       a non-centrality parameters (ncp) that corresponds to each of the sampled effect sizes (corrected for the square root
    #       of the number of data points in the replication data). This gives us the numerator for the last line of equation 6 on 
    #       p. 1461 (already averaging across all samples, i.e., aready also doing 1/M Sum_1toM)
    average.Likelihood.H1.delta <- mean(dt(trep,df.rep,sample.prior*sqrt.n.rep))
    options(warn=0)
    
    # TFJ:  Calculate the replication BF as 1) the average likelihood of the t-value and DFs observed in the replication
    #       (t_rep, df.rep) given the hypothesis of the original effect sizes divided by 2) the likelihood of a null 
    #       effect given t_rep.
    BF <- average.Likelihood.H1.delta/Likelihood.Y.H0
  }
  
  # For two sample t-test: 
  
  if (sample==2)  {
    df.rep <- n2 + m2 -2 
    sqrt.n.rep <- sqrt(1/(1/n2+1/m2))
    Likelihood.Y.H0 <- dt(trep,df.rep)
    
    sample.prior <- rnorm(M,prior.mudelta,prior.sdelta)
    options(warn=-1)
    # TFJ:  Before we asses the density of trep for each of the different samples drawn from the prior, we multiply the samples
    #       by the square root of the number of data points in the replication.
    average.Likelihood.H1.delta <- mean(dt(trep,df.rep,sample.prior*sqrt.n.rep))
    options(warn=0)
    BF <- average.Likelihood.H1.delta/Likelihood.Y.H0
  }
  
  #################################################################
  #STEP 3: Posterior distribution 
  #################################################################
  
  
  if( post == 1) { 
    
    options(warn=-1)
    likelihood <- dt(trep,df.rep,sample.prior*sqrt.n.rep)
    prior.density <- dnorm(sample.prior,prior.mudelta,prior.sdelta)
    likelihood.x.prior <- likelihood * prior.density
    
    LikelihoodXPrior <- function(x) {dnorm(x,prior.mudelta,prior.sdelta) * dt(trep,df.rep,x*sqrt.n.rep) }
    fact  <- integrate(LikelihoodXPrior,-Inf,Inf) 
    posterior.density <- likelihood.x.prior/fact$value
    PosteriorDensityFunction <- function(x) {(dnorm(x,prior.mudelta,prior.sdelta) * dt(trep,df.rep,x*sqrt.n.rep))/fact$value }
    options(warn=0)
    
    mean <- prior.mudelta
    sdh  <- prior.mudelta + .5*(max(sample.prior)[1] - prior.mudelta)
    sdl  <- prior.mudelta - .5*(prior.mudelta - min(sample.prior)[1])
    dev  <- 2
    
    # TFJ:  iteratively refined integration (the outer loops sets the mean, sdh, sdl and dev) each time
    #       it is run. In particular dev is decreased by an order of magnitude at each iteration of the
    #       outer loop.
    for ( j in 1:10) {
      rangem <- seq((mean-dev)[1],(mean+dev)[1],dev/10)
      rangesdh <- seq((sdh-dev)[1],(sdh+dev)[1],dev/10)
      rangesdl <- seq((sdl-dev)[1],(sdl+dev)[1],dev/10)
      perc <- matrix(0,length(rangem),3)
      
      
      I<-min(length(rangem), length(rangesdh),length(rangesdl) )
      for ( i in 1:I) { 
        options(warn=-1)
        vpercm <-  integrate(PosteriorDensityFunction, -Inf,  rangem[i])
        perc[i,1]<- vpercm$value
        vpercsh <-  integrate(PosteriorDensityFunction, -Inf,  rangesdh[i])
        perc[i,2]<- vpercsh$value
        vpercsl <-  integrate(PosteriorDensityFunction, -Inf,  rangesdl[i])
        perc[i,3]<- vpercsl$value
        options(warn=0)    
      }
      mean <- rangem[which(abs(perc[,1]-.5)== min(abs(perc[,1]-.5)))]
      sdh <-  rangesdh[which(abs(perc[,2]-pnorm(1))== min(abs(perc[,2]-pnorm(1))))]
      sdl <-  rangesdl[which(abs(perc[,3]-pnorm(-1))== min(abs(perc[,3]-pnorm(-1))))]
      dev <- dev/10
    }
    
    posterior.mean <- mean
    posterior.sd <- mean(c(abs(sdh- mean),abs(sdl - mean)))
  }
  
  if (post != 1) {
    posterior.mean <- 0
    posterior.sd <- 0
  } 
  
  ###########OUT 
  
  dat.SD=new.env()
  dat.SD$BF      <- BF
  dat.SD$prior.mean= round(prior.mudelta,2)
  dat.SD$prior.sd= round(prior.sdelta,2)
  dat.SD$post.mean= round(posterior.mean,2)
  dat.SD$post.sd= round(posterior.sd,2)
  dat.SD=as.list(dat.SD)
  
  
  ###########################################  
  #PLOT
  ########################################### 
  if (plot==1)
  {
    
    if (post != 1) {
      options(warn =-1)
      rp <- ReplicationPosterior(trep,prior.mudelta,prior.sdelta,n2,m2=1,sample=sample)
      options(warn =0)
      
      posterior.mean <- rp[[1]] 
      posterior.sd  <- rp[[2]]
    }
    
    # Set x limites to max of density over mean of prior or posterior (whichever is larger)
    max.prior.density = dnorm(prior.mudelta,prior.mudelta,prior.sdelta)
    max.post.density = dnorm(posterior.mean,posterior.mean,posterior.sd)
    
    if (max.prior.density > max.post.density) {
      high = max.prior.density
      x.text = prior.mudelta
    } else {
      high = max.post.density
      x.text = posterior.mean
    }
    
    if(yhigh==0) { 
      yhigh <- high + high/5 
    }
    
    scale <- 3
    
    # Set x limites to max of 3 STDEVs of prior or posterior (whichever is more extreme)
    if (is.null(limits.x)) {
      min.x <-  min((posterior.mean - scale*posterior.sd),(prior.mudelta - scale*prior.sdelta))
      max.x <-  max((posterior.mean + scale*posterior.sd),(prior.mudelta + scale*prior.sdelta))
    } else {
      min.x <-  limits.x[1]
      max.x <-  limits.x[2]
    }
    
    ggplot(
      data = data.frame(x = 0),
      mapping = aes(
        x = x
      )
    ) + 
      scale_x_continuous(expression(Effect ~ size ~ delta),
                         limits = c(min.x, max.x)
      ) +
      scale_y_continuous("Density",
                         limits = c(0, yhigh)
      ) +
      stat_function(
        fun = function(x) dnorm(x, prior.mudelta, prior.sdelta),
        aes(linetype = "Before replication"),
        size = .75
      ) +
      stat_function(
        fun = function(x) dnorm(x, posterior.mean, posterior.sd),
        aes(linetype = "After replication"),
        size = 1
      ) +
      scale_linetype_manual("",
                            values = c("Before replication" = 2,
                                       "After replication" = 1)
      ) +
      geom_point(
        x = 0, y = dnorm(0,prior.mudelta,prior.sdelta),
        color = "gray", size = 3
      ) +
      geom_point(
        x = 0, y = dnorm(0,posterior.mean,posterior.sd),
        color = "gray", size = 3
      ) + 
      geom_segment(
        x = 0, xend = 0,
        y = dnorm(0,prior.mudelta,prior.sdelta), yend= dnorm(0,posterior.mean,posterior.sd),
        color = "gray"
      ) +
      annotate(
        geom = "text",
        x = x.text, y = yhigh - (1/10) * yhigh,
        label = as.expression(bquote(BF[r0] == .(round(BF, digits=2))))
      ) +
      theme(
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm")
      )
    
  }
  
  return(dat.SD)
}

emplog = function(p, n) {
  log( (p * n + .5) / (n - p * n + .5) )
}

se = function(x) return(sd(x)/sqrt(sum(!is.na(x))))

inv.gamma <- function(x, a, b){
  ((b^a)/gamma(a))*((1/x)^(a-1))*exp(-b/x) #PDF Inv. Gamma
}
