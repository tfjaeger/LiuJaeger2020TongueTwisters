library(ggplot2)

## ----------------------------------------------------------------------------
#
# FUNCTIONS
#
## ----------------------------------------------------------------------------
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
    axis.title.y = element_text(angle=90, size=base_size+1, hjust= 0.4, vjust=0, face = "bold"), 
    legend.title = element_text(size=base_size+1, face = "bold", hjust= 0), 
    legend.text = element_text(size=base_size)
    #  1) top, 2) right, 3) bottom, 4) left
    # plot.margin = unit(margin, "lines")
  )
}


make.data.generator <- function(
                                estimates = cbind("Intercept" = 0, "Label" = 0, "Context" = 0, "LabelContext" = 0), 
                                n.subj.perCondition = 40,
                                n.subj = n.subj.perCondition * 4, # Total number of subjects (for between subject design)
                                n.obs = 7,                        # Number of items per subject
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

  n = n.obs * n.subj 
  d = expand.grid(Subject = 1:n.subj,
              Item = 1:n.obs,
              Intercept = 1)
  d = d[order(d$Subject),]
  d$Label = c(rep(-1, n / 2), # between subject, exactly balanced.
                  rep(1, n / 2)
  )
  d$Context = rep(c(rep(-1, n / 4), # between subject, exactly balanced.
              rep(1, n / 4)), 2
  )
  d$LabelContext = d$Label * d$Context
    
  generate.data <- function() {
    if (all(dim(item.ranef.covar) == 1)) {
      item.adjustment <- data.frame(sapply(1:n.obs, function(x) { rnorm(n = 1, mean = 0, item.ranef.covar**.5) } ))
    } else {
      item.adjustment <- as.data.frame(t(sapply(1:n.obs, function(x) t(rmvnorm(n = 1, sigma = item.ranef.covar)))))
    }
    names(item.adjustment) <- colnames(item.ranef.covar)
    item.adjustment$Item <- 1:n.obs
    
    if (all(dim(subject.ranef.covar) == 1)) {
      subject.adjustment <- data.frame(sapply(1:n.subj, function(x) { rnorm(n = 1, mean = 0, subject.ranef.covar**.5) } ))
    } else{
      subject.adjustment <- as.data.frame(t(sapply(1:n.subj, function(x) t(rmvnorm(n = 1, sigma = subject.ranef.covar)))))
    }
    names(subject.adjustment) <- colnames(subject.ranef.covar)
    subject.adjustment$Subject <- 1:n.subj
    
    # ********************* TO DO: actually, I do not understand what you are doing here. Why is the Std. Error added 
    # *********************        Why is the Std. Error added to beta, and why the multiplication?
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
}

fit.models <- function(d.sim) {
  require(lme4)
  
  m.tmp <- glmer(Outcome ~ 1 + Label*Context + 
                     (1 | Subject), # only by-subject intercepts; no item random effects 
                   d.sim, 
                   family = "binomial")
  m.tmp.coef <- coef(summary(m.tmp))
  
  simulation <- data.frame(
    Int.estimate = m.tmp.coef[1,1],
    Int.z = m.tmp.coef[1,3],
    Label.estimate = m.tmp.coef[2,1],
    Label.z = m.tmp.coef[2,3],
    Context.estimate = m.tmp.coef[3,1],
    Context.z = m.tmp.coef[3,3],
    LabelContext.estimate = m.tmp.coef[4,1],
    LabelContext.z = m.tmp.coef[4,3],
    Subj.var = VarCorr(m.tmp)[[1]][1],
    ConvergenceFailure = any(grepl("failed to converge", m.tmp@optinfo$conv$lme4$messages)),
    n.subj = length(unique(d.sim$Subject)),
    n.item = length(unique(d.sim$Item))
  )
  
  return(simulation)
}

## ----------------------------------------------------------------------------
#
# SIMULATION
#
## ----------------------------------------------------------------------------

# From Exp 1 in TT paper:
#
# Random effects:
#   Groups   Name        Variance Std.Dev.
# RandomID (Intercept) 0.4483   0.6696  
# Number of obs: 5600, groups:  RandomID, 160
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                               -0.09767    0.06016  -1.623    0.104    
# Label.sumSH vs S                           0.30070    0.06020   4.995 5.88e-07 ***
# Condition.sumnotTT vs TT                  -0.04943    0.06015  -0.822    0.411    
# Label.sumSH vs S:Condition.sumnotTT vs TT  0.02375    0.06014   0.395    0.693    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) L.SHvS C.TTvT
# Lbl.smSHvsS -0.007              
# Cndtn.TTvTT  0.006 -0.003       
# L.SHvS:C.vT -0.002  0.006 -0.006
#
#
# The effect size estimate for Label in Liu and Jaeger was 1.14 log-odds with an SE of .26, 
# under a .5 vs. -.5 sum-coding scheme.

# Using intercept and context effect from Experiment 1 of TT paper. Same for by-subject variance
# in intercept. Label effect is estimated as 50% of Liu and Jaeger (2018, Experiment 1). Label:Context
# interaction is estimated as 50% of the Label effect.
# (Keep in mind that we coded effects -1 and 1 in the TT paper but as -.5 and .5 in Liu and 
# Jaeger. The simulation assumes a coding of -1 and 1)
e = cbind("Intercept" = -.25, "Label" = .56, "Context" = -.07, "LabelContext" = .28)
my.simulator <- make.data.generator(
  estimates = e, 
  n.subj.perCondition = 40, 
  n.obs = 7, # Include all five blocks of 7 stims during Test block
  subject.ranef.covar = matrix(
    c(.9, rep(0, length(e)**2-1)), 
    nrow = length(e),
    dimnames = list(colnames(e), colnames(e))
  )
)


num.sims = 10000
d.sim = plyr::rdply(.n = num.sims, 
            fit.models(my.simulator()),
            .progress = "text")
summary(d.sim)
saveRDS(d.sim, "data/powersims.RDS", compress = T)

# Power = proportion of "true" (significant effects)
round(prop.table(table(d.sim$Label.z > 1.96)), 2)
round(prop.table(table(d.sim$LabelContext.z > 1.96)), 2)

## ----------------------------------------------------------------------------
#
# PLOTS
#
## ----------------------------------------------------------------------------
myGplot.defaults(type = "paper")
ggplot(d.sim, 
       aes(x = Label.z, y = LabelContext.z)) +
  geom_vline(xintercept = c(-1.96, 1.96), color = "red", linetype = 2) +
  geom_hline(yintercept = c(-1.96, 1.96), color = "red", linetype = 2) +
  geom_rect(xmin = -1.96,
            xmax = 1.96,
            ymin = -1.96,
            ymax = 1.96, fill = "red", alpha = .01) +
  geom_point(alpha = .5) +
  geom_density_2d() +
  theme_bw() +
  scale_x_continuous("Perceptual recalibration\n(z-value of Label effect)",
                     limits = c(-3,15)) +
  scale_y_continuous("Blocking of perceptual recalibration\n(z-value of interaction with Context)",
                     limits = c(-3,10))
ggsave("figures/power.pdf", height = 4, width = 4)


ggplot(d.sim, 
       aes(x = Label.estimate, y = LabelContext.estimate)) +
  geom_point(aes(shape = ifelse(abs(LabelContext.z) > 1.96, "yes", "no"),
                 alpha = ifelse(abs(LabelContext.z) > 1.96, 1, .7))) +
  geom_point(x = mean(d.sim$Label.estimate), 
             y = mean(d.sim$LabelContext.estimate),
             color = "blue",
             size = 2
  ) +
  geom_density_2d() +
  geom_vline(xintercept = 0, color = "black", linetype = 2) +
  geom_hline(yintercept = 0, color = "black", linetype = 2) + 
  theme_bw() +
  scale_shape_discrete("Condition significant?") +
  guides(alpha = "none")


d.sim.old = d.sim

