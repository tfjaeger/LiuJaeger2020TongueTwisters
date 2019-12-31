# Different intercepts
set.seed(1234)
num.sims = 10000
d <- NULL

# # Use multiple cores
# library(doParallel}
# registerDoParallel(cores=3)
# foreach(i=1:3) %dopar% sqrt(i)
# # to end parallel computing
# registerDoSEQ()

for (i in c(qlogis(1/2), qlogis(3/4), qlogis(7/8), qlogis(15/16), qlogis(31/32))) {
  for (j in c(0, .1, .2, .4)) {
    my.simulator <- make.data.generator( 
      estimates = cbind("Intercept" = i, "Condition" = j), 
      n.subj = 64, 
      n.obs = 64, # for now we're ignoring that items were crossed across participants (not all participants saw all items)
      subject.ranef.covar = matrix(
        c(.8, 0, 0, 0), 
        nrow = 2,
        dimnames = list(c("Intercept", "Condition"),
                        c("Intercept", "Condition"))
      ),
      item.ranef.covar = matrix(
        c(1.7, 0, 0, 0), 
        nrow = 2,
        dimnames = list(c("Intercept", "Condition"),
                        c("Intercept", "Condition"))
      )
    )
    
    d.sim = plyr::rdply(
      .n = num.sims, 
      fit.models(my.simulator()),
      .progress = "text")
    
    d.sim$Int.true = i
    d.sim$Cond.true = j
    
    if (is.null(d)) {
      d = d.sim
    } else {
      d = rbind(d, d.sim)
    }
  }
}

saveRDS(d, 
        file = paste("simulations", 
                     Sys.time(), 
                     sep = "-"),
        compress = T
)



# plot the distribution of z-value
# This creates an "NAs introduced through coersion warning" during plotting. That's ok, as far as I can tell
Int_labeller <- function(string) {
  ifelse(string == "(all)", 
         "(all)",
         paste("True mean\n", 
               round(plogis(as.numeric(as.character(string)))*100, 1), 
               "%",
               sep = "") 
  )
}

Cond_labeller <- function(string) {
  ifelse(string == "(all)", 
         "(all)",
         paste("True effect\n", 
               string, 
               " logits", 
               sep = "") 
  )
}

signed_log_trans <- function() {
  require(scales)
  
  trans_new("signed_log", 
            function(x) sign(x) * log(abs(x), base = 10), # transformation
            function(x) sign(x) * log(abs(x), base = 10), # inverse
            breaks = function(x) 10^x
  )
} 

breaks = trans_breaks("signed_log", function(x) 10^x),
labels = trans_format("signed_log", math_format(10^.x))

ggplot(d %>% filter(ConvergenceFailure == FALSE), 
       aes(x = Int.z, y = Cond.z)) +
  annotate("rect", xmin = -1.96, xmax = 1.96, ymin = -1.96, ymax = 1.96, fill = "red", alpha = .1) +
  geom_vline(xintercept = c(-1.96, 1.96), color = "red", linetype = 2) +
  geom_hline(yintercept = c(-1.96, 1.96), color = "red", linetype = 2) +
  geom_point(alpha = .1) +
  geom_density_2d() +
  scale_x_continuous("z-value of intercept estimate", trans = "signed_log") +
  scale_y_continuous("z-value of condition estimate", trans = "signed_log") +
  facet_grid(Cond.true ~ Int.true,
             margins = T,
             scales = "fixed", # change to "fixed" to have same scale throughout
             shrink = F,      # change if plotting determined only by summary statistics
             as.table = F,
             switch = "y",
             labeller = labeller(
               Int.true = Int_labeller,
               Cond.true = Cond_labeller
             )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(xlim = c(-50, 50), ylim = c(-50,50))

ggplot(d %>% filter(ConvergenceFailure == FALSE), 
       aes(x = Int.z, y = Cond.z)) +
  annotate("rect", xmin = -1.96, xmax = 1.96, ymin = -1.96, ymax = 1.96, fill = "red", alpha = .1) +
  geom_vline(xintercept = c(-1.96, 1.96), color = "red", linetype = 2) +
  geom_hline(yintercept = c(-1.96, 1.96), color = "red", linetype = 2) +
  geom_point(alpha = .1) +
  geom_density_2d() +
  scale_x_continuous("z-value of intercept estimate") +
  scale_y_continuous("z-value of condition estimate") +
  facet_grid(Cond.true ~ Int.true,
             margins = T,
             scales = "fixed", # change to "fixed" to have same scale throughout
             shrink = F,      # change if plotting determined only by summary statistics
             as.table = F,
             switch = "y",
             labeller = labeller(
               Int.true = Int_labeller,
               Cond.true = Cond_labeller
             )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# The above plot(s) should really be done by coordinate transform, but something goes wrong when I try that:
# coord_trans(x = "signed_log", y = "signed_log")
ggplot(d, 
       aes(x = Int.z, y = Cond.z)) +
  annotate("rect", xmin = -1.96, xmax = 1.96, ymin = -1.96, ymax = 1.96, fill = "red", alpha = .1) +
  geom_vline(xintercept = c(-1.96, 1.96), color = "red", linetype = 2) +
  geom_hline(yintercept = c(-1.96, 1.96), color = "red", linetype = 2) +
  geom_point(alpha = .5) +
  geom_density_2d() +
  facet_grid(Cond.true ~ Int.true,
             margins = T,
             scales = "free", # change to "fixed" to have same scale throughout
             shrink = F,      # change if plotting determined only by summary statistics
             as.table = F,
             switch = "y",
             labeller = labeller(
               Int.true = Int_labeller,
               Cond.true = Cond_labeller
             )) +
  scale_x_continuous("z-value of intercept estimate", trans = "signed_log") +
  scale_y_continuous("z-value of condition estimate", trans = "signed_log") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Get z distribution for true null effects         
d %>% 
  filter(ConvergenceFailure == FALSE & Cond.true == 0) %>%
  group_by(Int.true) %>%
  dplyr::summarise(
    zero.n = length(Cond.z),
    zero.z.mean = mean(Cond.z),
    zero.z.sd = sd(Cond.z)
  ) -> d.zero

# Now let's look at only those simulations that contained true non-zero effect.
# Let's join mean and sd of z-distribution for true null effect to rest of data frame.
# And then let's get rate of signififance with and without correction for Type I error rate
d %>%
  filter(ConvergenceFailure == FALSE & Cond.true != 0) %>%
  left_join(d.zero) %>%
  group_by(Int.true, Cond.true) %>%
  dplyr::summarise(
    n = length(Cond.z),
    corrected.CI.lower = mean(zero.z.mean - 1.96 * zero.z.sd), # CHECK (intent: get z distribution for null effect)
    corrected.CI.upper = mean(zero.z.mean + 1.96 * zero.z.sd), # CHECK (intent: get z distribution for null effect)
    significant.05 = sum(ifelse(abs(Cond.z) > 1.96, 1, 0)),
    prop_signif.05 = significant.05 / n,
    significant.corrected = sum(ifelse(Cond.z < corrected.CI.lower | Cond.z > corrected.CI.upper, 1, 0)),
    prop_signif.corrected = significant.corrected / n,
    significant.corrected.positive = sum(ifelse(Cond.z > corrected.CI.upper, 1, 0)),
    prop_signif.corrected.positive = significant.corrected.positive / n
  ) -> d.summary

d.summary %>%
  select(Int.true, 
         Cond.true, 
         prop_signif.05, 
         prop_signif.corrected, 
         prop_signif.corrected.positive)

# Plotting power as function of true intercept and effect size
# Uncorrected for Type I error
ggplot(d.summary,
       aes(x = round(plogis(Int.true)*100,1), y = prop_signif.05, color = as.factor(Cond.true))) +
  geom_point() +
  geom_line() +
  scale_x_continuous("True mean (percent)") +
  scale_y_continuous("Power") +
  scale_color_discrete("True effect (logits)")

# Corrected for Type I error
ggplot(d.summary,
       aes(x = round(plogis(Int.true)*100,1), y = prop_signif.corrected, color = as.factor(Cond.true))) +
  geom_point() +
  geom_line() +
  scale_x_continuous("True mean (percent)") +
  scale_y_continuous("Power") +
  scale_color_discrete("True effect (logits)")

# Corrected for Type I error and counting only correct direction of effect
ggplot(d.summary,
       aes(x = round(plogis(Int.true)*100,1), y = prop_signif.corrected.positive, color = as.factor(Cond.true))) +
  geom_point() +
  geom_line() +
  scale_x_continuous("True mean (percent)") +
  scale_y_continuous("Power") +
  scale_color_discrete("True effect (logits)")

# I wonder how to interpret the little upticks. could be due to differences in failure to converge?