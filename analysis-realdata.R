library(nlme)
library(tidyverse)
library(Deriv)
library(stringr)
library(LaplacesDemon)
library(purrr)
library(MASS)
library(xtable)
library(simsurv)
library(survival)
library(dplyr)
library(ggplot2)

## Source all functions  
(file.sources = list.files(path=here::here("src"),pattern="*.R$"))
(file.sources <- paste0(here::here("src"), "/", file.sources))
sapply(file.sources,source)

## Read in the data
data_longit <- read.csv("data1.txt", sep = "") %>%
  dplyr::select(-c(No, `No.1`)) %>%
  mutate(patid = as.integer(factor(patid))) %>%
  mutate(Day = Day/max(Day))

# Nonlinear function used for modeling
nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)

# Fit preliminary NLME model to longitudinal response to obtain starting values
# for fixed effects, random effect SDs, and residual SD
start0 <- c(p1=3,p2=1,p3=0.5) # random starting values
# Nonlinear least squares to obtain starting values for NLME
nls.fit  <- nls(RNA~nf1(p1,p2,p3, Day), data=data_longit, start=start0)
start <- coef(nls.fit)
nlme.fit <- nlme(RNA ~ nf1(p1, p2, p3, Day),
                 data = data_longit,
                 fixed = p1 + p2 + p3 ~ 1,
                 random = p1 + p3 ~ 1 | patid,
                 start = start,
                 control = nlmeControl(maxIter = 200, 
                                       msMaxIter = 200, 
                                       pnlsMaxIter = 100))

# Variance model
sigma2 <- list(model=~1+cd4+(1|patid),
               link='log',
               ran.dist="normal",
               str.fixed=c(2*log(nlme.fit$sigma), 0),
               lower.fixed=NULL,
               upper.fixed=NULL,
               fixName="alpha",
               ranName="v",
               dispName="sigma",
               str.disp=0.5
)

# Mean model (using residual variance modeled above)
nlmeObject_JM <- list(nf = "nf1",
                      model= RNA ~ nf(p1,p2,p3,Day),
                      var=c("Day"),
                      fixed = p1+p2+p3 ~1,
                      random = p1+p3 ~1,
                      family='normal', 
                      ran.dist='normal',
                      fixName="beta",
                      ranName="u",
                      dispName="d",
                      sigma=sigma2,    # residual dispersion model (include residual random eff)
                      ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff)
                      str.fixed=fixef(nlme.fit),  # starting value for fixed effect
                      str.disp=as.numeric(VarCorr(nlme.fit)[,"StdDev"][c("p1", "p3")]),  # starting value for fixed dispersion of random eff
                      lower.fixed=NULL, # lower bounds for fixed eff
                      upper.fixed=rep(100,3), # upper bounds for fixed eff
                      lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
                      upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

nlmeObjects_JM <- list(nlmeObject_JM)

# Jointly fit mean and variance models
JM <- try(Rnlme(nlmeObjects=nlmeObjects_JM, long.data=data_longit,
                idVar="patid", sd.method="HL", dispersion.SD = TRUE,
                independent.raneff="byModel"))

AIC_JM <- JM$AIC

## Survival data
data_dropout <- data_longit %>%
  group_by(patid) %>%
  summarise(event_time = max(Day)) %>%
  mutate(status = ifelse(event_time > 0.75, 0, 1)) %>%
  mutate(patid = as.character(patid))

data_rebound <- data_longit %>%
  group_by(patid) %>%
  mutate(suppressed = (RNA == 1.698970),
         first_suppressed = ifelse(any(suppressed),
                                   min(Day[suppressed]),
                                   1),
         rebound = (Day > first_suppressed) & (RNA > 1.698970)) %>%
  summarise(event_time = ifelse(any(rebound),
                           min(Day[rebound]),
                           max(Day)),
         status = ifelse(any(rebound), 1, 0)
  ) %>%
  ungroup() %>%
  mutate(patid = as.character(patid))

# Extract random effects from mean and variance models
pats <- row.names(ranef(nlme.fit))
est_u3 <- data.frame(pats, scale(ranef(nlme.fit)$p3, center=T,scale=T))
names(est_u3) <- c("patid", "est_u3")

est_v <- data.frame(pats, scale(JM$Bi$v0, center=T,scale=T))
names(est_v) <- c("patid", "est_v")

# Only use dropout for report analysis
data_surv_re <- data_dropout %>%
  merge(est_u3, by='patid', all=TRUE) %>%
  merge(est_v, by='patid', all=TRUE)
# data_surv_re <- data_rebound %>%
#   merge(est_u3, by='patid', all=TRUE) %>%
#   merge(est_v, by='patid', all=TRUE)

obj_surv <- Surv(data_dropout$event_time, data_dropout$status)
#obj_surv <- Surv(data_rebound$event_time, data_rebound$status)
fit_coxph <- coxph(obj_surv~est_u3+est_v , data=data_surv_re)
summary(fit_coxph)
AIC_coxph <- extractAIC(fit_coxph)
AIC_JM <- JM$AIC

## Analysis results
# Longitudinal fixed effect parameters
est_fixed <- JM$fixedest
SD_fixed <- JM$fixedSD
z_fixed <- abs(est_fixed / SD_fixed)
p_fixed <- 2 * pnorm(z_fixed, lower.tail = FALSE)
p_fixed

# Survival model parameters
est_surv <- fit_coxph$coefficients
expest_surv <- exp(est_surv)
SD_surv <- sqrt(diag(fit_coxph$var))
z_surv <- abs(est_surv / SD_surv)
p_surv <- 2 * pnorm(z_surv, lower.tail = FALSE)
results_all <- data.frame(Estimate = c(est_fixed, est_surv),
                          SD = c(SD_fixed, SD_surv),
                          pvalue = c(p_fixed, p_surv),
                          row.names = c("Beta1", "Beta2", "Beta3",
                                       "Alpha0", "Alpha1",
                                       "Gamma1", "Gamma2"))

## Summary statistics
data_longit <- data_longit %>%
  group_by(patid) %>%
  mutate(ni = n())

rna_mean <- mean(data_longit$RNA)
rna_sd <- sd(data_longit$RNA)
cd4_mean <- mean(data_longit$cd4)
cd4_sd <- sd(data_longit$cd4)

n <- max(data_longit$patid)
ni_min <- min(data_longit$ni)
ni_max <- max(data_longit$ni)

data_longit2 <- read.csv("data2.txt", sep = "")
data_longit2 <- na.exclude(data_longit2[,c(2,4,5,10)]) %>%
  mutate(patid = as.integer(factor(patid)))

# Survival event prevalence
1 - mean(data_dropout$status)
1 - mean(data_rebound$status)

## Visualizations
set.seed(1005)
plot.samp <- ggplot(data_longit[data_longit$patid%in%c(4:9),], 
                    aes(x = Day, y = RNA, color = patid, group = patid)) +
  ylim(min_y = 0, max_y = max(data_longit$RNA)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 1.698970, lty = 2, colour = "darkgrey") +
  geom_label(aes(x = 0.1, y = 1.5), label = "Detection limit", colour = "dimgrey") +
  guides(color = "none") +
  labs(x = "Day (standardized)", y = "Viral load (log-transformed)") + theme_pubclean()
plot.samp

plot.all <- ggplot(data_longit, 
                   aes(x = Day, y = RNA, color = patid, group = patid)) +
  ylim(min_y = 0, max_y = max(data_longit$RNA)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 1.698970, lty = 2, colour = "darkgrey") +
  geom_label(aes(x = 0.1, y = 1.5), label = "Detection limit", colour = "dimgrey") +
  guides(color = "none") +
  labs(x = "Day (standardized)", y = "Viral load (log-transformed)") + theme_pubclean()
plot.all

## Dropout plots
data_dropout_ordered <- data_dropout %>%
  mutate(patid = factor(patid, levels = rev(unique(patid))))   # top-to-bottom ordering
plot.dropout <- ggplot(data_dropout_ordered, aes(y = patid)) +
  # Horizontal lines
  geom_segment(aes(x = 0, xend = event_time, yend = patid)) +
  # Event / censor marks
  geom_point(aes(x = event_time,
                 shape = factor(status),
                 fill  = factor(status),
                 colour = factor(status)),
             size = 3) +
  scale_shape_manual(values = c("0" = 21, "1" = 21)) +   # both circles
  scale_fill_manual(name = "Event status",
                    labels = c("0" = "Censored", "1" = "Observed"),
                    values = c("0" = "white", "1" = "navy")) +
  scale_color_manual(values = c("0" = "navy", "1" = "navy")) +
  guides(shape = "none",
         color = "none",
         fill = guide_legend(override.aes = list(
           shape = 21,
           colour = "navy", 
           fill = c("white", "navy")))) +
  xlab("Event time") + ylab("Patient") +
  theme_pubclean()
plot.dropout

plot.dropout.5 <- ggplot(data_dropout_ordered[c(4:9),], aes(y = patid)) +
  # Horizontal lines
  geom_segment(aes(x = 0, xend = event_time, yend = patid)) +
  # Event / censor marks
  geom_point(aes(x = event_time,
                 shape = factor(status),
                 fill  = factor(status),
                 colour = factor(status)),
             size = 3) +
  scale_shape_manual(values = c("0" = 21, "1" = 21)) +   # both circles
  scale_fill_manual(name = "Event status",
                    labels = c("0" = "Censored", "1" = "Observed"),
                    values = c("0" = "white", "1" = "navy")) +
  scale_color_manual(values = c("0" = "navy", "1" = "navy")) +
  guides(shape = "none",
         color = "none",
         fill = guide_legend(override.aes = list(
           shape = 21,
           colour = "navy", 
           fill = c("white", "navy")))) +
  xlab("Event time") + ylab("Patient") +
  theme_pubclean()
plot.dropout.5

## Rebound plots
plot.rebound.longit <- ggplot(data_longit[data_longit$patid%in%c(17,19,24,29),], 
                   aes(x = Day, y = RNA, color = patid, group = patid)) +
  ylim(min_y = 0, max_y = max(data_longit$RNA)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 1.698970, lty = 2, colour = "darkgrey") +
  geom_label(aes(x = 0.1, y = 1.5), label = "Detection limit", colour = "dimgrey") +
  guides(color = "none") +
  labs(x = "Day (standardized)", y = "Viral load (log-transformed)") + theme_pubclean()
plot.rebound.longit

data_rebound_ordered <- data_rebound %>%
  mutate(patid = factor(patid, levels = rev(unique(patid))))

plot.rebound <- ggplot(data_rebound_ordered, aes(y = patid)) +
  # Horizontal lines
  geom_segment(aes(x = 0, xend = event_time, yend = patid), colour = "navy") +
  # Event / censor marks
  geom_point(aes(x = event_time,
                 shape = factor(status),
                 fill  = factor(status),
                 colour = factor(status)),
             size = 3) +
  scale_shape_manual(values = c("0" = 21, "1" = 21)) +   # both circles
  scale_fill_manual(name = "Event status",
                    labels = c("0" = "Censored", "1" = "Observed"),
                    values = c("0" = "white", "1" = "navy")) +
  scale_color_manual(values = c("0" = "navy", "1" = "navy")) +
  guides(shape = "none",
         color = "none",
         fill = guide_legend(override.aes = list(
           shape = 21,
           colour = "navy", 
           fill = c("white", "navy")))) +
  xlab("Event time") + ylab("Patient") +
  theme_pubclean()
plot.rebound

plot.rebound.4 <- ggplot(data_rebound[data_rebound$patid%in%c(17,19,24,29),], aes(y = patid)) +
  # Horizontal lines
  geom_segment(aes(x = 0, xend = event_time, yend = patid), colour = "navy") +
  # Event / censor marks
  geom_point(aes(x = event_time,
                 shape = factor(status),
                 fill  = factor(status),
                 colour = factor(status)),
             size = 3) +
  scale_shape_manual(values = c("0" = 21, "1" = 21)) +   # both circles
  scale_fill_manual(name = "Event status",
                    labels = c("0" = "Censored", "1" = "Observed"),
                    values = c("0" = "white", "1" = "navy")) +
  scale_color_manual(values = c("0" = "navy", "1" = "navy")) +
  guides(shape = "none",
         color = "none",
         fill = guide_legend(override.aes = list(
           shape = 21,
           colour = "navy", 
           fill = c("white", "navy")))) +
  xlab("Event time") + ylab("Patient") +
  theme_pubclean()
plot.rebound
