library(survival)
library(survminer)
library(ggplot2)
library(caret)
library(rpart)
library(rpart.plot)
library(dplyr)
library(tidyr)
library(DataExplorer)
library(corrplot)


cardio <- read.csv("/Users/yuedu/STA 592 Survival/Final project/echocardiogram.csv")
sum(is.na(cardio))
plot_missing(cardio)
cardio <- cardio[-which(cardio$pericardialeffusion==77),]
# attach(cardio)

## missing data imputation age
age.filled <- ifelse(is.na(cardio$age), mean(na.omit(cardio$age)), cardio$age)
cardio$age <- age.filled
sum(is.na(cardio$age))

fractionalshortening.filled <- ifelse(is.na(cardio$fractionalshortening), 
                                      mean(na.omit(cardio$fractionalshortening)), 
                               cardio$fractionalshortening)
cardio$fractionalshortening <- fractionalshortening.filled

epss.filled <- ifelse(is.na(cardio$epss), mean(na.omit(cardio$epss)), cardio$epss)
cardio$epss <- epss.filled

lvdd.filled <- ifelse(is.na(cardio$lvdd), mean(na.omit(cardio$lvdd)), cardio$lvdd)
cardio$lvdd <- lvdd.filled

wallmotion.score.filled <- ifelse(is.na(cardio$wallmotion.score), mean(na.omit(cardio$wallmotion.score)), 
                                  cardio$wallmotion.score)
cardio$wallmotion.score <- wallmotion.score.filled

mult.filled <- ifelse(is.na(cardio$mult), mean(na.omit(cardio$mult)), cardio$mult)
cardio$mult <- mult.filled

wallmotion.index <- ifelse(is.na(cardio$wallmotion.index), (cardio$wallmotion.score)/12, cardio$wallmotion.index)
cardio$mult <- mult.filled

attach(cardio)
pericardialeffusion <- as.factor(pericardialeffusion)

## Maximum likelihood estimation; weibull is good

fit_wei <- survreg(Surv(survival, 1-alive) ~ 1,dist = "weibull", data=cardio)
summary(fit_wei)

#fit_exp <- survreg(Surv(survival, 1-alive) ~ 1,dist = "exponential", data=cardio)
#summary(fit_exp)

anova(fit_exp, fit_wei)

# MLE for beta
coef(fit_wei)

## Fit Kaplan Miere
fit_KM <- survfit(Surv(survival, 1-alive) ~ 1,
               data=cardio, conf.type = "log-log")

summary(fit_KM)
plot(fit_KM)

## check goodness of fit
plot(fit_KM$surv,
     pweibull(fit_KM$time, 1 / fit_wei$scale, exp(coef(fit_wei)), lower.tail = FALSE),
     xlab = "KM", ylab = "Weibull")
title(main="Plot of Weibull fit againts KM")
abline(a = 0, b = 1, col = "red", lty = "dashed")

## test group difference
fit <- survfit(Surv(survival, 1-alive) ~ pericardialeffusion, data = cardio)
plot(fit, col = c("red", "blue"))
legend(32,1.1, legend = names(fit$strata),
       col = c("red", "blue"), lty = "solid")

## two group test whether lamda is equal
fit3 <- survreg(Surv(survival, 1-alive) ~ pericardialeffusion, data = cardio, dist = "weibull")
summary(fit3)

anova(fit_wei,fit3)
confint(fit3)

## assume that the shape parameters could be different
fit4 <- survreg(
  Surv(survival, 1-alive) ~ pericardialeffusion + strata(pericardialeffusion),
  data = cardio, dist = "weibull")

fit_null <- survreg(
  Surv(survival, 1-alive) ~ strata(pericardialeffusion),
  data = cardio, dist = "weibull")

anova(fit_null, fit4)

## test whether the shape parameters are equal or not
fit5 <- survreg(
  Surv(survival, 1-alive) ~ pericardialeffusion + strata(pericardialeffusion),
  data = cardio, dist = "weibull")
fit_null1 <- survreg(
  Surv(survival, 1-alive) ~ pericardialeffusion,
  data = cardio, dist = "weibull")
anova(fit_null1, fit5)

fit6 <- survfit(Surv(survival, 1-alive) ~pericardialeffusion, data = cardio)
summary(fit6)
plot(fit6, col = c("red", "blue"))
legend(35,1, legend = names(fit6$strata),
       col = c("red", "blue"), lty = "solid")
title(main="Check Goodness of fit")

## log-rank test
fit7 <- survdiff(Surv(survival, 1-alive)~pericardialeffusion, data=cardio)
fit7


## Log-location-scale regression for two groups

## Log-location-scale regression
# Weibull; model 1 is better
fit_1 <- survreg(Surv(survival, 1 - alive) ~age+pericardialeffusion+fractionalshortening+epss+
                 wallmotion.score+lvdd,data = cardio, dist = "weibull")
step(fit_1)

fit_3 <- survreg(Surv(survival, 1 - alive) ~1,data = cardio, dist = "weibull")
summary(fit_3)
AIC(fit_3)

fit_2 <- survreg(Surv(survival, 1 - alive) ~age+fractionalshortening+epss+
                   wallmotion.score,data = cardio, dist = "weibull")
AIC(fit_1)
step(fit_1)
summary(fit_1)
fit_8 <- survreg(Surv(survival, 1 - alive) ~age+I(pericardialeffusion==0)+fractionalshortening+epss+
                   wallmotion.score,data = cardio, dist = "weibull")
AIC(fit_8)

summary(fit_8)
step(fit_1)

anova(fit_wei, fit_1)
anova(fit_wei, fit_8)

AIC(fit_1)
BIC(fit_1)

fit_2 <- survreg(Surv(survival, 1 - alive) ~age+pericardialeffusion+fractionalshortening+epss+
                   wallmotion.score+lvdd,
                 data = cardio, dist = "weibull")
summary(fit_2)
AIC(fit_2)
BIC(fit_2)
anova(fit_wei,fit_2) #p-value is almost 0
anova(fit_1, fit_2) 

summary(fit_2)

## model selection
## lognormal
fit_ln <- survreg(Surv(survival, 1 - alive) ~ age+pericardialeffusion+fractionalshortening+epss+
                    wallmotion.score+lvdd, data = cardio, dist = "lognormal")
summary(fit_ln)
step(fit_ln)


## fit1 has the smallest AIC


fit_ln1 <- survreg(Surv(survival, 1 - alive) ~ wallmotion.score, data = cardio, dist = "lognormal")
summary(fit_ln1) 
AIC(fit_ln1)

fit_ln2 <- survreg(Surv(survival, 1 - alive) ~ 1, data = cardio, dist = "lognormal")
summary(fit_ln2)
AIC(fit_ln2)

AIC(fit_ln)

## loglogistic
fit_ll1 <- survreg(Surv(survival, 1 - alive) ~ age+pericardialeffusion+fractionalshortening+epss+
                   wallmotion.score+lvdd, data = cardio, dist = "loglogistic")
summary(fit_ll1)

step(fit_ll1)
fit_ll2 <- survreg(Surv(survival, 1 - alive) ~ wallmotion.score, data = cardio, dist = "loglogistic")
summary(fit_ll2)
AIC(fit_ll2)

fit_ll <- survreg(Surv(survival, 1 - alive) ~ age+pericardialeffusion+fractionalshortening+epss+
                    wallmotion.score, data = cardio, dist = "loglogistic")

anova(fit_ll, fit_ll2)
AIC(fit_ll)

summary(fit_ll)
fit_ll1 <- survreg(Surv(survival, 1 - alive) ~ 1, data = cardio, dist = "loglogistic")

anova(fit_ll1, fit_ll) #suggest the simple model is enough
AIC(fit_ll1)
summary(fit_ll)
AIC(fit_ln)


##compare AIC 
c(wei=AIC(fit_wei), ln=AIC(fit_ln1), ll=AIC(fit_ll2))
c(wei=BIC(fit_wei), ln=BIC(fit_ln1), ll=BIC(fit_ll2))

## the simplest weibull fit is the model with the smallest AIC and BIC


## Distribution of the error term
## weibull
epsilonhat <- (
  log(cardio$survival) - predict(fit_3, type = "linear")
) / fit_3$scale
km <- survfit(Surv(exp(epsilonhat), 1 - alive) ~ 1)
plot(km)
curve(exp(-x), col = "red", add = TRUE)
title(main="Plot of Error term (only with intercept)")

epsilonhat <- (
  log(cardio$survival) - predict(fit_ln1, type = "linear")
) / fit_1$scale
km <- survfit(Surv(exp(epsilonhat), 1 - alive) ~ 1)
plot(km)
curve(exp(-x), col = "red", add = TRUE)
title(main="Plot of Error term (all variable)")

epsilonhat <- (
  log(cardio$survival) - predict(fit_ln1, type = "linear")
) / fit_ln1$scale
km <- survfit(Surv(exp(epsilonhat), 1 - alive) ~ 1)
plot(km)
curve(exp(-x), col = "red", add = TRUE)

title(main="Plot of Error term (drop pericardial effusion)")

## lognormal
epsilonhat <- (
  log(cardio$survival) - predict(fit_ln1, type = "linear")
) / fit_ln1$scale
km <- survfit(Surv(exp(epsilonhat), 1 - alive) ~ 1)
plot(km)
curve(pnorm(log(x), lower.tail = FALSE), col = "red", add = TRUE)
title(main="Plot of error term (with only wall motion score)")


## loglogistic
epsilonhat <- (
  log(cardio$survival) - predict(fit_ll2, type = "linear")
) / fit_ll2$scale
km <- survfit(Surv(exp(epsilonhat), 1 - alive) ~ 1)

plot(km)
curve(1 / (1 + x), col = "red", add = TRUE)
title(main="Plot of error term (with only wall motion score)")

## fit cox
fit_reduced <- coxph(Surv(survival, 1-alive) ~ 1, data = cardio)
summary(fit_reduced)
fit_1 <- coxph(Surv(survival, 1-alive) ~ wallmotion.score, data = cardio)
anova(fit_reduced, fit1)
summary(fit_1)

Hhat <- basehaz(fit_reduced, centered = FALSE)
plot(Hhat$time, Hhat$hazard)




