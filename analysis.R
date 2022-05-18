library(tidyverse)
library(cowplot)
library(pROC)
library(lme4)
library(MASS)


#Code written by Greg Courtice, last updated September 16, 2021

#This multilevel logistic regression dose-response model using data from salmonid dose-response observations
#follows the work of courtice et al. (2022) (see the code for this work here https://github.com/gregcourtice/SS_Effects) 
#The model developed here separates ln(SSC) and ln(DoE), as opposed to modeling them together in a single parameter, ln(SSD)
#This model is used to quantify prediction uncertainty and evaluate the use of suspended sediment dose (SSD) as a management limit



#read in dose response dataset and filter only adult and juvenile salmonids, remove habitat-related observations
data <- read_csv("data.csv") %>%
  filter(Group %in% c("AS", "JS")) %>%
  filter(Habitat_Related == "No")

#more data tidying for multilevel modeling
data <- mutate(data, Study_Year = ifelse(Group == "AS", paste0(Study, Year), paste0(Study, Year, "J"))) %>%
  group_by(Study_Year) %>%
  mutate(NStudy = n()) %>%
  ungroup()

data <- mutate(data,
               LS = Group,
               Log_Conc = log(Concentration),
               Log_Dur = log(Duration),
               Dose = Log_Conc + Log_Dur,
               Effect = ifelse(SEV < 5, 1, ifelse(SEV > 4 & SEV < 8, 2, 3)),
               Effect1 = ifelse(SEV < 5, 0, 1),
               Effect2 = ifelse(SEV < 8, 0, 1),
               Species_NA = ifelse(is.na(Species), "Not Available", Species),
               Species_group = ifelse(Group == "AS", Species_NA, paste0(Species_NA, "J"))
)

#remove minor physiological effects due to ambiguous effect outcomes for logistic regression modeling
#see Courtice et al. (2022) and this paper for discussion on this matter.
m.data <- data%>%filter(Effect != 2)
#redeclare parameters for consistency
m.data <- m.data %>% mutate(ln.SSC = Log_Conc,
                            ln.DoE = Log_Dur,
                            ln.SSD = ln.SSC + ln.DoE,
                            #define low, medium, and high DoE bins to be used in figure 2 (ln(DoE)=0 is approx 1 hour, while ln(DoE)=5 is approx 6 days)
                            ln.DoE_bin = case_when(
                              ln.DoE <= 0 ~ "Low",
                              ln.DoE >0 & ln.DoE <= 5 ~ "Medium",
                              ln.DoE > 5 ~ "High"
                            ))

#Create axis tick labels
dur_labs <- c(
  "1.2 minutes",
  "10 minutes",
  "1 hour",
  "6 hours",
  "1 day",
  "1 week",
  "1 month",
  "6 months",
  "2 years"
)

conc_labs <- c("0.7",
               "5",
               "10",
               "25",
               "50",
               "100",
               "250",
               "1 000",
               "5 000",
               "50 000",
               "207 000")

#create color blind friendly color palette
palette  <- c(
  "#000000",
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7"
)

######Fit Models to be used in study ############



#Courtice et al. (2022) Model with SSD as the sole predictor:
fit <- lme4::glmer(as.factor(Effect1) ~
                     ln.SSD +
                     (1 | Study_Year),
                   family = "binomial",
                   data = m.data)

#New model, separating ln(SSC) and ln(DoE), and adding interaction:
fit.int <- lme4::glmer(
  as.factor(Effect1) ~
    ln.SSC + ln.DoE + ln.SSC:ln.DoE +
    (1 | Study_Year),
  family = "binomial",
  data = m.data
)

#Similar model to new one, but without an interaction term:
fit.no_int <- lme4::glmer(
  as.factor(Effect1) ~
    ln.SSC + ln.DoE +
    (1 | Study_Year),
  family = "binomial",
  data = m.data
)

#Defining the correlation matrices of new model and model without interaction

#new model
corMat.int <- vcov(fit.int)
mu.int <- c(fixef(fit.int)[1],
            fixef(fit.int)[2],
            fixef(fit.int)[3],
            fixef(fit.int)[4])

#model without interaction
corMat.no_int <- vcov(fit.no_int)
mu.no_int <- c(fixef(fit.no_int)[1],
               fixef(fit.no_int)[2],
               fixef(fit.no_int)[3])


###### Goodness of Fit Analysis #########
#ROC analysis: predict probabilities using new model with interaction, and then classify using roc()
pred.m.int <- m.data%>%mutate(logit = 
                                mu.int[1] + 
                                mu.int[2]*Log_Conc + 
                                mu.int[3]*Log_Dur + 
                                mu.int[4]*Log_Conc*Log_Dur,
                              prob = 1/(1+exp(-logit)))

ROC.m.int <- roc(pred.m.int,
                 Effect1,
                 prob
                 )
auc(ROC.m.int)



############# Monte Carlo simulation###########
#Simulate parameters from multivariate distribution
#estimate probability associated with CCME (2003) benchmark limit over range of SSC and DoE
#draw 10 000 parameter sets from multivariate normal distribution
#estimate probability based on parameters at SSC=25mg/L and DoE=24h (the benchmark limit)

set.seed(1)
dat1.int <-               
  as_tibble(mvrnorm(
    n = 10000,
    mu = mu.int,
    Sigma = corMat.int,
    empirical = FALSE
  )) %>%
  rename(
    alpha = `(Intercept)`,
    beta.SSC = ln.SSC,
    beta.DoE = ln.DoE,
    beta.int = `ln.SSC:ln.DoE`
  ) %>%
  mutate(
    draw = sample(1:10000, size = 10000, replace = FALSE),
    logit.CCME = alpha + beta.SSC * log(25) + beta.DoE * log(24) + beta.int * log(24) * log(25),
    prob.CCME = 1 / (1 + exp(-logit.CCME))
  )

#create data frame with range of ln(DoE), each ln(DoE) value has 10 000 draws
#Join the above parameter sets and probability predictions at the benchmark limit to this data frame by the draw numbers
#Each ln(DoE) value now has 10 000 parameter sets
pred.SSC.int <-
  crossing(ln_DoE = seq(from = -4, to = 10, by = 0.01),
           draw = 1:10000) %>%
  left_join(dat1.int, by = "draw")

#Using a ln(DoE) value, predict the corresponding ln(SSC) value that would be required to achieve the estimated probability at the benchmark limit
#based on the parameters used to estimate the benchmark probability
#the estimated probability is unique for each draw due to different parameters
#Also calculate the lower and upper 5% quantiles to get the 90% management band
pred.SSC.int <- pred.SSC.int %>%
  mutate(pred.ln.SSC = ifelse((logit.CCME - alpha - beta.DoE * ln_DoE) /
                                (beta.SSC + beta.int * ln_DoE) > 15,
                              NA,
                              (logit.CCME - alpha - beta.DoE * ln_DoE) / (beta.SSC + beta.int * ln_DoE)
  )) 

pred.SSC.int.quant <- pred.SSC.int%>%group_by(ln_DoE)%>%summarise(lower5 = quantile(pred.ln.SSC, probs = 0.05, na.rm = TRUE),
                                                                  upper5 = quantile(pred.ln.SSC, probs = 0.95, na.rm = TRUE))





#undertake same methods for model with no interaction
set.seed(1)
dat1.no_int <-
  as_tibble(mvrnorm(
    n = 10000,
    mu = mu.no_int,
    Sigma = corMat.no_int,
    empirical = FALSE
  )) %>%
  rename(alpha = `(Intercept)`,
         beta.SSC = ln.SSC,
         beta.DoE = ln.DoE) %>%
  mutate(
    draw = sample(1:10000, size = 10000, replace = FALSE),
    logit.CCME = alpha + beta.SSC * log(25) + beta.DoE * log(24),
    prob.CCME = 1 / (1 + exp(-logit.CCME))
  )

pred.SSC.no_int <-
  crossing(ln_DoE = seq(from = -4, to = 10, by = 0.01),
           draw = 1:10000) %>%
  left_join(dat1.no_int, by = "draw")


pred.SSC.no_int <- pred.SSC.no_int %>%
  mutate(pred.ln.SSC = ifelse((logit.CCME - alpha - beta.DoE * ln_DoE) /
                                (beta.SSC) > 15,
                              NA,
                              (logit.CCME - alpha - beta.DoE * ln_DoE) / (beta.SSC)
  ))

pred.SSC.no_int.quant <- pred.SSC.no_int%>%group_by(ln_DoE)%>%summarise(lower5 = quantile(pred.ln.SSC, probs = 0.05, na.rm = TRUE),
                                                                        upper5 = quantile(pred.ln.SSC, probs = 0.95, na.rm = TRUE))






ndraw <- 16:30  #select draws from simulation to present in figures



#Figure3: refined in Inkscape for presentation in the manuscript

p.ccme.iso.int <- ggplot() +
  theme_classic() +
  geom_path(
    data = pred.SSC.int.quant,
    mapping =aes(
    x = exp(ln_DoE),
    y = exp(lower5)
  ),
  color = "red",
  size = 2) +
  geom_path(
    data = pred.SSC.int.quant,
    mapping =aes(
    x = exp(ln_DoE),
    y = exp(upper5)
  ),
  color = "red",
  size = 2) +
  geom_path(
    data = filter(pred.SSC.int, draw == ndraw),
    mapping = aes(
      x = exp(ln_DoE),
      y = exp(pred.ln.SSC),
      group = draw
    )
  ) +
  geom_point(aes(x = 24,
                 y = 25),
             size = 5,
             color = "red") +
  xlab("Exposure Duration") +
  ylab(expression(atop(italic(SSC), (mg %.% L ^ -1)))) +
  scale_x_continuous(
    trans = "log",
    breaks = c(1.2 / 60, 10 / 60, 1, 6, 24, 168, 720, 24 *
                 30 * 6, 17520),
    labels = dur_labs,
    limits = c(1.2 / 60, 24 * 30)
  ) +
  scale_y_continuous(
    trans = "log",
    breaks = c(0.66, 5, 10, 25, 50, 100, 250, 1000, 5000, 50000, 207000),
    labels = conc_labs,
    limits = c(0.66, 10000)
  ) +
  geom_point(
    data = m.data,
    mapping = aes(
      x = Duration,
      y = Concentration,
      shape = as.factor(Effect1)
    ),
    size = 5
  ) +
  scale_shape_manual(
    labels = c("Behavioural   ",
               "Major Physiological and Lethal   "),
    values = c(21, 23)
  ) +
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(angle = 0,
                                    vjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))


#The following figure was not used in the paper, but is the same figure
#to the above Figure 3 though using the model without the interaction term.
p.ccme.iso.no_int <- ggplot() +
  theme_classic() +
  geom_path(
    data = pred.SSC.no_int.quant,
    mapping =aes(
      x = exp(ln_DoE),
      y = exp(lower5)
    ),
    color = "red",
    size = 2) +
  geom_path(
    data = pred.SSC.no_int.quant,
    mapping =aes(
      x = exp(ln_DoE),
      y = exp(upper5)
    ),
    color = "red",
    size = 2) +
  geom_path(
    data = filter(pred.SSC.no_int, draw == ndraw),
    mapping = aes(
      x = exp(ln_DoE),
      y = exp(pred.ln.SSC),
      group = draw
    )
  ) +
  geom_point(aes(x = 24,
                 y = 25),
             size = 5,
             color = "red") +
  xlab("Exposure Duration") +
  ylab(expression(atop(italic(SSC), (mg %.% L ^ -1)))) +
  scale_x_continuous(
    trans = "log",
    breaks = c(1.2 / 60, 10 / 60, 1, 6, 24, 168, 720, 24 *
                 30 * 6, 17520),
    labels = dur_labs,
    limits = c(1.2 / 60, 24 * 30)
  ) +
  scale_y_continuous(
    trans = "log",
    breaks = c(0.66, 5, 10, 25, 50, 100, 250, 1000, 5000, 50000, 207000),
    labels = conc_labs,
    limits = c(0.66, 10000)
  )+
  geom_point(
    data = m.data,
    mapping = aes(
      x = Duration,
      y = Concentration,
      shape = as.factor(Effect1)
    ),
    size = 5
  ) +
  scale_shape_manual(
    labels = c("Behavioural   ",
               "Major Physiological and Lethal   "),
    values = c(21, 23)
  ) +
  theme(legend.position = "none",
        axis.title.y = element_text(angle = 0,
                                    vjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))



#triptych plots (Figure 2)
#create data frame with range of ln(SSC) and 10 000 draws
#pull in parameters from simulation corresponding to draw numbers
pred.prob.int <-
  crossing(ln_SSC = seq(from = -1, to = 12, by = 0.01),
           draw = 1:10000) %>%
  left_join(dat1.int, by = "draw")

#Estimate the logit and probability for ranges of ln(SSC) for the three bins of DoE
#Use mean ln(DoE) value based on data within each bin
pred.prob.int <- pred.prob.int %>%
  mutate(pred.logit.low = alpha + ln_SSC * beta.SSC + mean(filter(m.data, ln.DoE_bin == "Low")$ln.DoE) * beta.DoE + ln_SSC * mean(filter(m.data, ln.DoE_bin == "Low")$ln.DoE) * beta.int,
         pred.prob.low = 1/(1+exp(-pred.logit.low)),
         pred.logit.med = alpha + ln_SSC * beta.SSC + mean(filter(m.data, ln.DoE_bin == "Medium")$ln.DoE) * beta.DoE + ln_SSC * mean(filter(m.data, ln.DoE_bin == "Medium")$ln.DoE) * beta.int,
         pred.prob.med = 1/(1+exp(-pred.logit.med)),
         pred.logit.high = alpha + ln_SSC * beta.SSC + mean(filter(m.data, ln.DoE_bin == "High")$ln.DoE) * beta.DoE + ln_SSC * mean(filter(m.data, ln.DoE_bin == "High")$ln.DoE) * beta.int,
         pred.prob.high = 1/(1+exp(-pred.logit.high))
         
  )

#upper and lower quantiles
pred.prob.int.q95 <- pred.prob.int%>%
  group_by(ln_SSC)%>%summarise(pred.logit.low.95lower = quantile(pred.logit.low, probs = 0.025, na.rm = TRUE),
                            pred.logit.low.95upper = quantile(pred.logit.low, probs = 0.975, na.rm = TRUE),
                            pred.logit.med.95lower = quantile(pred.logit.med, probs = 0.025, na.rm = TRUE),
                            pred.logit.med.95upper = quantile(pred.logit.med, probs = 0.975, na.rm = TRUE),
                            pred.logit.high.95lower = quantile(pred.logit.high, probs = 0.025, na.rm = TRUE),
                            pred.logit.high.95upper = quantile(pred.logit.high, probs = 0.975, na.rm = TRUE))


#mean prediction curves
pred.mean.int <- 
  crossing(ln.SSC = seq(from = -1, to = 12, by = 0.01),
           ln.DoE = c(mean(filter(m.data, ln.DoE_bin == "Low")$ln.DoE), 
                      mean(filter(m.data, ln.DoE_bin == "Medium")$ln.DoE), 
                      mean(filter(m.data, ln.DoE_bin == "High")$ln.DoE)
                      )
           )%>%
  mutate(mean.logit = fixef(fit.int)[1] + ln.SSC * fixef(fit.int)[2] + ln.DoE * fixef(fit.int)[3] + ln.SSC * ln.DoE * fixef(fit.int)[4])


#logit curves


#new model with interaction
#short duration
set.seed(3)
p.logit.DoE.low.int <-ggplot()+
  theme_classic()+
  geom_jitter(data = filter(m.data, ln.DoE_bin == "Low"),
              mapping = aes(
                x = ln.SSC,
                y = (Effect1-0.5)*100,
                shape = as.factor(Effect1)
              ),
              size = 5,
              # shape = 21,
              height = 1,
              width = 0.05
  ) +
  geom_path(
    data = pred.prob.int.q95,
    mapping = aes(
      x = ln_SSC,
      y = pred.logit.low.95lower
    ),
    size = 1.5,
    color = "darkgreen"
  )+
  geom_path(
    data = pred.prob.int.q95,
    mapping = aes(
      x = ln_SSC,
      y = pred.logit.low.95upper
    ),
    size = 1.5,
    color = "darkgreen"
  )+
  geom_path(
    data = filter(pred.prob.int, draw == ndraw),
    mapping = aes(
      x = ln_SSC,
      y = pred.logit.low,
      group = draw
    ),
    alpha = 0.2
  )+
  geom_path(
    data = filter(pred.mean.int, ln.DoE == min(ln.DoE)),
    mapping = aes(
      x = ln.SSC,
      y = mean.logit
    ),
    size = 1.5,
    color = "blue"
  )+
  ylim(-51,51)+
  scale_shape_manual(
    labels = c("Behavioural   ",
               "Major Physiological and Lethal   "),
    values = c(21, 23)
  ) +
  ggtitle(expression(""<="1h"))+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size = 16),
        legend.position = "none")+
  ylab("logit")

#medium duration
set.seed(3)
p.logit.DoE.med.int <-ggplot()+
  theme_classic()+
  geom_jitter(data = filter(m.data, ln.DoE_bin == "Medium"),
              mapping = aes(
                x = ln.SSC,
                y = (Effect1-0.5)*100,
                shape = as.factor(Effect1)
              ),
              size = 5,
              # shape = 21,
              height = 1,
              width = 0.05
  ) +
  geom_path(
    data = pred.prob.int.q95,
    mapping = aes(
      x = ln_SSC,
      y = pred.logit.med.95lower
    ),
    size = 1.5,
    color = "darkgreen"
  )+
  geom_path(
    data = pred.prob.int.q95,
    mapping = aes(
      x = ln_SSC,
      y = pred.logit.med.95upper
    ),
    size = 1.5,
    color = "darkgreen"
  )+
  geom_path(
    data = filter(pred.prob.int, draw == ndraw),
    mapping = aes(
      x = ln_SSC,
      y = pred.logit.med,
      group = draw
    ),
    alpha = 0.2
  )+
  geom_path(
    data = filter(pred.mean.int, ln.DoE == median(ln.DoE)),
    mapping = aes(
      x = ln.SSC,
      y = mean.logit
    ),
    size = 1.5,
    color = "blue"
  )+
  scale_shape_manual(
    labels = c("Behavioural   ",
               "Major Physiological and Lethal   "),
    values = c(21, 23)
  ) +
  ggtitle("1h to 6d") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 20),
        axis.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.text.y.left = element_blank(),
        legend.position = "none")

#long duration
set.seed(3)
p.logit.DoE.high.int <- ggplot()+
  theme_classic()+
  geom_jitter(data = filter(m.data, ln.DoE_bin == "High"),
              mapping = aes(
                x = ln.SSC,
                y =(Effect1-0.5)*100,
                shape = as.factor(Effect1)
              ),
              size = 5,
              # shape = 21,
              height = 1,
              width = 0.05
  ) +
  geom_path(
    data = pred.prob.int.q95,
    mapping = aes(
      x = ln_SSC,
      y = pred.logit.high.95lower
    ),
    size = 1.5,
    color = "darkgreen"
  )+
  geom_path(
    data = pred.prob.int.q95,
    mapping = aes(
      x = ln_SSC,
      y = pred.logit.high.95upper
    ),
    size = 1.5,
    color = "darkgreen"
  )+
  geom_path(
    data = filter(pred.prob.int, draw == ndraw),
    mapping = aes(
      x = ln_SSC,
      y = pred.logit.high,
      group = draw
    ),
    alpha = 0.2
  )+
  geom_path(
    data = filter(pred.mean.int, ln.DoE == max(ln.DoE)),
    mapping = aes(
      x = ln.SSC,
      y = mean.logit
    ),
    size = 1.5,
    color = "blue"
  )+
  scale_shape_manual(
    labels = c("Behavioural   ",
               "Major Physiological and Lethal   "),
    values = c(21, 23)
  ) +
  ggtitle(expression("">"6d"))  +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 20),
        axis.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.text.y = element_blank(),
        axis.line.y.right = element_blank(),
        legend.position = "none")+
  scale_y_continuous(sec.axis = dup_axis(breaks = c(-50,50),
                                       labels = c("BHV",
                                                  "MaPL"))
                   
)


#model without interaction
#same methods as model with interaction
pred.prob.no_int <-
  crossing(ln_SSC = seq(from = -1, to = 12, by = 0.01),
           draw = 1:10000) %>%
  left_join(dat1.no_int, by = "draw")


pred.prob.no_int <- pred.prob.no_int %>%
  
  
  mutate(pred.logit.low = alpha + ln_SSC * beta.SSC + mean(filter(m.data, ln.DoE_bin == "Low")$ln.DoE) * beta.DoE,
         pred.prob.low = 1/(1+exp(-pred.logit.low)),
         pred.logit.med = alpha + ln_SSC * beta.SSC + mean(filter(m.data, ln.DoE_bin == "Medium")$ln.DoE) * beta.DoE,
         pred.prob.med = 1/(1+exp(-pred.logit.med)),
         pred.logit.high = alpha + ln_SSC * beta.SSC + mean(filter(m.data, ln.DoE_bin == "High")$ln.DoE) * beta.DoE,
         pred.prob.high = 1/(1+exp(-pred.logit.high))
                            
  )


pred.prob.no_int.q95 <- pred.prob.no_int%>%
  group_by(ln_SSC)%>%summarise(pred.logit.low.95lower = quantile(pred.logit.low, probs = 0.025, na.rm = TRUE),
                               pred.logit.low.95upper = quantile(pred.logit.low, probs = 0.975, na.rm = TRUE),
                               pred.logit.med.95lower = quantile(pred.logit.med, probs = 0.025, na.rm = TRUE),
                               pred.logit.med.95upper = quantile(pred.logit.med, probs = 0.975, na.rm = TRUE),
                               pred.logit.high.95lower = quantile(pred.logit.high, probs = 0.025, na.rm = TRUE),
                               pred.logit.high.95upper = quantile(pred.logit.high, probs = 0.975, na.rm = TRUE))

#mean prediction curves
pred.mean.no_int <- 
  crossing(ln.SSC = seq(from = -1, to = 12, by = 0.01),
           ln.DoE = c(mean(filter(m.data, ln.DoE_bin == "Low")$ln.DoE), 
                      mean(filter(m.data, ln.DoE_bin == "Medium")$ln.DoE), 
                      mean(filter(m.data, ln.DoE_bin == "High")$ln.DoE)
           )
  )%>%
  mutate(mean.logit = fixef(fit.no_int)[1] + ln.SSC * fixef(fit.no_int)[2] + ln.DoE * fixef(fit.no_int)[3])



#logit curves for model without interaction
#short duration
  set.seed(3)
  p.logit.DoE.low.no_int <- ggplot()+
      theme_classic()+
    geom_jitter(data = filter(m.data, ln.DoE_bin == "Low"),
                mapping = aes(
                  x = ln.SSC,
                  y = (Effect1-0.5)*100,
                  shape = as.factor(Effect1)
                ),
                size = 5,
                height = 1,
                width = 0.05
    ) +
    geom_path(
      data = pred.prob.no_int.q95,
      mapping = aes(
        x = ln_SSC,
        y = pred.logit.low.95lower
      ),
      size = 1.5,
      color = "darkgreen"
    )+
    geom_path(
      data = pred.prob.no_int.q95,
      mapping = aes(
        x = ln_SSC,
        y = pred.logit.low.95upper
      ),
      size = 1.5,
      color = "darkgreen"
    )+
    geom_path(
      data = filter(pred.prob.no_int, draw == ndraw),
      mapping = aes(
        x = ln_SSC,
        y = pred.logit.low,
        group = draw
      ),
      alpha = 0.2
    )+
    geom_path(
      data = filter(pred.mean.no_int, ln.DoE == min(ln.DoE)),
      mapping = aes(
        x = ln.SSC,
        y = mean.logit
      ),
      size = 1.5,
      color = "blue"
    )+
    scale_shape_manual(
      labels = c("Behavioural   ",
                 "Major Physiological and Lethal   "),
      values = c(21, 23)
    ) +
    xlab("ln(SSC)")+
    ylab("logit")+
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 20),
          axis.text.y = element_blank(),
          axis.text = element_text(size = 16),
          legend.position = "none")
  
  
  #medium duration
  set.seed(3)
  p.logit.DoE.med.no_int <-ggplot()+
    theme_classic()+
    geom_jitter(data = filter(m.data, ln.DoE_bin == "Medium"),
                mapping = aes(
                  x = ln.SSC,
                  y = (Effect1-0.5)*100,
                  shape = as.factor(Effect1)
                ),
                size = 5,
                height = 1,
                width = 0.05
    ) +
    geom_path(
      data = pred.prob.no_int.q95,
      mapping = aes(
        x = ln_SSC,
        y = pred.logit.med.95lower
      ),
      size = 1.5,
      color = "darkgreen"
    )+
    geom_path(
      data = pred.prob.no_int.q95,
      mapping = aes(
        x = ln_SSC,
        y = pred.logit.med.95upper
      ),
      size = 1.5,
      color = "darkgreen"
    )+
    geom_path(
      data = filter(pred.prob.no_int, draw == ndraw),
      mapping = aes(
        x = ln_SSC,
        y = pred.logit.med,
        group = draw
      ),
      alpha = 0.2
    )+
    geom_path(
      data = filter(pred.mean.no_int, ln.DoE == median(ln.DoE)),
      mapping = aes(
        x = ln.SSC,
        y = mean.logit
      ),
      size = 1.5,
      color = "blue"
    )+
    scale_shape_manual(
      labels = c("Behavioural   ",
                 "Major Physiological and Lethal   "),
      values = c(21, 23)
    ) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 20),
          axis.text = element_text(size = 16),
          axis.text.y = element_blank(),
          legend.position = "none")+
    xlab("ln(SSC)")
  
  #long duration
  set.seed(3)
  p.logit.DoE.high.no_int <- ggplot()+
      theme_classic()+
    geom_jitter(data = filter(m.data, ln.DoE_bin == "High"),
                mapping = aes(
                  x = ln.SSC,
                  y =(Effect1-0.5)*100,
                  shape = as.factor(Effect1)
                ),
                size = 5,
                height = 1,
                width = 0.05
    ) +
    geom_path(
      data = pred.prob.no_int.q95,
      mapping = aes(
        x = ln_SSC,
        y = pred.logit.high.95lower
      ),
      size = 1.5,
      color = "darkgreen"
    )+
    geom_path(
      data = pred.prob.no_int.q95,
      mapping = aes(
        x = ln_SSC,
        y = pred.logit.high.95upper
      ),
      size = 1.5,
      color = "darkgreen"
    )+
    geom_path(
      data = filter(pred.prob.no_int, draw == ndraw),
      mapping = aes(
        x = ln_SSC,
        y = pred.logit.high,
        group = draw
      ),
      alpha = 0.2
    )+
    geom_path(
      data = filter(pred.mean.no_int, ln.DoE == max(ln.DoE)),
      mapping = aes(
        x = ln.SSC,
        y = mean.logit
      ),
      size = 1.5,
      color = "blue"
    )+
    scale_shape_manual(
      labels = c("Behavioural   ",
                 "Major Physiological and Lethal   "),
      values = c(21, 23)
    ) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 20),
          axis.text = element_text(size = 16),
          axis.text.y = element_blank(),
          axis.line.y.right = element_blank(),
          legend.position = "none")+
    scale_y_continuous(sec.axis = dup_axis(breaks = c(-50,50),
                                           labels = c("BHV",
                                                      "MaPL"))
                       
                                           )+
    xlab("ln(SSC)")
    
    
#plot figure 2 - Inkscape was used to refine figure
plot_grid(p.logit.DoE.low.int,    p.logit.DoE.med.int,    p.logit.DoE.high.int,  
          p.logit.DoE.low.no_int, p.logit.DoE.med.no_int, p.logit.DoE.high.no_int, 
          ncol = 3)
  


#calculate SSC-DoE pairs for dose values
dose <- tibble(DoE = seq(from = 0.1, to = 48, by = 0.001),
               SSC.10 = 24*10 / DoE,
               SSC.25 = 24*25 / DoE,
               SSC.50 = 24*50 / DoE,
               SSC.100 = 24*100 / DoE,
               SSC.175 = 24*175 / DoE,
               SSC.250 = 24*250 / DoE
)%>%
  pivot_longer(!DoE, 
               names_to = "SSCGrp",
               values_to = "SSC")
    

#Figure 4 - Inkscape was used to create figure from this base graph
ggplot(filter(pred.SSC.int.quant, lower5 <=550 &upper5 <=550&exp(ln_DoE)<=25)) +
  theme_classic() +
  geom_path(aes(
    x = exp(ln_DoE),
    y = exp(lower5)
  ),
  color = "red",
  size = 2) +
  geom_path(aes(
    x = exp(ln_DoE),
    y = exp(upper5)
  ),
  color = "red",
  size = 2) +
  geom_path(
    data = filter(dose, SSC<=550&DoE<=25),
    mapping = aes(
      x = DoE,
      y = SSC,
      group = SSCGrp
    ),
    linetype = "dotted"
  )+
  geom_point(aes(x = 24,
                 y = 25),
             size = 5,
             color = "red") +
  xlab("Exposure Duration") +
  ylab(expression(atop(italic(SSC), (mg %.% L ^ -1)))) +
  coord_cartesian(xlim = c(0, 24),
                  ylim = c(0, 500))+
  theme(axis.title.y = element_text(angle = 0,
                                    vjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))









