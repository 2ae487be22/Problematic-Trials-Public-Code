#### Analysis of retrospective cohort

# R version - analysed on 4.3.1
R.Version()$version.string

# Remove previous variables
rm(list = ls())

#### Import libraries ----
library(here)
library(data.table)
library(ggplot2)
library(dplyr)
library(lubridate)
library(table1)
library(survival)
library(ggsurvfit)
library(VGAM)


#### Import data ----

here::i_am("RetrospectiveCohortAnalysis.R")
Data1 <- read.csv(here("CombinedRetrospectiveDatabaseV2_anon.csv"))
Data1 <- data.table(Data1)

#### Data Cleaning ----

Data1[is.na(IPD.x) , IPD.x := 0]

Data1[, Zombie := factor(Dodgy.1, labels = c("No", "Yes"))]
Data1[, Concern := 0]
Data1[!is.na(IPD.y) , Concern := 1]
Data1[, Concern := factor(Concern, labels = c("No", "Yes"))]
Data1[X1y_pval < 0.05, Significant := "Yes"]
Data1[X1y_pval >= 0.05, Significant := "No"]
Data1[, Significant := factor(Significant)]

Data1[, Num_centres := factor(Num_centres)]

Data1[, DelaySubmission := dmy(Date_Submission) - dmy(End_recruitment)]
Data1[DelaySubmission < 0, DelaySubmission := NA]

Data1[Num_centres == "1", Num_centres_factor := "1"]
Data1[Num_centres != "1", Num_centres_factor := ">=2"]

Data1[is.na(Funding), Funding := "None"]

Data1[is.na(Conflict_of_interest), Conflict_of_interest := "No"]

Data1[is.na(Docs_available), Docs_available := "No"]

Data1[is.na(Code_available), Code_available := "No"]

Data1[(X1y_direction == "Expected" | X1y_direction == "Expected") & Significant == "No",]
Data1[X1y_direction == "Expected" & X1y_pval == 0.05,]

Data1[is.na(Analysis_match), Analysis_match := "NED"]

Data1$RoB_Randomisation <- factor(Data1$RoB_Randomisation, levels = c("High", "Some", "Low"))
Data1$RoB_Deviations <- factor(Data1$RoB_Deviations, levels = c("High", "Some", "Low"))
Data1$RoB_Missing <- factor(Data1$RoB_Missing, levels = c("High", "Some", "Low"))
Data1$RoB_Measurement <- factor(Data1$RoB_Measurement, levels = c("High", "Some", "Low"))
Data1$RoB_Selection <- factor(Data1$RoB_Selection, levels = c("High", "Some", "Low"))
Data1$RoB_Overall <- factor(Data1$RoB_Overall, levels = c("High", "Some", "Low"))

Data1[Outcome != "Accept", Accept := "No"]
Data1[Outcome == "Accept", Accept := "Yes"]
Data1$Accept <- as.factor(Data1$Accept)

Data1[is.na(X1y_match), X1y_match := "Missing"]

#### Summary ----

dim(Data1[Concern == "Yes"])

table1(~ Pre_reg + X1y_match + X2y_match + Analysis_match + Data_available + Code_available + Docs_available + Funding + Conflict_of_interest +
         Study_size + X1y_pval + Significant + X1y_direction + Num_centres_factor + Num_authors + as.numeric(DelaySubmission)| Zombie, data = Data1,
       render.continuous=c(.="Mean (SD)", .="Median (Q1 - Q3 [Min - Max])"))


# IPD available
dim(Data1[IPD.x == 1,])
dim(Data1[Concern == "Yes" & IPD.x != 1,])
dim(Data1[Zombie == "Yes" & IPD.x != 1,])


#### Publication bias ----

table1(~ Outcome | Significant, data = Data1[!is.na(Significant)])

chisq.test(table(Data1[!is.na(Significant)]$Accept, Data1[!is.na(Significant)]$Significant))

mylogit <- glm(Accept ~ Significant, data = Data1, family = "binomial")
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

# p-value distribution

Data1[!is.na(X1y_pval), .(mean = mean(X1y_pval), Quantiles = quantile(X1y_pval, prob=c(.25,.5,.75))), by = Dodgy.1]
wilcox.test(X1y_pval~Zombie, data = Data1)
wilcox.test(X1y_pval~Concern, data = Data1)

# p-value against Zombie
table(Data1$Significant,Data1$Zombie)

mylogit <- glm(Significant ~ Zombie, data = Data1, family = "binomial")
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

# Delay for significant/Zombie

table1(~ X1y_pval + Accept + as.numeric(DelaySubmission) |Significant, data = Data1[!is.na(Significant)])

wilcox.test(as.numeric(Data1$DelaySubmission) ~ Data1$Significant)
wilcox.test(as.numeric(Data1$DelaySubmission) ~ Data1$Zombie)

mysurv <- survfit(Surv(DelaySubmission, rep(1, dim(Data1[!is.na(DelaySubmission)])[[1]] )) ~ Significant, data = Data1[!is.na(DelaySubmission)]) 
ggsurvfit(mysurv)
mysurv2 <- coxph(Surv(DelaySubmission, rep(1, dim(Data1[!is.na(DelaySubmission)])[[1]] )) ~ Significant, data = Data1[!is.na(DelaySubmission)]) 
summary(mysurv2)

hist(as.numeric(Data1$DelaySubmission))

mysurv <- survfit(Surv(DelaySubmission, rep(1, dim(Data1[!is.na(DelaySubmission)])[[1]] )) ~ Zombie, data = Data1[!is.na(DelaySubmission)]) 
ggsurvfit(mysurv)
mysurv2 <- coxph(Surv(DelaySubmission, rep(1, dim(Data1[!is.na(DelaySubmission)])[[1]] )) ~ Zombie, data = Data1[!is.na(DelaySubmission)]) 
summary(mysurv2)

# Direction of effect

chisq.test(table(Data1[!is.na(X1y_direction)]$X1y_direction, Data1[!is.na(X1y_direction)]$Zombie))

#### Study quality ----

# sample size
sample_size_boxplot <- ggplot(Data1, aes(x = Zombie, y = Study_size)) + geom_boxplot()

mylm <- glm(Study_size ~ Zombie, data = Data1)
summary(mylm)
cbind(OR = coef(mylm), confint(mylm))

# 1y registration

chisq.test(table(Data1$X1y_match, Data1$Zombie))

# sensitivity test to excluding missing values

Data1[, Successful_1y_match := "No"]
Data1[X1y_match == "Yes", Successful_1y_match := "Yes"]

chisq.test(table(Data1$Successful_1y_match, Data1$Zombie))

mylogit <- glm(as.factor(Successful_1y_match) ~ Zombie, data = Data1, family = "binomial")
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))


# 2y registration

Data1[is.na(X2y_match), X2y_match := "Missing"]
chisq.test(table(Data1$X2y_match, Data1$Zombie))

# sensitivity test to excluding missing values
Data1[, Successful_2y_match := "No"]
Data1[X2y_match == "Yes", Successful_2y_match := "Yes"]

chisq.test(table(Data1$Successful_2y_match, Data1$Zombie))

mylogit <- glm(as.factor(Successful_2y_match) ~ Zombie, data = Data1, family = "binomial")
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

# Successful primary and secondary registration

Data1[, .N, by = .(Successful_1y_match, Successful_2y_match, Zombie)]
Data1[, .N, by = .(Successful_1y_match, Successful_2y_match, Accept)]


Data1[,Successful_full_reg := "No"]
Data1[Successful_1y_match == "Yes" & Successful_2y_match == "Yes",Successful_full_reg := "Yes"]

mylogit <- glm(as.factor(Successful_full_reg) ~ Zombie, data = Data1, family = "binomial")
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

# Analysis
chisq.test(table(Data1$Analysis_match, Data1$Zombie))

# sensitivity test to excluding missing values
Data1[, Successful_Analysis_match := "No"]
Data1[Analysis_match == "Yes", Successful_Analysis_match := "Yes"]

chisq.test(table(Data1$Successful_Analysis_match, Data1$Zombie))

mylogit <- glm(as.factor(Successful_Analysis_match) ~ Zombie, data = Data1, family = "binomial")
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

# Authors

wilcox.test(Data1$Num_authors ~ Data1$Zombie)

# Zero-truncated poisson
m1 <- vglm(Num_authors ~ Zombie, family = pospoisson(), data = Data1)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

# Results similar to standard poisson
mypois <- glm(Num_authors ~ Zombie, data = Data1, family = "poisson")
summary(mypois)
exp(cbind(OR = coef(mypois), confint(mypois)))



# Number of centres

wilcox.test(as.numeric(Data1$Num_centres) ~ Data1$Zombie)

# Conflict of interest

chisq.test(table(Data1$Conflict_of_interest, Data1$Zombie))

# funding


chisq.test(table(Data1$Funding, Data1$Zombie))



##### Histograms ----

# pvals_histogram <- ggplot(data= Data1, aes(x = X1y_pval, colour = Zombie)) + 
#   geom_histogram(fill = "white") + geom_vline(aes(xintercept = 0.05), linetype="dashed", linewidth=1) +
#   xlab("Primary outcome p-value") + ylab("Frequency") + labs(colour = "Zombie") + 
#   theme(legend.position="none") +
#   theme(axis.line = element_line(colour = "black"),
#     #panel.grid.minor = element_blank(),
#         #panel.border = element_blank(),
#         panel.background = element_blank())

pvals_histogram2 <- ggplot(data= Data1, aes(x = X1y_pval, fill = Zombie)) + 
  geom_histogram(colour = "black") + geom_vline(aes(xintercept = 0.05), linetype="dashed", linewidth=1) +
  xlab("Primary outcome p-value") + ylab("Frequency") + labs(colour = "Zombie") + 
  theme(legend.position="none") +
  theme(axis.line = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank())

ggsave("P-value_histogram.pdf", pvals_histogram2, dpi = 600, device = "pdf")

# pvals_histogram_dens <- ggplot(data= Data1, aes(x = X1y_pval, colour = Zombie)) + 
#   geom_histogram(aes(y = ..density..), fill = "white") + geom_vline(aes(xintercept = 0.05), linetype="dashed", linewidth=1)
# 
# ggplot(data= Data1, aes(x = X1y_pval, colour = Zombie)) + geom_density()

pvals_histogram3 <- ggplot(data = Data1, aes(x = X1y_pval, colour = Zombie)) + 
  stat_ecdf(pad = FALSE, linewidth = 1.5) + geom_vline(aes(xintercept = 0.05), linetype="dashed", linewidth=1) +
  xlab("Primary outcome p value") + ylab("Cumulative frequency") + labs(colour = "Zombie") + 
  theme(legend.position="none") +
  theme(axis.line.x = element_line(colour = "grey30", linewidth = 1, linetype='solid'),
        axis.line.y = element_line(colour = "grey30", linewidth = 1, linetype='solid'),
        panel.background = element_blank()
  )

ggsave("P-value_cumulative.pdf", pvals_histogram3, dpi = 600, device = "pdf")


##### Problems in trials -----

table1( ~ Discrepancies + False_results + False_data + Changed.resubmission + Plagiarism | Zombie, data = Data1[Concern == "Yes"], render.categorical=my.render.cat)

dim(Data1[Concern == "Yes" & False_data == "Y" & False_results == "N" & Discrepancies == "N" & (Changed.resubmission == "N" | is.na(Changed.resubmission)) & Plagiarism == "N",])
dim(Data1[Zombie == "Yes" & False_data == "Y" & False_results == "N" & Discrepancies == "N" & (Changed.resubmission == "N" | is.na(Changed.resubmission)) & Plagiarism == "N",])

dim(Data1[Concern == "Yes" & False_data == "Y" & False_results == "Y" & Discrepancies == "N" & (Changed.resubmission == "N" | is.na(Changed.resubmission)) & Plagiarism == "N",])
dim(Data1[Zombie == "Yes" & False_data == "Y" & False_results == "Y" & Discrepancies == "N" & (Changed.resubmission == "N" | is.na(Changed.resubmission)) & Plagiarism == "N",])


Data1[Concern == "Yes", .N, by = .(False_data, False_results, Plagiarism, Changed.resubmission, Discrepancies)]

##### Risk of bias ----

table1( ~ RoB_Randomisation + RoB_Deviations + RoB_Missing + RoB_Measurement + RoB_Selection + RoB_Overall | Zombie, data = Data1[Concern == "Yes"], render.categorical=my.render.cat)

##### Combined table for publication

table1( ~ Discrepancies + False_results + False_data + Changed.resubmission + Plagiarism + RoB_Overall | Zombie, data = Data1[Concern == "Yes"], render.categorical=my.render.cat)

##### Risk of bias by publication status

table1(~ Pre_reg + X1y_match + X2y_match + Analysis_match + Data_available + Code_available + Docs_available + Funding + Conflict_of_interest +
         Study_size + X1y_pval + Significant + X1y_direction + Num_centres_factor + Num_authors + as.numeric(DelaySubmission)| Accept, data = Data1,
       render.continuous=c(.="Mean (SD)", .="Median (Q1 - Q3 [Min - Max])"), render.categorical=my.render.cat)