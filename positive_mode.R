library(car)
library(readr)
library(lattice)
library(nlme)
library(ggplot2)
library(GGally)
library(nnet)
library(foreign)
library(biotools)
library(glmmML)
library(MASS)
library(lme4)
library(multcomp)
library(dplyr)
library(knitr)
library(xtable)
library(kableExtra)
library(DT)
library(glmnet)
library(corrplot)
library(ggpubr)
library(lmerTest)
library("merTools")
library(reshape2)
library(ggplot2)
library(GGally)
library(mgcv)
library(gplots)
library(tidyr)
library(bkmr)
library(factoextra) 
library(spatstat)
library(Hmisc)
library(gtsummary)
library(blme)
library(grpreg)
library(robustHD)
library(gWQS)
library(gridExtra)
library(ggcorrplot)
library(BWQS)
library(qwraps2)
library(MatchIt)
library(data.table)
library(DescTools)
library(mice)
library(ggrepel)

##########################################################################################

setwd("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot")
d <- as.data.frame(read.csv("All_mached_pairs_and_epi_data.csv"))
d <- d[!(colnames(d) %in% "X"),]
selected_case_control_final <- read.csv("Selected_Cases_Controls_Final List 8-13-2020.csv")
d <- d[d$casemasked_mrn %in% selected_case_control_final$casemasked_mrn,]
d <- merge(d, selected_case_control_final, by = "casemasked_mrn")
d <- d[,!(colnames(d) %in% c("contmasked_mrn.y","PLATE","X"))]
d <- d %>% rename(contmasked_mrn = contmasked_mrn.x)


data_case <- as.data.frame(d[,colnames(d)[c(1:97,194)]])
new_colnames <- rep(NA_character_,length(colnames(d)[c(1:97,194)]))
for(i in 1:length(colnames(d)[c(1:97,194)])){new_colnames[i] = strsplit(colnames(d)[c(1:97,194)][1:length(colnames(d)[c(1:97,194)])],"case")[[i]][2]}
colnames(data_case) <- new_colnames

#################################################################################

## Epi data

data_control <- as.data.frame(d[,c("contmasked_mrn",colnames(d)[c(99:193)],"caseDID","contDID")])
colnames(data_control) <- new_colnames

data_all <- read.csv("merged_pfas_epi_liver_data.csv")
data_all <- data_all[,!(colnames(data_all) %in% c("NAFLD__status"))]
data_all$SAMPLEID <- as.character(data_all$SAMPLEID) 


ff1 <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/ff1_matching.csv")
dim(ff1)

ff1_prob <- ff1[!is.na(ff1$td2_case_all),]
dim(ff1_prob)

ff1_prob <- ff1_prob[ff1_prob$self_reported_race == "African American" |ff1_prob$self_reported_race == "European American" 
                     |ff1_prob$self_reported_race == "Hispanic", ]
dim(ff1_prob)

ff1_prob <- ff1_prob[ff1_prob$gender == "Female"| ff1_prob$gender == "Male", ]
dim(ff1_prob)

ff1_prob <- ff1_prob[ff1_prob$td2_case_incident == 0| ff1_prob$td2_case_incident == 2, ]
dim(ff1_prob)

ff1_prob$participant <- ff1_prob$masked_mrn %in% data_all$masked_mrn + 0


# estimation of denominator of ip weights

denom.fit <- glm(participant ~ as.character(self_reported_race) + as.factor(td2_case_all) + age_at_enrollment + gender, 
                 family = binomial(), data = ff1_prob)
summary(denom.fit)

pd.qsmk <- predict(denom.fit, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(participant~1, family = binomial(), data = ff1_prob)
summary(numer.fit)

pn.qsmk <- predict(numer.fit, type = "response")

ff1_prob$sw <- ifelse(ff1_prob$participant == 0, ((1-pn.qsmk)/(1-pd.qsmk)),
                      (pn.qsmk/pd.qsmk))

summary(ff1_prob$sw)

rrt <- ff1_prob[ff1_prob$masked_mrn %in% data_all$masked_mrn, c("sw","masked_mrn")]

# Final Inverse probability weights adjusted for case-control design

data_all$sw <- rrt[order(rrt$masked_mrn,data_all$masked_mrn),"sw"]

##################################################################################
## Positive Mode

mapfile <- read.delim("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/HRM Data/HRE0013_BioMe_AllModes_Mapfile_14May21.txt")
refmet_pos <- read.delim("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/HRM Data/RefMet_C18-HILICpos_Targets_7.5ppm_10percent-Studies_14May21.txt")

# All metabolites were reported in at least 14 
min(refmet_pos$Max_number_studies_reported)
refmet_pos <- refmet_pos[,colnames(refmet_pos) %in% c(colnames(refmet_pos)[1:17],mapfile$Sample.ID[mapfile$Sample_Class == "Study_Sample"])]
# rownames(refmet_pos) <- refmet_pos$refmet_name_mode_id

y <- apply(refmet_pos[,c(18:376)], 1, function(x){sum(as.numeric(x) > 0)/length(as.numeric(x))})
refmet_pos <- refmet_pos[which(y>0.75),] # metabolites with at least 75% non-zero values

d_pos <- refmet_pos[,c(3,18:376)]
uniq <- unique(d_pos$refmet_name)
ccv <- NA_character_
for(i in 1:length(uniq)){
  ccv <- c(ccv, (names(which.max(apply(d_pos[which(d_pos$refmet_name == uniq[i]), -c(1)], 1, function(x){sd(as.numeric(x), na.rm = T)/mean(as.numeric(x), na.rm = T)})))))
}
ccv <- ccv[-1]
refmet_pos <- refmet_pos[ccv,]
refmet_pos$met <- paste0("pos",seq(1,dim(refmet_pos)[1]))
rownames(refmet_pos) <- paste0("pos",seq(1,dim(refmet_pos)[1]))

# transpose
t_refmet_pos <- data.table::transpose(refmet_pos)

# get row and colnames in order
colnames(t_refmet_pos) <- rownames(refmet_pos)
rownames(t_refmet_pos) <- colnames(refmet_pos)
t_refmet_pos <- t_refmet_pos[!(rownames(t_refmet_pos) %in% c(rownames(t_refmet_pos)[1:17],"met")),]

t_refmet_pos$SAMPLEID <- rownames(t_refmet_pos)
for(i in 1:length(t_refmet_pos$SAMPLEID)){
  if(length(strsplit(t_refmet_pos$SAMPLEID,"AB")[i][[1]]) == 2){
    t_refmet_pos$SAMPLEID[i] = (strsplit(t_refmet_pos$SAMPLEID,"AB")[i][[1]])[2]
  }
  else if (length(strsplit(t_refmet_pos$SAMPLEID,"X")[i][[1]]) == 2){
    t_refmet_pos$SAMPLEID[i] = (strsplit(t_refmet_pos$SAMPLEID,"X")[i][[1]])[2]
  }
  
}

t_refmet_pos <- apply(t_refmet_pos, 2, as.numeric)
t_refmet_pos <- as.data.frame(t_refmet_pos)
t_refmet_pos$SAMPLEID <- as.character(t_refmet_pos$SAMPLEID)


chunk_18 <- apply(t_refmet_pos[,!(colnames(t_refmet_pos) %in% c("SAMPLEID"))], 2, function(x)(log(x + 1, base = 2)))
chunk_18 <- as.data.frame(chunk_18)
chunk_18 <- as.data.frame(scale(chunk_18))
chunk_18$SAMPLEID <-  t_refmet_pos$SAMPLEID
t_refmet_pos <- as.data.frame(chunk_18)
metabolites_refmet_pos <- colnames(t_refmet_pos)[colnames(t_refmet_pos) != "SAMPLEID"]

data_all__pos <- merge(data_all,t_refmet_pos, by = "SAMPLEID", all.x	= T, all.y = T)

## Merge with more updated New PFAS - August 2021

aug_pfas <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/New PFAS data Aug 21.csv")
aug_pfas <- aug_pfas[(aug_pfas$Sample_Class %in% c("Study_Sample")),]

for(i in 1:length(aug_pfas$SAMPLEID)){
  if(length(strsplit(aug_pfas$SAMPLEID,"AB")[i][[1]]) == 2){
    aug_pfas$SAMPLEID[i] = (strsplit(aug_pfas$SAMPLEID,"AB")[i][[1]])[2]
  }
  else if (length(strsplit(aug_pfas$SAMPLEID,"X")[i][[1]]) == 2){
    aug_pfas$SAMPLEID[i] = (strsplit(aug_pfas$SAMPLEID,"X")[i][[1]])[2]
  }
  
}

data_all__pos <- merge(data_all__pos,aug_pfas, by = "SAMPLEID")


df <- data_all__pos[,colnames(data_all__pos) %in% c("PRS_score","age_at_enrollment","self_reported_race",
                                                    "smoking_at_enrollment","gender","bmi_at_enrollment","glucose_enrl","age_at_t2d",
                                                    "PFDA_Aug21","PFHpA_Aug21","PFHxS_Aug21","PFHpS_Aug21","PFNA_Aug21","PFOA_Aug21",
                                                    "PFOS_Aug21","status","selection_prop","month_yr_enrl")]
init = mice(df, maxit=0) 
meth = init$method
predM = init$predictorMatrix
# predM[, c("redcap_survey_identifier")]=0

set.seed(1234)
imputed = mice(df, method=meth, predictorMatrix=predM, m=10)
imputed <- complete(imputed, action = 10)
data_all__pos$glucose_enrl_imputed <- imputed$glucose_enrl
data_all__pos$bmi_at_enrollment_imputed <- imputed$bmi_at_enrollment
dim(data_all__pos)

data_pred <- data_all__pos[,colnames(data_all__pos) %in% c("PRS_score","age_at_enrollment","self_reported_race",
                                                           "smoking_at_enrollment","gender","bmi_at_enrollment_imputed","glucose_enrl_imputed","age_at_t2d",
                                                           "PFDA_Aug21","PFHpA_Aug21","PFHxS_Aug21","PFHpS_Aug21","PFNA_Aug21","PFOA_Aug21",
                                                           "PFOS_Aug21","status","selection_prop","month_yr_enrl","sw",metabolites_refmet_pos)]
data_pred$status <- ifelse(data_pred$status == "case", 1, 0)
data_pred$date_enrl <- rep(NA_real_, nrow(data_pred))
for(i in 1:nrow(data_pred)){
  
  x <- anytime::anydate(paste((strsplit(data_pred$month_yr_enrl[i],"-")[[1]][2]), " 1,", 2000 + as.numeric(strsplit(data_pred$month_yr_enrl[i],"-")[[1]][1])))
  mydates <- as.Date(c("2011-01-01"))
  data_pred$date_enrl[i] <- as.numeric((x - mydates[1])/365 )
  
}



##############################

d_logis_dummy <- as.data.frame(dummify(data_pred[,!(colnames(data_pred) %in% c("month_yr_enrl"))]))

new_quantile <- function(x, cuts ){
  
  y <- x[x!= 0 & !is.na(x)]
  qi <- unique(quantile(y, probs = seq(0, 1, by = 1/cuts), na.rm = TRUE))
  
  if(length(qi) == 1){ 
    qi = c(-Inf, qi)
  } else{ 
    qi[1] <- -Inf
    qi[length(qi)] <- Inf
  }
  
  x[which(x!= 0 & !is.na(x))] = cut(x[x!= 0 & !is.na(x)], breaks = qi, labels = FALSE, include.lowest = TRUE)
  
  return(x)
  
}




# Change PFAS into Tertiles


qqt <- as.data.frame(apply(d_logis_dummy[,c("PFDA_Aug21","PFHpA_Aug21","PFHxS_Aug21","PFHpS_Aug21","PFNA_Aug21","PFOA_Aug21",
                                            "PFOS_Aug21")], 2, function(x) new_quantile(x, cuts = 3)))
colnames(qqt) <- paste0(c("PFDA_Aug21","PFHpA_Aug21","PFHxS_Aug21","PFHpS_Aug21","PFNA_Aug21","PFOA_Aug21",
                          "PFOS_Aug21"), "_q")
d_logis_dummy <- cbind(d_logis_dummy, qqt)


# Change Metabolites into Tertiles

qqt <- as.data.frame(apply(d_logis_dummy[,colnames(d_logis_dummy) %in% paste0("pos",seq(1,dim(refmet_pos)[1]))], 2, function(x) new_quantile(x, cuts = 3)))
colnames(qqt) <- paste0(paste0("pos",seq(1,dim(refmet_pos)[1])), "_q")
d_logis_dummy <- cbind(d_logis_dummy, qqt)

d_logis_dummy$c_date_enrl <- ifelse(d_logis_dummy$date_enrl > 0, 1,0)

###################################################################################################################################################

# T2D vs. PFAS

## WQS

name_data_wqs <- c("PFDA_Aug21_q","PFHpA_Aug21_q","PFHxS_Aug21_q","PFHpS_Aug21_q","PFNA_Aug21_q","PFOA_Aug21_q",
                   "PFOS_Aug21_q")

# name_data_wqs <- c("PFDA_Aug21","PFHpA_Aug21","PFHxS_Aug21","PFHpS_Aug21","PFNA_Aug21","PFOA_Aug21",
#                    "PFOS_Aug21")


pos_pfas_t2d_wqs <- gwqs(status ~ wqs + self_reported_race.African.American
                         + self_reported_race.European.American + age_at_enrollment
                         + smoking_at_enrollment.No + gender.Female
                         + bmi_at_enrollment_imputed + c_date_enrl,
                         mix_name = name_data_wqs, data = d_logis_dummy, q = NULL, signal = "t2",
                         validation = 0, b = 1000, rs = T, seed= 123123123,
                         b1_pos = T,  family = "binomial",plan_strategy = "multicore")

summary(pos_pfas_t2d_wqs)

d_logis_dummy$wqs_pfas <- pos_pfas_t2d_wqs$wqs

# hist(data_pred$date_enrl[data_pred$status==1])
# hist(data_pred$date_enrl[data_pred$status==0])
# hist(data_pred$date_enrl)


#                                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                          -0.6154517  0.9237141  -0.666  0.50523    
# wqs                                   0.2704159  0.1803414   1.499  0.13375    
# self_reported_race.African.American  -0.0241991  0.2864808  -0.084  0.93268    
# self_reported_race.European.American  0.4592939  0.3090220   1.486  0.13720    
# age_at_enrollment                    -0.0250714  0.0105841  -2.369  0.01785 *  
# smoking_at_enrollment.No             -0.7964932  0.3126820  -2.547  0.01086 *  
# gender.Female                         0.0007211  0.2481268   0.003  0.99768    
# bmi_at_enrollment_imputed             0.0812953  0.0187186   4.343 1.41e-05 ***
# c_date_enrl                          -0.6782607  0.2344903  -2.892  0.00382 ** 

# write.csv(pos_pfas_t2d_wqs$wqs, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/pfas_index_t2d_positive_mode.csv",
#           row.names = F)


#                    mix_name  mean_weight
# PFOS_Aug21_q   PFOS_Aug21_q 6.204153e-01
# PFHxS_Aug21_q PFHxS_Aug21_q 2.182434e-01
# PFNA_Aug21_q   PFNA_Aug21_q 5.670706e-02
# PFOA_Aug21_q   PFOA_Aug21_q 5.367972e-02
# PFHpA_Aug21_q PFHpA_Aug21_q 4.955160e-02
# PFDA_Aug21_q   PFDA_Aug21_q 1.360001e-03
# PFHpS_Aug21_q PFHpS_Aug21_q 4.296387e-05



start.time <- Sys.time()

pos_pfas_t2d_wqs <- gwqsrh(status ~ wqs + self_reported_race.African.American
                           + self_reported_race.European.American + age_at_enrollment
                           + smoking_at_enrollment.No + gender.Female
                           + bmi_at_enrollment_imputed,
                           mix_name = name_data_wqs, data = d_logis_dummy, q = NULL, signal = "t2",
                           validation = 0.25, b = 500, rs = F, rh = 20,
                           b1_pos = T,  family = "binomial",plan_strategy = "multicore")

summary(pos_pfas_t2d_wqs)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

pos_pfas_t2d_wqs$final_weights

#                    mix_name    Estimate        2.5 %      97.5%
# PFOS_Aug21_q   PFOS_Aug21_q 0.439852378 2.297773e-01 0.70862230
# PFHpA_Aug21_q PFHpA_Aug21_q 0.207809587 4.512701e-02 0.38486808
# PFOA_Aug21_q   PFOA_Aug21_q 0.128671308 4.358734e-02 0.35455778
# PFHxS_Aug21_q PFHxS_Aug21_q 0.099206354 1.945692e-02 0.22839569
# PFNA_Aug21_q   PFNA_Aug21_q 0.067676970 8.573174e-03 0.12640295
# PFDA_Aug21_q   PFDA_Aug21_q 0.047951243 5.560138e-03 0.11362042
# PFHpS_Aug21_q PFHpS_Aug21_q 0.008832161 1.210982e-05 0.02714881


# PFAS WQS Index

wqs_pfas <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/pfas_index_t2d_positive_mode.csv")
wqs_pfas <- wqs_pfas$x
d_logis_dummy$PFAS_index <- wqs_pfas

#######################################################################################################

# Metabolites vs. PFAS index


scaled_HRM <- as.data.frame((d_logis_dummy[,(colnames(d_logis_dummy) %in% metabolites_refmet_pos)]))
data_all_HRM_refmet_lm <- d_logis_dummy[,(colnames(d_logis_dummy) %in% c(colnames(scaled_HRM), "PFDA_Aug21_q","PFHpA_Aug21_q","PFHxS_Aug21_q","PFHpS_Aug21_q","PFNA_Aug21_q","PFOA_Aug21_q",
                                                                         "PFOS_Aug21_q"))]

data_all_HRM_refmet_lm$ID <- as.character(seq(1,nrow(data_all_HRM_refmet_lm)))
data_all_HRM_refmet_lm$PFAS_index <- d_logis_dummy$PFAS_index
data_all_HRM_refmet_lm$self_reported_race.African.American <- d_logis_dummy$self_reported_race.African.American
data_all_HRM_refmet_lm$smoking_at_enrollment.No <- d_logis_dummy$smoking_at_enrollment.No
data_all_HRM_refmet_lm$self_reported_race.European.American <- d_logis_dummy$self_reported_race.European.American
data_all_HRM_refmet_lm$age_at_enrollment <- d_logis_dummy$age_at_enrollment
data_all_HRM_refmet_lm$smoking_at_enrollment.No <- d_logis_dummy$smoking_at_enrollment.No
data_all_HRM_refmet_lm$gender.Female <- d_logis_dummy$gender.Female
data_all_HRM_refmet_lm$bmi_at_enrollment_imputed <- d_logis_dummy$bmi_at_enrollment_imputed
data_all_HRM_refmet_lm$status <- d_logis_dummy$status
data_all_HRM_refmet_lm$sw <- d_logis_dummy$sw
data_all_HRM_refmet_lm$c_date_enrl <- d_logis_dummy$c_date_enrl

name_data_wqs <- c("PFDA_Aug21_q","PFHpA_Aug21_q","PFHxS_Aug21_q","PFHpS_Aug21_q","PFNA_Aug21_q","PFOA_Aug21_q",
                   "PFOS_Aug21_q")

model_bwqs_gaussian_lasso <- "data {

int<lower=0> N;              // number of individual
int<lower=0> C1;             // number of element in the mix
int<lower=0> K;              // number of covariates
matrix[N,C1] XC1;            // matrix of first mix
matrix[N,K] KV;	             // matrix of covariates
vector[C1] DalpC1;           // vector of the Dirichlet coefficients for first mix
vector[N] sw;                // IPW weights
real y[N];                   // outcome gaussian variable
}

parameters {

real <lower=0> sigma;
real mu;                              // intercepts
real beta;                            // coeffs by group
vector[K] delta;                      // covariates coefficients
real<lower=0> lambda_squared;         // penalization factor
simplex[C1] WC1;                      // weights of first mix

}
transformed parameters {

vector[N] Xb;
Xb = mu + (XC1*WC1)*beta  + KV*delta;
}
model {

mu ~ normal(0, 10);
sigma ~ inv_gamma(0.01,0.01);
lambda_squared ~ gamma(2,0.5);
beta ~ normal(0,lambda_squared);
for(j in 1:K) delta[j] ~ normal(0,K);
WC1 ~ dirichlet(DalpC1);
for(n in 1:N){
  target +=  normal_lpdf(y[n]| Xb[n], sigma) * sw[n];
}
}

"
m_lasso_data_challenge <- rstan::stan_model(model_code =  model_bwqs_gaussian_lasso)

data <- data_all_HRM_refmet_lm
bwqs_pos_pfas_met <- data.frame(mean = NA_real_, se_mean = NA_real_, sd = NA_real_,
                                lower = NA_real_, upper = NA_real_, n_eff = NA_real_,
                                Rhat = NA_real_)

start.time <- Sys.time()

for(i in 1:372){
  
  y_name  <- colnames(data_all_HRM_refmet_lm)[i]
  formula = as.formula( ~ self_reported_race.African.American
                        + self_reported_race.European.American + age_at_enrollment
                        + smoking_at_enrollment.No + gender.Female
                        + bmi_at_enrollment_imputed + c_date_enrl)
  
  KV_name <- all.vars(formula)
  mix_name_1 <- c("PFDA_Aug21_q","PFHpA_Aug21_q","PFHxS_Aug21_q","PFHpS_Aug21_q","PFNA_Aug21_q","PFOA_Aug21_q",
                  "PFOS_Aug21_q")
  
  X1 = data_all_HRM_refmet_lm[,c("PFDA_Aug21_q","PFHpA_Aug21_q","PFHxS_Aug21_q","PFHpS_Aug21_q","PFNA_Aug21_q","PFOA_Aug21_q",
                                 "PFOS_Aug21_q")]
  
  data_reg <- list(
    
    N   = nrow(data),
    C1  = length(mix_name_1),
    XC1 = cbind(X1),
    DalpC1 = rep(1, length(mix_name_1)),
    KV = data[,KV_name],
    K   = length(KV_name),
    sw = as.vector(data[,"sw"]),
    y = as.vector(data[,y_name])
  )
  
  
  
  fit_lasso <- rstan::sampling(m_lasso_data_challenge,
                               data = data_reg,
                               chains = 1,
                               iter = 1e3,
                               thin = 1,
                               refresh = 0, verbose = T,
                               control=list(max_treedepth = 20,
                                            adapt_delta = 0.999999999999999))
  
  
  
  sum_fit_lasso <- (summary(fit_lasso,
                            probs = c(0.025, 0.975))$summary)
  
  
  bwqs_pos_pfas_met <-  rbind(bwqs_pos_pfas_met,as.numeric(sum_fit_lasso[3,])) 
  
  write.csv(bwqs_pos_pfas_met, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/posmode_meta_vs_pfas_bwqs.csv",
            row.names = F)
  
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


# bwqs_pos_pfas_met <- bwqs_pos_pfas_met[-1,]
# write.csv(bwqs_pos_pfas_met, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/posmode_meta_vs_pfas_bwqs.csv",
#           row.names = F)



bwqs_pos_pfas_met <-   read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/posmode_meta_vs_pfas_bwqs.csv")
bwqs_pos_pfas_met$refmet_name <- refmet_pos$refmet_name
bwqs_pos_pfas_met$super_class <- refmet_pos$super_class
bwqs_pos_pfas_met$sub_class <- refmet_pos$sub_class
bwqs_pos_pfas_met$time <- refmet_pos$time
bwqs_pos_pfas_met$mz <- refmet_pos$mz
bwqs_pos_pfas_met$mode <- rep("positive",nrow(bwqs_pos_pfas_met))

se = (bwqs_pos_pfas_met$upper - bwqs_pos_pfas_met$mean)/1.96
pval <- pnorm(abs(bwqs_pos_pfas_met$mean/se), lower.tail = F) 
bwqs_pos_pfas_met$p.value <- pval
q <- qvalue::qvalue(pval, lambda = 0)
bwqs_pos_pfas_met$q.value <-  q$qvalues

write.csv(bwqs_pos_pfas_met, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/posmode_meta_vs_pfas_bwqs.csv",
          row.names = F)
# 

# sort(bwqs_pos_pfas_met$q.value)[38]    
# 0.1534536

write.table(bwqs_pos_pfas_met[,c("mz","time","p.value","mean")], sep = "\t",
            "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/bwqs_exwas_pfas_meta_mummichog.txt",row.names = F, col.names = F)



#### T2D status vs. metabolites

require(foreign)
require(sandwich)

data_pos_mode_status_met <- d_logis_dummy[,(colnames(d_logis_dummy) %in% paste0(paste0("pos",seq(1,dim(refmet_pos)[1])), "_q"))]
data_pos_mode_status_met$ID <- as.character(seq(1,nrow(data_pos_mode_status_met)))
data_pos_mode_status_met$self_reported_race.African.American <- d_logis_dummy$self_reported_race.African.American
data_pos_mode_status_met$smoking_at_enrollment.No <- d_logis_dummy$smoking_at_enrollment.No
data_pos_mode_status_met$self_reported_race.European.American <- d_logis_dummy$self_reported_race.European.American
data_pos_mode_status_met$age_at_enrollment <- d_logis_dummy$age_at_enrollment
data_pos_mode_status_met$smoking_at_enrollment.No <- d_logis_dummy$smoking_at_enrollment.No
data_pos_mode_status_met$gender.Female <- d_logis_dummy$gender.Female
data_pos_mode_status_met$bmi_at_enrollment_imputed <- d_logis_dummy$bmi_at_enrollment_imputed
data_pos_mode_status_met$status <- d_logis_dummy$status
data_pos_mode_status_met$c_date_enrl <- d_logis_dummy$c_date_enrl


d_lm_status <- data.frame(met = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)
for(i in 1:372){
  s_lm <- (glm(status ~ data_pos_mode_status_met[,i] + self_reported_race.African.American
               + self_reported_race.European.American + age_at_enrollment + gender.Female
               + bmi_at_enrollment_imputed + c_date_enrl , data = data_pos_mode_status_met, family=binomial))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  d_lm_status <- rbind(d_lm_status,c(colnames(data_pos_mode_status_met)[i], as.numeric(r.est[2,c(1,2,3,4)])))
}

d_lm_status <- d_lm_status[-1,]
d_lm_status$z.value <- as.numeric(d_lm_status$z.value)
d_lm_status$p.value <- as.numeric(d_lm_status$p.value)


d_lm_status$refmet_name <- refmet_pos$refmet_name
d_lm_status$super_class <- refmet_pos$super_class
d_lm_status$sub_class <- refmet_pos$sub_class
d_lm_status$time <- refmet_pos$time
d_lm_status$mz <- refmet_pos$mz

q <-qvalue::qvalue(as.numeric(d_lm_status$p.value), lambda=0)
d_lm_status$q.value <-  q$qvalues

write.table(d_lm_status[,c("mz","time","p.value","z.value")], sep = "\t",
            "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/t2d_meta_mummichog.txt",
            row.names = F, col.names = F)


write.csv(d_lm_status,
          "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/mixedmodel_exwas_t2d_meta_all.csv",
          row.names = F)


##############################################################################################################################

# Common Metabolites

# Positive mode

bwqs_pos_pfas_met <-   read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/posmode_meta_vs_pfas_bwqs.csv")
d_lm_status <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/mixedmodel_exwas_t2d_meta_all.csv")


bb <- bwqs_pos_pfas_met[bwqs_pos_pfas_met$p.value < 0.05,c("refmet_name")]
tt <- d_lm_status[d_lm_status$p.value < 0.05, c("refmet_name")]
cc <- intersect(bb,tt) 

ft <-bwqs_pos_pfas_met[bwqs_pos_pfas_met$refmet_name %in% cc,c("refmet_name","mean","p.value","q.value")]
ft <- ft[order(ft$refmet_name),]

dt <- d_lm_status[d_lm_status$refmet_name %in% cc,c("refmet_name","Value","p.value","q.value")]
dt <- dt[order(dt$refmet_name),]

common_tab <- data.frame(refmet_name = dt[,"refmet_name"], mean.pfas.met = ft$mean, pval.pfas.met = ft$p.value, 
                         qval.pfas.met = ft$q.value, mean.met.t2d = dt$Value, pval.met.t2d = dt$p.value,
                         qval.met.t2d = dt$q.value)



# Negative mode

bwqs_neg_pfas_met <-   read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/neg/negmode_meta_vs_pfas_bwqs.csv")
d_lm_status <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/neg/mixedmodel_exwas_t2d_meta_all.csv")

bb <- bwqs_neg_pfas_met[bwqs_neg_pfas_met$p.value < 0.05,c("refmet_name")]
tt <- d_lm_status[d_lm_status$p.value < 0.05, c("refmet_name")]
cc <- intersect(bb,tt) 

ft <-bwqs_neg_pfas_met[bwqs_neg_pfas_met$refmet_name %in% cc,c("refmet_name","mean","p.value","q.value")]
ft <- ft[order(ft$refmet_name),]

dt <- d_lm_status[d_lm_status$refmet_name %in% cc,c("refmet_name","Value","p.value","q.value")]
dt <- dt[order(dt$refmet_name),]

common_tab <- data.frame(refmet_name = dt[,"refmet_name"], mean.pfas.met = ft$mean, pval.pfas.met = ft$p.value, 
                         qval.pfas.met = ft$q.value, mean.met.t2d = dt$Value, pval.met.t2d = dt$p.value,
                         qval.met.t2d = dt$q.value)

common_tab

#               refmet_name mean.pfas.met pval.pfas.met qval.pfas.met mean.met.t2d pval.met.t2d qval.met.t2d
# 1    5-Hydroxy-tryptophan     0.3394450  7.232089e-04  1.062669e-02   -0.3158309  0.031736747    0.7303621
# 2          Glucoheptulose     0.2048479  1.256501e-02  8.180441e-02   -0.3400293  0.025262938    0.7303621
# 3 Sulfolithocholylglycine     0.5170315  4.024093e-09  1.142842e-06   -0.4648498  0.001309185    0.3718085

#################################################################################################################################

# Pathway analysis
# Positive

em_comp_pfas_meta <- read.delim("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/pos/result_pfas_meta/tables/ListOfEmpiricalCompounds.tsv")
colnames(em_comp_pfas_meta)[4:5] <- c("compounds_mummichog","compound_names_mummichog")


em_comp_meta_t2d <- read.delim("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\BIOME-PFAS pilot\\T2D\\Paper T2D\\pos\\result_meta_t2d\\tables\\ListOfEmpiricalCompounds.tsv")
colnames(em_comp_meta_t2d)[4:5] <- c("compounds_mummichog","compound_names_mummichog")


em_pathways_pfas_meta  <- read.delim("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\BIOME-PFAS pilot\\T2D\\Paper T2D\\pos\\result_pfas_meta\\tables\\mcg_pathwayanalysis_.tsv")
colnames(em_pathways_pfas_meta)[5:7] <- c("overlap_EID","overlap_features_ID","overlap_features_names")

em_pathways_meta_t2d  <- read.delim("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\BIOME-PFAS pilot\\T2D\\Paper T2D\\pos\\result_meta_t2d\\tables\\mcg_pathwayanalysis_.tsv")
colnames(em_pathways_meta_t2d)[5:7] <- c("overlap_EID","overlap_features_ID","overlap_features_names")

cutoff <- 0.05
common_paths_pos <- unique(c(em_pathways_pfas_meta$pathway[em_pathways_pfas_meta$p.value < cutoff], 
                             em_pathways_meta_t2d$pathway[em_pathways_meta_t2d$p.value < cutoff]))
common_paths_pos

# [1] "Tyrosine metabolism"                                       "Glutamate metabolism"                                     
# [3] "Ascorbate (Vitamin C) and Aldarate Metabolism"             "Bile acid biosynthesis"                                   
# [5] "Vitamin B3 (nicotinate and nicotinamide) metabolism"       "Alanine and Aspartate Metabolism"                         
# [7] "Putative anti-Inflammatory metabolites formation from EPA" "Arginine and Proline Metabolism"                          
# [9] "Vitamin B9 (folate) metabolism"                            "Beta-Alanine metabolism"                                  
# [11] "Glycine, serine, alanine and threonine metabolism"         "Butanoate metabolism"                                     
# [13] "Leukotriene metabolism"                                    "Nitrogen metabolism"  



# Negative

em_comp_pfas_meta <- read.delim("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/neg/result_pfas_meta/tables/ListOfEmpiricalCompounds.tsv")
colnames(em_comp_pfas_meta)[4:5] <- c("compounds_mummichog","compound_names_mummichog")


em_comp_meta_t2d <- read.delim("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\BIOME-PFAS pilot\\T2D\\Paper T2D\\neg\\result_meta_t2d\\tables\\ListOfEmpiricalCompounds.tsv")
colnames(em_comp_meta_t2d)[4:5] <- c("compounds_mummichog","compound_names_mummichog")


em_pathways_pfas_meta  <- read.delim("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\BIOME-PFAS pilot\\T2D\\Paper T2D\\neg\\result_pfas_meta\\tables\\mcg_pathwayanalysis_.tsv")
colnames(em_pathways_pfas_meta)[5:7] <- c("overlap_EID","overlap_features_ID","overlap_features_names")

em_pathways_meta_t2d  <- read.delim("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\BIOME-PFAS pilot\\T2D\\Paper T2D\\neg\\result_meta_t2d\\tables\\mcg_pathwayanalysis_.tsv")
colnames(em_pathways_meta_t2d)[5:7] <- c("overlap_EID","overlap_features_ID","overlap_features_names")


cutoff <- 0.05

common_paths_neg <- unique(c(em_pathways_pfas_meta$pathway[em_pathways_pfas_meta$p.value < cutoff], 
                             em_pathways_meta_t2d$pathway[em_pathways_meta_t2d$p.value < cutoff]))

common_paths_neg

intersect(common_paths_pos,common_paths_neg)



