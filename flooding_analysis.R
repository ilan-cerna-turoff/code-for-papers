# File: Imputation and analysis of flooding and mental health in Peru
# Date: 08/15/2023
# Author: Ilan Cerna-Turoff

#------------------------------------------------------------------------------
####    NOTES   ####
#------------------------------------------------------------------------------
# This script sets up analysis, imputes covariates for missingness, and conducts
# sensitivity analyses. The associated data can be requested from the UK data service at
# https://beta.ukdataservice.ac.uk/datacatalogue/series/series?id=2000060#!/access-data.
# Data preparation can be found at https://github.com/ilan-cerna-turoff/young_lives_data_preparation.
# For additional information, contact the corresponding author directly.

# Citation: Cerna-Turoff I, Kang H, Keyes K. El Niño-driven flooding 
# and mental health symptomology among adolescents and young adults in Peru. 
# [INSERT FINAL JOURNAL LINK HERE].

#------------------------------------------------------------------------------
#### 1. LOAD DATA  ####
#------------------------------------------------------------------------------
pacman::p_load("tidyverse","naniar","here","mice","MatchThem","cobalt","survey",
               "gtsummary","webshot2","kableExtra","EValue")

df1 <- read_rds(here("data","final","merged_final","pe.anydis.v2.rds")) 

#------------------------------------------------------------------------------
#### 2. CLEAN INDIVIDUAL AND COMMUNITY DATA ####
#------------------------------------------------------------------------------
# 2a. Removes missing label
df1 <- df1 %>%
  replace_with_na_all(condition = ~.x == "missing") 

# 2b. Remove clusters with sample size 1
clusters <- df1 %>%
  group_by(commid) %>%
  count() %>%
  ungroup() %>%
  filter(n > 1)

df1 <- left_join(clusters,df1, by = "commid")

rm(clusters)

# 2c. Create dataframe with community location of respondents
commid <- df1 %>%
  select(c(childid,commid)) %>%
  distinct()

# 2d. Remove other types of natural hazard exposures
noexposure <- df1 %>%
  filter(flood != "other") %>%
  mutate(flood = ifelse(flood == "none", 0, 1))

# 2e. Clean covariates 
covariates <- df1 %>%   
  dplyr::select(childid,chsex,chethnic,age,smoking,wellbeing,hhsize,ownhouse,
                jobloss,caredu,divorce,sochel) %>%
  distinct() %>%
  mutate(across(where(is.factor), ~ fct_drop(.)), #get rid of missing factor levels
         hhsize = ifelse(hhsize == "living outside the household", NA, hhsize),
         across(c(hhsize,wellbeing), as.character),
         across(c(hhsize,wellbeing), as.numeric),
         chethnic = as.factor(chethnic),
         chethnic = fct_collapse(chethnic, other = c("asiatic","native of the amazon",
                                                     "negro"))) 

noexposure <- noexposure %>%
  dplyr::select(c(childid,flood,finalscore)) %>% #add outcome here
  distinct() 

noexposure1 <- left_join(noexposure, covariates, by = "childid") #join covariates with exposure

rm(covariates,noexposure)

#------------------------------------------------------------------------------
#### 3. IMPUTATION - FULL SAMPLE ####
#------------------------------------------------------------------------------
# 3a. Explore missingness of covariates 
p_missing <- unlist(lapply(noexposure1, function(x) sum(is.na(x))))/nrow(noexposure1)
sort(p_missing[p_missing > 0], decreasing = TRUE) 

# NOTE: 
# Some guidance recommends taking out missing variables if over 25%, but 
# simulation study has shown that up to 90% missingness can be fine if missing at
# random (Madley-Dowd et al., 2019). This analysis is run without a percentage 
# threshold for inclusion.

# 3b. Retain childid but exclude from imputation
ini <- mice(noexposure1,maxit=0) # imputation without iterations

pred1 <- ini$predictorMatrix # this is your predictor matrix
pred1[,'childid'] <- 0 # set all id column values to zero to exclude it as a predictor
pred1[,'finalscore'] <- 0 # set all outcome column values to zero to exclude it as a predictor
pred1[,'flood'] <- 0 # set all treatment column values to zero to exclude it as a predictor

# 3c. Run imputation
imp <- mice(noexposure1, m = 5, seed = 9999, pred = pred1)

# 3d. Examine imputed results
# head(imp$imp$ownhouse)
# tail(imp$imp$ownhouse)
# colSums(is.na(imp$imp$ownhouse)) 
# test <- complete(imp, 1)
# densityplot(imp)
# head(imp$loggedEvents, 1) 

#------------------------------------------------------------------------------
#### 4. MATCHING - FULL SAMPLE ####
#------------------------------------------------------------------------------
# 4a. Mahalanobis distance with propensity score matching on imputed data 
matched <- matchthem(flood ~ chsex + chethnic + age + 
                       smoking + wellbeing + hhsize + ownhouse + jobloss + 
                       caredu + divorce + sochel,
                     datasets = imp,
                     approach = 'within',
                     method = 'full',
                     distance = 'glm',
                     link = binomial("logit"),
                     mahvars = ~ chsex + chethnic + age + 
                       smoking + wellbeing + hhsize + ownhouse + jobloss + 
                       caredu + divorce + sochel,
                     min.controls = 1,
                     max.controls = 3,
                     omit.fraction = .6)

# 4b. Examine balance
bal.tab(matched, un = T)

new.names <- c(chsex_male = "Sex",
               age = "Age",
               chethnic_other = "Race: Other",
               chethnic_mestizo = "Race: Mestizo",
               chethnic_white = "Race: White",
               smoking_yes = "Tobacco Use",
               wellbeing = "Wellbeing Score",
               hhsize = "Household Size",
               ownhouse_yes = "Household Ownership",
               jobloss_yes = "Recent Caregiver Jobloss",
               divorce_yes = "Recent Caregiver Divorce",
               caredu_none = "Caregiver Education: None",
               caredu_basic = "Caregiver Education: Basic",
               `caredu_above basic` = "Caregiver Education: Above Basic",
               sochel_yes = "Community Social Worker Access")

love.plot(matched, stars = "raw", drop.distance = TRUE, var.names = new.names, thresholds = c(m = .1))
#ggsave("loveplot1.png", width = 6, height = 6, dpi = 300)

#------------------------------------------------------------------------------
#### 4. ANALYSIS - FULL SAMPLE ####
#------------------------------------------------------------------------------
# 4a. Run the models on the imputed datasets
matched.models <- with(matched,
                       svyglm(finalscore ~ flood, family = quasipoisson(link = "log")))

# Note:
# svydesign not needed with MatchThem. Weights for clustering automatically constructed. 

# 4b. Pool across imputed datasets
summary(pool(matched.models), conf.int = TRUE)

# 4c. Exponentiate
inverse_logit = function(x){ 
  exp(x)
}

result.tibble <-  tibble(summary(pool(matched.models), conf.int = TRUE)) %>%
  filter(row_number()==2)

expected.result = inverse_logit(result.tibble$estimate)
low.result = inverse_logit(result.tibble$`2.5 %`)
high.result = inverse_logit(result.tibble$`97.5 %`)

# 4d. Sensitivity analysis
summary(evalue(RR(1.02), true = 2))  #if true effect was 2
 
#------------------------------------------------------------------------------
#### 5. IMPUTATION - STRATIFIED SAMPLE BY GENDER ####
#------------------------------------------------------------------------------
# 5a. Stratify
female <- noexposure1 %>%
  filter(chsex == "female") %>%
  select(-chsex) %>%
  distinct()

male <- noexposure1 %>%
  filter(chsex == "male") %>%
  select(-chsex) %>%
  distinct()

# 5b. Retain childid but exclude from imputation
# Female
ini.female <- mice(female,maxit=0) 

pred1.female <- ini.female$predictorMatrix 
pred1.female[,'childid'] <- 0 
pred1.female[,'finalscore'] <- 0 
pred1.female[,'flood'] <- 0 

# Male
ini.male <- mice(male,maxit=0) 

pred1.male <- ini.male$predictorMatrix
pred1.male[,'childid'] <- 0 
pred1.male[,'finalscore'] <- 0 
pred1.male[,'flood'] <- 0 

# 5c. Run imputations
imp1 <- mice(female, m = 5, seed = 9999, pred = pred1.female)
imp2 <- mice(male, m = 5, seed = 9999, pred = pred1.male)

#------------------------------------------------------------------------------
#### 6. MATCHING - STRATIFIED SAMPLE BY GENDER ####
#------------------------------------------------------------------------------
# 6a. Mahalanobis distance with propensity score matching on imputed data 
matched1 <- matchthem(flood ~ chethnic + age + 
                        smoking + wellbeing + hhsize + ownhouse + jobloss + 
                        caredu + divorce + sochel,
                      datasets = imp1,
                      approach = 'within',
                      method = 'full',
                      distance = 'glm',
                      link = binomial("logit"),
                      mahvars = ~ chethnic + age + smoking + wellbeing + hhsize + 
                        ownhouse + jobloss + caredu + divorce + sochel,
                      min.controls = 1,
                      max.controls = 3,
                      omit.fraction = .6)

matched2 <- matchthem(flood ~ chethnic + age + 
                        smoking + wellbeing + hhsize + ownhouse + jobloss + 
                        caredu + divorce + sochel,
                      datasets = imp2,
                      approach = 'within',
                      method = 'full',
                      distance = 'glm',
                      link = binomial("logit"),
                      mahvars = ~ chethnic + age + smoking + wellbeing + hhsize + 
                        ownhouse + jobloss + caredu + divorce + sochel,
                      min.controls = 1,
                      max.controls = 3,
                      omit.fraction = .6)

# 6b. Examine balance
bal.tab(matched1, un = T)
bal.tab(matched2, un = T)

new.names <- c(age = "Age",
               chethnic_other = "Race: Other",
               chethnic_mestizo = "Race: Mestizo",
               chethnic_white = "Race: White",
               smoking_yes = "Tobacco Use",
               wellbeing = "Wellbeing Score",
               hhsize = "Household Size",
               ownhouse_yes = "Household Ownership",
               jobloss_yes = "Recent Caregiver Jobloss",
               divorce_yes = "Recent Caregiver Divorce",
               caredu_none = "Caregiver Education: None",
               caredu_basic = "Caregiver Education: Basic",
               `caredu_above basic` = "Caregiver Education: Above Basic",
               sochel_yes = "Community Social Worker Access")

love.plot(matched1, stars = "raw", drop.distance = TRUE, var.names = new.names, thresholds = c(m = .1))
#ggsave("loveplot.female.png", width = 6, height = 6, dpi = 300)

love.plot(matched2, stars = "raw", drop.distance = TRUE, var.names = new.names, thresholds = c(m = .1))
#ggsave("loveplot.male.png", width = 6, height = 6, dpi = 300)

#------------------------------------------------------------------------------
#### 7. ANALYSIS - STRATIFIED SAMPLE BY GENDER ####
#------------------------------------------------------------------------------
# 7a. Run the models on the imputed datasets
matched.models1 <- with(matched1,
                        svyglm(finalscore ~ flood, family = quasipoisson(link = "log")))

matched.models2 <- with(matched2,
                        svyglm(finalscore ~ flood, family = quasipoisson(link = "log")))

# 7b. Pool across imputed datasets
summary(pool(matched.models1), conf.int = TRUE)
summary(pool(matched.models2), conf.int = TRUE)

# 7c. Exponentiate
# Female
tibble1 <-  tibble(summary(pool(matched.models1), conf.int = TRUE)) %>%
  filter(row_number()==2)

expected1 = inverse_logit(tibble1$estimate)
low1 = inverse_logit(tibble1$`2.5 %`)
high1 = inverse_logit(tibble1$`97.5 %`)

# Male
tibble2 <-  tibble(summary(pool(matched.models2), conf.int = TRUE)) %>%
  filter(row_number()==2)

expected2 = inverse_logit(tibble2$estimate)
low2 = inverse_logit(tibble2$`2.5 %`)
high2 = inverse_logit(tibble2$`97.5 %`)

# 7d. Sensitivity analyses
#Female
summary(evalue(RR(0.97), true = 2)) #if true effect was 2

#Male
summary(evalue(RR(1.06), true = 2)) 