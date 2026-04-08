# =====================================================================================================================================================================
# Explore the potential mediating effect of prenatal TMAO on the association between contaminants measured in pregnancy and child neurodevelopmental outcomes
# Purpose: 
# - Linear mixed effect models between between As, Hg, PFOS during pregnancy and each neurodevelopmental outcome, 
# with or without interaction an term between the contaminant and tertiles of TMAO. 
# - Likelihood Ratio Test between the model with or wihtout interaction term
# Requirements: dplyr, lmerTest, broom, lmtest
# =====================================================================================================================================================================

library(dplyr)
library(lmerTest)
library(broom)
library(lmtest)

# Data files
ADHD_df <- read.csv("data/longitudinal_merged_ADHD_df.csv")
CBCL_df <- read.csv("data/longitudinal_merged_CBCL_df.csv")
dnumeros_df <- read.csv("data/longitudinal_merged_dnumeros_df.csv")
general_cognition_df <- read.csv("data/longitudinal_merged_general_cognition_df.csv")
hitrtse_df <- read.csv("data/longitudinal_merged_hitrtse_df.csv")
raven_df <- read.csv("data/longitudinal_merged_raven_df.csv")


compute_TMAO_tertiles <- function(df){
  df <- df %>% mutate(TMAO_tertiles_1st = ntile(TMAO_1st, 3), TMAO_tertiles_1st = factor(TMAO_tertiles_1st))
  df <- df %>% mutate(TMAO_tertiles_3rd = ntile(TMAO_3rd, 3), TMAO_tertiles_3rd = factor(TMAO_tertiles_3rd))
  return(df)
}


ADHD_df<- compute_TMAO_tertiles(ADHD_df)
CBCL_df <- compute_TMAO_tertiles(CBCL_df)
dnumeros_df <- compute_TMAO_tertiles(dnumeros_df)
general_cognition_df <- compute_TMAO_tertiles(general_cognition_df)
hitrtse_df <- compute_TMAO_tertiles(hitrtse_df)
raven_df <- compute_TMAO_tertiles(raven_df)


# Transform the CBCL outcomes + ADHD (z score of the square root):
CBCL_df <- CBCL_df %>% mutate(
  Gen_Tot_scaled = scale(sqrt(CBCL_df$Gen_Tot)),
  Gen_Ext_scaled = scale(sqrt(CBCL_df$Gen_Ext)),
  Gen_Int_scaled = scale(sqrt(CBCL_df$Gen_Int))
)

ADHD_df <-ADHD_df %>% mutate(ADHD_scaled = scale(sqrt(ADHD_df$ADHD)))


# CONTAMINANTS FROM 1st OR 3rd TRIMESTER
metabolites_1st <- c("As_1st", "PFOS_1st")
metabolites_3rd <- c("As_3rd", "Hg_3rd")


prepare_confounders <- function(df) {
  # obstetric complications: birth weight <2500gr or >4000gr or GA <37 weeks
  df$obs_complications<-factor(ifelse((df$peso<2500 | df$peso>4000 | df$sges <37),"yes","no"))
  
  # transform to factors
  df$paridad<-as.factor(df$paridad)
  df$estudios <- as.factor(df$estudios)
  df$season <- as.factor(df$season)
  df$cohort <- as.factor(df$cohort)
  df$id_inma <- as.factor(df$id_inma)
  return(df)
}


ADHD_df<- prepare_confounders(ADHD_df)
CBCL_df <- prepare_confounders(CBCL_df)
dnumeros_df <- prepare_confounders(dnumeros_df)
general_cognition_df <- prepare_confounders(general_cognition_df)
hitrtse_df <- prepare_confounders(hitrtse_df)
raven_df <- prepare_confounders(raven_df)



def_formula <- function(outcome, metabolite, trim) {
  
  TMAO_tertiles <- ifelse(trim==1, "TMAO_tertiles_1st", "TMAO_tertiles_3rd")
  
  formula_main <- paste0(outcome, " ~ log2(",metabolite,") +", TMAO_tertiles," +", paste(main_confounders, collapse = " + "), " + (1|id_inma)")
  formula_main_interaction <- paste0(outcome, " ~ log2(",metabolite,") +", TMAO_tertiles, " +", paste(main_confounders, collapse = " + "), " + (1|id_inma) + ", TMAO_tertiles,":log2(",metabolite,")")
  
  formula_FU1 <- paste0(outcome, " ~ log2(",metabolite,") +", TMAO_tertiles," +", paste(c(main_confounders, "obs_complications"), collapse = " + "), " + (1|id_inma)")
  formula_FU1_interaction <- paste0(outcome, " ~ log2(",metabolite,") +", TMAO_tertiles," +", paste(c(main_confounders, "obs_complications"), collapse = " + "), " + (1|id_inma) + ", TMAO_tertiles,":log2(",metabolite,")")
  
  formula_FU2_1st <- paste0(outcome, " ~ log2(",metabolite,") +", TMAO_tertiles," +", paste(c(main_confounders, "Pescados_gr_s12"), collapse = " + "), " + (1|id_inma)")
  formula_FU2_1st_interaction <- paste0(outcome, " ~ log2(",metabolite,") +", TMAO_tertiles," +", paste(c(main_confounders, "Pescados_gr_s12"), collapse = " + "), " + (1|id_inma) + ", TMAO_tertiles,":log2(",metabolite,")")
  
  
  formula_FU2_3rd <- paste0(outcome, "~ log2(",metabolite,") +", TMAO_tertiles," +", paste(c(main_confounders, "PescadoGr_gr_s32"), collapse = " + "), "+ (1|id_inma)")
  formula_FU2_3rd_interaction <- paste0(outcome, "~ log2(",metabolite,") +", TMAO_tertiles," +", paste(c(main_confounders, "PescadoGr_gr_s32"), collapse = " + "), "+ (1|id_inma) + ", TMAO_tertiles,":log2(",metabolite,")")

  return(list(formula_main = formula_main, formula_main_interaction = formula_main_interaction,
              formula_FU1 = formula_FU1, formula_FU1_interaction = formula_FU1_interaction, 
              formula_FU2_1st = formula_FU2_1st, formula_FU2_1st_interaction = formula_FU2_1st_interaction, 
              formula_FU2_3rd = formula_FU2_3rd, formula_FU2_3rd_interaction = formula_FU2_3rd_interaction))
}


results_MAIN_df <- data.frame()
results_FU1_df <- data.frame()
results_FU2_df <- data.frame()

all_tertiles <- c(1,2,3)

### FOR 1st SEMESTER METABOLITES  
for (metabolite in metabolites_1st) {
  
  print(metabolite)
  
  # If metabolite is Arsenic --> Remove cohort from the list of coufounders (only measured in INMA Sabadell, not Gipuzkoa)
  if (metabolite == "As_1st") {
    main_confounders <- c("estudios","season","edadm","paridad","age_neuro")
  }
  
  else {
    main_confounders <- c("estudios","season","edadm","paridad","cohort","age_neuro")
  }
  
  
  #######  DNUMEROS2 #####
  outcome <- "dnumeros2"
  data <- dnumeros_df
  data <- data %>% filter(TMAO_tertiles_1st %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=1)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_1st <- formulas$formula_FU2_1st
  formula_FU2_1st_interaction <- formulas$formula_FU2_1st_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_1st)
  print(formula_FU2_1st_interaction)

  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_1st, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_1st_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval

  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
#   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)

#   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
#   

  
  
#   #### DNUMEROS 3 ######
#   
  outcome <- "dnumeros3"
  data <- dnumeros_df
  #data <- data %>% filter(TMAO_tertiles_1st == t)
  data <- data %>% filter(TMAO_tertiles_1st %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=1)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_1st <- formulas$formula_FU2_1st
  formula_FU2_1st_interaction <- formulas$formula_FU2_1st_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_1st)
  print(formula_FU2_1st_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_1st, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_1st_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
#   
#   ####### HITRTSE #####
#   

  outcome <- "hitrtse"
  data <- hitrtse_df
  #data <- data %>% filter(TMAO_tertiles_1st == t)
  data <- data %>% filter(TMAO_tertiles_1st %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=1)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_1st <- formulas$formula_FU2_1st
  formula_FU2_1st_interaction <- formulas$formula_FU2_1st_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_1st)
  print(formula_FU2_1st_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_1st, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_1st_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)

#   
#   ####### GENERAL COGNITION - MOTOR #####
#   
#   
  outcome <- "motor"
  data <- general_cognition_df
  #data <- data %>% filter(TMAO_tertiles_1st == t)
  data <- data %>% filter(TMAO_tertiles_1st %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=1)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_1st <- formulas$formula_FU2_1st
  formula_FU2_1st_interaction <- formulas$formula_FU2_1st_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_1st)
  print(formula_FU2_1st_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_1st, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_1st_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
# 

   
#   ####### GENERAL COGNITION - MENTAL #####
#   
#   
  outcome <- "mental"
  data <- general_cognition_df
  #data <- data %>% filter(TMAO_tertiles_1st == t)
  data <- data %>% filter(TMAO_tertiles_1st %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=1)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_1st <- formulas$formula_FU2_1st
  formula_FU2_1st_interaction <- formulas$formula_FU2_1st_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_1st)
  print(formula_FU2_1st_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_1st, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_1st_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
#   
#   
#   ####### RAVEN PERCENTILE #####
#   
  outcome <- "percentile_correct_raven"
  data <- raven_df
  #data <- data %>% filter(TMAO_tertiles_1st == t)
  data <- data %>% filter(TMAO_tertiles_1st %in% all_tertiles)

  formulas <- def_formula(outcome, metabolite, trim=1)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_1st <- formulas$formula_FU2_1st
  formula_FU2_1st_interaction <- formulas$formula_FU2_1st_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_1st)
  print(formula_FU2_1st_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_1st, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_1st_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
#   
#   ####### ADHD #####
#   
  outcome <- "ADHD_scaled"
  data <- ADHD_df
  #data <- data %>% filter(TMAO_tertiles_1st == t)
  data <- data %>% filter(TMAO_tertiles_1st %in% all_tertiles)

  formulas <- def_formula(outcome, metabolite, trim=1)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_1st <- formulas$formula_FU2_1st
  formula_FU2_1st_interaction <- formulas$formula_FU2_1st_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_1st)
  print(formula_FU2_1st_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_1st, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_1st_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)

#   
#   ####### CBCL Gen Tot #####
#   
  outcome <- "Gen_Tot_scaled"
  data <- CBCL_df
  #data <- data %>% filter(TMAO_tertiles_1st == t)
  data <- data %>% filter(TMAO_tertiles_1st %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=1)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_1st <- formulas$formula_FU2_1st
  formula_FU2_1st_interaction <- formulas$formula_FU2_1st_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_1st)
  print(formula_FU2_1st_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_1st, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_1st_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)

#   
#   ####### CBCL Gen Ext #####
#   
#   
  outcome <- "Gen_Ext_scaled"
  data <- CBCL_df
  #data <- data %>% filter(TMAO_tertiles_1st == t)
  data <- data %>% filter(TMAO_tertiles_1st %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=1)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_1st <- formulas$formula_FU2_1st
  formula_FU2_1st_interaction <- formulas$formula_FU2_1st_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_1st)
  print(formula_FU2_1st_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_1st, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_1st_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  # print(paste0("FIRST PVAL IS ", lrt$'Pr(>Chisq)'[1]))
  # print(paste0("SECOND PVAL IS ", lrt$'Pr(>Chisq)'[2]))
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
#   
#   
#   
#   ####### CBCL Gen Int #####
#   
#   
  outcome <- "Gen_Int_scaled"
  data <- CBCL_df
  #data <- data %>% filter(TMAO_tertiles_1st == t)
  data <- data %>% filter(TMAO_tertiles_1st %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=1)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_1st <- formulas$formula_FU2_1st
  formula_FU2_1st_interaction <- formulas$formula_FU2_1st_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_1st)
  print(formula_FU2_1st_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_1st, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_1st_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "1st",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
  
}


# ###################################################################################


# ### FOR 3rd SEMESTER METABOLITES  
for (metabolite in metabolites_3rd) {

  print(metabolite)

  if (metabolite == "As_3rd") {
    main_confounders <- c("estudios","season","edadm","paridad","age_neuro")
  }

  else {
    main_confounders <- c("estudios","season","edadm","paridad","cohort","age_neuro")
  }


  #######  DNUMEROS2 #####
  outcome <- "dnumeros2"
  data <- dnumeros_df
  #data <- data %>% filter(TMAO_tertiles_3rd == t)
  data <- data %>% filter(TMAO_tertiles_3rd %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=3)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_3rd <- formulas$formula_FU2_3rd
  formula_FU2_3rd_interaction <- formulas$formula_FU2_3rd_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_3rd)
  print(formula_FU2_3rd_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_3rd, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_3rd_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
#   
#   
#   
#   #### DNUMEROS 3 ######
#   
  outcome <- "dnumeros3"
  data <- dnumeros_df
  #data <- data %>% filter(TMAO_tertiles_3rd == t)
  data <- data %>% filter(TMAO_tertiles_3rd %in% all_tertiles)

  formulas <- def_formula(outcome, metabolite, trim=3)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_3rd <- formulas$formula_FU2_3rd
  formula_FU2_3rd_interaction <- formulas$formula_FU2_3rd_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_3rd)
  print(formula_FU2_3rd_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_3rd, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_3rd_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)  
#   
#   ####### HITRTSE #####

  outcome <- "hitrtse"
  data <- hitrtse_df
  #data <- data %>% filter(TMAO_tertiles_3rd == t)
  data <- data %>% filter(TMAO_tertiles_3rd %in% all_tertiles)

  formulas <- def_formula(outcome, metabolite, trim=3)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_3rd <- formulas$formula_FU2_3rd
  formula_FU2_3rd_interaction <- formulas$formula_FU2_3rd_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_3rd)
  print(formula_FU2_3rd_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_3rd, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_3rd_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
   
#   ####### GENERAL COGNITION - MOTOR #####
#   

  outcome <- "motor"
  data <- general_cognition_df
  #data <- data %>% filter(TMAO_tertiles_3rd == t)
  data <- data %>% filter(TMAO_tertiles_3rd %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=3)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_3rd <- formulas$formula_FU2_3rd
  formula_FU2_3rd_interaction <- formulas$formula_FU2_3rd_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_3rd)
  print(formula_FU2_3rd_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_3rd, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_3rd_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)

#   ####### GENERAL COGNITION - MENTAL #####


  outcome <- "mental"
  data <- general_cognition_df
  #data <- data %>% filter(TMAO_tertiles_3rd == t)
  data <- data %>% filter(TMAO_tertiles_3rd %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=3)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_3rd <- formulas$formula_FU2_3rd
  formula_FU2_3rd_interaction <- formulas$formula_FU2_3rd_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_3rd)
  print(formula_FU2_3rd_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_3rd, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_3rd_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)

#   ####### RAVEN PERCENTILE #####

  outcome <- "percentile_correct_raven"
  data <- raven_df
  #data <- data %>% filter(TMAO_tertiles_3rd == t)
  data <- data %>% filter(TMAO_tertiles_3rd %in% all_tertiles)

  
  formulas <- def_formula(outcome, metabolite, trim=3)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_3rd <- formulas$formula_FU2_3rd
  formula_FU2_3rd_interaction <- formulas$formula_FU2_3rd_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_3rd)
  print(formula_FU2_3rd_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_3rd, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_3rd_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
  
  
#   ####### ADHD #####

  outcome <- "ADHD_scaled"
  data <- ADHD_df
  #data <- data %>% filter(TMAO_tertiles_3rd == t)
  data <- data %>% filter(TMAO_tertiles_3rd %in% all_tertiles)

  
  formulas <- def_formula(outcome, metabolite, trim=3)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_3rd <- formulas$formula_FU2_3rd
  formula_FU2_3rd_interaction <- formulas$formula_FU2_3rd_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_3rd)
  print(formula_FU2_3rd_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_3rd, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_3rd_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
  
  
#   ####### CBCL Gen Tot #####

  outcome <- "Gen_Tot_scaled"
  data <- CBCL_df
  #data <- data %>% filter(TMAO_tertiles_3rd == t)
  data <- data %>% filter(TMAO_tertiles_3rd %in% all_tertiles)
  
  formulas <- def_formula(outcome, metabolite, trim=3)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_3rd <- formulas$formula_FU2_3rd
  formula_FU2_3rd_interaction <- formulas$formula_FU2_3rd_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_3rd)
  print(formula_FU2_3rd_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_3rd, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_3rd_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)

  
#   ####### CBCL Gen Ext #####


  outcome <- "Gen_Ext_scaled"
  data <- CBCL_df
  #data <- data %>% filter(TMAO_tertiles_3rd == t)
  data <- data %>% filter(TMAO_tertiles_3rd %in% all_tertiles)

  formulas <- def_formula(outcome, metabolite, trim=3)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_3rd <- formulas$formula_FU2_3rd
  formula_FU2_3rd_interaction <- formulas$formula_FU2_3rd_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_3rd)
  print(formula_FU2_3rd_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_3rd, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_3rd_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
  
  
#   ####### CBCL Gen Int #####


  outcome <- "Gen_Int_scaled"
  data <- CBCL_df
  #data <- data %>% filter(TMAO_tertiles_3rd == t)
  data <- data %>% filter(TMAO_tertiles_3rd %in% all_tertiles)
  formulas <- def_formula(outcome, metabolite, trim=3)
  formula_main <- formulas$formula_main
  formula_main_interaction <- formulas$formula_main_interaction
  formula_FU1 <- formulas$formula_FU1
  formula_FU1_interaction <- formulas$formula_FU1_interaction
  formula_FU2_3rd <- formulas$formula_FU2_3rd
  formula_FU2_3rd_interaction <- formulas$formula_FU2_3rd_interaction
  
  
  print(formula_main)
  print(formula_main_interaction)
  print(formula_FU1)
  print(formula_FU1_interaction)
  print(formula_FU2_3rd)
  print(formula_FU2_3rd_interaction)
  
  
  
  print(outcome)
  
  # Remove participants with NA value in the outcome
  filtered_data <- data %>% filter(!is.na(.[[outcome]]))
  
  # Remove participants with NA value in the metabolite
  filtered_data <- filtered_data %>% filter(!is.na(.[[metabolite]]))
  
  # Number of participants with non NA values for that outcome (= sample size)
  n <- nrow(filtered_data)
  
  model_main <- lmer(formula_main, data = filtered_data)
  model_main_interaction <- lmer(formula_main_interaction, data = filtered_data)
  
  model_FU1 <- lmer(formula_FU1, data = filtered_data)
  model_FU1_interaction <- lmer(formula_FU1_interaction, data = filtered_data)
  
  model_FU2 <- lmer(formula_FU2_3rd, data = filtered_data)
  model_FU2_interaction <- lmer(formula_FU2_3rd_interaction, data = filtered_data)
  
  CI_main <- confint(model_main)
  CI_main_interaction <- confint(model_main_interaction)
  CI_FU1 <- confint(model_FU1)
  CI_FU1_interaction <- confint(model_FU1_interaction)
  CI_FU2 <-  confint(model_FU2)
  CI_FU2_interaction <-  confint(model_FU2_interaction)
  
  # MAIN
  results_MAIN<- as.data.frame(summary(model_main)$coefficients)[2,]
  results_MAIN <-results_MAIN %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main[4,1],
           high_CI = CI_main[4,2])
  
  results_MAIN_interaction<- as.data.frame(summary(model_main_interaction)$coefficients)[2,]
  results_MAIN_interaction <-results_MAIN_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_main_interaction[4,1],
           high_CI = CI_main_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_main, model_main_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_MAIN$pval_LRT <- pval
  results_MAIN_interaction$pval_LRT <- pval
  
  
  results_MAIN_df <- rbind(results_MAIN_df, results_MAIN, results_MAIN_interaction)
  
  #   # FU1
  results_FU1<- as.data.frame(summary(model_FU1)$coefficients)[2,]
  results_FU1 <-results_FU1 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1[4,1],
           high_CI = CI_FU1[4,2])
  
  results_FU1_interaction<- as.data.frame(summary(model_FU1_interaction)$coefficients)[2,]
  results_FU1_interaction <-results_FU1_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU1_interaction[4,1],
           high_CI = CI_FU1_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU1, model_FU1_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU1$pval_LRT <- pval
  results_FU1_interaction$pval_LRT <- pval
  
  
  results_FU1_df <- rbind(results_FU1_df, results_FU1, results_FU1_interaction)
  
  #   # FU2
  results_FU2<- as.data.frame(summary(model_FU2)$coefficients)[2,]
  results_FU2 <-results_FU2 %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = FALSE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2[4,1],
           high_CI = CI_FU2[4,2])
  
  results_FU2_interaction<- as.data.frame(summary(model_FU2_interaction)$coefficients)[2,]
  results_FU2_interaction <-results_FU2_interaction %>% # The second line of the model summary is the metabolite info results
    mutate(outcome = outcome,
           metabolite = metabolite,
           interaction_term = TRUE,
           trimester = "3rd",
           sample_size = n,
           low_CI = CI_FU2_interaction[4,1],
           high_CI = CI_FU2_interaction[4,2])
  
  # LIKELIHOOD RATIO TEST BETWEEN MODELS WITH/WITHOUT INTERACTION TERM
  lrt <- lrtest(model_FU2, model_FU2_interaction)
  pval <- lrt$'Pr(>Chisq)'[2]
  
  results_FU2$pval_LRT <- pval
  results_FU2_interaction$pval_LRT <- pval
  
  
  results_FU2_df <- rbind(results_FU2_df, results_FU2, results_FU2_interaction)
}

results_MAIN_df <- results_MAIN_df %>% select(outcome, metabolite, everything())
results_FU1_df <- results_FU1_df %>% select(outcome, metabolite, everything())
results_FU2_df <- results_FU2_df %>% select(outcome, metabolite, everything())



write.csv(results_MAIN_df, file = "results/results_contaminants_MAIN_interaction.csv", row.names=FALSE)
write.csv(results_FU1_df, file = "results/results_contaminants_FU1_interaction.csv", row.names=FALSE)
write.csv(results_FU2_df, file = "results/results_contaminants_FU2_interaction.csv", row.names=FALSE)
