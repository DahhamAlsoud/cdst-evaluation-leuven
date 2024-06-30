# define initial variables ----
regression_vrs <- c("gender", "age_yrs", "disease_duration_yrs", "smoking", "bmi_kg_m2", "albumin_g_L","crp_mg_L", 
                    "concomitant_cs", "concomitant_imm", "concomitant_5asa", "previous_antitnf", 
                    "intestinal_resections_before", "age_at_diagnosis_yrs", 
                    "disease_location", "upper_gi_disease", "disease_behaviour","peri_anal_disease",
                    "disease_extent", "mes",
                    "fistulizing_disease_prior", "fistulizing_disease_current", "clinical_remission_cdst") # compared to tblone vrs, we excluded ses_cd & fcal due to missings
# vdz-cd ----
## prepare regression df ----
vdz_cd_regression_df <- vdz_cd_included %>% select(any_of(regression_vrs))

vdz_cd_regression_df %<>% mutate(across(any_of(continuous_vrs), ~ round(as.numeric(.), 1)))

vdz_cd_regression_df %<>% mutate(disease_behaviour = if_else(str_detect(disease_behaviour, "B3"), "B3", disease_behaviour))

vdz_cd_regression_df %<>% mutate(prior_smoker = if_else(smoking == 0, 0, 1), .after = smoking) %>% select(-smoking)

vdz_cd_regression_df %<>% get_age_montreal() %>% select(-age_at_diagnosis_yrs)

vdz_cd_regression_df %<>% recode_steroids()

vdz_cd_regression_df %<>% mutate(concomitant_imm = if_else(concomitant_imm == "0", 0, 1),
                             concomitant_5asa = if_else(concomitant_5asa == "0", 0, 1))

vdz_cd_regression_df <- recode_vrs(vdz_cd_regression_df, data_dictionary = ibd_data_dict, vrs = categoricals_to_recode, factor = TRUE)

## define the formula ----
pred <- names(vdz_cd_regression_df)[-length(vdz_cd_regression_df)]
formula <- paste("clinical_remission_cdst~", paste(pred, collapse = "+"))

## estimate the full model ----
full_mod <- glm(formula, data = vdz_cd_regression_df, x = T, y = T, family = "binomial")
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
summary(full_mod)

vdz_cd_full_mod <- full_mod
vdz_cd_full_est <- full_est
vdz_cd_full_se <- full_se

## select and estimate the best model ----
sel_mod <- step(glm(formula, data = vdz_cd_regression_df, x = T, y = T, family = "binomial"), 
                direction = "backward",
                scope = list(upper = formula, 
                             lower = formula(clinical_remission_cdst~1)),
                trace = 0)

summary(sel_mod)

sel_est <- coef(sel_mod)[c("(Intercept)", names(coef(full_mod)))]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", names(coef(full_mod)))
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", names(coef(full_mod)))]
sel_se[is.na(sel_se)] <- 0
names(sel_se) <- c("(Intercept)", names(coef(full_mod)))

vdz_cd_sel_mod <- sel_mod
vdz_cd_sel_est <- sel_est
vdz_cd_sel_se <- sel_se

## perform bootstraps ----
bootnum <- 1000
boot_est <-  boot_se <- matrix(0, ncol = length(names(coef(full_mod))) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", names(coef(full_mod)))))

set.seed(5437854)
for (i in 1:bootnum) {
  data_id <- sample(1:dim(vdz_cd_regression_df)[1], replace = T)
  boot_mod <- step(glm(formula, data = vdz_cd_regression_df[data_id, ], 
                      x=T, y=T, family = "binomial"), 
                   scope = list(upper = formula, 
                                lower = 
                                  formula(clinical_remission_cdst ~ 1)),
                   direction = "backward", trace = 0)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## get frequencies ----
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))

vdz_cd_vrs_frq <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se, 
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
vdz_cd_vrs_frq <- vdz_cd_vrs_frq[order(vdz_cd_vrs_frq[,"boot_inclusion"], decreasing=T),]
vdz_cd_vrs_frq

# vdz-uc ----
## prepare regression df ----
vdz_uc_regression_df <- vdz_uc_included %>% select(any_of(regression_vrs))

vdz_uc_regression_df %<>% select(-intestinal_resections_before) # zero variance (all is 0)

vdz_uc_regression_df %<>% mutate(across(any_of(continuous_vrs), ~ round(as.numeric(.), 1)))

vdz_uc_regression_df %<>% mutate(prior_smoker = if_else(smoking == 0, 0, 1), .after = smoking) %>% select(-smoking)

vdz_uc_regression_df %<>% get_age_montreal() %>% select(-age_at_diagnosis_yrs)

vdz_uc_regression_df %<>% recode_steroids()

vdz_uc_regression_df %<>% mutate(concomitant_imm = if_else(concomitant_imm == "0", 0, 1),
                                 concomitant_5asa = if_else(concomitant_5asa == "0", 0, 1))


vdz_uc_regression_df <- recode_vrs(vdz_uc_regression_df, data_dictionary = ibd_data_dict, vrs = categoricals_to_recode, factor = TRUE)

## define the formula ----
pred <- names(vdz_uc_regression_df)[-length(vdz_uc_regression_df)]
formula <- paste("clinical_remission_cdst~", paste(pred, collapse = "+"))

## estimate the full model ----
full_mod <- glm(formula, data = vdz_uc_regression_df, x = T, y = T, family = "binomial")
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
summary(full_mod)

vdz_uc_full_mod <- full_mod
vdz_uc_full_est <- full_est
vdz_uc_full_se <- full_se

## select and estimate the best model ----
sel_mod <- step(glm(formula, data = vdz_uc_regression_df, x = T, y = T, family = "binomial"), 
                direction = "backward",
                scope = list(upper = formula, 
                             lower = formula(clinical_remission_cdst~1)),
                trace = 0)

summary(sel_mod)

sel_est <- coef(sel_mod)[c("(Intercept)", names(coef(full_mod)))]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", names(coef(full_mod)))
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", names(coef(full_mod)))]
sel_se[is.na(sel_se)] <- 0
names(sel_se) <- c("(Intercept)", names(coef(full_mod)))

vdz_uc_sel_mod <- sel_mod
vdz_uc_sel_est <- sel_est
vdz_uc_sel_se <- sel_se

## perform bootstraps ----
bootnum <- 1000
boot_est <-  boot_se <- matrix(0, ncol = length(names(coef(full_mod))) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", names(coef(full_mod)))))

set.seed(5437854)
for (i in 1:bootnum) {
  data_id <- sample(1:dim(vdz_uc_regression_df)[1], replace = T)
  boot_mod <- step(glm(formula, data = vdz_uc_regression_df[data_id, ], 
                       x=T, y=T, family = "binomial"), 
                   scope = list(upper = formula, 
                                lower = 
                                  formula(clinical_remission_cdst ~ 1)),
                   direction = "backward", trace = 0)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## get frequencies ----
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))

vdz_uc_vrs_frq <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se, 
                              rmsdratio, boot_relbias, boot_median, boot_025per, 
                              boot_975per), 4)
vdz_uc_vrs_frq <- vdz_uc_vrs_frq[order(vdz_uc_vrs_frq[,"boot_inclusion"], decreasing=T),]
vdz_uc_vrs_frq


# ust-cd ----
## prepare regression df ----
ust_cd_regression_df <- ust_cd_included %>% select(any_of(regression_vrs))

ust_cd_regression_df %<>% select(-c(concomitant_imm, concomitant_5asa)) # zero variance (all is 0)

ust_cd_regression_df %<>% mutate(across(any_of(continuous_vrs), ~ round(as.numeric(.), 1)))

ust_cd_regression_df %<>% mutate(disease_behaviour = if_else(str_detect(disease_behaviour, "B3"), "B3", disease_behaviour))

ust_cd_regression_df %<>% mutate(prior_smoker = if_else(smoking == 0, 0, 1), .after = smoking) %>% select(-smoking)

ust_cd_regression_df %<>% get_age_montreal() %>% select(-age_at_diagnosis_yrs)

ust_cd_regression_df %<>% recode_steroids()

ust_cd_regression_df <- recode_vrs(ust_cd_regression_df, data_dictionary = ibd_data_dict, vrs = categoricals_to_recode, factor = TRUE)

## define the formula ----
pred <- names(ust_cd_regression_df)[-length(ust_cd_regression_df)]
formula <- paste("clinical_remission_cdst~", paste(pred, collapse = "+"))

## estimate the full model ----
full_mod <- glm(formula, data = ust_cd_regression_df, x = T, y = T, family = "binomial")
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
summary(full_mod)

ust_cd_full_mod <- full_mod
ust_cd_full_est <- full_est
ust_cd_full_se <- full_se

## select and estimate the best model ----
sel_mod <- step(glm(formula, data = ust_cd_regression_df, x = T, y = T, family = "binomial"), 
                direction = "backward",
                scope = list(upper = formula, 
                             lower = formula(clinical_remission_cdst~1)),
                trace = 0)

summary(sel_mod)

sel_est <- coef(sel_mod)[c("(Intercept)", names(coef(full_mod)))]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", names(coef(full_mod)))
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", names(coef(full_mod)))]
sel_se[is.na(sel_se)] <- 0
names(sel_se) <- c("(Intercept)", names(coef(full_mod)))

ust_cd_sel_mod <- sel_mod
ust_cd_sel_est <- sel_est
ust_cd_sel_se <- sel_se

## perform bootstraps ----
bootnum <- 1000
boot_est <-  boot_se <- matrix(0, ncol = length(names(coef(full_mod))) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", names(coef(full_mod)))))

set.seed(5437854)
for (i in 1:bootnum) {
  data_id <- sample(1:dim(ust_cd_regression_df)[1], replace = T)
  boot_mod <- step(glm(formula, data = ust_cd_regression_df[data_id, ], 
                       x=T, y=T, family = "binomial"), 
                   scope = list(upper = formula, 
                                lower = 
                                  formula(clinical_remission_cdst ~ 1)),
                   direction = "backward", trace = 0)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## get frequencies ----
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))

ust_cd_vrs_frq <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se, 
                              rmsdratio, boot_relbias, boot_median, boot_025per, 
                              boot_975per), 4)
ust_cd_vrs_frq <- ust_cd_vrs_frq[order(ust_cd_vrs_frq[,"boot_inclusion"], decreasing=T),]
ust_cd_vrs_frq

# export results ----
vdz_cd_vrs_frq %<>% as.data.frame.matrix() %>% rownames_to_column() 
export(vdz_cd_vrs_frq, "vdz_cd_vrs_frq.xlsx")
vdz_uc_vrs_frq %<>% as.data.frame.matrix() %>% rownames_to_column()
export(vdz_uc_vrs_frq, "vdz_uc_vrs_frq.xlsx")
ust_cd_vrs_frq %<>% as.data.frame.matrix() %>% rownames_to_column()
export(ust_cd_vrs_frq, "ust_cd_vrs_frq.xlsx")
