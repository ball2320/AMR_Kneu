
# Data preparation --------------------------------------------------------

rm(list = ls())
library(sf)
library(sp)
library(dplyr)


# Maps section ------------------------------------------------------------
# Please skip this section to other step


setwd("C:/Users/naekaphirat/Desktop/KP/shp")


## Load global shapefiles (SF downloadable from IHME)

world <- st_read("GBD2023_analysis_final_loc_set_22.shp")
class(world)
head(world)


## Convert to data frame and remove geometry column
world2 <- world %>%
  st_drop_geometry() %>%  # Removes spatial data
  select(loc_id, level, loc_name, loc_nm_sh)  # Select only required columns

#write.csv(world2, "world2.csv", row.names = FALSE)




# Lasso section -----------------------------------------------------------


setwd("C:/Users/naekaphirat/Desktop/KP")
## Load pathogen with antimicrobial resistant profiles
library(openxlsx)
kp <- read.xlsx("KP_overall.xlsx", sheet = "analysis2") 

kp <- kp %>%
  filter(abx_class == "3gc") %>% # Keep only rows contain specific abx
  select(ihmelocid, locid, year_id, source, abx_name, abx_class, susceptible, resistant, sample_size) # Select only required columns

# Rename column names
kp <- kp %>%
  rename(location_id = locid, sum = sample_size)

# Load covariates
#covs <- read.csv('covariates_GBD_2023.csv', stringsAsFactors = F) #Obtained from Freddie
covs <- read.csv('Gram_covs_1990to2021.csv', stringsAsFactors = F) #Obtained from Gisela (Abx included)

#Forecasting missing information from 2021- 2024

forecast <- covs[covs$year_id==2021,]
for(i in c(2022,2023,2024)){
  forecast$year_id <- as.numeric(paste0(i))
  covs <- rbind(covs,forecast)
}
forecast <- covs[covs$location_id == 44855,]
forecast$location_id <- 60908L
covs <- rbind(covs,forecast)
forecast$location_id <- 94364L #95069
covs <- rbind(covs,forecast)
forecast$location_id <- 95069L
covs <- rbind(covs,forecast)


#specify the family you are modelling (currently can use binomial or gaussian)
family <- 'binomial'

#specify transformation of the data to do - 'logit', 'log' or NULL
transformation <- NULL

#specify what you columns are
p <- kp$prop                       #the proportion of your indicator successes
n <- kp$resistant            #the number of your indicator successes
d <- kp$sum                     #the denominator (sample size)
w <- NULL                      #the weights to use


#covs_to_include <- c("anc1_coverage_prop", "anc4_coverage_prop", "diabetes_prev_age_std_prop", "latitude", "mean_temperature", "pollution_outdoor_pm25", "pop_dens_300_500_psqkm_pct", "pop_dens_500_1000_psqkm_pct", "pop_dens_over_1000_psqkm_pct", "sanitation_prop", "water_prop", "alc_binge_prop", "maternal_educ_yrs_pc", "prop_urban", "sdi", "prop_pop_agg", "haqi", "tb_infection_prev", "sev_agestd_nutrition_stunting", "adult_hiv_death_rate_both_sexes", "universal_health_coverage")
covs_to_include <- colnames(covs)
#, "under_5_pop"
#covs_to_include <- c("anc1_coverage_prop", "anc4_coverage_prop")

# Rename columns dynamically
rename_if_exists <- function(df, old_name, new_name) {
  if (!is.null(old_name) && old_name %in% colnames(df)) {
    df <- df %>% rename(!!sym(new_name) := all_of(old_name))
  }
  return(df)
}


kp <- rename_if_exists(kp, "sum", "d")
kp <- rename_if_exists(kp, "prop", "p")
kp <- rename_if_exists(kp, "resistant", "n")

# Compute `n` if missing but `p` and `d` exist
if (is.null(n) && !is.null(p) && !is.null(d)) {
  kp <- kp %>% mutate(n = p * d)
}

# Assign default weight if `w` is missing
kp <- if (is.null(w)) {
  kp %>% mutate(w = 1)
} else {
  rename_if_exists(w, "w")
}

# Apply transformations if specified
if (!is.null(transformation)) {
  if (transformation == "log") {
    if (family == "binomial") {
      kp <- kp %>% mutate(across(c(n, d, p), log))
    } else if (family == "gaussian") {
      kp <- kp %>% mutate(n = log(n))
    }
  } else if (transformation == "logit" && family == "binomial") {
    kp <- kp %>% mutate(p = log(p / (1 - p)))
  } else if (transformation == "logit" && family == "gaussian") {
    message("Logit transformation should not be used with Gaussian data")
  }
}

# Shuffle data and assign 5-folds
set.seed(123)  # For reproducibility
kp <- kp %>%
  sample_frac(1) %>%
  mutate(fold_id = cut(seq_len(n()), breaks = 5, labels = FALSE),
         a_rowid = row_number())

library(data.table)
# Filter and merge covariates
covs <- covs %>%
  select(any_of(c(covs_to_include, "location_id", "year_id"))) %>%
  data.table()

# Merge covariates with mydata (keeping all rows from mydata)
#kp <- merge(data.table(kp), covs, by = intersect(names(kp), names(covs)), all.x = TRUE)

kp <- merge(data.table(kp), covs, by = c("location_id", "year_id"), all.x = TRUE)

# Remove NAs based on model type ; I don't need P column at the moment as of resistant number provided
required_cols <- if (family == "binomial") c("n", "d", names(covs)) else c("n", names(covs))
kp <- na.omit(kp, cols = required_cols)

#Add small nunber to prevent an issue occurred from zero value
kp$n <- pmax(kp$n, 1e-6)  # Replace 0 with small positive value
kp$d <- pmax(kp$d, 1e-6)


#2. Run the lasso 
#define what your response variable is in the data
library(glmnet)

if(family == 'binomial'){
  response <- cbind(failures   = pmax(kp$d - kp$n, 1e-6), 
                    successes = pmax(kp$n, 1e-6))
}else if(family == 'gaussian'){
  response <- kp$n
}

#define variables to include (as a matrix)
vars <- as.matrix(kp[, covs_to_include, with = F])
colnames(vars) <- covs_to_include

#fit cross validated lasso to select lambda
cv_lasso = cv.glmnet(x = vars , y= response, family = family, alpha = 1, weights = kp$w, nfolds = 5, foldid = kp$fold_id)

#If prefer accuracy using this option
#Print λ that gives the minimum cross-validated error (best fit).
cv_lasso$lambda.min

#If prefer simpler model
#print out the model coefficients.
cv_lasso$lambda.1se
coef(cv_lasso, s = "lambda.1se")

#trial other lambdas - investigate trying to get ~ the 10 most important covariates
coef(cv_lasso, s = 0.0005)


# Results in PDF ----------------------------------------------------------


# Load required libraries
library(glmnet)
library(gridExtra)

# Set the output PDF file
#pdf("Lasso_Results_3rd_cep.pdf", width = 8, height = 12)
#pdf("Lasso_Results_4th_cep.pdf", width = 8, height = 12)
pdf("Lasso_Results.pdf", width = 8, height = 12)

# -----------------------------------
# 1. Plot the Lasso Cross-Validation
# -----------------------------------
plot(cv_lasso)
title(main = "Lasso Cross-Validation Plot", line = 2)

# Extract optimal lambda values
lambda_min <- cv_lasso$lambda.min
lambda_1se <- cv_lasso$lambda.1se

# Extract number of selected variables
nvars_min <- sum(coef(cv_lasso, s = lambda_min) != 0) - 1
nvars_1se <- sum(coef(cv_lasso, s = lambda_1se) != 0) - 1

# Add a text summary below the plot
text(x = -6, y = 0.3,
     labels = paste("Optimal Lambda (min):", round(lambda_min, 5),
                    "\nSelected Variables:", nvars_min),
     pos = 4)

text(x = -5, y = 0.35,
     labels = paste("Optimal Lambda (1se):", round(lambda_1se, 5),
                    "\nSelected Variables:", nvars_1se),
     pos = 4)

# -----------------------------------
# 2. Page for Coefficient Matrix at lambda.min
# -----------------------------------
coef_min <- coef(cv_lasso, s = lambda_min)
coef_min_df <- as.data.frame(as.matrix(coef_min))
coef_min_df <- coef_min_df[coef_min_df[,1] != 0, , drop = FALSE] # Remove zero coefficients
coef_min_df$Variable <- rownames(coef_min_df)

# Create a clean page for lambda.min coefficients
plot.new()
title("Coefficient Estimates at Lambda.min")
grid.table(coef_min_df, rows = NULL)


# -----------------------------------
# 3. Page for Coefficient Matrix at lambda.1se
# -----------------------------------
coef_1se <- coef(cv_lasso, s = lambda_1se)
coef_1se_df <- as.data.frame(as.matrix(coef_1se))
coef_1se_df <- coef_1se_df[coef_1se_df[,1] != 0, , drop = FALSE] # Remove zero coefficients
coef_1se_df$Variable <- rownames(coef_1se_df)

# Create a clean page for lambda.1se coefficients
plot.new()
title("Coefficient Estimates at Lambda.1se")
grid.table(coef_1se_df, rows = NULL)


# -----------------------------------
# 4. Summary Page
# -----------------------------------
plot.new()
title("Lasso Regression Summary")

text(0.3, 0.9, paste("Optimal Lambda (minimize):", round(lambda_min, 5)))
text(0.3, 0.8, paste("Number of Selected Variables:", nvars_min))

text(0.3, 0.7, paste("Optimal Lambda (1se):", round(lambda_1se, 5)))
text(0.3, 0.6, paste("Number of Selected Variables:", nvars_1se))

text(0.3, 0.4, paste("Note: Variables with zero coefficients were excluded."))


# -----------------------------------
# 5. Close the PDF
# -----------------------------------
dev.off()

# Success message
cat("✅ PDF file 'Lasso_Results_ has been successfully generated!\n")




#~~~~~~~~~~~~~~~~~~~~~#
#       END           #
#~~~~~~~~~~~~~~~~~~~~~#



# -----------------------------------
# Collinearity checking 
# -----------------------------------



# Correlation matrix ------------------------------------------------------
# Problematic if |r| > 0.7 (R-squared > 0.49)

library(corrplot)

cor_matrix <- cor(covs[, c("location_id", "year_id", "contra_mod_prev_prop", "dtp3_coverage_prop",
                           "hib3_coverage_prop", "hospital_beds_per1000", "latitude", "pigs_pc",
                           "pollution_outdoor_pm25", "maternal_educ_yrs_pc",
                           "age_std_hiv_prev", "he_cap", "frac_oop_hexp", "tb_strains_transmission_risk",
                           "sev_agestd_nutrition_wasting", "sev_scalar_agestd_hiv", "bcg_vacc_cov_prop",
                           "physicians_pc", "pharmacists_pc", "sev_agestd_smoking_direct", "lpc",
                           "smok_daily_prev_agstd_both", "cce", "gee", "rle", "raw_ddd_per_1000",
                           "prop_j01e_final", "prop_j01f_final")], use = "complete.obs")


cor_matrix <- cor(covs[, c("location_id", "year_id", "contra_mod_prev_prop",
                           "pigs_pc",
                           "pollution_outdoor_pm25",
                           "age_std_hiv_prev", "tb_strains_transmission_risk",
                           "sev_scalar_agestd_hiv",
                           "physicians_pc",
                           "smok_daily_prev_agstd_both", "raw_ddd_per_1000")], use = "complete.obs")

cor_matrix_plot <- cor_matrix
cor_matrix_plot[abs(cor_matrix_plot) < 0.4] <- NA
corrplot.mixed(cor_matrix_plot, upper = "color", lower = "number",
               tl.cex = 0.7, number.cex = 0.6)

# R-squared < 0.4
cor_matrix <- cor(covs[, c("location_id", "year_id",
                           "pigs_pc",
                           "pollution_outdoor_pm25",
                           "age_std_hiv_prev", "tb_strains_transmission_risk",
                           "sev_scalar_agestd_hiv",
                           "physicians_pc",
                           "smok_daily_prev_agstd_both", "raw_ddd_per_1000")], use = "complete.obs")

print(cor_matrix)
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.7)


# Variance Inflation factor -----------------------------------------------
# < 5 is acceptable, > 5 is moderate, > 10 is problematic

# Reference
#Collinearity: a review of methods to deal with it and a simulation study evaluating their performance.
#Ecography, 36(1), 27–46.
#DOI: 10.1111/j.1600-0587.2012.07348.x
#In ecological and spatial modeling, they review simulation results and recommend |r| ≥ 0.7–0.9 as practical thresholds for collinearity screening.


library(car)
model <- lm(n ~ location_id + year_id + pigs_pc + pollution_outdoor_pm25 + age_std_hiv_prev +
              tb_strains_transmission_risk + sev_scalar_agestd_hiv + physicians_pc +
              smok_daily_prev_agstd_both + raw_ddd_per_1000, data = kp)

#Compute VIF values
vif_values <- vif(model)

# Create color categories based on thresholds
bar_colors <- ifelse(vif_values > 5, "darkred",               # severe
                     ifelse(vif_values > 2, "red",            # moderate
                            "lightgreen"))                    # acceptable

#Visualize VIF values and represent the possibility of problematic
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = bar_colors, xlim = c(0, 10))
abline(v = 2, lty = 2, col = "lightgreen")   # Acceptable
abline(v = 5, lty = 2, col = "red")   # Moderate threshold
abline(v = 10, lty = 2, col = "darkred")  # Severe threshold

# Add background text
text(x = 3.5, y = length(vif_values) / 2, labels = "Moderate (2–5)", col = "red", cex = 0.9, srt = 90)
text(x = 7.5, y = length(vif_values) / 2, labels = "Severe (>5)", col = "darkred", cex = 0.9, srt = 90)


# -----------------------------------
# Univariate checking 
# -----------------------------------

# Define outcome and predictor list
response <- "n"
predictors <- c("location_id", "year_id", "pigs_pc", "pollution_outdoor_pm25", "age_std_hiv_prev",
                "tb_strains_transmission_risk", "sev_scalar_agestd_hiv", "physicians_pc",
                "smok_daily_prev_agstd_both", "raw_ddd_per_1000")  # add all variable names here

# Run univariate models automatically
results <- lapply(predictors, function(var) {
  formula <- as.formula(paste(response, "~", var))
  model <- lm(formula, data = kp)
  
  data.frame(
    Variable = var,
    Estimate = coef(summary(model))[2, "Estimate"],
    StdError = coef(summary(model))[2, "Std. Error"],
    p_value = coef(summary(model))[2, "Pr(>|t|)"],
    Adj_R2 = summary(model)$adj.r.squared,
    AIC = AIC(model),
    logLik = as.numeric(logLik(model))
  )
})

# Combine all results into a data frame
univariate_results <- do.call(rbind, results)
# Rank predictors by raking the lowest AIC
univariate_results <- univariate_results[order(univariate_results$AIC), ]

# View the summary table
print(univariate_results)

#Model comparison
model_a <- lm(n ~ location_id + year_id + pigs_pc + pollution_outdoor_pm25 + age_std_hiv_prev +
              tb_strains_transmission_risk + sev_scalar_agestd_hiv + physicians_pc +
              smok_daily_prev_agstd_both + raw_ddd_per_1000, data = kp)

model_b <- lm(n ~ year_id + pigs_pc + age_std_hiv_prev +
              sev_scalar_agestd_hiv + physicians_pc +
              smok_daily_prev_agstd_both + raw_ddd_per_1000, data = kp)

model_c <- lm(n ~ year_id + pigs_pc + age_std_hiv_prev +
                sev_scalar_agestd_hiv + physicians_pc +
                smok_daily_prev_agstd_both, data = kp)

model_d <- lm(n ~ year_id + pigs_pc + age_std_hiv_prev +
                sev_scalar_agestd_hiv +
                smok_daily_prev_agstd_both, data = kp)


AIC(model_a, model_b, model_c, model_d)
BIC(model_a, model_b, model_c, model_d)

# Choosing between model_c and model_d

# https://easystats.github.io/bayestestR/reference/bayesfactor_models.html
# A Bayes factor greater than 1 can be interpreted as evidence against the null,
# at which one convention is that a Bayes factor greater than 3 can be considered as "substantial"
# evidence against the null (and vice versa, a Bayes factor smaller than 1/3 indicates substantial
# evidence in favor of the null-model) (Wetzels et al. 2011).

library(bayestestR)
bic_to_bf(BIC(model_c), denominator = BIC(model_a)) #*****
bic_to_bf(BIC(model_d), denominator = BIC(model_a))
