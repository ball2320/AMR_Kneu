
rm(list = ls())
setwd("C:/Users/naekaphirat/Desktop/kp")

# Load abx data
sum_dat <- readRDS("C:/Users/naekaphirat/Desktop/IHME_datasets/sum_dat")

# Load covariates
covs <- read.csv('covariates_GBD_2023.csv', stringsAsFactors = F)

kp <- merge(sum_dat, covs, by = c("location_id", "year_id"), all.x = TRUE)

# Filtering abx_class
kp <- kp %>%
  filter(res > 0, sum > 10, abx_class == "carbapenem")

library(dplyr)
kp2 <- kp %>% select(-loc_year, -source, -location_id, -year_id, -hosp_group,
                     -specimen_group, -abx_class, -sus, -unknown, -sum)


# List of abx interested
int_abx <- c("3gc", "carbapenem", "aminoglycoside", "fluoroquinolone")


library(mpath)
library(glmnet)
library(gridExtra)

# Remove predictors with only one unique value (to avoid contrast errors)
kp2 <- kp2[, sapply(kp2, function(x) length(unique(x)) > 1)]


# Negative Binomial LASSO selection ---------------------------------------------------------


# Shuffle data and assign 5-folds
set.seed(123)  # For reproducibility

cv_nb <- cv.glmregNB(res ~ ., data = kp2, nfolds = 5)  # performs cross-validation to choose lambda

# Cross-validation curve
plot(cv_nb)
title("Negative Binomial LASSO Cross-Validation", line = 2)

# Optimal lambda values
cv_nb$lambda.min
cv_nb$lambda.1se

# Coefficients at optimal lambda
coef_nb <- coef(cv_nb, s = "lambda.min")
print(coef_nb)

# List selected variables (non-zero coefficients)
selected_vars <- rownames(coef_nb)[which(coef_nb != 0)]
selected_vars

pred_nb <- predict(cv_nb, newx = kp2, s = "lambda.min", type = "response")
head(pred_nb)

coef_df <- as.data.frame(as.matrix(coef_nb))
coef_df <- subset(coef_df, V1 != 0)
coef_df$Variable <- rownames(coef_df)

grid.table(coef_df, rows = NULL)



# Filter and merge covariates
covs <- covs %>%
  select(any_of(c(covs_to_include, "location_id", "year_id"))) %>%
  data.table()

# Merge covariates with mydata (keeping all rows from mydata)
kp <- merge(data.table(kp), covs, by = c("location_id", "year_id"), all.x = TRUE)
# Remove all rows that contain any "NA"
kp <- na.omit(kp)


#Add small nunber to prevent an issue occurred from zero value
kp$res <- pmax(kp$res, 1e-6)  # Replace 0 with small positive value
kp$sum <- pmax(kp$sum, 1e-6)


#2. Run the lasso 
#define what your response variable is in the data
library(glmnet)

if(family == 'binomial'){
  response <- cbind(failures   = pmax(kp$sum - kp$res, 1e-6), 
                    successes = pmax(kp$res, 1e-6))
}else if(family == 'gaussian'){
  response <- kp$res
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

# Save Lasso Results to PDF
pdf("Prev_Lasso_carbapenem.pdf", width = 8, height = 12)

# --- 1. Cross-Validation Plot ---
plot(cv_lasso, main = "Lasso Cross-Validation")

lambda_min <- cv_lasso$lambda.min
lambda_1se <- cv_lasso$lambda.1se

nvars_min <- sum(coef(cv_lasso, s = lambda_min) != 0) - 1
nvars_1se <- sum(coef(cv_lasso, s = lambda_1se) != 0) - 1

text(-6, 0.3, sprintf("Lambda (min): %.5f\nSelected Vars: %d", lambda_min, nvars_min), pos = 4)
text(-5, 0.35, sprintf("Lambda (1se): %.5f\nSelected Vars: %d", lambda_1se, nvars_1se), pos = 4)

# --- Helper function to extract and plot coefficients ---
plot_coef_table <- function(cv_fit, lambda, label) {
  coefs <- coef(cv_fit, s = lambda)
  df <- subset(data.frame(Variable = rownames(as.matrix(coefs)),
                          Coefficient = as.numeric(coefs)),
               Coefficient != 0)
  plot.new()
  title(paste("Coefficient Estimates at", label))
  grid.table(df, rows = NULL)
}

# --- 2. Coefficients at lambda.min ---
plot_coef_table(cv_lasso, lambda_min, "Lambda.min")

# --- 3. Coefficients at lambda.1se ---
plot_coef_table(cv_lasso, lambda_1se, "Lambda.1se")

# --- 4. Summary Page ---
plot.new()
title("Lasso Regression Summary")
summary_text <- sprintf(
  "Lambda (min): %.5f\nSelected Vars: %d\n\nLambda (1se): %.5f\nSelected Vars: %d\n\nNote: Zero coefficients excluded.",
  lambda_min, nvars_min, lambda_1se, nvars_1se
)
text(0.3, 0.7, summary_text, pos = 4)

# --- 5. Close PDF ---
dev.off()




