
rm(list = ls())
#setwd("~/Library/CloudStorage/OneDrive-Nexus365/Desktop/GRAM/Data")
setwd("C:/Users/naekaphirat/Desktop/kp")


# Load required libraries
library(rstan)
library(parallel)
options(mc.cores = parallel::detectCores()) # set stan to use all avalable cores
rstan_options(auto_write = TRUE) # save chain results
library(ggplot2)
library(dplyr)
library(readr)
library(loo)
library(readxl)

# Load cleaned data
kp <- readRDS("C:/Users/naekaphirat/Desktop/IHME_datasets/sum_dat")

kp <- kp %>%
     group_by(location_id, year_id, source, abx_class) %>%
     summarize(
         resistant = sum(res, na.rm = TRUE),
         susceptible = sum(sus, na.rm = TRUE),
         unknown = sum(unknown, na.rm = TRUE),
         sample_size = sum(res + sus, na.rm = TRUE),
         .groups = 'drop',
     )

# Change region_id to location_id for some countries.

#kp$location_id[kp$location_id == 4625] <- 95 # South East England -> UK
#kp$location_id[kp$location_id == 53674] <- 51 # Wielkopolskie -> Poland
#kp$location_id[kp$location_id == 53667] <- 51 # Opolskie -> Poland
#kp$location_id[kp$location_id == 53663] <- 51 # Lubelskie -> Poland
#kp$location_id[kp$location_id == 53662] <- 51 # Łódzkie -> Poland
#kp$location_id[kp$location_id == 53660] <- 51 # DolnoÅ›lÄ…skie -> Poland
#kp$location_id[kp$location_id == 53621] <- 165 # Sindh -> Pakistan
#kp$location_id[kp$location_id == 53620] <- 163 # Punjab -> India
#kp$location_id[kp$location_id == 53618] <- 165 # Islamabad Capital Territory -> Pakistan
#kp$location_id[kp$location_id == 53614] <- 163 # National Capital Region -> India
#kp$location_id[kp$location_id == 44955] <- 62 # Nizhny Novgorod Oblast -> Russia
#kp$location_id[kp$location_id == 44870] <- 142 # Fars -> Iran
#kp$location_id[kp$location_id == 44861] <- 179 # Addis Ababa -> Ethiopia
#kp$location_id[kp$location_id == 35646] <- 180 # Nairobi -> Kenya
#kp$location_id[kp$location_id == 35630] <- 180 # Kilifi -> Kenya
#kp$location_id[kp$location_id == 25349] <- 214 # Plateau -> Nigeria ???
#kp$location_id[kp$location_id == 25337] <- 214 # Kano -> Nigeria
#kp$location_id[kp$location_id == 25332] <- 214 # FCT (Abuja) -> Nigeria
#kp$location_id[kp$location_id == 4875] <- 51 # West Bengal -> India
#kp$location_id[kp$location_id == 4874] <- 51 # Uttarakhand -> India
#kp$location_id[kp$location_id == 4873] <- 51 # Uttar Pradesh -> India
#kp$location_id[kp$location_id == 4870] <- 51 # Tamil Nadu -> India
#kp$location_id[kp$location_id == 4865] <- 51 # Orissa -> India
#kp$location_id[kp$location_id == 4860] <- 51 # Maharashtra -> India
#kp$location_id[kp$location_id == 4857] <- 51 # Kerala -> India
#kp$location_id[kp$location_id == 4856] <- 51 # Karnataka -> India
#kp$location_id[kp$location_id == 4849] <- 51 # Delhi -> India
#kp$location_id[kp$location_id == 4841] <- 51 # Andhra Pradesh -> India
#kp$location_id[kp$location_id == 4775] <- 135 # Sao Paulo -> Brazil
#kp$location_id[kp$location_id == 4755] <- 135 # CearÃ¡ -> Brazil
#kp$location_id[kp$location_id == 4751] <- 135 # Alagoas -> Brazil
#kp$location_id[kp$location_id == 4670] <- 130 # Tamaulipas -> Mexico
#kp$location_id[kp$location_id == 4656] <- 130 # Jalisco-> Mexico
#kp$location_id[kp$location_id == 4651] <- 130 # Distrito Federal-> Mexico
#kp$location_id[kp$location_id == 4643] <- 130 # Aguascalientes-> Mexico
#kp$location_id[kp$location_id == 572] <- 102 # Wisconsin -> United States
#kp$location_id[kp$location_id == 570] <- 102 # Washington -> United States
#kp$location_id[(kp$location_id <= 572) & (kp$location_id >= 500)] <- 102 # -> United States
#kp$location_id[kp$location_id == 354] <- 6 # Hong Kong -> China
#kp$location_id[kp$location_id == 490] <- 198 # Western Cape -> Zimbabwe
#kp$location_id[kp$location_id == 485] <- 198 # KwaZulu-Natal -> Zimbabwe






# Select only abx of interest , specimen_group == "blood"
kp <- kp %>%
  filter(abx_class == "3gc", sample_size >= 30) # Keep only rows contain specific abx
colnames(kp)[colnames(kp) =="location_id"] <- "loc_id"

# Add IHME countries id
shp_id <- read.csv("C:/Users/naekaphirat/Desktop/oraya/shp_id.csv", stringsAsFactors = F)
colnames(shp_id)[colnames(shp_id) =="loc_name"] <- "country"
kp <- left_join(kp, shp_id, by = "loc_id")

# Temporary. Need to check again
# Remove NA from main dataset
kp <- kp[!is.na(kp$country), ]

#kp <- kp %>%
#  filter(!country == "Uruguay")

# For model 06_all only
kp <- kp %>%
  filter(!country %in% c("Bhutan", "Libya", "Moldova", "Mongolia", "Rwanda", "Sri Lanka", "Uruguay"))  


library(sf)
library(sp)

# Load shapefiles (same source with Oraya)
locs <- st_read("C:/Users/naekaphirat/Desktop/oraya/geo_export_39720b2d-42d1-46a5-8d79-f5c21ec46ccf.shp")



#reformat country names
change_name <- function(df) {
  name_map <- c(
    "Czech Republic" = "Czechia",
    "D. P. R. of Korea" = "South Korea",
    "Islamic Republic of Iran" = "Iran",
    "Lao People's Democratic Republic" = "Laos",
    "United Kingdom" = "UK",
    "United States of America" = "USA",
    "United States" = "USA",
    "R. B. de Venezuela" = "Venezuela",
    "Puerto Rico (U.S.)" = "Puerto Rico",
    "Russian Federation" = "Russia",
    "Vietnam" = "Viet Nam",
    "Syrian Arab Republic" = "Syria",
    "Arab Republic of Egypt" = "Egypt",
    "Slovak Republic" = "Slovakia",
    "occupied Palestinian territory" = "Palestine",
    "Palestinian Territory" = "Palestine",
    "Syrian Arab Republic" = "Syria",
    "Turkey" = "Türkiye",
    "United Republic of Tanzania" = "Tanzania",
    "Swaziland" = "Eswatini",
    "Congo" = "Congo (Brazzaville)",
    "Congo DRC" = "DR Congo",
    "Brunei Darussalam" = "Brunei",
    "The Former Yugoslav Republic of Macedonia" = "North Macedonia",
    "Hong Kong" = "Hong Kong Special Administrative Region of China",
    "US Virgin Islands" = "Virgin Islands",
    "Micronesia" = "Federated States of Micronesia",
    "Hong Kong" = "Hong Kong Special Administrative Region of China",
    "Cape Verde" = "Cabo Verde",
    "Gambia" = "The Gambia",
    "Sao Tome and Principe" = "São Tomé and Príncipe",
    "Bahamas" = "The Bahamas",
    
    #for owid
    "Reunion" = "Réunion",
    "Democratic Republic of Congo" = "DR Congo",
    "East Timor" = "Timor-Leste",
    "Micronesia (country)" = "Federated States of Micronesia",
    "South Georgia and the South Sandwich Islands" = "South Georgia",
    "Cote d'Ivoire" = "Côte d'Ivoire",
    "Vatican" = "Vatican City"
  )
  
  df$country <- dplyr::recode(df$country, !!!name_map)
  return(df)
}

# Sort out empty geometry
locs <- locs %>% 
  filter(!st_is_empty(geometry))


# Join pathogen with shapefiles
kp_locs <- left_join(kp, locs, by = "country")

# Sort out empty geometry rows
kp_locs <- kp_locs %>% 
  filter(!st_is_empty(geometry))
# Sort out NA loc_id rows
kp_locs <- kp_locs[!is.na(kp_locs$loc_id), ]
# Select only unique country name from the dataset
kp_locs_countries <- unique(kp_locs$country)


# Select only countries used
kp_locs <- kp_locs[kp_locs$country %in% kp_locs_countries, ]
shp_id <- shp_id[shp_id$country %in% kp_locs_countries, ]

# Add index to country
shp_id$country_numeric <- 1:nrow(shp_id)
shp_id$row_id <- 1:nrow(shp_id)

# Join shapefiles with IDs
map_filtered <- left_join(shp_id, locs, by = "country")

# Spatial preparation -----------------------------------------------------

map_filtered <- st_as_sf(map_filtered)

library(spdep)
library(spData)

nb <- poly2nb(map_filtered)
adj <- nb2mat(nb,style="B", zero.policy=TRUE)

library(netdiffuseR)
edge_list <- data.frame(adjmat_to_edgelist(adj))
edge_list$node1_locId <- shp_id$country_numeric[match(edge_list[,"ego"], shp_id$row_id)]
edge_list$node2_locId <- shp_id$country_numeric[match(edge_list[,"alter"], shp_id$row_id, nomatch = NA)]

# calculate the centroids of each country
map_filtered_centroids <- st_centroid(st_geometry(map_filtered))

#map_filtered_centroids <- st_centroid(map_filtered)
# calculate the distance matrix between these centroids
#distance_mat <- cbind(st_distance(map_filtered_centroids))

D <- as.matrix(st_distance(map_filtered_centroids)) # Same format with distance_mat
D <- D / max(D) # standardize distance matrix
D2 <- D/(6378100*pi) # For regional data
D2 <- D2 / max(D2) # standardize distance matrix
D <- units::drop_units(D)
D2 <- units::drop_units(D2)

# #test distance-based neighbours -knearneigh ensures all areas have k neighbours
dist_nb <- knn2nb(knearneigh(st_centroid(map_filtered), k=1), row.names=map_filtered$loc_id)
dist_adj <- nb2mat(dist_nb,style="B", zero.policy=TRUE)
edge_list_dist <- data.frame(adjmat_to_edgelist(dist_adj))
edge_list_dist$node1_locId <- shp_id$country_numeric[match(edge_list_dist[,"ego"], shp_id$row_id)]
edge_list_dist$node2_locId <- shp_id$country_numeric[match(edge_list_dist[,"alter"], shp_id$row_id, nomatch = NA)]

N_edges <- nrow(edge_list_dist)
node1 <- edge_list_dist[,"node1_locId"]
node2 <- edge_list_dist[,"node2_locId"]


# Covariates preparation ---------------------------------------------------

covs_IHME <- read.csv("C:/Users/naekaphirat/Desktop/KP/Gram_covs_1990to2021.csv", stringsAsFactors = F)
colnames(covs_IHME)[colnames(covs_IHME) =="location_id"] <- "loc_id"

iden_vars <- c("loc_id", "year_id")
covs <- setdiff(names(covs_IHME), iden_vars)

covs_IHME[covs] <- scale(covs_IHME[covs]) #scale all covariates

#NAs in covs
covs_IHME$raw_ddd_per_1000 <- NULL # remove
covs_IHME$adj_id <- NULL # remove
covs_IHME$sex_id <- NULL # remove
covs_IHME$age_group_id <- NULL # remove

#covar_list <- c("prop_j01c_final","prop_j01d_final", "physicians_pc", "age_std_hiv_prev") # For 3GC
#covar_list <- c("prop_j01c_final","prop_j01d_final")
#covar_list <- c("ddd_per_1000_fitted")
covar_list <- c("hiv_prev_pct")

library(dplyr)

kp_locs_matched <- kp_locs %>%
  inner_join(
    covs_IHME %>% select(loc_id, year_id, all_of(covar_list)),
    by = c("loc_id", "year_id")
  )



cov_mat <- as.matrix(kp_locs_matched[covar_list])
K <- length(covar_list)
kp_locs <- kp_locs_matched

# Stan data preparation
data2 <- kp_locs %>% dplyr::rename(Country = country, Year = year_id)
data2 <- data2 %>% mutate(resistant = floor(resistant + 0.5)) # round up or down by cut-off value
data2 <- data2 %>% mutate(Country_id = as.numeric(as.factor(Country)))
data2 <- data2 %>% mutate(Year = as.numeric(Year), Year_id = Year - 1989) 
data2 <- data2 %>% mutate(Source_id = as.numeric(as.factor(source)))
data2$prop <- data2$resistant / data2$sample_size

#data2 <- data2 %>% mutate(new_source = as.integer(factor(paste(Country_id, Source_id, sep = "_"))))
# Temporary
# Remove country without country information
data2 <- data2[!is.na(data2$Country_id), ]

# Sort data
data2 <- data2[order(data2$Country_id, data2$Year_id), ]

data2$Country_index <- as.integer(factor(data2$Country_id))
data2$Source_index  <- as.integer(factor(data2$Source_id))

##stan-specific vars
datalist <- list(
  N = nrow(data2),
  C = length(unique(data2$Country_id)),
  S = length(unique(data2$Source_id)),
  y = as.integer(data2$resistant),
  t = data2$Year_id,
  cases = as.integer(data2$sample_size),
  #country = as.integer(as.factor(data2$Country_id)),
  country = data2$Country_index,
  #source = as.integer(as.factor(data2$Source_id)),
  source = data2$Source_index,
  N_edges = N_edges,
  node1 = node1,
  node2 = node2,
  D = D,
  D2 = D2,
  calc_likelihood = 1,
  K = K,
  cov_mat = cov_mat
)

# Rstan fit and extract  -----------------------------------------------------
# Fit and extraction function

fit_and_extract <- function(stan_file, model_name, datalist, data2) {
  
  sm <- stan_model(file = stan_file)
  
  fit <- sampling(
    sm,
    data = datalist,
    iter = 200,
    chains = 2,
    control = list(adapt_delta = 0.95, max_treedepth = 30),
    cores = 4
  )
  
  ysum <- as.data.frame(
    summary(fit, pars = "y_pred", probs = c(0.05, 0.5, 0.95))$summary
  )
  
  colnames(ysum)[c(1,4,5,6)] <- c("mean","q5","q50","q95")
  
  cbind(
    ysum,
    t = data2$Year_id,
    resis_prop = data2$prop,
    country = data2$Country_id,
    country_name = data2$Country,
    isolates = data2$sample_size,
    source_name = data2$source,
    year_id = data2$Year,
    model = model_name      
  )
}


model_pred_all <- bind_rows(
  
  fit_and_extract("01_betaGammaI0_country.stan", "01_country", datalist, data2)#,
  #fit_and_extract("02_betaGammaI0_covar.stan", "02_covars", datalist, data2)#,
  #fit_and_extract("03_betaGammaI0_source.stan", "03_source", datalist, data2)#,
  #fit_and_extract("04_spatial_cmdStan.stan", "04_sp", datalist, data2),
  #fit_and_extract("05_spatialcmdStanSource.stan", "05_sp+source", datalist, data2),
  #fit_and_extract("06_spatialcmdStanSourceCovar.stan", "06_sp+source+covars", datalist, data2)
  
)

library(ggforce)
p <- ggplot(model_pred_all, aes(x = year_id)) +
  
  geom_ribbon(aes(ymin = q5, ymax = q95, fill = model, group = model), alpha = 0.15) +
  
  geom_line(aes(y = q50, color = model), linewidth = 0.9) +
  
  geom_point(
    aes(
      y = resis_prop,
      size = isolates#,
      #color = isolates <= 100
    ),
    shape = 21,
    fill = NA,   # hollow circles
    stroke = 0.3, # control the thickness of circle
    color = ifelse(model_pred_all$isolates <= 100, "blue", "black")
  ) +
  
  facet_wrap_paginate(~country_name, ncol = 4, nrow = 3, page = 1) +
  
  labs(
    title = "Model Comparison of 3GC",
    x = "Year",
    y = "Proportion Resistant",
    color = "Model",
    fill  = "Model"
  ) +
  
  ylim(0,1) +
  
  theme_minimal()
#other resp urine blood
num_pages <- n_pages(p)

pdf("TEST01.pdf", width = 10, height = 8)

for (i in 1:num_pages) {
  print(p + facet_wrap_paginate(~country_name, ncol = 4, nrow = 3, page = i))
}

dev.off()

# End of Rstan fit and extract  -------------------------------------------

# CMDstan fit and extract  ------------------------------------------------
setwd("C:/Users/naekaphirat/Desktop/kp/local_output")
library(cmdstanr)

# Settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4
p_chains <- 4


# Compile model #
message("Compiling Stan model...") #
mod <- cmdstan_model("04_spatial_cmdStan.stan")

# Run sampling
message("Running sampling...") #
fit <- mod$sample(
  data = datalist,
  iter_warmup = n_iter/2,
  iter_sampling = n_iter/2,
  chains = n_chains,
  parallel_chains = p_chains,
  adapt_delta = 0.95,
  max_treedepth = 15
)

# Save result
fit$save_object("fit_results.rds")

message("Job finished successfully.") #
# End of CMDstan fit and extract  -----------------------------------------

# Correlation checked -----------------------------------------------------

kp_checked <- kp %>%
  group_by(loc_id, year_id) %>%
  summarize(
    sum_resistant = sum(resistant, na.rm = TRUE),
    sum_sample = sum(sample_size, na.rm = TRUE),
    .groups = 'drop',
  )

kp_checked$res_prop <- kp_checked$sum_resistant / kp_checked$sum_sample
kp_checked$sum_resistant <- NULL
kp_checked$sum_sample <- NULL

kp_checked2 <- kp_checked %>% inner_join(covs_IHME, by = c("loc_id", "year_id"))

summary(kp_checked2) # represent type of varaibles
colSums(is.na(kp_checked2)) # NAs count in each variables  


library(corrplot)

cor_matrix <- cor(kp_checked2, use = "complete.obs")

cor_matrix_sq <- (cor_matrix)^2 # calcualte R-square

cor_matrix_plot <- cor_matrix
cor_matrix_plot[abs(cor_matrix_plot) < 0.4] <- NA
corrplot.mixed(cor_matrix_plot, upper = "color", lower = "number",
               tl.cex = 0.7, number.cex = 0.6)

print(cor_matrix)
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.7)

# End of correlation checked ----------------------------------------------

