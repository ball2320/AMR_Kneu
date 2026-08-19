# Outlier detection ####
rm(list = ls())
setwd("C:/Users/naekaphirat/Desktop/KP")

library(lme4)
library(foreign)
library(ggplot2)
library(ggforce)
library(splines)
library(data.table)
# Setup ####

library(openxlsx)
mydata <- read.xlsx("KP_overall.xlsx", sheet = "analysis2")

mydata <- mydata[mydata$sample_size>=5,]
mydata <- mydata[mydata$year_id <2020,]

library(dplyr)
mydata <- mydata %>%
  #filter(abx_class == "3gc") %>%
  select(ihmelocid, locid, year_id, source, abx_name, abx_class, susceptible, resistant, sample_size)
#Rename column name
mydata <- mydata %>%
  rename(location_id = locid)

#Select a single study
#mydata <- mydata %>%
#  filter(source == "ecdc_ears")

#Screening out some countries that contain incomplete information
mydata <- mydata %>%
  filter(!ihmelocid %in% c("Slovak Republic", "Hong Kong", "Kosovo"))

#Create the resistance proportion (val)
mydata <- mydata %>%
  mutate(val = resistant / sample_size)

covs <- read.csv("Gram_covs_1990to2021.csv")
covs <- covs[covs$year_id <2020,]
covs <- covs[covs$national == 1,]
#These covs are selected by LASSO
#covs <- covs[names(covs) %in% c('tb_strains_transmission_risk',  
#                                'cce', 'he_cap', 'year_id',
#                                'smok_prev_agestd_both',
#                                'pollution_outdoor_pm25',
#                                'location_id')]

covs <- covs[names(covs) %in% c('year_id', 'pigs_pc', 'age_std_hiv_prev',
                                'sev_scalar_agestd_hiv', 'physicians_pc',
                                'smok_daily_prev_agstd_both', 'location_id')]

covs[,3:7] <- scale(covs[,3:7])

library(sf)
library(sp)
locs <- st_read("C:/Users/naekaphirat/Desktop/KP/shp/GBD2023_analysis_final_loc_set_22.shp")
locs <- locs[locs$level == 3,]
locs2 <- data.frame(locs$loc_id, locs$region_id, locs$spr_reg_id, locs$loc_name)
names(locs2)[1] <- "location_id"
names(locs2)[2] <- "region_id"
names(locs2)[3] <- "spr_reg_id"
names(locs2)[4] <- "country"
mydata <- merge(mydata,locs2)



#merge data and covs
mydata <- merge(mydata, covs, by = c('location_id', 'year_id'), all.x = T, all.y = F)


#or binomial model:
response = cbind(successes = round(mydata$val*mydata$sample_size,0),
                 failures = mydata$sample_size - round(mydata$val*mydata$sample_size,0))


model1 <- glmer(response ~ 1 + 
                  year_id + pigs_pc + age_std_hiv_prev +
                  sev_scalar_agestd_hiv + physicians_pc +
                  smok_daily_prev_agstd_both +
                  (year_id|country) +
                  (1 |spr_reg_id/region_id/country), data = mydata, family = 'binomial')

#Minimized model
#model1 <- glmer(response ~ 1 + 
#                  year_id +
#                  (year_id|country) +
#                  (1 |spr_reg_id/region_id/country), data = mydata, family = 'binomial')

summary(model1)

#Using glmmTNB instead

library(glmmTMB)

#model11 <- glmmTMB(response ~ 1 + 
#                  tb_strains_transmission_risk +
#                  cce + he_cap + year_id + 
#                  smok_prev_agestd_both +
#                  pollution_outdoor_pm25 +
#                  (year_id|country) +
#                  (1 |spr_reg_id/region_id/country), data = mydata, family = binomial(),
#                  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

#summary(model11)

covs <- merge(covs, locs, by.x = 'location_id', by.y = 'loc_id')
colnames(covs)[colnames(covs) == 'loc_name'] <- 'country'
colnames(covs)[colnames(covs) == 'region_id'] <- 'region'
colnames(covs)[colnames(covs) == 'spr_reg_id'] <- 'super_region'
covs$region <- as.character(covs$region)
covs$super_region <- as.character(covs$super_region)

covs$region[covs$region == 5] <- 'East Asia '
covs$region[covs$region == 9] <- 'Southeast Asia'
covs$region[covs$region == 65] <- 'High-income Asia Pacific'
covs$region[covs$region == 138] <- 'North Africa & Middle East'
covs$region[covs$region == 159] <- 'South Asia'
covs$region[covs$region == 167] <- 'Central Sub-Saharan Africa'
covs$region[covs$region == 174] <- 'Eastern Sub-Saharan Africa'
covs$region[covs$region == 192] <- 'Southern Sub-Saharan Africa'
covs$region[covs$region == 199] <- 'Western Sub-Saharan Africa'

covs$super_region[covs$super_region == 64] <- 'High Income'
covs$super_region[covs$super_region == 137] <- 'North Africa & Middle East'
covs$super_region[covs$super_region == 158] <- 'South Asia'
covs$super_region[covs$super_region == 4] <- 'Southeast Asia, East Asia & Oceania'
covs$super_region[covs$super_region == 166] <- 'Sub-Saharan Africa'

# covs <- covs[names(covs) %in% names(mydata)]
colnames(covs)[colnames(covs) == 'region'] <- 'region_id'
colnames(covs)[colnames(covs) == 'super_region'] <- 'spr_reg_id'
covs$country <- as.character(covs$country) 
covs$pred <- predict(model1, newdata = covs, type = 'response', allow.new.levels = TRUE)
summary(covs$pred)


covs <- merge(covs, mydata[,c('location_id', 'year_id', 'val', 'sample_size')], by = c('location_id', 'year_id'), all.x = T, all.y = T)


covs$sample_size_bins <- NA
covs$sample_size_bins[covs$sample_size<50] <- '10-49'
covs$sample_size_bins[covs$sample_size>=50 &covs$sample_size<100 ] <- '50-99'
covs$sample_size_bins[covs$sample_size>=100 &covs$sample_size<500 ] <- '100-499'
covs$sample_size_bins[covs$sample_size>=500 ] <- '500+'
covs$sample_size_bins <-  as.factor(covs$sample_size_bins)
covs$sample_size_bins <- factor(covs$sample_size_bins, levels = c("10-49", "50-99", "100-499", "500+"))


# Plot 
pdf("outliers_glmer2.pdf",
    height = 17,
    width = 12)

for(i in 1:(ceiling(length(unique(covs$spr_reg_id))))){
  print(ggplot(covs[covs$spr_reg_id==unique(covs$spr_reg_id)[i],])+
          geom_line(aes(x = year_id, y = pred))+
          geom_point(aes(x = year_id, y = val))+
          ylim(0,1)+
          facet_wrap_paginate(~country, page = i)
  )}

dev.off()

saveRDS(model1, 'outlier_glmer2')

# Define outliers ####


#Points are considered outlier if they are GLMM pred +/- n*MAD ####
# Calculate the MAD for each country and apply this as this differs greatly by country
# MAD is the median absolute deviation and in this case is preferable to MSE due to highly heterogenous data
# Try various values of n
# Want to outlier ~10% of the data

# Need to remove geometry out of dataset before executing
covs <- covs[, !names(covs) %in% "geometry"]

#Save covs, model1 and mydata
#saveRDS(covs, "covs.rds")
#saveRDS(model1, "model1.rds")
#saveRDS(mydata, "mydata.rds")

#Reload covs, model1 and mydata
#covs <- readRDS("covs.rds")
#model1 <- readRDS("model1.rds")
#mydata <- readRDS("mydata.rds")

library(stats)

covs <- data.table(covs)

MADs <-  covs[,.(upper_bound = pred + 2*mad(pred[!is.na(val)], val[!is.na(val)]),
                 lower_bound = pred - 2*mad(pred[!is.na(val)],val[!is.na(val)])),
              by = c('country')]

MADs <- MADs[,2:3]
covs <- cbind(covs, MADs)
covs$upper_bound[covs$upper_bound>1] <- 1
covs$lower_bound[covs$lower_bound<0] <- 0

covs <- covs[!is.na(covs$spr_reg_id),]
pdf('outlier_boundary2.pdf',
    height = 17,
    width = 12) 

for(i in 1:(ceiling(length(unique(covs$spr_reg_id))))){ 
  print(ggplot(covs[covs$spr_reg_id==unique(covs$spr_reg_id)[i],])+ 
          geom_line(aes(x = year_id, y = pred))+
          geom_ribbon(aes(x = year_id, ymin = lower_bound, ymax = upper_bound, colour = 'red', fill = 'red', alpha  =0.5))+
          geom_point(aes(x = year_id, y = val))+
          ylim(0,1)+
          facet_wrap_paginate(~country, page = i)
  )} 

dev.off()

#define outliers in the dataset
MADs <- covs[,.(country, year_id, lower_bound, upper_bound)]
MADs <-  unique(MADs)

mydata <- merge(mydata, MADs, by = c('country', 'year_id'))

mydata$is_outlier[mydata$val < mydata$lower_bound] <- 1
mydata$is_outlier[mydata$val > mydata$upper_bound] <- 2

mydata$is_outlier[mydata$val<mydata$lower_bound |mydata$val>mydata$upper_bound] <- 1
mydata$is_outlier[mydata$val>mydata$lower_bound & mydata$val<mydata$upper_bound] <- 0
outliers <- mydata[mydata$is_outlier == 1,]
#outliers <- unique(outliers[c('nid', 'year_id')])
write.csv(outliers, 'outlier2.csv', row.names = F)
write.csv(mydata, 'outlier_data2.csv')
pdf('outlier_identified2.pdf',
    height = 8.3, width = 11.7)

for(i in 1:length(unique(mydata$region_id))){
  subset <- mydata[mydata$region_id == unique(mydata$region_id)[i],]
  print(
    ggplot(subset, aes(x = year_id, y = val, color=factor(is_outlier))) +
      geom_point()+
      ylim(0,1)+
      facet_wrap(~country,nrow = ceiling(sqrt(length(unique(subset$location_id)))))
    
  )
}

dev.off()

library(ggrepel)
pdf('MAD_outliers2.pdf',
    height = 6, width = 15)
ggplot(mydata, aes(x = year_id, y = val)) +
  geom_point()+
  geom_label_repel(aes(label = is_outlier), size = 3,max.overlaps = 1000)+
  geom_point(color = ifelse(mydata$is_outlier == "", "grey50", "red"))

dev.off()
