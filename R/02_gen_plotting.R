
rm(list = ls())
setwd("C:/Users/naekaphirat/Desktop/IHME_datasets")

library(dplyr)

sum_dat <- readRDS("C:/Users/naekaphirat/Desktop/IHME_datasets/sum_dat")

library(sf)
library(sp)

locs <- st_read("C:/Users/naekaphirat/Desktop/KP/shp/GBD2023_analysis_final_loc_set_22.shp")
locs <- locs[locs$level == 3,]
locs2 <- data.frame(locs$loc_id, locs$region_id, locs$spr_reg_id, locs$loc_name)
names(locs2)[1] <- "location_id"
names(locs2)[2] <- "region_id"
names(locs2)[3] <- "spr_reg_id"
names(locs2)[4] <- "country"
data2 <- merge(sum_dat,locs2)

library(ggplot2)
library(ggforce)

int_abx <- c("3gc", "carbapenem", "aminoglycoside", "fluoroquinolone")

for (abx in int_abx) {

data3 <- data2 %>%
  filter(res > 0, sum > 10, abx_class == abx) %>% # Keep abx_class
  mutate(prop = r / sum)              # Keep only rows contain resistant value
                                     
if (nrow(data3) == 0) next            # Skip empty cell

# Your original plot, modified slightly
p <- #ggplot(data3, aes(x = year_id, y = prop, color = hosp_group, size = r)) +
     ggplot(data3, aes(x = year_id, y = prop, color = specimen_group, size = res)) +
  geom_point(alpha = 0.7) +
  facet_wrap_paginate(~country, ncol = 4, nrow = 3, page = 1) +  # key line changed
  scale_size_continuous(range = c(1, 6), guide = guide_legend(title = "cases")) +
  labs(
    title = "Observed Resistance by Country and Source",
    x = "Year",
    y = "Proportion of Resistant"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Check how many pages are needed
num_pages <- n_pages(p)

pdf_filename <- paste0(abx, ".pdf")
pdf(pdf_filename, width = 10, height = 8)


for (i in 1:num_pages) {
  print(
    p + facet_wrap_paginate(~country, ncol = 4, nrow = 3, page = i)
  )
}
dev.off()

}
