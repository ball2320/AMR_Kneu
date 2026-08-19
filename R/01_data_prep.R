# --- Step 1: Set your working directory or folder path ---
# Replace this with your folder path

rm(list = ls())
setwd("C:/Users/naekaphirat/Desktop/IHME_datasets")

# --- Step 2: List all CSV files in the folder ---
csv_files <- list.files(pattern = "\\.csv$", full.names = TRUE)

# --- Step 3: Read and combine all CSV files ---
# Load the necessary package
library(dplyr)

# Use lapply to read and tag each file, then bind them together
kp <- lapply(csv_files, function(file) {
  # Read each CSV
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # Keep only rows contain K.pneumoniae
  df <- df %>% 
    filter(grepl("klebsiella_pneumoniae", pathogen, ignore.case = TRUE))
  
  # Add a 'source' column using the filename (without the full path)
  df$source <- tools::file_path_sans_ext(basename(file))
  
  return(df)
}) %>%
  bind_rows()   # Combine all dataframes

# --- Step 4: View or save the combined dataframe ---
# View the first few rows
kp <- kp %>%
  #filter(pathogen == "klebsiella_pneumoniae") %>%
  select(nid, location_id, year_id, hosp, raw_antibiotic, abx_class, pathogen, resistance, specimen, source)

# Assign "No Data" to No data available cell
kp <- kp %>%
  mutate(hosp = ifelse(is.na(hosp), "ND", hosp))





# Transform many cases row into single case row ---------------------------
# in waiting list (indonesia, Relavra, DEU, Janis(japan),)

rm(list = ls())
setwd("C:/Users/naekaphirat/Desktop/IHME_datasets/waiting list")
library(dplyr)
library(tidyr)

# --- Step 1: Read the original CSV file ---
# Replace 'input.csv' with your actual file name
df <- read.csv("ECDC.csv", stringsAsFactors = FALSE)

#To ensure "cases" is numeric
df <- df %>%
  mutate(cases = round(as.numeric(cases)))

# --- Step 2: Expand rows based on the 'cases' column ---
expanded_df <- df %>%
  uncount(weights = cases) %>%   # repeats each row 'cases' times
  mutate(cases = 1)              # reset cases to 1

# --- Step 3: Write the expanded dataset to a new CSV file ---
write.csv(expanded_df, "sECDC.csv", row.names = FALSE)



# Bar plot of top five most frequent value --------------------------------

library(ggplot2)
library(scales)

# Example: we’ll use the column named 'species'
top5 <- kp %>%
  count(specimen, sort = TRUE) %>%       # Count how many times each species appears
  filter(!is.na(specimen)) %>%          # Remove NA from sorting
  slice_max(n, n = 5)                  # Keep only the top 5

# View top 5 values
print(top5)

# --- Plot ---
mid_y <- max(top5$n) / 2

pdf("top5_specimens.pdf", width = 10, height = 8)

ggplot(top5, aes(x = reorder(specimen, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(y = mid_y, aes(label = scales::comma(n)),  # show actual numbers on bars
            hjust = 0.5,
            vjust = 0.5,
            size = 5) +
  scale_y_continuous(
    labels = label_number(scale_cut = cut_short_scale())) +
  coord_flip() +  # Flip for better readability
  labs(
    title = "Top 5 Most Frequent Specimen",
    x = "Specimen",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
  plot.title = element_text(face = "bold", hjust = 0.5),
  axis.title.y = element_text(face = "bold"),
  axis.title.x = element_text(face = "bold"))

dev.off()

# Make it to desired format -----------------------------------------------


library(dplyr)
library(readr)
library(stringr)

setwd("C:/Users/naekaphirat/Desktop/IHME_datasets/waiting list")
# Read your dataset (replace with actual path)
#data <- read_csv("CHAIN.csv")

# Check data structure
# head(data)

#Grouping specimen types
kp <- kp %>%
  mutate(
    specimen_group = case_when(
      str_detect(tolower(specimen), "urine|urin|urinary") ~ "urine",
      str_detect(tolower(specimen), "blood") ~ "blood",
      str_detect(tolower(specimen), "resp|respiratory!") ~ "resp",
      TRUE ~ "other"
    )
  )

#Grouping health facilities
kp <- kp %>%
  mutate(
    hosp_group = case_when(
      str_detect(tolower(hosp), "community") ~ "community",
      str_detect(tolower(hosp), "hospital") ~ "hospital",
      TRUE ~ "unknown"
    )
  )

#Rename abx_class
kp <- kp %>%
  mutate(abx_class = case_when(
    abx_class == "first_gen_ceph" ~ "1gc",
    abx_class == "second_gen_ceph" ~ "2gc",
    abx_class == "third_gen_ceph" ~ "3gc",
    abx_class == "fourth_gen_ceph" ~ "4gc",
    abx_class == "fifth_gen_ceph" ~ "5gc",
    abx_class == "aminopenicilli" ~ "aminopenicillin",
    abx_class == "anti_pseudomonal_penicillin" ~ "AntiPseudoPen",
    abx_class == "beta_lactamase_inhibitor" ~ "bli",
    abx_class == "third_gen_ceph_beta_lactamase_inhibitor" ~ "3gc_bli",
    TRUE ~ abx_class  # keep all others unchanged
  ))


# Summarize counts by desired grouping
summary_data <- kp %>%
  group_by(location_id, year_id, source, hosp_group, specimen_group, abx_class) %>%
  summarise(
    number_of_resistant = sum(resistance == "resistant", na.rm = TRUE),
    number_of_susceptible = sum(resistance == "susceptible", na.rm = TRUE),
    unknown = sum(resistance == "unknown" | is.na(resistance)),
    total = n(),
    .groups = "drop"
  )

# Rename some column names
summary_data <- summary_data %>%
  rename(
    res = number_of_resistant,
    sus = number_of_susceptible,
    unk = unknown,
    sum = total
  )

# View the summarized data
print(summary_data)

# Optional: write to CSV
write_csv(summary_data, "aggregated_resistance_summary.csv")

setwd("C:/Users/naekaphirat/Desktop/IHME_datasets")
saveRDS(summary_data, 'sum_dat')

