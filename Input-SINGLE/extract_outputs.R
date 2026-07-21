# Install ggplot2 and tidyr if not already installed
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(tidyr)) install.packages("tidyr")
if (!require(dplyr)) install.packages("dplyr")

library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Map scenarios to their generated IDs and types
# This matches the loop order generated in run_scenarios.R

scenario_names <- c("control", "biomass", "SEbiomass", "survival", "establishment", "sterility", "seednumber")
scenarios <- data.frame(
  name = scenario_names,
  id = 1:7,
  type = c("control", rep("treatment", 6)), # 1 is control, 2-7 are treatments
  stringsAsFactors = FALSE
)

# Initialize empty dataframe
all_herbfact <- data.frame()
for (scenario in scenario_names) {
  file_name <- paste0("HerbFact_", scenario, ".txt")
  if (file.exists(file_name)) {
    # Read the text file
    temp_data <- read.table(file_name, header = TRUE, sep = "\t")
    # Add a column for the scenario name
    temp_data$Scenario <- scenario 
    # Optional: Add a 'Week' column since HerbFact applies per week
    temp_data$Week <- 1:nrow(temp_data) 
    # Bind to the master table
    all_herbfact <- rbind(all_herbfact, temp_data)
  } else {
    warning(paste("File not found:", file_name))
  }
}
# Save the combined table
write.csv(all_herbfact, "../output-SINGLE/Combined_HerbFact_All_Scenarios.csv", row.names = FALSE)

all_fieldedge <- data.frame()
for (scenario in scenario_names) {
  file_name <- paste0("Fieldedge_", scenario, ".txt")
  if (file.exists(file_name)) {
    # Read the text file
    temp_data <- read.table(file_name, header = TRUE, sep = "\t")
    # Add a column for the scenario name
    temp_data$Scenario <- scenario 
    # Bind to the master table
    all_fieldedge <- rbind(all_fieldedge, temp_data)
  } else {
    warning(paste("File not found:", file_name))
  }
}
# Save the combined table
write.csv(all_fieldedge, "../output-SINGLE/Combined_Fieldedge_All_Scenarios.csv", row.names = FALSE)

all_data <- data.frame()
# 2. Read all scenario outputs into a single dataframe
for (i in 1:nrow(scenarios)) {
  scenario_name <- scenarios$name[i]
  scenario_id <- scenarios$id[i]
  scenario_type <- scenarios$type[i]
  # Construct the exact filename based on the Model-files directory listing
  file_name <- paste0("../Model-files/Pt__type_", scenario_type, "_Scenario_", scenario_id, "_MCrun_1.txt")
  if (file.exists(file_name)) {
    # The output is a tab-separated text file, so we use read.table
    temp_data <- read.table(file_name, header = TRUE, sep = "\t")
    temp_data$Scenario <- scenario_name # Add a column to identify the scenario
    # Standardize 'Time' column name to lowercase 'time' for consistency in ggplot
    if("Time" %in% colnames(temp_data)) names(temp_data)[names(temp_data) == "Time"] <- "time"
    all_data <- rbind(all_data, temp_data)
  } else {
    warning(paste("File not found:", file_name))
  }
}
# 3. Reshape the data from wide to long format for ggplot faceting
long_data <- all_data %>%
  pivot_longer(
    cols = -any_of(c("time", "Scenario", "PFT")), # Added "PFT" here so it isn't pivoted into the 'Value' column
    names_to = "OutputVariable",
    values_to = "Value"
  )
# 4. Calcul des métriques relatives au contrôle
# a. Isoler les valeurs du contrôle
control_data <- long_data %>%
  filter(Scenario == "control") %>%
  select(time, PFT, OutputVariable, Control_Value = Value)

# b. Joindre avec le data long et calculer les nouvelles colonnes
long_data <- long_data %>%
  left_join(control_data, by = c("time", "PFT", "OutputVariable")) %>%
  mutate(
    Ratio_to_Control = Value / Control_Value,              # Division par le contrôle
    Diff_from_Control = Control_Value - Value              # Contrôle moins scénario
  )

# Save RDS
saveRDS(long_data, "../output-SINGLE/long_data.rds")
cat("Clean output file saved as long_data.rds\n")
