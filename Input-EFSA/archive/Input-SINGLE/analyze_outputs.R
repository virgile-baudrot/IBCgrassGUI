library(ggplot2)
# Read RDS
long_data <- readRDS("output-SINGLE/long_data.rds")

#  Generate the faceted graphic
ggplot(long_data, aes(x = time, y = Value, color = Scenario)) +
  geom_line(linewidth = 0.2) +
  # scale_y_sqrt() +
  facet_wrap(~ OutputVariable, scales = "free_y", ncol=1) +
  theme_minimal() +
  labs(
    title = "Model Outputs Over Time Across Scenarios",
    x = "Time",
    y = "Output Value",
    color = "Scenarios"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )
ggplot(long_data, aes(x = time, y = Ratio_to_Control, color = Scenario)) +
  geom_line(linewidth = 0.2) +
  facet_wrap(~ OutputVariable, scales = "free_y", ncol=1) +
  theme_minimal() +
  labs(
    title = "Model Outputs Over Time Across Scenarios",
    x = "Time",
    y = "Output Value",
    color = "Scenarios"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )
ggplot(long_data, aes(x = time, y = Diff_from_Control, color = Scenario)) +
  geom_line(linewidth = 0.2) +
  facet_wrap(~ OutputVariable, scales = "free_y", ncol=1) +
  theme_minimal() +
  labs(
    title = "Model Outputs Over Time Across Scenarios",
    x = "Time",
    y = "Output Value",
    color = "Scenarios"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )

###

filt_data <- long_data[long_data[["OutputVariable"]] == "cover", ]
ggplot(filt_data, aes(x = time, y = Value, color = Scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Scenario) +
  theme_minimal() +
  labs(
    title = "Model Outputs Over Time Across Scenarios",
    x = "Time",
    y = "Output Value",
    color = "Scenarios"
  )
filt_data <- long_data[long_data[["OutputVariable"]] == "Inds", ]
ggplot(filt_data, aes(x = time, y = Value, color = Scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Scenario) +
  theme_minimal() +
  labs(
    title = "Model Outputs Over Time Across Scenarios",
    x = "Time",
    y = "Output Value",
    color = "Scenarios"
  )
filt_data <- long_data[long_data[["OutputVariable"]] == "repromass", ]
ggplot(filt_data, aes(x = time, y = Value, color = Scenario)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Scenario) +
  theme_minimal() +
  labs(
    title = "Model Outputs Over Time Across Scenarios",
    x = "Time",
    y = "Output Value",
    color = "Scenarios"
  )
filt_data <- long_data[long_data[["OutputVariable"]] == "seedlings", ]
ggplot(filt_data, aes(x = time, y = Value, color = Scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Scenario) +
  theme_minimal() +
  labs(
    title = "Model Outputs Over Time Across Scenarios",
    x = "Time",
    y = "Output Value",
    color = "Scenarios"
  )
filt_data <- long_data[long_data[["OutputVariable"]] == "seeds", ]
ggplot(filt_data, aes(x = time, y = Value, color = Scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Scenario) +
  theme_minimal() +
  labs(
    title = "Model Outputs Over Time Across Scenarios",
    x = "Time",
    y = "Output Value",
    color = "Scenarios"
  )
filt_data <- long_data[long_data[["OutputVariable"]] == "shootmass", ]
ggplot(filt_data, aes(x = time, y = Value, color = Scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Scenario, scale="free_y") +
  theme_minimal() +
  labs(
    title = "Model Outputs Over Time Across Scenarios",
    x = "Time",
    y = "Output Value",
    color = "Scenarios"
  )
