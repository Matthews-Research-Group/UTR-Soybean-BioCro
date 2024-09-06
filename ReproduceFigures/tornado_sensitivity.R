# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)

# Read the CSV data
data <- read_csv("../ParameterOptimization/prcc_ranking_1e+05.csv")

# Process the data
tornado_data <- data %>%
  select(parameter = "", value = original) %>%
  mutate(abs_value = abs(value)) %>%
  arrange(desc(abs_value)) %>%
  mutate(parameter = factor(parameter, levels = parameter))

# Create the tornado graph
ggplot(tornado_data, aes(x = value, y = parameter, fill = value > 0)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(value, 3), 
                x = ifelse(value > 0, 0.01, -0.01)),
            hjust = ifelse(tornado_data$value > 0, 0, 1),
            size = 3) +
  scale_fill_manual(values = c("TRUE" = "#82ca9d", "FALSE" = "#8884d8")) +
  labs(title = "Sensitivity Analysis Tornado Graph",
       x = "Impact",
       y = "Parameter") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 8)) +
  scale_x_continuous(limits = c(-0.4, 0.4)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

# Save the plot
ggsave("tornado_graph.png", width = 10, height = 8, dpi = 300)

