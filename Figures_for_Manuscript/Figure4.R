############### Figure 4a ##################


library(ggplot2)

############### Number of significant DAs from each tools #########
cold_hot_Office <-  read.table(file = "cold_hot_Office_data.txt", sep = "\t", header = TRUE)

cold_hot_Office$Tools
cold_hot_Office$Tools_Values

fig4a <- ggplot(spike10, aes(x = Tools_Values, y = Tools, color = Tools, fill = Tools)) +
  geom_dotplot(binwidth = 2, stackdir = "center", dotsize = 0.7) +  # Set binwidth
  scale_color_manual(values = c("forestgreen", "firebrick", "gold", "purple", "darkorange", 
                                "dodgerblue",  "cyan", 
                                "magenta", "black")) +
  scale_fill_manual(values = c("forestgreen", "firebrick", "gold", "purple", "darkorange", 
                               "dodgerblue",  "cyan", 
                               "magenta", "black")) +
  scale_x_continuous(limits = c(-10, 50)) +
  theme(
    panel.background = element_blank(),  # Remove background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    text = element_text(size = 14, colour = "black", face = "bold"),
    axis.line = element_line(color = "black")  # Black x and y axis lines
  ) +
  labs(
    title = "False Positive Analyses",
    x = "Percentage of Significant after multiple test corrections",
    y = "Consensus Tools"
  ) +
  facet_grid(~ Group)

fig4a

ggsave("fig4a.pdf",  fig4a = fig4a, width = 15, height = 10)


############### Figure 4b ##################


############### Number of significant DAs from each tools output of simulation data #########

spike10_simulation <-  read.table(file = "spike10_simulation.txt", sep = "\t", header = TRUE)


# Define group-specific max values and x-axis limits
group_specs <- data.frame(
  Group = unique(spike10_simulation$Group),  # Adjust to match your group names
  max_value = c(40, 60, 90, 140, 270, 500)
  #  x_limit = c(45, 65, 95, 145, 275, 505)  # Slightly higher than max for better visualization
)

unique(spike10_simulation$Group)

# Define dotted lines for each group based on their scale
group_dotted_lines <- data.frame(
  Group = rep(unique(spike10_simulation$Group), times = c(1, 1, 1, 1, 1, 1)),  # Different number of lines per group
  line_pos = c(
    # Group 1 (max 40): lines at 10, 25, 40
    25,
    # Group 2 (max 60): lines at 15, 30, 50
    50,
    # Group 3 (max 90): lines at 25, 50, 75
    75,
    # Group 4 (max 140): lines at 25, 50, 75, 125
    125,
    # Group 5 (max 270): lines at 50, 100, 150, 200, 250
    250,
    # Group 6 (max 500): lines at 100, 200, 300, 400, 500
    500
  )
)

fig4b <- ggplot(spike10_simulation, aes(x = Tools_Values, y = Tools, color = Tools, fill = Tools)) +
  # Add dotted vertical lines (different per group)
  geom_vline(data = group_dotted_lines, aes(xintercept = line_pos), 
             linetype = "dotted", color = "navyblue", size = 0.8) +
  
  # Add max value vertical lines
  # geom_vline(data = group_specs, aes(xintercept = max_value), 
  #           linetype = "dashed", color = "red", size = 1.2) +
  
  geom_boxplot() +
  
  
  scale_color_manual(values = c("forestgreen", "firebrick", "gold", "purple", "darkorange", 
                                "dodgerblue", "cyan", 
                                "magenta", "black")) +
  scale_fill_manual(values = c("forestgreen", "firebrick", "gold", "purple", "darkorange", 
                               "dodgerblue",  "cyan", 
                               "magenta", "black")) +
  
  # No global scale_x_continuous - let each facet scale independently
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(size = 14, colour = "black", face = "bold"),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "True Positive Analyses",
    x = "Significant FDR counts out of 500",
    y = "Consensus Tools"
  ) +
  # Key change: scales = "free_x" allows different x-axis scales per group
  facet_grid(~ Group, scales = "free_x")

fig4b

ggsave("fig4b.pdf",  plot = fig4b, width = 15, height = 10)



