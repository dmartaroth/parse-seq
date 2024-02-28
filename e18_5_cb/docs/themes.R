# ## ################################# ## #
#                  THEMES                 #
# ## ################################# ## #
# Custom colors
my_colors <- c("#FFB6C1", "#ADD8E6", "#FFD700", "#98FB98", "#FFA07A")
error <- red $ bold
note1 <- red $ bold
note2 <- black $ bold
note3 <- blue $ bold

# Custom theme for barplot
custom_theme_bar <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(color = "white", fill = "white"),
      panel.grid.major = element_blank(),           # No major gridlines
      panel.grid.minor = element_blank(),           # No minor gridlines
      legend.position = "none",                     # No legend
      axis.ticks = element_line(color = "black"),
      axis.ticks.x = element_line(color = NA),
      axis.title = element_text(color = "black",size = 8),
      axis.text = element_text(size = 8),           # Size of axis text
      axis.line.y = element_line(colour = NA), # Set y-axis line color
      axis.line.x = element_line(colour = NA),  
      axis.text.y = element_text(colour = "black"), # Set y-axis text color
      axis.ticks.length = unit(0.2, "cm"),          # Shorten tick length
      panel.border = element_blank(),   
      panel.grid.major.y = element_line(colour = "gray", linetype = 3),  # Dashed horizontal gridlines
      panel.grid.major.x = element_blank(),         # No vertical gridlines
      text = element_text(colour = "black", size = 8),                        # Regular text
      plot.title = element_text(face = "bold", size = 8,hjust = 0),     # Bold plot title
      # plot.margin = margin(20, 20, 20, 20),         # Adjust plot margins
      plot.caption = element_text(color = "black")  # Caption color and position
    )
}

# Custom theme for density plot
custom_theme_density <- function() {
  theme_minimal() +
    theme(
      # Background and Gridlines
      plot.background = element_rect(fill = "white", color = NA),  # White background
      panel.background = element_rect(fill = "white"),            # Panel background
      panel.grid.major = element_blank(),                          # No major gridlines
      panel.grid.minor = element_blank(),                          # No minor gridlines
      panel.grid.major.y = element_line(colour = "gray", linetype = 3),  # Dashed horizontal gridlines
      panel.grid.major.x = element_blank(),                        # No vertical gridlines
      
      # Axis
      axis.ticks = element_line(color = "black"),                  # Color of axis ticks
      axis.line.y = element_line(colour = "black"),                # Color of y-axis line
      axis.line.x = element_line(colour = NA),                     # Remove x-axis line
      axis.title = element_text(color = "black"),                  # Color of axis titles
      axis.text = element_text(size = 10),                         # Size of axis text
      axis.text.y = element_text(colour = "black"),                # Color of y-axis text
      axis.ticks.length = unit(0.2, "cm"),                         # Shorten tick length
      
      # Panel
      panel.border = element_rect(color = NA, fill = NA),     # Panel border
      
      # Text
      plot.title = element_text(face = "bold", vjust = 0.5),       # Bold plot title, centered vertically
      text = element_text(),                                       # Regular text
      
      # Legend and Caption
      legend.position = "right",                                    # No legend
      plot.caption = element_text(color = "black", hjust = 0)      # Caption color and position
    )
}

# Custom theme for scatter plots
custom_theme_scatter <- function() {
  theme_minimal() +
    theme(
      # Background and Gridlines
      plot.background = element_rect(fill = "white", color = NA),  # White background
      panel.background = element_rect(fill = "white"),            # Panel background
      panel.grid.major = element_blank(),                          # No major gridlines
      panel.grid.minor = element_blank(),                          # No minor gridlines
      
      # Axis
      axis.line = element_line(color = "black"),                   # Color of axis lines
      axis.title = element_text(color = "black"),                  # Color of axis titles
      axis.text = element_text(color = "black", size = 8),         # Size and color of axis text
      
      # Legend and Caption
      legend.position = "none",                                    # No legend
      plot.caption = element_text(color = "black", hjust = 0.5),    # Caption color and position
      
      # Remove axis ticks
      axis.ticks = element_blank()                                 # No axis ticks
    )
}

# Custom violin plot theme
custom_theme_violin <- function() {
  theme_minimal() +
    theme(
      legend.position = "none",  # Remove legend
      panel.background = element_rect(fill = "white"),  # White background for panels
      panel.grid.major.y = element_line(colour = "grey98"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),  # Color of axis lines
      axis.title = element_text(color = "black"),  # Color of axis titles
      axis.text = element_text(color = "black"),  # Color of axis text
      axis.ticks = element_line(color = "black"),  # Color of axis ticks
      text = element_text(color = "black"),  # Color of text
      plot.title = element_text(face = "bold", hjust = 0.5),  # Bold plot title, centered
      plot.caption = element_text(hjust = 0)  # Caption position
    )
}