library(ggplot2)

# 1. Create summary data frame from your Kentucky regression results
# Converting log-duration back to weeks using exp() for easy reading
did_data <- data.frame(
    Period = factor(c("Before Policy Change", "After Policy Change", 
                      "Before Policy Change", "After Policy Change"),
                    levels = c("Before Policy Change", "After Policy Change")),
    Group = c("Low Earners (Control)", "Low Earners (Control)",
              "High Earners (Treatment)", "High Earners (Treatment)"),
    Weeks = c(exp(1.1),                # 3.00 weeks
              exp(1.1 + 0.01),         # 3.03 weeks
              exp(1.1 + 0.26),         # 3.90 weeks
              exp(1.1 + 0.26 + 0.01 + 0.19)) # 4.76 weeks
)

# 2. Create the counterfactual line (what would have happened without the cap increase)
counterfactual <- data.frame(
    Period = factor(c("Before Policy Change", "After Policy Change"),
                    levels = c("Before Policy Change", "After Policy Change")),
    Weeks = c(exp(1.1 + 0.26),         # 3.90 weeks
              exp(1.1 + 0.26 + 0.01))  # 3.94 weeks
)

# 3. Build the LinkedIn-ready plot
ggplot(did_data, aes(x = Period, y = Weeks, group = Group, color = Group)) +
    # Counterfactual dashed line
    geom_line(data = counterfactual, aes(x = Period, y = Weeks, group = 1),
              linetype = "dashed", color = "#E63946", linewidth = 1, alpha = 0.6) +
    # Main solid trend lines
    geom_line(linewidth = 1.5) +
    geom_point(size = 4.5) +
    # Direct data labels on points
    geom_text(aes(label = paste0(round(Weeks, 2), " wks")), 
              vjust = -1.2, size = 4, fontface = "bold", show.legend = FALSE) +
    # Annotation highlighting the +21% causal jump
    annotate("segment", x = 2.04, xend = 2.04, y = 3.94, yend = 4.76,
             arrow = arrow(ends = "both", length = unit(0.2, "cm")), 
             color = "#1D3557", linewidth = 0.8) +
    annotate("text", x = 2.07, y = 4.35, label = "+21% Causal\nJump (p = 0.005)", 
             hjust = 0, size = 3.8, fontface = "bold", color = "#1D3557") +
    # Clean styling & colors
    scale_color_manual(values = c("Low Earners (Control)" = "#457B9D", 
                                  "High Earners (Treatment)" = "#E63946")) +
    scale_y_continuous(limits = c(2.5, 5.4), breaks = seq(2.5, 5.0, 0.5)) +
    expand_limits(x = c(1, 2.45)) + # Extra right margin for the annotation
    labs(
        title = "Higher Payouts Led to Longer Recovery Times",
        subtitle = "Average workers' compensation claim duration in Kentucky (1980 Policy Change)",
        x = NULL,
        y = "Average Claim Duration (Weeks)",
        color = NULL,
        caption = "Source: Meyer, Viscusi, & Durbin (1995) Kentucky Dataset | Dashed line = expected trend without policy change"
    ) +
    theme_minimal(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold", size = 16, color = "#1D3557"),
        plot.subtitle = element_text(size = 11.5, color = "#4A5568", margin = margin(b = 15)),
        plot.caption = element_text(size = 9, color = "#718096", margin = margin(t = 15)),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 11.5, color = "#2D3748"),
        plot.margin = margin(20, 25, 15, 20)
    )

# 4. Save as a high-resolution image optimized for LinkedIn
ggsave("kentucky_did_linkedin.png", width = 8, height = 5.5, dpi = 300, bg = "white")
