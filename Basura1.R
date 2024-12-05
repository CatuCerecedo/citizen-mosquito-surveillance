cor_predictions <- cor_predictions %>% filter(comparison != "count_suit", year == "2022")

# Define una paleta de colores personalizada
custom_colors <- c("count_ma" = "#d53e4f", "suit_ma" = "#3288bd")

# # the mean annual values
# mean_values <- cor_predictions %>%
#   group_by(year) %>%
#   summarise(mean_value = mean(value, na.rm = TRUE))

# Crear el gráfico con ggplot2
ggplot(cor_predictions, aes(x = month, y = value, color = comparison, shape = comparison, linetype = comparison)) +
  geom_point(size = 7) +
  geom_line(size = 2) +
  geom_hline(yintercept = mean(cor_predictions$value), linetype = "dashed") +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c("count_ma" = 16, "suit_ma" = 17)) +
  scale_linetype_manual(values = c("count_ma" = 1, "suit_ma" = 10)) +
  labs(
    y = "Spearman Correlation (S)",
    x = "Month",
    color = "Comparison",      # Aquí se especifica "Comparison" para la leyenda unificada
    shape = "Comparison",      # Se especifica también para shape
    linetype = "Comparison",   # Y para linetype
    title = "2022"
  ) +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  scale_y_continuous(breaks = seq(-0.5, 0.9, 0.2)) +
  theme_classic() +
  # facet_wrap(~year) +
  theme(
    legend.position = "bottom",    # La leyenda se posiciona abajo
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(hjust = 2, size = 20),
  ) +
  guides(                       # Unificación de la leyenda
    color = guide_legend(title = "Comparison"),
    shape = guide_legend(title = "Comparison"),
    linetype = guide_legend(title = "Comparison")
  )

df <- df %>% filter(comparison != "count_suit")
ggplot(df, aes(x = females, y = value, group = comparison, color = comparison)) +
  geom_point(size = 7, alpha = 0.9) +
  geom_smooth(method = lm, aes(fill = comparison)) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "\nAverage Counts of Tiger Mosquitoes from Traps",
    y = "Spearman Spatial Correlation (S)\n",
    color = "Comparison",
    fill = "Comparison"
  ) +
  # facet_wrap(~year) +
  theme_classic() +
  theme(
    legend.position = "right",    # La leyenda se posiciona abajo
    legend.title = element_text(size = 22, face = "bold", margin = margin(b=20)),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(hjust = 2, size = 20),
    legend.key.spacing.y = unit(0.5, 'cm')
  ) 


ggsave(file = paste0(loc.fig, "plot5.png"), 
       units = "cm", height = 20, width = 35, bg = "white")


df <- df %>% filter(comparison != "count_suit")
ggplot(df, aes(x = n_reports, y = value, group = comparison, color = comparison)) +
  geom_point(size = 7, alpha = 0.9) +
  geom_smooth(method = lm, aes(fill = comparison)) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "\nNumber of Validated Reports of Mosquito Alert",
    y = "Spearman Spatial Correlation (S)\n",
    color = "Comparison",
    fill = "Comparison"
  ) +
  # theme_classic()
  theme_classic() +
  theme(
    legend.position = "right",    # La leyenda se posiciona abajo
    legend.title = element_text(size = 22, face = "bold", margin = margin(b=20)),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(hjust = 2, size = 20),
    legend.key.spacing.y = unit(0.5, 'cm')
  ) 

ggsave(file = paste0(loc.fig, "plot6.png"), 
       units = "cm", height = 20, width = 35, bg = "white")
