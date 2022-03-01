# plotly and movement along a great circle
library(tidyverse)
library(plotly)

set.seed(10)
df <- tibble(
  x = rnorm(2),
  y = rnorm(2),
  z = rnorm(2)
) %>%
  mutate(normalize = sqrt(x^2 + y^2 + z^2)) %>%
  mutate(normalize = if_else(normalize == 0, 1, normalize)) %>%
  mutate(x = x / normalize,
         y = y / normalize,
         z = z / normalize) %>%
  select(-normalize) %>%
  rowid_to_column("id")

# fig <- plot_ly(df, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines',
#                opacity = 1)
# fig

max_angle <- df %>%
  pivot_wider(id_cols = id,
              names_from = id,
              values_from = c(x, y, z)) %>%
  mutate(interior = x_1 * x_2 + 
  y_1 * y_2 + 
  z_1 * z_2) %>%
  mutate(angle = acos(interior)) %>%
  pull(angle)

df2 <- df %>%
  pivot_wider(id_cols = id,
              names_from = id,
              values_from = c(x, y, z)) %>%
  expand_grid(angle = seq(0, max_angle, length.out = 100)) %>%
  mutate(normal_x = y_1 * z_2 - z_1 * y_2,
         normal_y = z_1 * x_2 - x_1 * z_2,
         normal_z = x_1 * y_2 - y_1 * x_2) %>%
  mutate(normal_norm = sqrt(normal_x^2 + normal_y^2 + normal_z^2)) %>%
  mutate(normal_x = normal_x / normal_norm,
         normal_y = normal_y / normal_norm,
         normal_z = normal_z / normal_norm) %>%
  mutate(perp_x = normal_y * z_1 - normal_z * y_1,
         perp_y = normal_z * x_1 - normal_x * z_1,
         perp_z = normal_x * y_1 - normal_y * x_1,
         normal_norm_check = sqrt(normal_x^2 + normal_y^2 + normal_z^2)) %>%
  mutate(perp_norm = sqrt(perp_x^2 + perp_y^2 + perp_z^2)) %>%
  mutate(perp_x = perp_x / perp_norm,
         perp_y = perp_y / perp_norm,
         perp_z = perp_z / perp_norm) %>%
  mutate(vector_check = x_1 * perp_x + y_1 * perp_y + z_1 * perp_z,
         x_1 = cos(angle) * x_1 + sin(angle) * perp_x,
         y_1 = cos(angle) * y_1 + sin(angle) * perp_y,
         z_1 = cos(angle) * z_1 + sin(angle) * perp_z) %>%
  mutate(norm_check = sqrt(x_1^2 + y_1^2 + z_1^2))

df3 <- df2 %>%
  select(x_1, y_1, z_1) %>%
  rename(x = x_1,
         y = y_1,
         z = z_1) %>%
  rowid_to_column("id") %>%
  mutate(id = id + 2) %>%
  rbind(df) %>%
  mutate(color = if_else(id < 3, 0, 1))

df4 <- df2 %>%
  filter(angle < pi) %>%
  select(x_1, y_1, z_1) %>%
  rename(x = x_1,
         y = y_1,
         z = z_1) %>%
  rowid_to_column("id") %>%
  mutate(id = id + 2) %>%
  rbind(df) %>%
  mutate(color = if_else(id < 3, 0, 1))

fig <- plot_ly(df3, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines',
               opacity = 1, line = list(width = 6, color = "blue"))
fig <- layout(fig, scene = list(xaxis = list(range = c(-1,1), color = "red"), 
             yaxis = list(range = c(-1,1)), 
             zaxis = list(range = c(-1,1))))
fig

              