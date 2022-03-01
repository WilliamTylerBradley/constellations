# https://en.wikipedia.org/wiki/Stereographic_map_projection
size <- 1000
projection <- tibble(
  theta = 2 * pi * seq(0, size - 1) / ((1 + sqrt(5)) / 2),
  phi = acos(1 - 2 * (seq(0, size - 1) + .5) / size)) %>%
  mutate(x = cos(theta) * sin(phi),
         y = sin(theta) * sin(phi),
         z = cos(phi)) %>%
  rowid_to_column("id")

ggplot(data = projection, 
       aes(x, y)) +
  geom_point()

projection <- projection %>% 
  mutate(R = 1) %>%
  mutate(lat = asin(z / R),
         long = atan2(y, x)) %>%
  mutate(new_x = if_else(z < 0, x / (1 - z), x / (1 + z)),
         new_y = if_else(z < 0, y / (1 - z), y / (1 + z)))
         
ggplot(data = projection) +
  geom_point(aes(long, lat, color = z))

ggplot(data = projection) +
  geom_point(aes(x = new_x, y = new_y, color = lat)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))

ggplot(data = projection[projection$z < 0, ]) +
  geom_point(aes(x = new_x, y = new_y, color = lat)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))
ggplot(data = projection[projection$z >= 0, ]) +
  geom_point(aes(x = new_x, y = new_y, color = lat)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))

projection <- projection %>% 
  mutate(R = 1) %>%
  mutate(lat = asin(z / R),
         long = atan2(y, x)) %>%
  mutate(new_x = if_else(z < 0, x / (1 - z), x / (1 + z) + 2),
         new_y = if_else(z < 0, y / (1 - z), y / (1 + z)))

ggplot(data = projection) +
  geom_point(aes(x = new_x, y = new_y, color = lat)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))

projection <- projection %>%
  mutate(z_check = if_else(z < 0, 0, 1))
ggplot(data = projection) +
  geom_point(aes(x = new_x, y = new_y, color = z_check)) +
  coord_equal()
