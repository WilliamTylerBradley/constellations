# edges look terrible

# https://mathworld.wolfram.com/GnomonicProjection.html
library(tidyverse)

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

# this one
projection <- projection %>% 
  mutate(R = 1) %>%
  mutate(lat = asin(z / R),
         long = atan2(y, x),
         central_long = 0 * pi/180,
         central_lat = if_else(z < 0, 90 * pi/180, -90 * pi/180)) %>%
  mutate(cos_c = sin(central_lat) * sin(lat) + cos(central_lat) * cos(lat) * cos(long - central_long)) %>%
  mutate(new_x = cos(lat) * sin(long - central_long) / cos_c,
         new_y = (cos(central_lat) * sin(lat) - sin(central_lat) * cos(lat) * cos(long - central_long)) / cos_c) %>%
  mutate(rho = sqrt(new_x^2 + new_y^2)) %>%
  mutate(c = atan(rho)) %>%
  mutate(lat_check = asin(cos(c) * sin(central_lat) + (new_y * sin(c) * cos(central_lat)) / rho),
         long_check = central_long + atan2(new_x * sin(c), (rho * cos(central_lat) * cos(c) - new_y * sin(central_lat) * sin(c))) %% pi)


ggplot(data = projection) +
  geom_point(aes(long, lat, color = z))

ggplot(data = projection) +
  geom_point(aes(x = new_x, y = new_y, color = lat)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))

ggplot(data = projection[projection$lat > 5*pi/180, ]) +
  geom_point(aes(x = new_x, y = new_y, color = lat)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))

ggplot(data = projection[projection$lat < -5*pi/180, ]) +
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
