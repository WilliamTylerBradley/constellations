#https://en.wikipedia.org/wiki/Nicolosi_globular_projection
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

long_0_new <- pi/2
projection <- projection %>% 
  mutate(R = 1,
         long_0 = long_0_new) %>%
  mutate(lat = asin(z / R),
         long = atan2(y, x)) %>%
  mutate(b = (pi / (2 * (long - long_0))) - ((2 * (long - long_0)) / pi),
         c = (2 * lat) / pi) %>%
  mutate(d = (1 - c^2) / (sin(lat) - c)) %>%
  mutate(M = (((b * sin(lat)) / d) - (b / 2)) / (1 + (b^2/d^2)),
         N = (((d^2 * sin(lat)) / b^2) + (d / 2)) / (1 + (d^2 / b^2))) %>%
  mutate(new_x = (pi / 2) * R * (M + sign(long - long_0) * sqrt(M^2 + (cos(lat)^2 / (1 + (b^2 / d^2)) ))),
         new_y = (pi / 2) * R * (N + -1 * sign(lat) * sqrt(N^2 - ((d^2 / b^2) * sin(lat)^2 + d * sin(lat) - 1) / (1 + d^2 / b^2)))) %>%
  mutate(new_x = if_else(long - long_0 == 0, 0, new_x),
         new_y = if_else(long - long_0 == 0, R * lat, new_y)) %>%
  mutate(new_x = if_else(lat == 0, R * (long - long_0), new_x),
         new_y = if_else(lat == 0, 0, new_y)) %>%
  mutate(new_x = if_else(abs(long - long_0) == (pi / 2), R * (long - long_0) * cos(lat), new_x),
         new_y = if_else(abs(long - long_0) == (pi / 2), pi / 2 * R * sin(lat), new_y)) %>%
  mutate(new_x = if_else(abs(lat) == (pi / 2), 0, new_x),
         new_y = if_else(abs(lat) == (pi / 2), R * lat, new_y)) %>%
  mutate(graph = ifelse(long - long_0 < (pi / 2) &
                          long - long_0 > -(pi / 2), 1, 0))

ggplot(data = projection) +
  geom_point(aes(long, lat, color = long))

ggplot(data = projection) +
  geom_point(aes(x = new_x, y = new_y, color = long)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))
ggplot(data = projection[projection$graph == 1, ]) +
  geom_point(aes(x = new_x, y = new_y, color = long)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))
ggplot(data = projection[projection$graph == 0, ]) +
  geom_point(aes(x = new_x, y = new_y, color = long)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))

# Flip to other side
long_0_new <- -pi/2
projection <- projection %>% 
  mutate(R = 1,
         long_0 = long_0_new) %>%
  mutate(lat = asin(z / R),
         long = atan2(y, x)) %>%
  mutate(b = (pi / (2 * (long - long_0))) - ((2 * (long - long_0)) / pi),
         c = (2 * lat) / pi) %>%
  mutate(d = (1 - c^2) / (sin(lat) - c)) %>%
  mutate(M = (((b * sin(lat)) / d) - (b / 2)) / (1 + (b^2/d^2)),
         N = (((d^2 * sin(lat)) / b^2) + (d / 2)) / (1 + (d^2 / b^2))) %>%
  mutate(new_x = (pi / 2) * R * (M + sign(long - long_0) * sqrt(M^2 + (cos(lat)^2 / (1 + (b^2 / d^2)) ))),
         new_y = (pi / 2) * R * (N + -1 * sign(lat) * sqrt(N^2 - ((d^2 / b^2) * sin(lat)^2 + d * sin(lat) - 1) / (1 + d^2 / b^2)))) %>%
  mutate(new_x = if_else(long - long_0 == 0, 0, new_x),
         new_y = if_else(long - long_0 == 0, R * lat, new_y)) %>%
  mutate(new_x = if_else(lat == 0, R * (long - long_0), new_x),
         new_y = if_else(lat == 0, 0, new_y)) %>%
  mutate(new_x = if_else(abs(long - long_0) == (pi / 2), R * (long - long_0) * cos(lat), new_x),
         new_y = if_else(abs(long - long_0) == (pi / 2), pi / 2 * R * sin(lat), new_y)) %>%
  mutate(new_x = if_else(abs(lat) == (pi / 2), 0, new_x),
         new_y = if_else(abs(lat) == (pi / 2), R * lat, new_y)) %>%
  mutate(graph = ifelse(long - long_0 < (pi / 2) &
                          long - long_0 > -(pi / 2), 1, 0))

ggplot(data = projection) +
  geom_point(aes(long, lat, color = long))

ggplot(data = projection) +
  geom_point(aes(x = new_x, y = new_y, color = long)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))
ggplot(data = projection[projection$graph == 1, ]) +
  geom_point(aes(x = new_x, y = new_y, color = long)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))
ggplot(data = projection[projection$graph == 0, ]) +
  geom_point(aes(x = new_x, y = new_y, color = long)) +
  coord_equal() +
  scale_color_gradient(limits = c(-pi, pi))

## reset here ##

projection <- projection %>% 
  mutate(R = 1,
         long_0 = -(pi / 2)) %>%
  mutate(lat = asin(z / R),
         long = atan2(y, x)) %>%
  mutate(b = (pi / (2 * (long - long_0))) - ((2 * (long - long_0)) / pi),
         c = (2 * lat) / pi) %>%
  mutate(d = (1 - c^2) / (sin(lat) - c)) %>%
  mutate(M = (((b * sin(lat)) / d) - (b / 2)) / (1 + (b^2/d^2)),
         N = (((d^2 * sin(lat)) / b^2) + (d / 2)) / (1 + (d^2 / b^2))) %>%
  mutate(new_x = (pi / 2) * R * (M + sign(long - long_0) * sqrt(M^2 + (cos(lat)^2 / (1 + (b^2 / d^2)) ))),
         new_y = (pi / 2) * R * (N + -1 * sign(lat) * sqrt(N^2 - ((d^2 / b^2) * sin(lat)^2 + d * sin(lat) - 1) / (1 + d^2 / b^2)))) %>%
  mutate(new_x = if_else(long - long_0 == 0, 0, new_x),
         new_y = if_else(long - long_0 == 0, R * lat, new_y)) %>%
  mutate(new_x = if_else(lat == 0, R * (long - long_0), new_x),
         new_y = if_else(lat == 0, 0, new_y)) %>%
  mutate(new_x = if_else(abs(long - long_0) == (pi / 2), R * (long - long_0) * cos(lat), new_x),
         new_y = if_else(abs(long - long_0) == (pi / 2), pi / 2 * R * sin(lat), new_y)) %>%
  mutate(new_x = if_else(abs(lat) == (pi / 2), 0, new_x),
         new_y = if_else(abs(lat) == (pi / 2), R * lat, new_y)) %>%
  mutate(graph = ifelse(long - long_0 < (pi / 2) &
                          long - long_0 > -(pi / 2), 1, 2),
         save_new_x = new_x,
         save_new_y = new_y) %>% 
  mutate(R = 1,
         long_0 = (pi / 2)) %>%
  mutate(lat = asin(z / R),
         long = atan2(y, x)) %>%
  mutate(b = (pi / (2 * (long - long_0))) - ((2 * (long - long_0)) / pi),
         c = (2 * lat) / pi) %>%
  mutate(d = (1 - c^2) / (sin(lat) - c)) %>%
  mutate(M = (((b * sin(lat)) / d) - (b / 2)) / (1 + (b^2/d^2)),
         N = (((d^2 * sin(lat)) / b^2) + (d / 2)) / (1 + (d^2 / b^2))) %>%
  mutate(new_x = (pi / 2) * R * (M + sign(long - long_0) * sqrt(M^2 + (cos(lat)^2 / (1 + (b^2 / d^2)) ))),
         new_y = (pi / 2) * R * (N + -1 * sign(lat) * sqrt(N^2 - ((d^2 / b^2) * sin(lat)^2 + d * sin(lat) - 1) / (1 + d^2 / b^2)))) %>%
  mutate(new_x = if_else(long - long_0 == 0, 0, new_x),
         new_y = if_else(long - long_0 == 0, R * lat, new_y)) %>%
  mutate(new_x = if_else(lat == 0, R * (long - long_0), new_x),
         new_y = if_else(lat == 0, 0, new_y)) %>%
  mutate(new_x = if_else(abs(long - long_0) == (pi / 2), R * (long - long_0) * cos(lat), new_x),
         new_y = if_else(abs(long - long_0) == (pi / 2), pi / 2 * R * sin(lat), new_y)) %>%
  mutate(new_x = if_else(abs(lat) == (pi / 2), 0, new_x),
         new_y = if_else(abs(lat) == (pi / 2), R * lat, new_y)) %>%
  mutate(new_x = if_else(graph == 2, new_x, save_new_x),
         new_y = if_else(graph == 2, new_y, save_new_y)) %>%
  select(-save_new_x, -save_new_y) 

ggplot(data = projection[projection$graph == 1, ]) +
  geom_point(aes(x = new_x, y = new_y)) +
  coord_equal() +
  scale_x_continuous(limits = c(-(pi / 2), (pi / 2))) +
  scale_y_continuous(limits = c(-(pi / 2), (pi / 2)))
ggplot(data = projection[projection$graph == 2, ]) +
  geom_point(aes(x = new_x, y = new_y)) +
  coord_equal() +
  scale_x_continuous(limits = c(-(pi / 2), (pi / 2))) +
  scale_y_continuous(limits = c(-(pi / 2), (pi / 2)))
