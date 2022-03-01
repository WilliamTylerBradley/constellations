# try Sample Elimination for Generating Poisson Disk Sample Sets
# http://www.cemyuksel.com/research/sampleelimination/sampleelimination.pdf
library(tidyverse)

starting_n <- 500
df <- tibble(
  x = rnorm(starting_n),
  y = rnorm(starting_n),
  z = rnorm(starting_n)
) %>%
  mutate(normalize = sqrt(x^2 + y^2 + z^2)) %>%
  mutate(x = x / normalize,
         y = y / normalize,
         z = z / normalize) %>%
  select(-normalize) %>%
  rowid_to_column("id") %>%
  mutate(max_distance = .25)

# https://math.stackexchange.com/questions/1304169/distance-between-two-points-on-a-sphere
df <- df %>%
  rename(id2 = id,
         x2 = x,
         y2 = y,
         z2 = z,
         max_distance2 = max_distance) %>%
  expand_grid(df) %>%
  filter(id != id2) %>%
  mutate(distance = acos(x*x2 + y*y2 + z*z2)) 


while(any(df$max_distance >= df$distance) | any(df$max_distance2 >= df$distance)) {

  drop_id <- df %>%
    filter(distance == min(distance)) %>%
    slice(1) %>%
    pull(id)
  
  df <- df %>%
    filter(id != drop_id & id2 != drop_id)

}

df <- df %>%
  select(id, x, y, z, max_distance) %>%
  unique()

# https://en.wikipedia.org/wiki/Nicolosi_globular_projection
## Use these two to graph the halfs
projection <- df %>% 
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

ggplot(data = projection) +
  geom_point(aes(x = long, y = lat, color = graph))
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





