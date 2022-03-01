# try Mitchell's best-candidate algorithm
# https://observablehq.com/@mbostock/best-candidate-circles
library(tidyverse)

points_n <- 30
distance_limit <- .01

df <- tibble(
  x = c(rnorm(1), rep(0, points_n - 1)),
  y = c(rnorm(1), rep(0, points_n - 1)),
  z = c(rnorm(1), rep(0, points_n - 1))
) %>%
  mutate(normalize = sqrt(x^2 + y^2 + z^2)) %>%
  mutate(normalize = if_else(normalize == 0, 1, normalize)) %>%
  mutate(x = x / normalize,
         y = y / normalize,
         z = z / normalize) %>%
  select(-normalize) %>%
  rowid_to_column("id") 

candidates_n <- 25
for(i in seq(2, points_n)) {
  
  candidates <- tibble(
    can_x = rnorm(candidates_n),
    can_y = rnorm(candidates_n),
    can_z = rnorm(candidates_n)
  ) %>%
    mutate(normalize = sqrt(can_x^2 + can_y^2 + can_z^2)) %>%
    mutate(can_x = can_x / normalize,
           can_y = can_y / normalize,
           can_z = can_z / normalize) %>%
    select(-normalize) %>%
    rowid_to_column("can_id") 
  
  best_candidate <- df %>%
    filter(x != 0 | y != 0 | z != 0) %>%
    expand_grid(candidates) %>%
    mutate(distance = acos(x * can_x + 
                             y * can_y + 
                             z * can_z)) %>%
    group_by(can_id) %>%
    summarise(distance = min(distance)) %>%
    filter(distance == max(distance)) %>%
    filter(distance > distance_limit) %>%
    slice(1) %>%
    pull(can_id)
  
  if(length(best_candidate) > 0) {
    best_candidate <- candidates %>%
      filter(can_id == best_candidate) %>%
      mutate(can_id = i) %>%
      rename_with(~ gsub("can_", "", .x))
    
    df <- df %>%
      rows_update(best_candidate, by = "id")
  }

}

df <- df %>%
  filter(x != 0 | y != 0 | z != 0)

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
  geom_point(aes(x = long, y = lat, color = graph)) +
  scale_x_continuous(limits = c(-pi, pi)) +
  scale_y_continuous(limits = c(-(pi / 2), (pi / 2)))
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





