library(tidyverse)

disk_sampling <- function(points_n, candidates_n, distance_limit) {
  
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
  
  return(df)
}

# https://en.wikipedia.org/wiki/Nicolosi_globular_projection
## Use these two to graph the halfs
projection <- function(df) {
  
  df <- df %>% 
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
    select(-R, -long_0, -b, -c, -d, -M, -N, 
           -save_new_x, -save_new_y)
  
  return(df)
}

constellations <- disk_sampling(20, 25, .01) %>%
  select(-id) %>%
  rowid_to_column("id") %>%
  projection() %>%
  rename_with(~ paste0("constellation_", .x))

stars <- rbind(disk_sampling(15, 25, .01),
               disk_sampling(15, 25, .01),
               disk_sampling(15, 25, .01),
               disk_sampling(15, 25, .01),
               disk_sampling(15, 25, .01),
               disk_sampling(15, 25, .01)) %>%
  select(-id) %>%
  rowid_to_column("id") %>%
  projection() %>%
  rename_with(~ paste0("star_", .x))

sky <- constellations %>%
  expand_grid(stars) %>%
  mutate(distance = acos(constellation_x * star_x + 
                           constellation_y * star_y + 
                           constellation_z * star_z)) %>%
  filter(star_graph == constellation_graph) %>% # make sure they cross over line
  group_by(star_id) %>%
  filter(distance == min(distance))

# move to the constellation so the boundaries are cleaner
move_percentage <- .33
star_movement <- sky %>%
  mutate(star_lat = star_lat + (constellation_lat - star_lat) * move_percentage,
         star_long = star_long + (constellation_long - star_long) * move_percentage) %>%
  mutate(star_x = sin(pi/2 - star_lat) * cos(star_long),
         star_y = sin(pi/2 - star_lat) * sin(star_long),
         star_z = cos(pi/2 - star_lat))

star_movement <- star_movement %>%
  select(star_id, star_x, star_y, star_z) %>%
  rename_with(~ gsub("star_", "", .x)) %>%
  projection() %>%
  rename_with(~ paste0("star_", .x))

sky <- sky %>%
  rows_update(star_movement, by = "star_id")

# Next connecting the star
# need to get intersecting line segments

max_stars_constellation <- sky %>%
  group_by(constellation_id) %>%
  summarise(cnt = n_distinct(star_id)) %>%
  filter(cnt == max(cnt)) %>%
  slice(1) %>%
  pull(constellation_id)

max_stars_constellation <- sky %>%
  group_by(constellation_id) %>%
  summarise(cnt = n_distinct(star_id))

max_stars_constellation <- 2

sky <- sky %>%
  filter(constellation_id == max_stars_constellation) %>%
  ungroup()

connections <- sky %>%
  select(starts_with("star")) %>%
  rename_with(~ paste0(gsub("star_", "", .x), "_1")) %>%
  select(id_1, x_1, y_1, z_1)

connections <- connections %>%
  rename_with(~ paste0(gsub("1", "", .x), "2")) %>%
  expand_grid(connections) %>%
  filter(id_1 < id_2) %>%
  mutate(a_1_2 = y_1 * z_2 - z_1 * y_2,
         b_1_2 = z_1 * x_2 - x_1 * z_2,
         c_1_2 = x_1 * y_2 - y_1 * x_2) %>%
  rowid_to_column("line_id_1_2")

intersections <- connections %>%
  rename_with(~ gsub("1", "3", .x)) %>%
  rename_with(~ gsub("2", "4", .x)) %>% 
  expand_grid(connections) %>%
  filter(line_id_1_2 < line_id_3_4) %>%
  filter((id_1 != id_3 & id_1 != id_4) &
           (id_2 != id_3 & id_2 != id_4)) %>% # remove ones that share a star
  mutate(l_1 = b_1_2 * c_3_4 - c_1_2 * b_3_4,
         m_1 = c_1_2 * a_3_4 - a_1_2 * c_3_4,
         n_1 = a_1_2 * b_3_4 - b_1_2 * a_3_4) %>%
  mutate(normal = sqrt(l_1^2 + m_1^2 + n_1^2)) %>%
  mutate(l_1 = l_1 / normal,
         m_1 = m_1 / normal,
         n_1 = n_1 / normal) %>%
  mutate(l_2 = -l_1,
         m_2 = -m_1,
         n_2 = -n_1) %>%
  mutate(normal_1 = sqrt(x_1^2 + y_1^2 + z_1^2),
         normal_2 = sqrt(x_2^2 + y_2^2 + z_2^2),
         normal_3 = sqrt(x_3^2 + y_3^2 + z_3^2),
         normal_4 = sqrt(x_4^2 + y_4^2 + z_4^2),
         normal_l_1 = sqrt(l_1^2 + m_1^2 + n_1^2),
         normal_l_2 = sqrt(l_2^2 + m_2^2 + n_2^2)) %>%  # technically, all 1's
  mutate(angle_1_l_1 = acos(round((x_1 * l_1 + y_1 * m_1 + z_1 * n_1) / 
                                    (normal_1 * normal_l_1), 7)) * 180/pi,
         angle_l_1_2 = acos(round((x_2 * l_1 + y_2 * m_1 + z_2 * n_1) / 
                                    (normal_2 * normal_l_1), 7)) * 180/pi,
         angle_1_l_2 = acos(round((x_1 * l_2 + y_1 * m_2 + z_1 * n_2) / 
                                    (normal_1 * normal_l_2), 7)) * 180/pi,
         angle_l_2_2 = acos(round((x_2 * l_2 + y_2 * m_2 + z_2 * n_2) / 
                                    (normal_2 * normal_l_2), 7)) * 180/pi,
         angle_1_2 = acos(round((x_1 * x_2 + y_1 * y_2 + z_1 * z_2) / 
                                  (normal_1 * normal_2), 7)) * 180/pi,
         angle_3_l_1 = acos(round((x_3 * l_1 + y_3 * m_1 + z_3 * n_1) / 
                                    (normal_3 * normal_l_1), 7)) * 180/pi,
         angle_l_1_4 = acos(round((x_4 * l_1 + y_4 * m_1 + z_4 * n_1) / 
                                    (normal_4 * normal_l_1), 7)) * 180/pi,
         angle_3_l_2 = acos(round((x_3 * l_2 + y_3 * m_2 + z_3 * n_2) / 
                                    (normal_3 * normal_l_2), 7)) * 180/pi,
         angle_l_2_4 = acos(round((x_4 * l_2 + y_4 * m_2 + z_4 * n_2) / 
                                    (normal_4 * normal_l_2), 7)) * 180/pi,
         angle_3_4 = acos(round((x_3 * x_4 + y_3 * y_4 + z_3 * z_4) / 
                                  (normal_3 * normal_4), 7)) * 180/pi) %>%
  mutate(sum_angle_1_l_1_2 = angle_1_l_1 + angle_l_1_2,
         sum_angle_1_l_2_2 = angle_1_l_2 + angle_l_2_2,
         sum_angle_3_l_1_4 = angle_3_l_1 + angle_l_1_4,
         sum_angle_3_l_2_4 = angle_3_l_2 + angle_l_2_4) %>%
  mutate(on_segment_1_l_1_2 = 
           if_else(abs(sum_angle_1_l_1_2 - angle_1_2) < .001, 1, 0),
         on_segment_1_l_2_2 = 
           if_else(abs(sum_angle_1_l_2_2 - angle_1_2) < .001, 1, 0),
         on_segment_3_l_1_4 = 
           if_else(abs(sum_angle_3_l_1_4 - angle_3_4) < .001, 1, 0),
         on_segment_3_l_2_4 = 
           if_else(abs(sum_angle_3_l_2_4 - angle_3_4) < .001, 1, 0)) %>%
  mutate(intersects = if_else((on_segment_1_l_1_2 == 1 & on_segment_3_l_1_4) |
                                (on_segment_1_l_2_2 == 1 & on_segment_3_l_2_4),
                              1, 0)) %>%
  filter(intersects == 1)

#write.csv(intersections, "testing.csv")

connections <- connections %>%
  select(-a_1_2, -b_1_2, -c_1_2) %>%
  pivot_longer(-line_id_1_2,
               names_to = c(".value", "star"),
               names_sep =  "_") %>%
  projection() %>%
  select(-graph) %>%
  pivot_wider(id_cols = line_id_1_2,
              names_from = star,
              values_from = c(id, x, y, z, lat, long,
                              new_x, new_y))

connections <- connections %>%
  mutate(intersects = if_else(line_id_1_2 %in% intersections$line_id_1_2 |
                               line_id_1_2 %in% intersections$line_id_3_4, 
         1, 0))


ggplot(data = sky[sky$star_graph == 1, ]) +
  geom_segment(data = connections, aes(x = long_1, y = lat_1, 
                                       xend = long_2, yend = lat_2,
                                       color = intersects)) +
  geom_text(aes(x = star_long, y = star_lat, label = star_id)) +
  geom_text(data = connections,
            aes(x = (long_1 + long_2) / 2, 
                y = (lat_1 + lat_2) / 2, label = line_id_1_2)) +
  coord_equal()

ggplot(data = sky) +
  geom_segment(data = connections, aes(x = new_x_1, y = new_y_1, 
                                       xend = new_x_2, yend = new_y_2,
                                       color = intersects)) +
  geom_text(aes(x = star_new_x, y = star_new_y, label = star_id)) +
  geom_text(data = connections,
            aes(x = (new_x_1 + new_x_2) / 2, 
                y = (new_y_1 + new_y_2) / 2, label = line_id_1_2)) +
  coord_equal() #+
  # scale_x_continuous(limits = c(-(pi / 2), (pi / 2))) +
  # scale_y_continuous(limits = c(-(pi / 2), (pi / 2))) +
  # theme(legend.position = "none") 
ggplot(data = sky[sky$star_graph == 2, ]) +
  geom_point(aes(x = star_new_x, y = star_new_y, color = as.factor(constellation_id),
                 size = constellation_id)) +
  geom_point(aes(x = constellation_new_x, y = constellation_new_y, color = as.factor(constellation_id)),
             shape = 4) +
  coord_equal() +
  scale_x_continuous(limits = c(-(pi / 2), (pi / 2))) +
  scale_y_continuous(limits = c(-(pi / 2), (pi / 2))) +
  theme(legend.position = "none") 

## Try using plotly to see if they really cross or not
# maybe change projection to conformal?
## Peirce quincuncial projection
# sterographic
# https://astronomy.stackexchange.com/questions/18219/map-projections-used-for-star-maps