# https://blog.mbedded.ninja/mathematics/geometry/spherical-geometry/finding-the-intersection-of-two-arcs-that-lie-on-a-sphere/
# https://www.dirkbertels.net/computing/greatCircles_files/great_circles_070206.pdf
library(tidyverse)

df <- tibble(
  x_1 = c(5896, 0, 4033, 1102),
  y_1 = c(2146, 3186, 711, 6250),
  z_1 = c(1106, 5517, 4880, 555)
) %>%
  rowid_to_column("id_1")

df <- df %>%
  rename_with(~ paste0(gsub("1", "", .x), "2")) %>%
  expand_grid(df) %>%
  filter(id_1 < id_2) %>%
  mutate(a_1_2 = y_1 * z_2 - z_1 * y_2,
         b_1_2 = z_1 * x_2 - x_1 * z_2,
         c_1_2 = x_1 * y_2 - y_1 * x_2) %>%
  rowid_to_column("line_id_1_2")
  
df <- df %>%
  rename_with(~ gsub("1", "3", .x)) %>%
  rename_with(~ gsub("2", "4", .x)) %>% 
  expand_grid(df) %>%
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
         n_2 = -n_1) # l and l_2 are the intersections

# Check if intersections are in segments
df2 <- df %>% 
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
           if_else(abs(sum_angle_1_l_1_2 - angle_1_2) < .00005, 1, 0),
         on_segment_1_l_2_2 = 
           if_else(abs(sum_angle_1_l_2_2 - angle_1_2) < .00005, 1, 0),
         on_segment_3_l_1_4 = 
           if_else(abs(sum_angle_3_l_1_4 - angle_3_4) < .00005, 1, 0),
         on_segment_3_l_2_4 = 
           if_else(abs(sum_angle_3_l_2_4 - angle_3_4) < .00001, 1, 0)) %>%
  mutate(intersects = if_else((on_segment_1_l_1_2 == 1 & on_segment_3_l_1_4) |
                                (on_segment_1_l_2_2 == 1 & on_segment_3_l_2_4),
                              1, 0)) %>%
  filter(intersects == 1)


# Check which intersections are in the bounding box of just the four points?




