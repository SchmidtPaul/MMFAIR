pacman::p_load(dplyr, purrr, readr, tibble, tidyr)   

# heterogeneous_error_variance
agridat::mcconway.turnip %>%
  as_tibble() %>%
  mutate(densf = density %>% as.factor) %>% 
  write_delim(path="SAS/data_heterogeneous_error_variance.txt")
