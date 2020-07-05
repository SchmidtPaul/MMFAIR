pacman::p_load(dplyr, purrr, readr, tibble, tidyr)   

# heterogeneous_error_variance
agridat::mcconway.turnip %>%
  as_tibble() %>%
  mutate(densf = density %>% as.factor) %>% 
  write_delim(path="SAS/data_heterogeneous_error_variance.txt")

 # autoregressive_time_series
agriTutorial::sorghum %>%
  rename(block = Replicate, plot = factplot) %>% 
  dplyr::select(y, variety, block, plot, factweek, varweek) %>% 
  write_delim(path="SAS/data_autoregressive_time_series.txt")
