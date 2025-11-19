library(tidyverse)

historic <- read_csv2("Data/climate_sailerhausen.csv")
historic_b1 <- historic %>% 
  filter(
      site == "b1"
  ) %>% 
  select(year, month, tmean, prec)

historic_b2 <- historic %>% 
  filter(
    site == "b2"
  ) %>% 
  select(year, month, tmean, prec)

projection_b1 <- read_csv2("lfu_cordex_ensemble_monthly_b1.csv")
projection_b2 <- read_csv2("lfu_cordex_ensemble_monthly_b2.csv")

add_historic <- function(x, historic, cutoff = 1990) {
  x <- filter(x, year > cutoff)
  x <- rbind(historic, x) %>% 
    arrange(year, month)
  return(x)
}

b1_padded <- projection_b1 %>% 
  rename(tmean = temp) %>% 
  group_by(rcp, gcm, rcm) %>% 
  nest() %>% 
  mutate(padded = purrr::map(data, add_historic, historic_b1, 1990)) %>% 
  select(-data) %>% 
  unnest(padded)

b2_padded <- projection_b2 %>% 
  rename(tmean = temp) %>% 
  group_by(rcp, gcm, rcm) %>% 
  nest() %>% 
  mutate(padded = purrr::map(data, add_historic, historic_b2, 1990)) %>% 
  select(-data) %>% 
  unnest(padded)

write_excel_csv2(b1_padded, "Data/cordex_b1_padded.csv")
write_excel_csv2(b2_padded, "Data/cordex_b2_padded.csv")

b2_padded %>% 
  group_by(rcp, year) %>% 
  summarise(temp = mean(tmean)) %>% 
  ggplot(aes(year, temp, colour = rcp)) +
  geom_line()


