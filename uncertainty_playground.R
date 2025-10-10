source("External/Vs.R", local = FALSE)
library(VSLiteR)
library(dplR) # for reading in dendro data
library(ggplot2) # for plotting later
library(tidyverse)
library(DEoptim)
library(purrr)

beech <- read.rwl("Data/buche_chrono.rwl")
climate <- read.csv2("Data/climate_bausenberg.csv")[, c(1, 2, 3, 6)]
spruce <- read.rwl("Data/spruce_alpine.rwl")
climate_alpine <- read.csv2("Data/climate_alpine.csv")
cordex_alpine <- read.csv2("Data/lfu_cordex_ensemble_monthly_alpine.csv")
cordex_beech <- read.csv2("Data/lfu_cordex_ensemble_monthly_beech.csv")

input_historic <- make_vsinput_historic(beech, climate)

beech_params <- vs_params(input_historic$trw,
                          input_historic$tmean,
                          input_historic$prec,
                          input_historic$syear,
                          input_historic$eyear,
                          .phi = 50) # approx. latitude in degrees

beech_params2 <- beech_params

input_transient <- make_vsinput_transient(climate)

beech_forward <- vs_run_forward(beech_params2,
                                input_transient$tmean,
                                input_transient$prec,
                                input_transient$syear,
                                input_transient$eyear,
                                .phi = 50)


make_forward_model <- function(x) {
  input <- make_vsinput_transient(x)
  vs_run_forward(beech_params2,
                 input$temp,
                 input$prec,
                 input$syear,
                 input$eyear,
                 .phi = 50)
}

forward_beech <- cordex_beech %>% 
  filter(year > 2020) %>% 
  group_by(rcp, gcm, rcm) %>% 
  filter(length(year) == 80 * 12) %>% 
  nest() %>% 
  mutate(vs_forward = purrr::map(data, make_forward_model)) %>% 
  unnest(vs_forward)

# Option 1: show all combinations per RCP
forward_beech %>% 
  ggplot(aes(year, trw, colour = gcm, linetype = rcm)) +
  geom_line() +
  facet_wrap(~ rcp)

# Option 2: compute uncertainty over GCM/RCM combinations
forward_beech %>% 
  group_by(rcp, year) %>% 
  summarise(
    trw_sd = sd(trw),
    trw = mean(trw)
  ) %>% 
  ggplot(aes(year, trw)) +
  geom_ribbon(aes(ymin = trw - trw_sd, ymax = trw + trw_sd)) +
  # geom_line(colour = "red") +
  geom_smooth(colour = "red", se = FALSE) +
  facet_wrap(~ rcp)

# Option 3: compute uncertainty over GCM/RCM combinations
smooth_trw <- function(x) {
  trw <- predict(loess(trw ~ year, data = x))
  data.frame(
    year = x$year,
    trw = trw
  )
}

forward_beech %>% 
  group_by(rcp, gcm, rcm) %>% 
  nest() %>% 
  mutate(trw = map(data, smooth_trw)) %>%
  unnest(trw) %>% 
  group_by(rcp, year) %>% 
  summarise(
    trw_sd = sd(trw),
    trw = mean(trw)
  ) %>% 
  ggplot(aes(year, trw)) +
  geom_ribbon(aes(ymin = trw - trw_sd, ymax = trw + trw_sd)) +
  geom_line(colour = "red") +
  facet_wrap(~ rcp)
