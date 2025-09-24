devtools::install_github("suztolwinskiward/VSLiteR")
source("https://raw.githubusercontent.com/cszang/runvslite/refs/heads/main/runvslite.R")
install.packages("DEoptim")
pkgbuild::find_rtools()


#Prepare data
library(VSLiteR)
library(dplR) # for reading in dendro data
library(ggplot2) # for plotting later
library(tidyverse)
library(DEoptim)

#Prepare Input data
beech <- read.rwl("buche_chrono.rwl")
climate <- read.csv2("climate_bausenberg.csv")[, c(1, 2, 3, 6)]
beech
climate

#Calibration Run
input_historic <- make_vsinput_historic(beech, climate)

beech_params <- vs_params(input_historic$trw,
                          input_historic$tmean,
                          input_historic$prec,
                          input_historic$syear,
                          input_historic$eyear,
                          .phi = 50) # approx. latitude in degrees

#Run model forward
input_transient <- make_vsinput_transient(climate)

beech_forward <- vs_run_forward(beech_params,
                                input_transient$tmean,
                                input_transient$prec,
                                input_transient$syear,
                                input_transient$eyear,
                                .phi = 50)

#Compare the model and observation
beech_observed <- data.frame(
  year = as.numeric(rownames(beech)),
  observed = scale(beech$bee)
)

beech_combined <- merge(beech_observed, beech_forward)
beech_combined$modelled <- scale(beech_combined$trw)
beech_combined$trw <- NULL
beech_combined <- tidyr::pivot_longer(beech_combined, -1, names_to = "variant",
                                      values_to = "rwi")

ggplot(beech_combined, aes(year, rwi)) +
  geom_line(aes(colour = variant))


# 1st sensitivity analysis
bp1 <- c(T1 = 2, T2 = 9, M1 = 0.01, M2 = 0.5)


bf1 <- vs_run_forward(bp1,
                                input_transient$tmean,
                                input_transient$prec,
                                input_transient$syear,
                                input_transient$eyear,
                                .phi = 50)

#Compare the model and observation

bc1 <- merge(beech_observed, bf1)
bc1$modelled <- scale(bc1$trw)
bc1$trw <- NULL
bc1 <- tidyr::pivot_longer(bc1, -1, names_to = "variant",
                                      values_to = "rwi")

ggplot(bc1, aes(year, rwi)) +
  geom_line(aes(colour = variant))



# 2nd Sensitivity analysis

bp2 <- c(T1 = 5, T2 = 13, M1 = 0.01, M2 = 0.2)

bf2 <- vs_run_forward(bp2,
                      input_transient$tmean,
                      input_transient$prec,
                      input_transient$syear,
                      input_transient$eyear,
                      .phi = 50)

#Compare the model and observation

bc2 <- merge(beech_observed, bf2)
bc2$modelled <- scale(bc2$trw)
bc2$trw <- NULL
bc2 <- tidyr::pivot_longer(bc2, -1, names_to = "variant",
                           values_to = "rwi")

ggplot(bc2, aes(year, rwi)) +
  geom_line(aes(colour = variant))


# 3rd Sensitivity analysis

bp3 <- c(T1 = 7, T2 = 18, M1 = 0.03, M2 = 0.3)

bf3 <- vs_run_forward(bp3,
                      input_transient$tmean,
                      input_transient$prec,
                      input_transient$syear,
                      input_transient$eyear,
                      .phi = 50)

#Compare the model and observation

bc3 <- merge(beech_observed, bf3)
bc3$modelled <- scale(bc3$trw)
bc3$trw <- NULL
bc3 <- tidyr::pivot_longer(bc3, -1, names_to = "variant",
                           values_to = "rwi")

ggplot(bc3, aes(year, rwi)) +
  geom_line(aes(colour = variant))

# 4th Sensitivity analysis

bp4 <- c(T1 = 3.5, T2 = 11, M1 = 0.02, M2 = 0.1)

bf4 <- vs_run_forward(bp4,
                      input_transient$tmean,
                      input_transient$prec,
                      input_transient$syear,
                      input_transient$eyear,
                      .phi = 50)

#Compare the model and observation

bc4 <- merge(beech_observed, bf4)
bc4$modelled <- scale(bc4$trw)
bc4$trw <- NULL
bc4 <- tidyr::pivot_longer(bc4, -1, names_to = "variant",
                           values_to = "rwi")

ggplot(bc4, aes(year, rwi)) +
  geom_line(aes(colour = variant))

##########################################################
# Correlation of various parameters

Tp <- c(T1 = 0, T2 = 9, M1 = 0.01, M2 = 0.1)

Tf <- vs_run_forward(Tp,
                      input_transient$tmean,
                      input_transient$prec,
                      input_transient$syear,
                      input_transient$eyear,
                      .phi = 50,)

Tc <- merge(beech_observed, Tf)
Tc$modelled <- scale(Tc$trw)
Tc$trw <- NULL
Tc <- tidyr::pivot_longer(Tc, -1, names_to = "variant",
                           values_to = "rwi")

ggplot(bc3, aes(year, rwi)) +
  geom_line(aes(colour = variant))

cor (Tc$rwi [Tc$variant == "modelled"], 
     Tc$rwi [Tc$variant == "observed"])


cor (beech_combined$rwi [beech_combined$variant == "modelled"], 
     beech_combined$rwi [beech_combined$variant == "observed"])

#####################################################################
# Data for spruce

spruce <- read.rwl("spruce_alpine.rwl")
climate_alpine <- read.csv2("climate_alpine.csv")


#Calibration Run
input_historic_s <- make_vsinput_historic(spruce, climate_alpine)

spruce_params <- vs_params(input_historic_s$trw,
                          input_historic_s$temp,
                          input_historic_s$prec,
                          input_historic_s$syear,
                          input_historic_s$eyear,
                          .phi = 50) # approx. latitude in degrees

#Run model forward
input_transient_s <- make_vsinput_transient(climate_alpine)

spruce_forward <- vs_run_forward(spruce_params,
                                input_transient_s$temp,
                                input_transient_s$prec,
                                input_transient_s$syear,
                                input_transient_s$eyear,
                                .phi = 50)

#Compare the model and observation
spruce_observed <- data.frame(
  year = as.numeric(rownames(spruce)),
  observed = scale(spruce$std)
)

spruce_combined <- merge(spruce_observed, spruce_forward)
spruce_combined$modelled <- scale(spruce_combined$trw)
spruce_combined$trw <- NULL
spruce_combined <- tidyr::pivot_longer(spruce_combined, -1, names_to = "variant",
                                      values_to = "rwi")

ggplot(spruce_combined, aes(year, rwi)) +
  geom_line(aes(colour = variant))

cor (spruce_combined$rwi [spruce_combined$variant == "modelled"], 
     spruce_combined$rwi [spruce_combined$variant == "observed"])


#sensitivity analysis section

sp <- c(T1 = 6.5, T2 = 11.9, M1 = 0.01, M2 = 0.1)

sf <- vs_run_forward(sp,
                     input_transient_s$temp,
                     input_transient_s$prec,
                     input_transient_s$syear,
                     input_transient_s$eyear,
                     .phi = 50)

sc <- merge(spruce_observed, sf)
sc$modelled <- scale(sc$trw)
sc$trw <- NULL
sc <- tidyr::pivot_longer(sc, -1, names_to = "variant",
                          values_to = "rwi")

ggplot(sc, aes(year, rwi)) +
  geom_line(aes(colour = variant))

cor (sc$rwi [sc$variant == "modelled"], 
     sc$rwi [sc$variant == "observed"])
