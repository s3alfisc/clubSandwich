library(clubSandwich)
library(plm)

rm(list=ls())
data(MortalityRates)

# subset for deaths in motor vehicle accidents, 1970-1983
MV_deaths <- subset(MortalityRates, cause=="Motor Vehicle" & 
                      year <= 1983 & !is.na(beertaxa))

MV_deaths <- within(MV_deaths, {
  year <- factor(year)
  legal_cent <- legal - tapply(legal, state, mean)[factor(state)]
  beer_cent <- beertaxa - tapply(beertaxa, state, mean)[factor(state)]
})

# within estimation

plm_within <- plm(mrate ~ legal + beertaxa, data = MV_deaths, 
                      effect = "individual", index = c("state","year"))


# between-series estimation

plm_between <- plm(mrate ~ 0 + legal + beertaxa, data = MV_deaths, 
                  effect = "individual", index = c("state","year"),
                  model = "between")

# random effects estimation

plm_random <- plm(mrate ~ 0 + legal + beertaxa, data = MV_deaths, 
                  effect = "individual", index = c("state","year"),
                  model = "random")


# artificial Hausman test

plm_Hausman <- plm(mrate ~ 0 + legal + beertaxa + legal_cent + beer_cent, 
                   data = MV_deaths,
                   effect = "individual", index = c("state","year"),
                   model = "random")

coef(plm_within)
coef(plm_between)
coef(plm_random)
coef(plm_Hausman)

all.equal(coef(plm_Hausman)[1:2], coef(plm_between))
all.equal(coef(plm_Hausman)[3:4], coef(plm_within) - coef(plm_between), check.attributes = FALSE)
