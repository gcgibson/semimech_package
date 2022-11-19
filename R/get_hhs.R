
source("R/av7fun.R")
library(dplyr)
get_hhs <- function() {

  require(RSocrata)
  require(readr)
  require(dplyr)

  # Socrata API token for this project (do not share)
  token <- "llBMOEMWQQKKIg7XEAwuNVri4"

  # Use Socrata API
  hhs <- read.socrata("https://healthdata.gov/resource/g62h-syeh.csv",
                      app_token = token)

  ## Read in data for states
  load("data/fips.rda")
  load("data/uspop.rda")

  colnames(fips) <- c("fips", "state_full", "state", "alphacount")
  colnames(uspop) <- c("state_full", "pop")

  # Region of

  fips <- fips %>%
    mutate(fips = str_pad(fips, 2, pad = "0"))

  # Clean data
  # 1. Attach state fips codes and population data
  # 2. Calculate 7-day average of flu and covid, also on per-cap
  hhs <- hhs %>%
    mutate(date = as.Date(date)) %>%
    left_join(fips) %>%
    left_join(uspop) %>%
    filter(!is.na(state_full)) %>% # For now, only 50 states + DC
    # Create more convenient names for admissions
    mutate(flu = previous_day_admission_influenza_confirmed,
           covid = previous_day_admission_adult_covid_confirmed,
           flu_coverage = previous_day_admission_influenza_confirmed_coverage,
           covid_coverage = previous_day_admission_adult_covid_confirmed_coverage) %>%
    # Per-capita
    mutate(flu_per_cap = flu / pop * 1e5) %>%
    mutate(covid_per_cap = covid / pop * 1e5) %>%
    # 7-day average for admissions
    arrange(state, date) %>%
    group_by(state) %>%
    mutate(flu_coverage_av7 = av7fun_center(flu_coverage, date),
           covid_coverage_av7 = av7fun_center(covid_coverage, date),
           flu_av7 = av7fun(flu, date),
           covid_av7 = av7fun(covid, date),
           flu_per_cap_av7 = av7fun(flu_per_cap, date),
           covid_per_cap_av7 = av7fun(covid_per_cap, date),
           flu_av7_lag = av7fun_lag(flu, date),
           flu_av7_center = av7fun_center(flu, date),
           flu_per_cap_av7_lag = av7fun_lag(flu_per_cap, date),
           flu_per_cap_av7_center = av7fun_center(flu_per_cap, date),
           covid_av7_lag = av7fun_lag(covid, date),
           covid_av7_center = av7fun_center(covid, date),
           covid_per_cap_av7_lag = av7fun_lag(covid_per_cap, date),
           covid_per_cap_av7_center = av7fun_center(covid_per_cap, date)) %>%
    # mutate(ratio_log = log(flu_per_cap_av7) - log(covid_per_cap_av7)) %>%
    ungroup() %>%
    mutate(date = date - 1)

  ## Return
  # hhs

  #Remove US Virgin Islands and American Samoa
  hhs %>% filter(!(state %in% c("VI", "AS")))

}


