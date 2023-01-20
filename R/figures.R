

generate_figs_2 <- function(){
  library(stringr)
  source('../R/get_hhs.R')
  hhs <- get_hhs()

  ### compute US

  hhs_us <- hhs %>% dplyr::group_by(date) %>% dplyr::summarize(state="US",previous_day_admission_adult_covid_confirmed = sum(previous_day_admission_adult_covid_confirmed,na.rm=T))
  hhs_us$covid_per_cap <- hhs_us$previous_day_admission_adult_covid_confirmed/300000000

  hhs <- rbind(hhs_us,hhs %>% select(state,date,previous_day_admission_adult_covid_confirmed,covid_per_cap))
}


