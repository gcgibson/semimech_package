extract_waves <- function(hhs,last_date){

  hhs <- hhs[hhs$date > "2021-02-01",]
  hhs <- hhs[hhs$date <= last_date,]
  wave_mat <- matrix(NA,nrow=10000,ncol=1000)

  wave_mat_idx <- 1
  for (state in unique(hhs$state)){
    tmp <- hhs[hhs$state == state,]$covid_per_cap
    inc_dat_sub_sub_tmp_dens <- sample(1:length(tmp),prob=tmp/sum(tmp),size=100000,replace = T)
    bw_d <- density(inc_dat_sub_sub_tmp_dens,from = 1, to=length(tmp),n = length(tmp))


    tmp_lowess_diff <- diff(bw_d$y)
    wave_starts <- c()
    for (t in 2:length(tmp_lowess_diff)){
      if(tmp_lowess_diff[t-1] < 0 & tmp_lowess_diff[t] >0){
        wave_starts <- c(wave_starts,t)
      }
    }
    if (length(wave_starts) > 2){
      for (wave_start_idx in 2:length(wave_starts)){
          tmptmp <- tmp[wave_starts[wave_start_idx-1]:wave_starts[wave_start_idx]]
          wave_mat[wave_mat_idx,1:length(tmptmp)] <- tmptmp
          wave_mat_idx <- wave_mat_idx +1
      }
    }

  }
  first_non_na <-which(is.na(wave_mat[,1]))[1]

  return (wave_mat[1:first_non_na,])
}


#' @export
#'
#'

generate_samples <- function(start,dates,region,inc_dat){
  library(dplyr)
  library(gam)
  library(glmnet)

  weights_mat <- matrix(NA, nrow=1000000,ncol=2)
  weights_mat_idx <- 1
  for (date_idx in 1:length(dates)){
    # set local date
    last_date <- dates[date_idx]
    # extract historical waves
    wave_mat <- extract_waves(inc_dat,last_date)

    ## hardcoded manipulation for now
    inc_dat_sub <- inc_dat %>%
      dplyr::select(state, date, previous_day_admission_adult_covid_confirmed)
    inc_dat_sub <- inc_dat_sub %>% dplyr::select(state,date,previous_day_admission_adult_covid_confirmed)
    inc_dat_sub$covid <- inc_dat_sub$previous_day_admission_adult_covid_confirmed
    inc_dat_sub <- inc_dat_sub %>% dplyr::select(state,date,covid)


    inc_dat_sub <- inc_dat_sub[inc_dat_sub$date <= last_date,]
    inc_dat_sub <- inc_dat_sub[inc_dat_sub$date > start,]

    ### create the list which will house fitted objects by state
    df_to_submit_list <- list()

    for (state in unique(inc_dat_sub$state)){
      local_pop <-uspop[uspop$location== state.name[which(state.abb == state )],]$Pop_Est_2019
      inc_dat_sub_sub <- inc_dat_sub[inc_dat_sub$state == state,]


      ### wave mat stuff
      wave_starts <- c(0)
      inc_dat_sub_sub_tmp <- inc_dat_sub_sub$covid
      inc_dat_sub_sub_tmp_dens <- sample(1:length(inc_dat_sub_sub_tmp),prob=inc_dat_sub_sub_tmp/sum(inc_dat_sub_sub_tmp),size=100000,replace = T)
      bw_d <- density(inc_dat_sub_sub_tmp_dens,from = 1, to=length(inc_dat_sub_sub_tmp),n = length(inc_dat_sub_sub_tmp))

      tmp_lowess_diff <- diff(bw_d$y)
      wave_starts <- c()
      for (t in 2:length(tmp_lowess_diff)){
        if(tmp_lowess_diff[t-1] < 0 & tmp_lowess_diff[t] >0){
          wave_starts <- c(wave_starts,t)
        }
      }

      inc_dat_sub_sub_tmp_train <- inc_dat_sub_sub_tmp[tail(wave_starts,1):length(inc_dat_sub_sub_tmp)]




      wave_mat_test <- wave_mat[,1:(length(inc_dat_sub_sub_tmp_train)+30)]
      wave_mat_test <-wave_mat_test[complete.cases(wave_mat_test),]
      if (is.null(dim(wave_mat_test))){
        wave_mat_test <- matrix(wave_mat_test,nrow=1)
      }
      wave_mat_test <- t(wave_mat_test)
      wave_mat_train <- wave_mat_test[1:length(inc_dat_sub_sub_tmp_train),]


      ###GAM stuff
      fit_spline <- gam(covid ~ s(t,df=6)  ,data=data.frame(covid=inc_dat_sub_sub_tmp_train,
                                                            t= 1:length(inc_dat_sub_sub_tmp_train),
                                                            wave_mat=wave_mat_train))






      #fit <- solve(t(basis)%*%basis)%*%t(basis)%*%inc_dat_sub_sub_tmp_train
      blorp <- data.frame(t=1:(length(inc_dat_sub_sub_tmp_train) +30),
                          wave_mat = wave_mat_test)#,
      #  wave_mat = wave_mat_test)
      preds_spline <- predict(fit_spline,newdata = blorp)
      preds_spline<- preds_spline[(length(preds_spline)-29):length(preds_spline)]

      if (!is.null(dim(wave_mat_train))){
        if (dim(wave_mat_train)[2] > 3){
          fit_historical <- SuperLearner(Y = inc_dat_sub_sub_tmp_train,X=data.frame(wave_mat_train),SL.library = c("SL.earth","SL.glmnet"))
          mse_historical <- sum(abs(inc_dat_sub_sub_tmp_train-predict(fit_historical,data.frame(wave_mat_train))$pred))

          preds_historical <- predict(fit_historical,newdata = data.frame(wave_mat_test))$pred
          preds_historical<- preds_historical[(length(preds_historical)-29):length(preds_historical)]
          weights <- c(1/sum(abs(fit_spline$residuals)),1/mse_historical)
          weights  <- weights/sum(weights)
        } else{
          weights <- c(1,1e-10)
          preds_historical <- rep(0,30)
        }
      } else{
        weights <- c(1,1e-10)
        preds_historical <- rep(0,30)
      }
      print (weights)
      weights_mat[weights_mat_idx,] <- weights
      weights_mat_idx <- weights_mat_idx + 1
      ###combine

      previous_diffs <- sd(diff(inc_dat_sub_sub_tmp))
      #### KCDE style differences
      wave_mat_sds <- rep(0,30)
      diff_by_horzizon <-rep(0,30)
      over_dispersion <-rep(0,30)

      for (h in 1:30){
        diff_by_horzizon[h] <- sd(diff(inc_dat_sub_sub_tmp,lag=h))
        wave_mat_diffs <- apply(wave_mat,2,function(x){
          return (diff(x,lag=h,na.rm=T))
        })

        wave_mat_sds[h] <- sd(wave_mat_diffs*local_pop/1e5,na.rm = T)
      }
      mean_done <- weights[1]*preds_spline + weights[2]*preds_historical
      plot(inc_dat_sub_sub_tmp_train,type='l',xlim=c(1,length(inc_dat_sub_sub_tmp_train)+30),ylim=c(0,max(mean_done,inc_dat_sub_sub_tmp_train)))
      lines(seq(length(inc_dat_sub_sub_tmp_train)+1,length(inc_dat_sub_sub_tmp_train)+30),preds_historical)
      preds_mat <- matrix(nrow=1000,ncol=30)
      for (row_idx in 1:1000){

        if (all(!is.na(mean_done))){
          preds_mat[row_idx,] <- rpois(length(mean_done),pmax(mean_done,1e-10)) + rnorm(length(mean_done),0,sd=diff_by_horzizon)
        }
      }

      preds_mat <- apply(preds_mat,2,function(x){pmax(0,x)})
      df_to_submit_list[[state]] <- preds_mat
    }


  }
  return (df_to_submit_list)
}




#' @export
#'
#'

generate_submission_file <- function(date,samples,model_name,sub_dir){
  library(dplyr)
  quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

  last_date <- date
  df_to_submit_list <- list()
  for (state in names(samples)){
    preds_mat <- samples[[state]]
    ncol_preds <- ncol(preds_mat)

    df_to_submit <- data.frame(location = rep(state,30*23),forecast_date = rep(last_date,30*23))



    horizon_quantile_mat <- matrix(NA,nrow=length(quantiles),ncol=30)

    for (h in 1:30){
      q_idx <- 1
      for (quant_lvl in quantiles){
        horizon_quantile_mat[q_idx,h] <- quantile(preds_mat[,h],probs=quant_lvl,na.rm = T)
        q_idx <- q_idx + 1
      }

    }
    df_to_submit$value <- c(horizon_quantile_mat)
    df_to_submit$type <- "quantile"
    df_to_submit$target_end_date <-rep(last_date + 1:30,each=23)
    df_to_submit$quantile <- rep(quantiles,30)
    df_to_submit$target <-rep(paste(1:30,"day ahead inc hosp") ,each=23)

    df_to_submit_list[[state]] <-df_to_submit
  }
  df_to_submit <-do.call(rbind.data.frame, df_to_submit_list)

  load('../data/fips.rda')
  colnames(fips) <- c("fips", "state_full", "state", "alphacount")
  glimpse(fips)
  library(stringr)
  fips <- fips %>%
    mutate(fips = str_pad(fips, 2, pad = "0"))
  fips$location <- fips$state
  df_to_submit <- df_to_submit %>% left_join(fips,by="location")
  df_to_submit$location <- df_to_submit$fips
  df_to_submit <- df_to_submit %>% dplyr::select(target,location,forecast_date,target_end_date,quantile,value,type)
  write.csv(file = paste0("../",sub_dir,"/",last_date,"-",model_name,".csv"),x=df_to_submit)
}

