

generate_figs <- function(){
  source('../R/get_hhs.R')

  extract_waves_local_flu <- function(hhs,last_date,start_date){

    hhs <- hhs[hhs$date <= last_date,]
    hhs <- hhs[hhs$date >= start_date,]
    wave_mat <- matrix(NA,nrow=10000,ncol=1000)

    wave_mat_idx <- 1
    for (state in unique(hhs$state)){
      tmp <- hhs[hhs$state == state,]$previous_day_admission_influenza_confirmed
      tmp <- tmp[!is.na(tmp)]
      if (length(tmp) > 100){
        if (!all(tmp==0)){

          tmp <- tmp[!is.na(tmp)]
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
      }
    }
    first_non_na <-which(is.na(wave_mat[,1]))[1]
    ret_list <- list()
    ret_list[[1]] <- wave_mat[1:first_non_na,]
    ret_list[[2]] <- wave_starts
    return (ret_list)
  }









  library(stringr)
  hhs <- get_hhs()

  ### compute US

  hhs_us <- hhs %>% dplyr::group_by(date) %>% dplyr::summarize(state="US",previous_day_admission_influenza_confirmed = sum(previous_day_admission_influenza_confirmed,na.rm=T))
  hhs_us$flu_per_cap <- hhs_us$previous_day_admission_influenza_confirmed/300000000



  ###WAVE MAT FIG

  hhs_us <- hhs_us[hhs_us$date >= "2020-06-01",]
  library(ggplot2)
  wave_mat_ret <- extract_waves_local_flu(hhs_us,last_date = as.Date("2023-01-01"),start_date = as.Date("2020-06-01"))

  wave_mat <- wave_mat_ret[[1]]
  wave_starts <- wave_mat_ret[[2]]
  p <- ggplot(data=data.frame(y=wave_mat[1,1:300],x=1:300),aes(x=x,y=y)) + geom_line(col='#F8766D')
  p <- p + geom_line(data=data.frame(y=wave_mat[2,1:300],x=1:300),aes(x=x,y=y),col='#AFA100')
  p <- p + geom_line(data=data.frame(y=wave_mat[3,1:300],x=1:300),aes(x=x,y=y),col='#00BF7D')
  p <- p + geom_line(data=data.frame(y=wave_mat[4,1:300],x=1:300),aes(x=x,y=y),col='#00B0F6')
  p <- p + geom_line(data=data.frame(y=wave_mat[5,1:300],x=1:300),aes(x=x,y=y),col='#C77CFF')
  #p <- p + geom_line(data=data.frame(y=wave_mat[6,1:300],x=1:300),aes(x=x,y=y,col='Wave 6'))

  p <- p + theme_minimal() + ylab("US Incident Hospitalizations Per Capita")
  p



  hhs_us$wave <- unlist(lapply(1:length(hhs_us$date) , function(x){
    if(x < wave_starts[1]){
      return (NA)
    } else if(x >= wave_starts[1] & x < wave_starts[2]){
      return (1)
    } else if(x >= wave_starts[2] & x < wave_starts[3]){
      return (2)
    }else if(x >= wave_starts[3] & x < wave_starts[4]){
      return (3)
    }else if(x >= wave_starts[4] & x < wave_starts[5] ){
      return (4)
    } else{
      return ("Current Wave")
    }
  }))

  hhs_us$wave <- as.factor(hhs_us$wave)
  p2 <- ggplot(hhs_us,aes(x=date,y=flu_per_cap)) + geom_line(aes(col=wave))
  p2 <- p2 + theme_bw() + ylab("US Incident Hospitalizations Per Capita")

  p2
  library(cowplot)
  wave_det <- cowplot::plot_grid(p2  ,p+theme_bw(),nrow=2)
  ggsave("wave_det_flu.png",plot = wave_det,height = 8,width = 10)
}
