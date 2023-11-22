# source("/home/ben/research/NOC/projects/EcoWind/analysis/buoy_validation.R")

   flag_regional <- TRUE

   flag_single_panel <- FALSE

# Buoys.
   #buoy_list <- c("41001","41010","41025","42001")
   buoy_list <- c("46001","46002","46005","46006","46035","46059","46066","46070","46072","46073","46075","46078")
   buoy_list <- c("46075")

# Observation per hour.
   buoy_obs_freq <- c(1,2,1,1)

#-----------------------------------------------------------------------#
# Function to load WW3 1-D spectra.
# i_ncols     = number of data colums in the file (usually 6)
# i_spec_bins = number of spectral bins
# i_steps     = number of time steps
   func_load_spectra <- function(data_path,i_ncols,i_spec_bins,i_steps) {
      mat_spectra <- t(matrix(scan(data_path,quiet=TRUE),i_ncols))
      arr_spectra <- array(0,dim=c(i_spec_bins,2,i_steps))
      for ( i in 1:i_steps ) { arr_spectra[,,i] <- mat_spectra[( ( ( i-1 ) * i_spec_bins ) + 1):( i * i_spec_bins ),c(-2,-4,-5,-6)] }
      arr_spectra
   }

   spectral_bins <- 30

# Integration bounds.
   lower_bound <- 0.01
   upper_bound <- 0.70
   rel.tol <- .Machine$double.eps^0.25

#-----------------------------------------------------------------------#
# Load observed data.
# NDBC/NCEI data from USACE.
#-----------------------------------------------------------------------#
   #array_list_hs <-
   array_list_tm <- array(list(),dim=c(12,1,length(buoy_list),2))
   array_list_hs <- array(list(),dim=c(12,1,length(buoy_list),2))

   for (b.idx in 1:length(buoy_list)) {

      file_name_buoy <- paste("/backup/datasets/buoys/USACE/NDBC_complete_records/", list.files(path=paste("/backup/datasets/buoys/USACE/NDBC_complete_records/",sep=""), pattern=buoy_list[b.idx]), sep="")
      df_buoy_csv <- read.csv(file_name_buoy)
      vec_buoy_time1 <- strptime(df_buoy_csv[,1], "%Y-%m-%d %H:%M:%S")
      vec_buoy_hs1 <- df_buoy_csv$waveHs
      vec_buoy_tm1 <- df_buoy_csv$waveTm
# Remove missing data.
      vec_buoy_time <- vec_buoy_time1[!is.na(vec_buoy_hs1)]
      vec_buoy_hs <- vec_buoy_hs1[!is.na(vec_buoy_hs1)]
      vec_buoy_tm <- vec_buoy_tm1[!is.na(vec_buoy_hs1)]

#   wind.nc = nc_open("/home/ben/research/waves/laptop_large_files/ecmwf_2011-TC_025_uv.nc")
#   nc_time_len <- wind.nc$dim$time$len
#   wind_time <- ncvar_get(wind.nc, varid="time", start=c(1), count=c(nc_time_len))
#
#   mat_xlim <- c(1,124)
#   print(paste("Plot start time:",as.POSIXct(wind_time[mat_xlim[1]]*3600, origin = '1900-01-01', tz='GMT')))
#   time_labels <- wind_time[seq(1,124,4)]
#
#   date_start_41001 <- 1
##   date_range_41001 <- c(date_start_41001:(date_start_41001 + (mat_xlim[2]-mat_xlim[1]) ))
#   date_range_41001 <- seq(1,,6,124)
#   date_start_41010 <- 1
##   date_range_41010 <- c(date_start_41010:(date_start_41010 + (mat_xlim[2]-mat_xlim[1]) ))
#   date_range_41010 <- seq(1,,12,124)
#   date_start_42001 <- 1
##   date_range_42001 <- c(date_start_42001:(date_start_42001 + (mat_xlim[2]-mat_xlim[1]) ))
#   date_range_42001 <- seq(1,,6,124)
#
#-----------------------------------------------------------------------#
# Load WW3 data from "tab" files.
#-----------------------------------------------------------------------#
# Get length of WW3 output file.
   #file_name1 <- paste("/home/ben/research/waves/experiments/global/cam5_global/",exp_name,"/output/11Nov173816/tab_",buoy_list[1],"_1",sep="")
      vec_WW3_months <- c("01","02","03","04","05","06","07","08","09","10","11","12")

      vec_buoy_time_YM <- format( vec_buoy_time, "%Y%m")
      vec_file_WW3_tab <- vec_file_WW3_1D <- c()
      list_WW3_time <- list()
      list_WW3_hs <- list_WW3_tm2 <- list_WW3_tm1 <- list_WW3_tz <- list_WW3_1D_params <- list()
      vec_buoy_month_offset <- c()

      y_idx <- 1
      WW3_year <- 2016

      n_months <- 3
      for (m_idx in 1:3) {

# Read in WW3 param table files.
         vec_file_WW3_tab[m_idx] <- paste("/home/ben/research/NOC/projects/EcoWind/output/",WW3_year,"/output_",vec_WW3_months[m_idx],"/outp_params/tab_",vec_WW3_months[m_idx],"_",buoy_list[b.idx],sep="")
         mat_WW3_tab_raw <- t(matrix(scan(vec_file_WW3_tab[m_idx],skip=3),nrow=12))
         i_time_steps <- dim(mat_WW3_tab_raw)[1]

         list_WW3_time[[m_idx]] <- strptime( apply(X=mat_WW3_tab_raw[,1:4],MAR=1,FUN=paste,collapse=" "), "%Y%m%d %H %M %S")
         list_WW3_hs[[m_idx]] <- mat_WW3_tab_raw[,5]

# Read in WW3 1D spectra files.
         vec_file_WW3_1D[m_idx] <- paste("/home/ben/research/NOC/projects/EcoWind/output/",WW3_year,"/output_",vec_WW3_months[m_idx],"/outp_1D/1D_",vec_WW3_months[m_idx],"_",buoy_list[b.idx],sep="")
# Create a "CLEAN" file containing only the spectral data.
         system(command=paste("cat ",vec_file_WW3_1D[m_idx]," | grep '^  0\\.' > ",vec_file_WW3_1D[m_idx],"_CLEAN",sep=""),intern=T)
         array_WW3_1D <- func_load_spectra(paste(vec_file_WW3_1D[m_idx],"_CLEAN",sep=""),6,30,i_time_steps)
# Matrix to store itegrated (1D) wave parameters.
         mat_1D_params <- matrix(NA,nrow=i_time_steps,ncol=8)
# Loop over time steps.
         for ( t_idx in 1:i_time_steps ) {
# Fit spline function to spectra.
            func_spline <- splinefun(rbind(c(0,0),array_WW3_1D[,,t_idx]),method = c("monoH.FC"))
            func_spline_moment <- function(x,n) { x^n * func_spline(x) }
# m^2
            mat_1D_params[t_idx,1] <- integrate(func_spline_moment, lower_bound, upper_bound, 2, rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol, stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)[[1]]
# m^1
#            mat_1D_params[t_idx,2] <- integrate(func_spline_moment, lower_bound, upper_bound, 1, rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol, stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)[[1]]
# m^0
            mat_1D_params[t_idx,3] <- integrate(func_spline_moment, lower_bound, upper_bound, 0, rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol, stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)[[1]]
# m^-1
            mat_1D_params[t_idx,4] <- integrate(func_spline_moment, lower_bound, upper_bound, -1, rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol, stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)[[1]]
# Significant wave height: Hs = 4 * sqrt( m^0 )
            mat_1D_params[t_idx,5] <- 4 * sqrt( mat_1D_params[t_idx,3] )
# Energy period:           Te =  m^-1 / m^0
            mat_1D_params[t_idx,6] <- mat_1D_params[t_idx,4] / mat_1D_params[t_idx,3]
# Zero crossing period:    Tz =  sqrt( m^0 / m^2 )
            mat_1D_params[t_idx,7] <- sqrt( mat_1D_params[t_idx,3] / mat_1D_params[t_idx,1] )
# Tm01:                   Tm1 =  m^0 / m^1
#         mat_1D_params[t_idx,8] <- mat_1D_params[t_idx,3] / mat_1D_params[t_idx,2]
         }

         list_WW3_tz[[m_idx]] <- mat_1D_params[,7]
         #list_WW3_tm1[[m_idx]] <- mat_1D_params[,8]
         #list_WW3_tm2[[m_idx]] <- mat_1D_params[,1]
         list_WW3_1D_params[[m_idx]] <- mat_1D_params
      
         vec_match_months <- c(paste(WW3_year,vec_WW3_months[m_idx],sep=""))
         Lvec_buoy_time_month_idx <- vec_buoy_time_YM == vec_match_months
         vec_buoy_time_month_idx <- which(Lvec_buoy_time_month_idx)
         Lvec_buoy_time_month_idx[which(vec_buoy_time_YM == vec_match_months)[1]-1] <- TRUE
         vec_buoy_time_month <- vec_buoy_time[ Lvec_buoy_time_month_idx ]
         #vec_buoy_hs_month <- vec_buoy_hs[ Lvec_buoy_time_month_idx ]
         vec_buoy_month_offset[m_idx] <- which(Lvec_buoy_time_month_idx)[1]

# Find matching start date (usually beginning of monthly time series).
         list_match_buoy_idx <- list()
         list_match_WW3_idx <- list()
         vec_buoy_min <- min( abs(as.numeric(vec_buoy_time_month) - as.numeric(list_WW3_time[[m_idx]][1])) )
         if ( vec_buoy_min < 60*30+1 ) {
            buoy_start_idx <- which( abs(as.numeric(vec_buoy_time_month) - as.numeric(list_WW3_time[[m_idx]][1])) == vec_buoy_min )
            list_match_buoy_idx[[1]] <- buoy_start_idx
            list_match_WW3_idx[[1]] <- 1
# Find time step differences in buoy time series.
            vec_buoy_time_diff <- c()
            for (t_idx in 2:length(as.numeric(vec_buoy_time_month))) { vec_buoy_time_diff[t_idx-1] <- as.numeric(vec_buoy_time_month[t_idx])-as.numeric(vec_buoy_time_month[t_idx-1]) }
 
## Correction for sub-hourly time series (e.g. 46059).
#         if ( median(vec_buoy_time_diff) == 600 ) {
#            vec_buoy_time_month_idx_SAMP <- c(buoy_start_idx,buoy_start_idx+seq(6,,6,floor(length(vec_buoy_time_month) / 6)-1))
#            vec_buoy_time_month_idx <- which(Lvec_buoy_time_month_idx)[vec_buoy_time_month_idx_SAMP]
#            vec_buoy_time_month <- vec_buoy_time_month[vec_buoy_time_month_idx_SAMP]
#         }
## Find time step differences in SAMPLED buoy time series.
## Then locate any missing time steps.
#         vec_buoy_time_diff <- c()
#         for (t_idx in 2:length(as.numeric(vec_buoy_time_month))) { vec_buoy_time_diff[t_idx-1] <- as.numeric(vec_buoy_time_month[t_idx])-as.numeric(vec_buoy_time_month[t_idx-1]) }

            vec_buoy_miss <- which( vec_buoy_time_diff > 1.01*range(vec_buoy_time_diff)[1] )
            if ( length(vec_buoy_miss) >= 1 ) {
               for ( miss_idx in 1:length(vec_buoy_miss) ) {
                  list_match_buoy_idx[[miss_idx]] <- list_match_buoy_idx[[miss_idx]][1]:vec_buoy_miss[miss_idx]
                  list_match_WW3_idx[[miss_idx]] <- list_match_WW3_idx[[miss_idx]][1]:(list_match_WW3_idx[[miss_idx]][1]+length(list_match_buoy_idx[[miss_idx]])-1)
# To-do:
# Add test for > 30 minute difference.
                  WW3_min_loc <- min( abs(as.numeric(list_WW3_time[[m_idx]]) - as.numeric(vec_buoy_time_month[vec_buoy_miss[miss_idx]+1])) )
                  WW3_start_idx_loc <- which( abs(as.numeric(list_WW3_time[[m_idx]]) - as.numeric(vec_buoy_time_month[vec_buoy_miss[miss_idx]+1])) == WW3_min_loc )
                  list_match_buoy_idx[[miss_idx+1]] <- vec_buoy_miss[miss_idx]+1
                  list_match_WW3_idx[[miss_idx+1]] <- WW3_start_idx_loc
               }
               list_match_WW3_idx[[miss_idx+1]] <- list_match_WW3_idx[[miss_idx+1]][1]:length(list_WW3_time[[m_idx]])
               list_match_buoy_idx[[miss_idx+1]] <- list_match_buoy_idx[[miss_idx+1]][1]:length(unlist(list_match_WW3_idx))
            } else {
               list_match_buoy_idx[[1]] <- list_match_buoy_idx[[1]][1]:length(vec_buoy_time_month)
               list_match_WW3_idx[[1]] <- list_match_WW3_idx[[1]][1]:length(list_match_buoy_idx[[1]])
            }
# Capture buoy and WW3 paired parameters.
            array_list_hs[[m_idx,y_idx,b.idx,1]] <- vec_buoy_hs[Lvec_buoy_time_month_idx][unlist(list_match_buoy_idx)]
            array_list_hs[[m_idx,y_idx,b.idx,2]] <- list_WW3_hs[[m_idx]][unlist(list_match_WW3_idx)]
            #array_list_tm[[m_idx,y_idx,b.idx]] <- cbind(buoy=vec_buoy_tm[Lvec_buoy_time_month_idx][unlist(list_match_buoy_idx)],WW3_tz=list_WW3_tz[[m_idx]][unlist(list_match_WW3_idx)])
            array_list_tm[[m_idx,y_idx,b.idx,1]] <- vec_buoy_tm[Lvec_buoy_time_month_idx][unlist(list_match_buoy_idx)]
            array_list_tm[[m_idx,y_idx,b.idx,2]] <- list_WW3_tz[[m_idx]][unlist(list_match_WW3_idx)]
         } else {
            print(paste(" Date:",vec_match_months,"No match within 30 mins"))
         }
      }
   }

#=======================================================================#
# Plotting.
#-----------------------------------------------------------------------#
# Set first 3 days to NA.
   for (b_idx in 1:length(buoy_list)) {
      array_list_hs[[1,1,b_idx,1]][1:168] <- NA
      array_list_hs[[1,1,b_idx,2]][1:168] <- NA
      array_list_tm[[1,1,b_idx,1]][1:168] <- NA
      array_list_tm[[1,1,b_idx,2]][1:168] <- NA
   }

# Loop over months and years to get all data (or subsets, e.g. seasonal).
   df_reg <- data.frame(
                        buoy_hs=unlist(sapply(X=1:length(buoy_list),FUN=function(x) { unlist(array_list_hs[,,x,1]) })),
                        ww3_hs=unlist(sapply(X=1:length(buoy_list),FUN=function(x) { unlist(array_list_hs[,,x,2]) })),
                        buoy_tm=unlist(sapply(X=1:length(buoy_list),FUN=function(x) { unlist(array_list_tm[,,x,1]) })),
                        ww3_tm=unlist(sapply(X=1:length(buoy_list),FUN=function(x) { unlist(array_list_tm[,,x,2]) })),
                        buoy_lab=unlist( sapply(X=1:length(buoy_list),FUN=function(x) { rep(buoy_list[x],length(unlist(array_list_hs[,,x,1]))) }) ) )

   list_mat_plot_hs <- list_mat_plot_tm <- list()
   list_hs_rmse <- list_hs_rmsd <- list_hs_bias <- list_hs_cor <- list()

#----------------------------------------------------------------------------------#
# Single-panel plotting.
   if ( flag_single_panel ) {
# Plotting parameters.
      pl_mfrow <- c(1,1)
      pl_oma <- c(6,7,2,6)
      pl_mar <- c(18,18,16,7)
      pl_mgp <- c(12,6,0)
      pl_cex <- 6; pl_cex_main <- 11; pl_cex_lab <- 7; pl_cex_axis <- 9; pl_cex_leg <- 8

      vec_plot_idx <- 2
   
      if ( flag_OS ) {
         plot_title <- c("J-3 Offshore","S6-MF LR Offshore","S6-MF HR Offshore")[Sidx]
      } else {
         plot_title <- c("J-3 Nearshore","S6-MF LR Nearshore","S6-MF HR Nearshore")[Sidx]
      }

      cex_mtext <- 8
   
      if ( length(buoy_list) <= 2 ) {
         fig_scatter_filename <- paste0("./figures/",paste(buoy_list[b_idx]),"/scatter_",paste(buoy_list[b_idx],collapse='_'),".png")
      } else {
         fig_scatter_filename <- paste0("./figures/scatter_ALL.png")
         #fig_scatter_filename <- paste0("./figures/",fig_lab_region,"/scatter_adapt_",fig_lab_region,"_","M",m_limit,"_",vec_tandem_labs[S_idx],"_",buoy_radius,"km_BUOYS_",gsub('[.]','',cor_thresh),".png")
#----------------------------------------------------------------------------------#
# Multi-panel plotting.
      }
   } else {
# Plotting parameters.
      pl_mfrow <- c(length(buoy_list),2)
      pl_oma <- c(2,2,2,2)
      pl_mar <- c(12,14,8,5)
      pl_mgp <- c(9,3,0)
      pl_cex <- 4; pl_cex_main <- 6; pl_cex_lab <- 6; pl_cex_axis <- 5; pl_cex_leg <- 5

      vec_plot_idx <- 1:length(buoy_list)

      plot_title <- buoy_list[b_idx]

      cex_mtext <- 4

      if ( length(buoy_list) <= 2 ) {
         fig_scatter_filename <- paste0("./figures/",paste(buoy_list[b_idx]),"/scatter_",paste(buoy_list[b_idx],collapse='_'),".png")
         #fig_scatter_adapt_nearest_file_name <- paste0("./figures/",paste(buoy_list[b_idx_list[Bidx]]),"/scatter_adapt_",fig_lab_region,"_","M",m_limit,"_",vec_tandem_labs[S_idx],"_",buoy_radius,"km_",paste(buoy_list[b_idx_list[Bidx]],collapse='_'),"_",gsub('[.]','',cor_thresh),"_MULTI.png")
      } else {
         fig_scatter_filename <- paste0("./figures/scatter_ALL.png")
      }
   }


   #for (b_idx in 1:length(buoy_list)) {
   #   plot(list_mat_plot_hs[[b_idx]],xlim=c(0,15),ylim=c(0,15),main=paste(WW3_year,"NDBC",buoy_list[b_idx]))
   #   mtext(side=3, line=-9, adj=0.03, cex=cex_mtext, outer=FALSE, text=paste("Correlation: ",format(list_hs_cor[[b_idx]],digits=3),sep=''))
   #   mtext(side=3, line=-12, adj=0.03, cex=cex_mtext, outer=FALSE, text=paste("RMSE: ",format(round(list_hs_rmse[[b_idx]],2),nsmall=2),sep=''))
   #
   #   abline(0,1)
   #
   #   plot(list_mat_plot_tm[[b_idx]],xlim=c(0,15),ylim=c(0,15),main=paste(WW3_year,"NDBC",buoy_list[b_idx]))
   #   abline(0,1)
   #}

#----------------------------------------------------------------------------------#
# Open file.
   if ( flag_regional ) {
# https://r-graph-gallery.com/2d-density-plot-with-ggplot2
# https://stats.stackexchange.com/questions/12392/how-to-compare-two-datasets-with-q-q-plot-using-ggplot2
      require(ggplot2)
      fig_scatter_filename <- "./figures/ggplot_scatter_test"
# Data for Q-Q plot.
      qq_data <- as.data.frame(qqplot(df_reg$buoy_hs, df_reg$ww3_hs, plot.it=FALSE))

      p1 <- ggplot(df_reg, aes(x=buoy_hs, y=ww3_hs) ) +
      geom_hex(bins = 70) +
      xlim(0,10) + ylim(0,10) +
      ggtitle(paste0("Region",buoy_list[1])) + xlab("Buoys") + ylab("WW3") +
      geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 5) +
      scale_fill_continuous(type = "viridis") +
      theme(
            plot.title = element_text(size = 110, hjust = 0.5, margin = margin(t = 50, r = 0, b = 50, l = 0)),
            axis.title.x=element_text(size = 100),
            axis.title.y=element_text(size = 100, margin = margin(t = 0, r = -20, b = 0, l = 0)),
            #panel.grid.minor = element_blank(),
            #panel.grid.major = element_blank(),
            #panel.background = element_rect(fill = "black"),

            strip.text = element_text(size = 50, margin = margin(25,25,25,25)),
            strip.background = element_rect(fill = "white"),
            panel.spacing.x = unit(0, "lines"),
            panel.spacing.y = unit(0, "lines"),
            axis.text.y = element_text(size = 80, margin = unit(c(5, 5, 5, 5), "mm")),
            axis.text.x = element_text(size = 80, margin = unit(c(5, 5, 5, 5), "mm")),
            axis.ticks.x = element_blank(),

            legend.position = "right",
            legend.margin = margin(0,75,0,0),
            legend.key.width = unit(1.5, "inch"),
            legend.key.height = unit(2, "inch"),
            legend.title = element_text(size = 75, margin = margin(25,0,0,0)),
            legend.title.align = 0.5,
            legend.text = element_text(size = 60, margin = margin(0,0,0,25))
         ) +
      geom_point(data = qq_data, aes(x=x, y=y), color = "orange", size = 6)

      png(fig_scatter_filename, width = 3000, height = 3000)
      plot(p1)
      dev.off()
      system(paste("okular",fig_scatter_filename,"&> /dev/null &"))

   } else {
#----------------------------------------------------------------------------------#
# Open file.
   png(fig_scatter_filename, width = 3000, height = 1500*length(buoy_list))
   par(mfrow=pl_mfrow,oma=pl_oma,mar=pl_mar,mgp=pl_mgp)
# Plot multiple panels, for different buoys.
   for (p_idx in 1:length(vec_plot_idx)) {
      plot_idx <- vec_plot_idx[p_idx]
# Data for buoy.
      df_reg_buoy <- df_reg[df_reg$buoy_lab==buoy_list[p_idx],]
# Find paired data (in case of asymettrically missing data).
      Lvec_hs_paired <- !is.na(df_reg_buoy$buoy_hs) & !is.na(df_reg_buoy$ww3_hs)
      vec_buoy_hs_paired <- df_reg_buoy$buoy_hs[Lvec_hs_paired]
      vec_ww3_hs_paired <- df_reg_buoy$ww3_hs[Lvec_hs_paired]
      Lvec_tm_paired <- !is.na(df_reg_buoy$buoy_tm) & !is.na(df_reg_buoy$ww3_tm)
      vec_buoy_tm_paired <- df_reg_buoy$buoy_hs[Lvec_tm_paired]
      vec_ww3_tm_paired <- df_reg_buoy$ww3_hs[Lvec_tm_paired]
# Regression.
      #lm_hs <- eval(parse(text=paste("lm(",vec_data_sets[plot_idx]," ~ buoy_hs,data=df_reg)",sep="")))
      lm_hs <- lm(ww3_hs ~ buoy_hs,data=df_reg_buoy)
# RMSE
      hs_rmse <- sqrt(mean(lm_hs$residuals^2))
      #hs_rmse1 <- vec_sat-df_reg$buoy_hs
# RMSD
      hs_rmsd <- sqrt( mean( ( df_reg_buoy$ww3_hs - df_reg_buoy$buoy_hs )^2 ,na.rm=T ) )
# Correlation.
      hs_cor <- cor(df_reg_buoy$buoy_hs,df_reg_buoy$ww3_hs,use="pairwise.complete.obs")
# Bias.
      hs_bias <- mean(vec_ww3_hs_paired) - mean(vec_buoy_hs_paired)
      #hs_high_bias <- quantile(vec_ww3_hs_paired,probs=high_bias) - quantile(vec_buoy_hs_paired,probs=high_bias)

# Plot Hs data.
      plot(df_reg_buoy$buoy_hs, df_reg_buoy$ww3_hs, xlim=c(0,10), ylim=c(0,10), main=paste0(buoy_list[p_idx],": Hs"), xlab="Buoy", ylab="WW3", pch=19, cex=pl_cex, cex.main=pl_cex_main, cex.lab=pl_cex_lab, cex.axis=pl_cex_axis)
      #points(df_reg$buoy_hs[df_reg$track_dist > 25], vec_sat[df_reg$track_dist > 25], pch=19, cex=(pl_cex-1), col="red")

      #if ( any( plot_idx == c(2,4) ) ) {
      #   #text(df_reg$buoy_hs[hs_rmse1 > 1], vec_sat[hs_rmse1 > 1]+0.4, labels=df_reg$buoy_lab[hs_rmse1 > 1],cex=4)
      #   text(df_reg$buoy_hs[hs_rmse1 < -2], vec_sat[hs_rmse1 < -2]+0.4, labels=df_reg$buoy_lab[hs_rmse1 < -2],cex=4)
      #}

      #segments(x0=mean(vec_buoy_hs_paired), y0=0, y1=mean(vec_ww3_hs_paired), col="red", lwd=4)
      #segments(x0=0, x1=mean(vec_buoy_hs_paired), y0=mean(vec_ww3_hs_paired), col="red", lwd=4)
      #points(mean(vec_buoy_hs_paired), mean(vec_ww3_hs_paired,na.rm=T), col="red", pch=19)
      par(new=T)
      qqplot(vec_buoy_hs_paired,vec_ww3_hs_paired,xlim=c(0,10),ylim=c(0,10),pch=19,cex=3,col="orange",axes=F,xlab="",ylab="")
      abline(0, 1, col="blue", lwd=5.0)

      sat_mean <- mean(vec_ww3_hs_paired)
      mtext(side=3, line=-6, adj=0.03, cex=cex_mtext, outer=FALSE, text=paste("Mean bias (Sat mean = ",format(sat_mean,digits=2),"): ",format(hs_bias,digits=2),sep=''))
      #sat_Q95 <- quantile(vec_ww3_hs_paired,probs=high_bias)
      #mtext(side=3, line=-6, adj=0.03, cex=cex_mtext, outer=FALSE, text=paste("Q95 bias (Sat Q95 = ",format(sat_Q95,digits=2),"): ",format(hs_high_bias,digits=2),sep=''))
      mtext(side=3, line=-12, adj=0.03, cex=cex_mtext, outer=FALSE, text=paste("Correlation: ",format(hs_cor,digits=3),sep=''))
      mtext(side=3, line=-18, adj=0.03, cex=cex_mtext, outer=FALSE, text=paste("RMSD: ",format(round(hs_rmsd,3),nsmall=2),sep=''))
      mtext(side=3, line=-24, adj=0.03, cex=cex_mtext, outer=FALSE, text=paste("RMSE: ",format(round(hs_rmse,3),nsmall=2),sep=''))
      mtext(side=3, line=-30, adj=0.03, cex=cex_mtext, outer=FALSE, text=paste("N:",sum(!is.na(Lvec_hs_paired))))
      #mtext(side=3, line=-36, adj=0.03, cex=cex_mtext, outer=FALSE, text=paste("Total 1 Hz:",vec_1Hz_count[plot_idx]))
# Legend for data and Q-Q plot.
      if ( plot_idx == 1 ) { legend(x=5.5,y=2.0,legend=c("Hs","Q-Q plot"),col=c("black","orange"),pch=c(19,19),cex=pl_cex_leg) }

# Plot Tm data.
      plot(df_reg_buoy$buoy_tm, df_reg_buoy$ww3_tm, xlim=c(0,15), ylim=c(0,15), main=paste0(buoy_list[p_idx],": Tm"), xlab="Buoy", ylab="WW3", pch=19, cex=pl_cex, cex.main=pl_cex_main, cex.lab=pl_cex_lab, cex.axis=pl_cex_axis)
      abline(0, 1, col="blue", lwd=5.0)
   }

   dev.off()
   system(paste("okular",fig_scatter_filename,"&> /dev/null &"))

   }

##=======================================================================#
### Plotting.
###-----------------------------------------------------------------------#
### 2.a Plot wave stats (Hs, Tp).
###-----------------------------------------------------------------------#
### Simualtion length (days).
##   sim_len <- 31
##
### Time steps per day in NetCDF.
##   ncdf_step <- 4
##
##   X11()
##   #pdf("/home/ben/research/code/R/ww3/global/exp_6_ecmwf_2011_JAN_025_loSpec/hs_ecmwf_025_loSpec_buoys.pdf", paper="USr",width=12, height=12)
##   par(mfrow=c(2,2),oma=c(0,0,0,0),mgp=c(3,1,0),mar=c(4,4,4,5))
##
###   for (b.idx in 1)
##   for (b.idx in 1:length(buoy_list))
##   {
##      plot_title <- paste("Simulated and observed wave statistics at NDBC", buoy_list[b.idx])
##      plot(NULL,xlim=wind_time[c(1,(sim_len*ncdf_step))],ylim=c(-9,7),xaxt='n',yaxt='n',xlab=NA,ylab=NA,main=plot_title)
##      axis(side=1,at=time_labels,labels=as.POSIXct(time_labels*3600, origin = '1900-01-01', tz='GMT'))
##      abline(v=time_labels,lty=4)
##      axis(side=2,at=seq(0,5),labels=seq(0,5))
##      mtext(expression(H[s]~(m)),side=2,line=2)
##      abline(h=0,lty=1)
##      abline(h=c(2.5,5),lty=2)
##
### Find the required date string (from WW3 output data.)
##      ww3_year <- substr(array_ww[1,1,b.idx],start=1,stop=4)
##      ww3_month <- substr(array_ww[1,1,b.idx],start=5,stop=6)
##      ww3_day <- substr(array_ww[1,1,b.idx],start=7,stop=8)
##      ww3_hour <- array_ww[1,2,b.idx]
##
### Start a counter.
##      test_ww3_date <- strptime(paste(ww3_year,"-",ww3_month,"-",ww3_day,sep=""),format="%Y-%m-%d",tz="GMT")
##print(paste("WW3 date:",test_ww3_date))
##      t_time <- 1
##      test_obs_date <- strptime(paste(array_buoy_obs[t_time,1,b.idx],"-",array_buoy_obs[t_time,2,b.idx],"-",array_buoy_obs[t_time,3,b.idx]," ",array_buoy_obs[t_time,4,b.idx],":00:00",sep=""),format="%Y-%m-%d %H:%M:%S",tz="GMT")
##print(paste("Obs date1:",test_obs_date))
##
##      while ( test_obs_date != test_ww3_date ) {
##         t_time <- t_time + 1
##         test_obs_date <- strptime(paste(array_buoy_obs[t_time,1,b.idx],"-",array_buoy_obs[t_time,2,b.idx],"-",array_buoy_obs[t_time,3,b.idx]," ",array_buoy_obs[t_time,4,b.idx],":00:00",sep=""),format="%Y-%m-%d %H:%M:%S",tz="GMT")
##            if ( b.idx == 3 ) {
##	          t_time <- 3201
##                  test_obs_date <- 1
##                  test_ww3_date <- 1
##               }
##          }
##print(paste("Obs date2:",test_obs_date))
##
### Observed data.
##      obs_plot_seq <- seq(wind_time[5],,(1/buoy_obs_freq[b.idx]),(sim_len*24*buoy_obs_freq[b.idx]))
##      obs_plot_data <- array_buoy_obs[t_time:(t_time - 1 + sim_len*24*buoy_obs_freq[b.idx]),9,b.idx]
##      print(paste("T_Time:",t_time))
##
##      points(obs_plot_seq, obs_plot_data, pch=1, cex=0.8, col="orange")
### Print diagnostics.
###      print(paste("Buoy:", buoy_list[b.idx],"Start time:",as.POSIXct(min(obs_plot_seq)*3600, origin = '1900-01-01', tz='GMT')))
###      print(paste("Buoy:", buoy_list[b.idx],"Stop time:",as.POSIXct(max(obs_plot_seq)*3600, origin = '1900-01-01', tz='GMT')))
###      print(paste("Data:", array_buoy_obs[1:2,,b.idx]))
###      print(paste("Data:", array_buoy_obs[1:2,,b.idx]))
### WW3 Hs data.
###   for (i in seq(4,,4,20))
###   {
###      lines(wind_time[sim_start_idx:sim_stop_idx],mat_ww_41010_sample[,5,i],col="red",lwd=0.5)
###   }
##      #lines(wind_time[sim_start_idx:sim_stop_idx], mat_ww_41010[seq(1,,2,53),5],col="darkred",lwd=2.5)
##      lines(seq(wind_time[5],,3,(sim_len*8)), array_ww[1:(sim_len*8),5,b.idx],col="black",lwd=2.5)
##   }
### Highlight points.
###   points(wind_time[770], mat_ww_41010[32,5],col="darkred",pch=4,cex=3.0)
###   points(wind_time[770], mat_ww_41010[32,5],col="darkred",pch=1,cex=2.0)
### Boxplot of MC data.
###   Hs_MC_output <- read.table("../R_Hs_max/Hs_MC_output")
###   boxplot(Hs_MC_output,add=TRUE,at=wind_time[770],pars = list(boxwex = 12000),yaxt='n',outline=TRUE)
##
#### Plot Tp data.
###   par(new=TRUE)
###   plot(NULL,xlim=wind_time[mat_xlim],ylim=c(2.0,20),xlab=NA,ylab=NA,xaxt='n',yaxt='n')
###   mtext(expression(T[p]~(s)),side=4,line=3)
###   axis(side=4,at=seq(3,9))
###   abline(h=5,lty=2)
###   points(seq(wind_time[1],,(1/buoy_obs_freq[b.idx]),(14*24*buoy_obs_freq[b.idx])),array_buoy_obs[1:(14*24*buoy_obs_freq[b.idx]),10,b.idx],pch=3,cex=0.8,col="darkblue")
#### WW3 Tp data.
####   for (i in seq(4,,4,20))
####   {
####      lines(wind_time[sim_start_idx:sim_stop_idx],1/mat_ww_41010_sample[,10,i],col="blue",lwd=0.5)
####   }
###   lines(wind_time[sim_start_idx:sim_stop_idx],1/array_ww[seq(1,,2,53),10,b.idx],col="darkblue",lwd=2.5)
#### Boxplot of MC data.
####   Tp_MC_output <- read.table("../R_Tp_max/Tp_MC_output")
####   boxplot(Tp_MC_output,add=TRUE,at=wind_time[770],pars = list(boxwex = 12000),yaxt='n')
###
### Legend.
###      legend(
###             x=wind_time[mat_xlim[1]],y=20,
###             c(expression(H[s]),expression(WW3~simulated~H[s]),expression(T[p]),expression(WW3~simulated~T[p])),
###             lty=c(-1,1,-1,1),lwd=c(-1,2,-1,2),pch=c(1,-1,3,-1),col=c("darkred","darkred","darkblue","darkblue"),bty="o",bg="white",
###             x.intersp=0.7,y.intersp=1.0,cex=1.1,xjust=0.1
###            )
##
###   }
###dev.off()
##
###-----------------------------------------------------------------------#
### 2.b Scatter plot of obs. vs simulated data.
###-----------------------------------------------------------------------#
##   X11()
##   #pdf("/home/ben/research/code/R/ww3/global/exp_6_ecmwf_2011_JAN_025_loSpec/scatter_ecmwf_025_loSpec_buoys.pdf", paper="USr",width=12, height=12)
##   par(mfrow=c(2,2),oma=c(0,0,0,0),mgp=c(3,1,0),mar=c(4,4,4,5)+0.1)
##
##   for (b.idx in 1:length(buoy_list)) {
##         plot(array_buoy_obs[seq(1,,(3*buoy_obs_freq[b.idx]),mat_ww_dims[1,1]),9,b.idx],array_ww[1:mat_ww_dims[1,1],5,b.idx],xlim=c(0,4),ylim=c(0,4),xlab=paste("Observation at NDBC",buoy_list[b.idx]),ylab="WW3 data")
##         abline(a=0,b=1)
##      }
###   dev.off()
##
###   X11()
###   plot(mat_ndbc_41010[scatter_range_41010,10],1/mat_ww_41010[,10],xlim=c(2,8),ylim=c(2,8),xlab="Observation at NDBC 41010",ylab="WW3 data")
###   abline(a=0,b=1)
##
