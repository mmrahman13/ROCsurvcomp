########################********************##########################
################ -------------------------- ##################
###################### Double censoring ######################
################ -------------------------- ##################
########################********************##########################

######################################################################
####### Functions for imputation with Double censoring ########
######################################################################

######################################################################
####### Find the number of LEFT censored obs before 1st event time OR 1st right censored

Calc_number_LEFTcen_before_1steventORrightcen <- function(cens_stat) {

  # First occurrence of event (0)
  min_index_event <- which(cens_stat == 0)[1]

  # First occurrence of right censoring (1)
  min_index_rightcen <- which(cens_stat == 1)[1]

  # Determine earliest of event or right censoring
  if (length(min_index_rightcen) == 0) {
    min_index_eventORrightcen <- min_index_event
  } else {
    min_index_eventORrightcen <- min(min_index_event, min_index_rightcen)
  }

  # Number of left-censored (-1) before 1st event time OR 1st right censored
  n_first_leftcen <- min_index_eventORrightcen - 1

  return(n_first_leftcen)
}


######################################################################
####### Find the number of RIGHT censored obs after last event time OR last right censored
Calc_number_RIGHTcen_after_lasteventORleftcen <- function(cens_stat) {

  # Last occurrence of event (0)
  max_index_event <- tail(which(cens_stat == 0), 1)

  # Last occurrence of left censoring (-1)
  max_index_leftcen <- tail(which(cens_stat == -1), 1)

  # Determine the last index between event or left censoring
  if (length(max_index_leftcen) == 0) {
    max_index_eventORleftcen <- max_index_event
  } else {
    max_index_eventORleftcen <- max(max_index_event, max_index_leftcen)
  }

  # Number of right-censored (1) after after last event time OR last right censored
  n <- length(cens_stat)
  n_last_rightcen <- n - max_index_eventORleftcen

  return(n_last_rightcen)
}



######################################################################
############## Our Method with double censoring ###############
####### ------------- Only ROC Length ------------- ##########
######################################################################

roc_DoubleCenSurvival_test <- function(time1, censor1, time2, censor2, boots, progress = TRUE) {

  ##### -------------------------------------------------------------------
  ##### censoring status (0 = event, 1 = right censored, -1 = left censored)

  # double censored data for trt and ctrl
  time_tt = time1
  cens_stat_t = censor1
  n_trt = length(time_tt)

  time_cc = time2
  cens_stat_c = censor2
  n_ctrl = length(time_cc)

  # Log transform
  time_t = log(1 + time_tt)
  time_c = log(1 + time_cc)

  # Sort by time
  idx_t = order(time_t)
  time_trt = time_t[idx_t]
  cens_stat_trt = cens_stat_t[idx_t]

  idx_c = order(time_c)
  time_ctrl = time_c[idx_c]
  cens_stat_ctrl = cens_stat_c[idx_c]

  # Box-Cox (with censoring) to estimate parameters and lambda
  res <- boxcox_censor(time_trt, time_ctrl, cens_stat_trt, cens_stat_ctrl)
  g = res$g
  exitflag = res$exitflag
  gval = res$gval

  lambda = g[5]
  mu_trt_hat  = g[1];  sd_trt_hat  = g[2]
  mu_ctrl_hat = g[3];  sd_ctrl_hat = g[4]

  # Transform data to Box-Cox scale
  time_trt_box  = ((time_trt^lambda)-1) / lambda
  time_ctrl_box = ((time_ctrl^lambda)-1) / lambda


  ###########################################################
  ################# Censored obs imputation #################

  ########### Trt arm
  # Find the number of RIGHT censored obs after last event time OR last right censored
  n_last_censor_trt = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_trt);
  # Find the number of LEFT censored obs before 1st event time OR 1st right censored
  n_first_leftcen_trt = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_trt);

  ########### Ctrl arm
  # Find the number of RIGHT censored obs after last event time OR last right censored
  n_last_censor_ctrl = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_ctrl);
  # Find the number of LEFT censored obs before 1st event time OR 1st right censored
  n_first_leftcen_ctrl = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_ctrl);


  #### Check if both arm's 1st & last obs are event

  if (n_last_censor_trt == 0 && n_first_leftcen_trt  == 0 && n_last_censor_ctrl == 0 && n_first_leftcen_ctrl == 0) {

    time_trt_box_imput = time_trt_box
    time_ctrl_box_imput = time_ctrl_box
    cens_stat_trt_imput = cens_stat_trt
    cens_stat_ctrl_imput = cens_stat_ctrl

    # Compute ROC length
    f_trt <- function(x) kde_double_censored(x, time_trt_box_imput, cens_stat_trt_imput)
    f_ctrl <- function(x) kde_double_censored(x, time_ctrl_box_imput, cens_stat_ctrl_imput)
    integrand <- function(x) sqrt(f_trt(x)^2 + f_ctrl(x)^2)

    length_roc = integrate(integrand, -Inf, Inf)$value
    length_roc_diff = length_roc - sqrt(2)

  } else {

    length_roc_sample = numeric(30)

    #### Otherwise do imputation
    for (imput in 1:30) {

      time_trt_box_imput = time_trt_box
      time_ctrl_box_imput = time_ctrl_box
      cens_stat_trt_imput = cens_stat_trt
      cens_stat_ctrl_imput = cens_stat_ctrl

      impute_sample1_right = numeric(n_last_censor_trt)
      impute_sample1_left = numeric(n_first_leftcen_trt)
      impute_sample2_right = numeric(n_last_censor_ctrl)
      impute_sample2_left = numeric(n_first_leftcen_ctrl)


      ########### Trt arm ##########
      ####### Censor in RIGHT tail
      if (n_last_censor_trt == 0) {
        time_trt_box_imput = time_trt_box_imput
        cens_stat_trt_imput = cens_stat_trt_imput
      } else {
        # Impute the censored obs after last event time/left censored
        for (i_t in 1:n_last_censor_trt) {
          sample1_R = time_trt_box[n_trt - n_last_censor_trt + i_t]   # Initial sample
          while (sample1_R <= time_trt_box[n_trt - n_last_censor_trt + i_t]) {
            sample1_R = rnorm(1, mean = mu_trt_hat, sd = sd_trt_hat)
          }
          impute_sample1_right[i_t] = sample1_R
        }
        # Substitute the obs in the main data
        time_trt_box_imput[(n_trt - n_last_censor_trt + 1):n_trt] = impute_sample1_right
        cens_stat_trt_imput[(n_trt - n_last_censor_trt + 1):n_trt] = 0
      }


      ####### Censor in LEFT tail
      if (n_first_leftcen_trt == 0) {
        time_trt_box_imput = time_trt_box_imput
        cens_stat_trt_imput = cens_stat_trt_imput
      } else {
        # Impute the censored obs before first event time/right censored
        for (i_t in 1:n_first_leftcen_trt) {
          sample1_L = time_trt_box[i_t]   # Initial sample
          while (sample1_L >= time_trt_box[i_t]) {
            sample1_L = rnorm(1, mean = mu_trt_hat, sd = sd_trt_hat)
          }
          impute_sample1_left[i_t] = sample1_L
        }

        # Substitute the obs in the main data
        time_trt_box_imput[1:n_first_leftcen_trt] = impute_sample1_left
        cens_stat_trt_imput[1:n_first_leftcen_trt] = 0
      }


      ########### Ctrl arm ##########
      ####### Censor in RIGHT tail
      if (n_last_censor_ctrl == 0) {
        time_ctrl_box_imput = time_ctrl_box_imput
        cens_stat_ctrl_imput = cens_stat_ctrl_imput
      } else {
        # Impute the censored obs after last event time/left censored
        for (i_c in 1:n_last_censor_ctrl) {
          sample2_R = time_ctrl_box[n_ctrl - n_last_censor_ctrl + i_c]   # Initial sample
          while (sample2_R <= time_ctrl_box[n_ctrl - n_last_censor_ctrl + i_c]) {
            sample2_R = rnorm(1, mean = mu_ctrl_hat, sd = sd_ctrl_hat)
          }
          impute_sample2_right[i_c] = sample2_R
        }

        # Substitute the obs in the main data
        time_ctrl_box_imput[(n_ctrl - n_last_censor_ctrl + 1):n_ctrl] = impute_sample2_right
        cens_stat_ctrl_imput[(n_ctrl - n_last_censor_ctrl + 1):n_ctrl] = 0
      }


      ####### Censor in LEFT tail
      if (n_first_leftcen_ctrl == 0) {
        time_ctrl_box_imput = time_ctrl_box_imput
        cens_stat_ctrl_imput = cens_stat_ctrl_imput
      } else {
        # Impute the censored obs before first event time/right censored
        for (i_c in 1:n_first_leftcen_ctrl) {
          sample2_L = time_ctrl_box[i_c]   # Initial sample
          while (sample2_L >= time_ctrl_box[i_c]) {
            sample2_L = rnorm(1, mean = mu_ctrl_hat, sd = sd_ctrl_hat)
          }
          impute_sample2_left[i_c] = sample2_L
        }

        # Substitute the obs in the main data
        time_ctrl_box_imput[1:n_first_leftcen_ctrl] = impute_sample2_left
        cens_stat_ctrl_imput[1:n_first_leftcen_ctrl] = 0
      }


      # Compute ROC length
      f_trt <- function(x) kde_double_censored(x, time_trt_box_imput, cens_stat_trt_imput)
      f_ctrl <- function(x) kde_double_censored(x, time_ctrl_box_imput, cens_stat_ctrl_imput)
      integrand <- function(x) sqrt(f_trt(x)^2 + f_ctrl(x)^2)

      length_roc_sample[imput] = integrate(integrand, -Inf, Inf)$value
    }

    length_roc = mean(length_roc_sample)
    length_roc_diff = length_roc - sqrt(2)
  }


  ########################################
  #################### Permutation
  ########################################

  time = c(time_t, time_c)    # Merge unsorted values
  cens = c(cens_stat_t, cens_stat_c)
  n = n_trt + n_ctrl

  length_roc_boot_diff = numeric(boots)

  show_progress <- progress && interactive()
  if (show_progress) {
    cat("Permutation test in progress...\n")
    bar_width <- 50
  }

  for (b in 1:boots) {

    # ----------- Progress bar update ----------

    if (show_progress && (b %% 10 == 0 || b == boots)) {
      prog <- b / boots
      n_stars <- floor(prog * bar_width)
      n_spaces <- bar_width - n_stars

      bar <- paste0(
        "|",
        paste(rep("*", n_stars), collapse = ""),
        paste(rep(" ", n_spaces), collapse = ""),
        "| ",
        sprintf("%3d%%", floor(prog * 100))
      )

      cat("\r", bar, sep = "")
      if (b == boots) cat("\n")  # move to next line at end
      flush.console()
    }

    # ------------------------------------------


    # Generate permutation sample from original sample
    at = sample(1:n, n, replace = FALSE)
    time_t_boot = time[at[1:n_trt]]
    cens_stat_t_boot = cens[at[1:n_trt]]
    time_c_boot = time[at[(n_trt + 1):n]]
    cens_stat_c_boot = cens[at[(n_trt + 1):n]]

    # Sort based on time values
    data1_boot = cbind(time_t_boot, cens_stat_t_boot)
    data_sorted1_boot = data1_boot[order(data1_boot[,1]),]
    time_trt_boot = data_sorted1_boot[, 1]
    cens_stat_trt_boot = data_sorted1_boot[, 2]

    data2_boot = cbind(time_c_boot, cens_stat_c_boot)
    data_sorted2_boot = data2_boot[order(data2_boot[,1]),]
    time_ctrl_boot = data_sorted2_boot[, 1]
    cens_stat_ctrl_boot = data_sorted2_boot[, 2]

    # Perform Boxcox
    res <- boxcox_censor(time_trt_boot,time_ctrl_boot,cens_stat_trt_boot,cens_stat_ctrl_boot)
    g = res$g
    exitflag = res$exitflag
    gval = res$gval

    lambda_boot = g[5]
    mu_trt_hat_boot = g[1];  sd_trt_hat_boot = g[2]
    mu_ctrl_hat_boot = g[3];  sd_ctrl_hat_boot = g[4]

    # Transform the dataset
    time_trt_box_boot = ((time_trt_boot^lambda_boot)-1)/lambda_boot
    time_ctrl_box_boot = ((time_ctrl_boot^lambda_boot)-1)/lambda_boot


    ###########################################################
    ################# Censored obs imputation #################

    ########### Trt arm
    # Find the number of RIGHT censored obs after last event time OR last right censored
    n_last_censor_trt_boot = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_trt_boot)
    # Find the number of LEFT censored obs before 1st event time OR 1st right censored
    n_first_leftcen_trt_boot = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_trt_boot)

    ########### Ctrl arm
    # Find the number of RIGHT censored obs after last event time OR last right censored
    n_last_censor_ctrl_boot = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_ctrl_boot)
    # Find the number of LEFT censored obs before 1st event time OR 1st right censored
    n_first_leftcen_ctrl_boot = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_ctrl_boot)


    #### Check if both arm's last obs is event
    if (n_last_censor_trt_boot == 0 && n_first_leftcen_trt_boot == 0 && n_last_censor_ctrl_boot == 0 && n_first_leftcen_ctrl_boot == 0) {

      time_trt_box_imput_boot = time_trt_box_boot
      time_ctrl_box_imput_boot = time_ctrl_box_boot
      cens_stat_trt_imput_boot = cens_stat_trt_boot
      cens_stat_ctrl_imput_boot = cens_stat_ctrl_boot

      # Compute ROC length
      f_trt_boot <- function(x) kde_double_censored(x, time_trt_box_imput_boot, cens_stat_trt_imput_boot)
      f_ctrl_boot <- function(x) kde_double_censored(x, time_ctrl_box_imput_boot, cens_stat_ctrl_imput_boot)
      integrand <- function(x) sqrt(f_trt_boot(x)^2 + f_ctrl_boot(x)^2)

      length_roc_boot = integrate(integrand, -Inf, Inf)$value
      length_roc_boot_diff[b] = length_roc_boot - sqrt(2)

    } else {

      length_roc_sample_boot = numeric(30)

      #### Otherwise do imputation
      for (imput_boot in 1:30) {

        time_trt_box_imput_boot = time_trt_box_boot
        time_ctrl_box_imput_boot = time_ctrl_box_boot
        cens_stat_trt_imput_boot = cens_stat_trt_boot
        cens_stat_ctrl_imput_boot = cens_stat_ctrl_boot

        impute_sample1_right_boot = numeric(n_last_censor_trt_boot)
        impute_sample1_left_boot = numeric(n_first_leftcen_trt_boot)
        impute_sample2_right_boot = numeric(n_last_censor_ctrl_boot)
        impute_sample2_left_boot = numeric(n_first_leftcen_ctrl_boot)


        ########### Trt arm ##########

        ####### Censor in RIGHT tail
        if (n_last_censor_trt_boot == 0) {
          time_trt_box_imput_boot = time_trt_box_imput_boot
          cens_stat_trt_imput_boot = cens_stat_trt_imput_boot
        } else {
          # Impute the censored obs after last event time/left censored
          for (i_t_boot in 1:n_last_censor_trt_boot) {
            sample1_R_boot = time_trt_box_boot[n_trt - n_last_censor_trt_boot + i_t_boot]   # Initial sample
            while (sample1_R_boot <= time_trt_box_boot[n_trt - n_last_censor_trt_boot + i_t_boot]) {
              sample1_R_boot = rnorm(1, mean = mu_trt_hat_boot, sd = sd_trt_hat_boot)
            }
            impute_sample1_right_boot[i_t_boot] = sample1_R_boot
          }
          # Substitute the obs in the main data
          time_trt_box_imput_boot[(n_trt - n_last_censor_trt_boot + 1):n_trt] = impute_sample1_right_boot
          cens_stat_trt_imput_boot[(n_trt - n_last_censor_trt_boot + 1):n_trt] = 0
        }



        ####### Censor in LEFT tail
        if (n_first_leftcen_trt_boot == 0) {
          time_trt_box_imput_boot = time_trt_box_imput_boot
          cens_stat_trt_imput_boot = cens_stat_trt_imput_boot
        } else {
          # Impute the censored obs before first event time/right censored
          for (i_t_boot in 1:n_first_leftcen_trt_boot) {
            sample1_L_boot = time_trt_box_boot[i_t_boot]   # Initial sample
            while (sample1_L_boot >= time_trt_box_boot[i_t_boot]) {
              sample1_L_boot = rnorm(1, mean = mu_trt_hat_boot, sd = sd_trt_hat_boot)
            }
            impute_sample1_left_boot[i_t_boot] = sample1_L_boot
          }

          # Substitute the obs in the main data
          time_trt_box_imput_boot[1:n_first_leftcen_trt_boot] = impute_sample1_left_boot
          cens_stat_trt_imput_boot[1:n_first_leftcen_trt_boot] = 0
        }


        ########### Ctrl arm ##########
        ####### Censor in RIGHT tail
        if (n_last_censor_ctrl_boot == 0) {
          time_ctrl_box_imput_boot = time_ctrl_box_imput_boot
          cens_stat_ctrl_imput_boot = cens_stat_ctrl_imput_boot
        } else {
          # Impute the censored obs after last event time/left censored
          for (i_c_boot in 1:n_last_censor_ctrl_boot) {
            sample2_R_boot = time_ctrl_box_boot[n_ctrl - n_last_censor_ctrl_boot + i_c_boot]   # Initial sample
            while (sample2_R_boot <= time_ctrl_box_boot[n_ctrl - n_last_censor_ctrl_boot + i_c_boot]) {
              sample2_R_boot = rnorm(1, mean = mu_ctrl_hat_boot, sd = sd_ctrl_hat_boot)
            }
            impute_sample2_right_boot[i_c_boot] = sample2_R_boot
          }
          # Substitute the obs in the main data
          time_ctrl_box_imput_boot[(n_ctrl - n_last_censor_ctrl_boot + 1):n_ctrl] = impute_sample2_right_boot
          cens_stat_ctrl_imput_boot[(n_ctrl - n_last_censor_ctrl_boot + 1):n_ctrl] = 0
        }


        ####### Censor in LEFT tail
        if (n_first_leftcen_ctrl_boot == 0) {
          time_ctrl_box_imput_boot = time_ctrl_box_imput_boot
          cens_stat_ctrl_imput_boot = cens_stat_ctrl_imput_boot
        } else {
          # Impute the censored obs before first event time/right censored
          for (i_c_boot in 1:n_first_leftcen_ctrl_boot) {
            sample2_L_boot = time_ctrl_box_boot[i_c_boot]   # Initial sample
            while (sample2_L_boot >= time_ctrl_box_boot[i_c_boot]) {
              sample2_L_boot = rnorm(1, mean = mu_ctrl_hat_boot, sd = sd_ctrl_hat_boot)
            }
            impute_sample2_left_boot[i_c_boot] = sample2_L_boot
          }

          # Substitute the obs in the main data
          time_ctrl_box_imput_boot[1:n_first_leftcen_ctrl_boot] = impute_sample2_left_boot
          cens_stat_ctrl_imput_boot[1:n_first_leftcen_ctrl_boot] = 0
        }


        # Compute ROC length
        f_trt_boot <- function(x) kde_double_censored(x, time_trt_box_imput_boot, cens_stat_trt_imput_boot)
        f_ctrl_boot <- function(x) kde_double_censored(x, time_ctrl_box_imput_boot, cens_stat_ctrl_imput_boot)
        integrand <- function(x) sqrt(f_trt_boot(x)^2 + f_ctrl_boot(x)^2)

        length_roc_sample_boot[imput_boot] = integrate(integrand, -Inf, Inf)$value
      }

      length_roc_boot = mean(length_roc_sample_boot)
      length_roc_boot_diff[b] = length_roc_boot - sqrt(2)
    }
  }


  # Calculate p-values
  prop_length = sum(length_roc_boot_diff > length_roc_diff)
  pval_length = prop_length / length(length_roc_boot_diff)

  ####### Result #######
  df_result <- data.frame(
    estimate = sprintf("%.4f", length_roc),
    p_value = pval_length
  )

  rownames(df_result) <- "ROC Length"

  # Message text
  message_text <- paste(
    "Hypothesis testing method: ROC length",
    "Significance level: 5% (two-sided)",
    paste0("Number of permutations used for inference: ", boots),
    sep = "\n"
  )

  # Create the list object & Assign a custom class name to this list
  out <- list(message = message_text, result = df_result)
  class(out) <- "survival_comp_roc"
  return(out)
}






######################################################################
############## Our Method with double censoring ###############
####### ------------ Only OVL ------------- #######
######################################################################

ovl_DoubleCenSurvival_test <- function(time1, censor1, time2, censor2, boots, progress = TRUE) {

  ##### -------------------------------------------------------------------
  ##### censoring status (0 = event, 1 = right censored, -1 = left censored)

  # double censored data for trt and ctrl
  time_tt = time1
  cens_stat_t = censor1
  n_trt = length(time_tt)

  time_cc = time2
  cens_stat_c = censor2
  n_ctrl = length(time_cc)

  # Log transform
  time_t = log(1 + time_tt)
  time_c = log(1 + time_cc)

  # Sort by time
  idx_t = order(time_t)
  time_trt = time_t[idx_t]
  cens_stat_trt = cens_stat_t[idx_t]

  idx_c = order(time_c)
  time_ctrl = time_c[idx_c]
  cens_stat_ctrl = cens_stat_c[idx_c]

  # Box-Cox (with censoring) to estimate parameters and lambda
  res <- boxcox_censor(time_trt, time_ctrl, cens_stat_trt, cens_stat_ctrl)
  g = res$g
  exitflag = res$exitflag
  gval = res$gval

  lambda = g[5]
  mu_trt_hat  = g[1];  sd_trt_hat  = g[2]
  mu_ctrl_hat = g[3];  sd_ctrl_hat = g[4]

  # Transform data to Box-Cox scale
  time_trt_box  = ((time_trt^lambda)-1) / lambda
  time_ctrl_box = ((time_ctrl^lambda)-1) / lambda


  ###########################################################
  ################# Censored obs imputation #################

  ########### Trt arm
  # Find the number of RIGHT censored obs after last event time OR last right censored
  n_last_censor_trt = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_trt);
  # Find the number of LEFT censored obs before 1st event time OR 1st right censored
  n_first_leftcen_trt = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_trt);

  ########### Ctrl arm
  # Find the number of RIGHT censored obs after last event time OR last right censored
  n_last_censor_ctrl = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_ctrl);
  # Find the number of LEFT censored obs before 1st event time OR 1st right censored
  n_first_leftcen_ctrl = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_ctrl);


  #### Check if both arm's last obs is event

  if (n_last_censor_trt == 0 && n_first_leftcen_trt  == 0 && n_last_censor_ctrl == 0 && n_first_leftcen_ctrl == 0) {

    time_trt_box_imput = time_trt_box
    time_ctrl_box_imput = time_ctrl_box
    cens_stat_trt_imput = cens_stat_trt
    cens_stat_ctrl_imput = cens_stat_ctrl

    # Compute ROC length
    f_trt <- function(x) kde_double_censored(x, time_trt_box_imput, cens_stat_trt_imput)
    f_ctrl <- function(x) kde_double_censored(x, time_ctrl_box_imput, cens_stat_ctrl_imput)

    # Compute Overlap
    integrand1 <- function(x) pmin(f_trt(x), f_ctrl(x))
    overlap = integrate(integrand1, -Inf, Inf)$value
    overlap_inv = 1 + (-overlap)

  } else {

    overlap_sample = numeric(30)

    for (imput in 1:30) {

      time_trt_box_imput = time_trt_box
      time_ctrl_box_imput = time_ctrl_box
      cens_stat_trt_imput = cens_stat_trt
      cens_stat_ctrl_imput = cens_stat_ctrl

      impute_sample1_right = numeric(n_last_censor_trt)
      impute_sample1_left = numeric(n_first_leftcen_trt)
      impute_sample2_right = numeric(n_last_censor_ctrl)
      impute_sample2_left = numeric(n_first_leftcen_ctrl)


      ########### Trt arm ##########
      ####### Censor in RIGHT tail
      if (n_last_censor_trt == 0) {
        time_trt_box_imput = time_trt_box_imput
        cens_stat_trt_imput = cens_stat_trt_imput
      } else {
        # Impute the censored obs after last event time/left censored
        for (i_t in 1:n_last_censor_trt) {
          sample1_R = time_trt_box[n_trt - n_last_censor_trt + i_t]   # Initial sample
          while (sample1_R <= time_trt_box[n_trt - n_last_censor_trt + i_t]) {
            sample1_R = rnorm(1, mean = mu_trt_hat, sd = sd_trt_hat)
          }
          impute_sample1_right[i_t] = sample1_R
        }
        # Substitute the obs in the main data
        time_trt_box_imput[(n_trt - n_last_censor_trt + 1):n_trt] = impute_sample1_right
        cens_stat_trt_imput[(n_trt - n_last_censor_trt + 1):n_trt] = 0
      }


      ####### Censor in LEFT tail
      if (n_first_leftcen_trt == 0) {
        time_trt_box_imput = time_trt_box_imput
        cens_stat_trt_imput = cens_stat_trt_imput
      } else {
        # Impute the censored obs before first event time/right censored
        for (i_t in 1:n_first_leftcen_trt) {
          sample1_L = time_trt_box[i_t]   # Initial sample
          while (sample1_L >= time_trt_box[i_t]) {
            sample1_L = rnorm(1, mean = mu_trt_hat, sd = sd_trt_hat)
          }
          impute_sample1_left[i_t] = sample1_L
        }

        # Substitute the obs in the main data
        time_trt_box_imput[1:n_first_leftcen_trt] = impute_sample1_left
        cens_stat_trt_imput[1:n_first_leftcen_trt] = 0
      }


      ########### Ctrl arm ##########
      ####### Censor in RIGHT tail
      if (n_last_censor_ctrl == 0) {
        time_ctrl_box_imput = time_ctrl_box_imput
        cens_stat_ctrl_imput = cens_stat_ctrl_imput
      } else {
        # Impute the censored obs after last event time/left censored
        for (i_c in 1:n_last_censor_ctrl) {
          sample2_R = time_ctrl_box[n_ctrl - n_last_censor_ctrl + i_c]   # Initial sample
          while (sample2_R <= time_ctrl_box[n_ctrl - n_last_censor_ctrl + i_c]) {
            sample2_R = rnorm(1, mean = mu_ctrl_hat, sd = sd_ctrl_hat)
          }
          impute_sample2_right[i_c] = sample2_R
        }

        # Substitute the obs in the main data
        time_ctrl_box_imput[(n_ctrl - n_last_censor_ctrl + 1):n_ctrl] = impute_sample2_right
        cens_stat_ctrl_imput[(n_ctrl - n_last_censor_ctrl + 1):n_ctrl] = 0
      }


      ####### Censor in LEFT tail
      if (n_first_leftcen_ctrl == 0) {
        time_ctrl_box_imput = time_ctrl_box_imput
        cens_stat_ctrl_imput = cens_stat_ctrl_imput
      } else {
        # Impute the censored obs before first event time/right censored
        for (i_c in 1:n_first_leftcen_ctrl) {
          sample2_L = time_ctrl_box[i_c]   # Initial sample
          while (sample2_L >= time_ctrl_box[i_c]) {
            sample2_L = rnorm(1, mean = mu_ctrl_hat, sd = sd_ctrl_hat)
          }
          impute_sample2_left[i_c] = sample2_L
        }

        # Substitute the obs in the main data
        time_ctrl_box_imput[1:n_first_leftcen_ctrl] = impute_sample2_left
        cens_stat_ctrl_imput[1:n_first_leftcen_ctrl] = 0
      }

      # Compute Overlap
      f_trt <- function(x) kde_double_censored(x, time_trt_box_imput, cens_stat_trt_imput)
      f_ctrl <- function(x) kde_double_censored(x, time_ctrl_box_imput, cens_stat_ctrl_imput)

      integrand1 <- function(x) pmin(f_trt(x), f_ctrl(x))
      overlap_sample[imput] = integrate(integrand1, -Inf, Inf)$value
    }

    overlap = mean(overlap_sample)
    overlap_inv = 1 - overlap
  }


  ########################################
  #################### Permutation
  ########################################

  time = c(time_t, time_c)    # Merge unsorted values
  cens = c(cens_stat_t, cens_stat_c)
  n = n_trt + n_ctrl

  overlap_inv_boot = numeric(boots)

  show_progress <- progress && interactive()
  if (show_progress) {
    cat("Permutation test in progress...\n")
    bar_width <- 50
  }

  for (b in 1:boots) {

    # ----------- Progress bar update ----------

    if (show_progress && (b %% 10 == 0 || b == boots)) {
      prog <- b / boots
      n_stars <- floor(prog * bar_width)
      n_spaces <- bar_width - n_stars

      bar <- paste0(
        "|",
        paste(rep("*", n_stars), collapse = ""),
        paste(rep(" ", n_spaces), collapse = ""),
        "| ",
        sprintf("%3d%%", floor(prog * 100))
      )

      cat("\r", bar, sep = "")
      if (b == boots) cat("\n")  # move to next line at end
      flush.console()
    }

    # ------------------------------------------


    # Generate permutation sample from original sample
    at = sample(1:n, n, replace = FALSE)
    time_t_boot = time[at[1:n_trt]]
    cens_stat_t_boot = cens[at[1:n_trt]]
    time_c_boot = time[at[(n_trt + 1):n]]
    cens_stat_c_boot = cens[at[(n_trt + 1):n]]

    # Sort based on time values
    data1_boot = cbind(time_t_boot, cens_stat_t_boot)
    data_sorted1_boot = data1_boot[order(data1_boot[,1]),]
    time_trt_boot = data_sorted1_boot[, 1]
    cens_stat_trt_boot = data_sorted1_boot[, 2]

    data2_boot = cbind(time_c_boot, cens_stat_c_boot)
    data_sorted2_boot = data2_boot[order(data2_boot[,1]),]
    time_ctrl_boot = data_sorted2_boot[, 1]
    cens_stat_ctrl_boot = data_sorted2_boot[, 2]

    # Perform Boxcox
    res <- boxcox_censor(time_trt_boot,time_ctrl_boot,cens_stat_trt_boot,cens_stat_ctrl_boot)
    g = res$g
    exitflag = res$exitflag
    gval = res$gval

    lambda_boot = g[5]
    mu_trt_hat_boot = g[1];  sd_trt_hat_boot = g[2]
    mu_ctrl_hat_boot = g[3];  sd_ctrl_hat_boot = g[4]

    # Transform the dataset
    time_trt_box_boot = ((time_trt_boot^lambda_boot)-1)/lambda_boot
    time_ctrl_box_boot = ((time_ctrl_boot^lambda_boot)-1)/lambda_boot


    ###########################################################
    ################# Censored obs imputation #################

    ########### Trt arm
    # Find the number of RIGHT censored obs after last event time OR last right censored
    n_last_censor_trt_boot = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_trt_boot)
    # Find the number of LEFT censored obs before 1st event time OR 1st right censored
    n_first_leftcen_trt_boot = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_trt_boot)

    ########### Ctrl arm
    # Find the number of RIGHT censored obs after last event time OR last right censored
    n_last_censor_ctrl_boot = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_ctrl_boot)
    # Find the number of LEFT censored obs before 1st event time OR 1st right censored
    n_first_leftcen_ctrl_boot = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_ctrl_boot)


    #### Check if both arm's last obs is event
    if (n_last_censor_trt_boot == 0 && n_first_leftcen_trt_boot == 0 && n_last_censor_ctrl_boot == 0 && n_first_leftcen_ctrl_boot == 0) {

      time_trt_box_imput_boot = time_trt_box_boot
      time_ctrl_box_imput_boot = time_ctrl_box_boot
      cens_stat_trt_imput_boot = cens_stat_trt_boot
      cens_stat_ctrl_imput_boot = cens_stat_ctrl_boot

      # Compute Overlap
      f_trt_boot <- function(x) kde_double_censored(x, time_trt_box_imput_boot, cens_stat_trt_imput_boot)
      f_ctrl_boot <- function(x) kde_double_censored(x, time_ctrl_box_imput_boot, cens_stat_ctrl_imput_boot)

      integrand1 <- function(x) pmin(f_trt_boot(x), f_ctrl_boot(x))
      overlap_boot = integrate(integrand1, -Inf, Inf)$value
      overlap_inv_boot[b] = 1 + (-overlap_boot)

    } else {

      overlap_inv_sample_boot = numeric(30)

      #### Otherwise do imputation
      for (imput_boot in 1:30) {

        time_trt_box_imput_boot = time_trt_box_boot
        time_ctrl_box_imput_boot = time_ctrl_box_boot
        cens_stat_trt_imput_boot = cens_stat_trt_boot
        cens_stat_ctrl_imput_boot = cens_stat_ctrl_boot

        impute_sample1_right_boot = numeric(n_last_censor_trt_boot)
        impute_sample1_left_boot = numeric(n_first_leftcen_trt_boot)
        impute_sample2_right_boot = numeric(n_last_censor_ctrl_boot)
        impute_sample2_left_boot = numeric(n_first_leftcen_ctrl_boot)


        ########### Trt arm ##########

        ####### Censor in RIGHT tail
        if (n_last_censor_trt_boot == 0) {
          time_trt_box_imput_boot = time_trt_box_imput_boot
          cens_stat_trt_imput_boot = cens_stat_trt_imput_boot
        } else {
          # Impute the censored obs after last event time/left censored
          for (i_t_boot in 1:n_last_censor_trt_boot) {
            sample1_R_boot = time_trt_box_boot[n_trt - n_last_censor_trt_boot + i_t_boot]   # Initial sample
            while (sample1_R_boot <= time_trt_box_boot[n_trt - n_last_censor_trt_boot + i_t_boot]) {
              sample1_R_boot = rnorm(1, mean = mu_trt_hat_boot, sd = sd_trt_hat_boot)
            }
            impute_sample1_right_boot[i_t_boot] = sample1_R_boot
          }
          # Substitute the obs in the main data
          time_trt_box_imput_boot[(n_trt - n_last_censor_trt_boot + 1):n_trt] = impute_sample1_right_boot
          cens_stat_trt_imput_boot[(n_trt - n_last_censor_trt_boot + 1):n_trt] = 0
        }



        ####### Censor in LEFT tail
        if (n_first_leftcen_trt_boot == 0) {
          time_trt_box_imput_boot = time_trt_box_imput_boot
          cens_stat_trt_imput_boot = cens_stat_trt_imput_boot
        } else {
          # Impute the censored obs before first event time/right censored
          for (i_t_boot in 1:n_first_leftcen_trt_boot) {
            sample1_L_boot = time_trt_box_boot[i_t_boot]   # Initial sample
            while (sample1_L_boot >= time_trt_box_boot[i_t_boot]) {
              sample1_L_boot = rnorm(1, mean = mu_trt_hat_boot, sd = sd_trt_hat_boot)
            }
            impute_sample1_left_boot[i_t_boot] = sample1_L_boot
          }

          # Substitute the obs in the main data
          time_trt_box_imput_boot[1:n_first_leftcen_trt_boot] = impute_sample1_left_boot
          cens_stat_trt_imput_boot[1:n_first_leftcen_trt_boot] = 0
        }


        ########### Ctrl arm ##########
        ####### Censor in RIGHT tail
        if (n_last_censor_ctrl_boot == 0) {
          time_ctrl_box_imput_boot = time_ctrl_box_imput_boot
          cens_stat_ctrl_imput_boot = cens_stat_ctrl_imput_boot
        } else {
          # Impute the censored obs after last event time/left censored
          for (i_c_boot in 1:n_last_censor_ctrl_boot) {
            sample2_R_boot = time_ctrl_box_boot[n_ctrl - n_last_censor_ctrl_boot + i_c_boot]   # Initial sample
            while (sample2_R_boot <= time_ctrl_box_boot[n_ctrl - n_last_censor_ctrl_boot + i_c_boot]) {
              sample2_R_boot = rnorm(1, mean = mu_ctrl_hat_boot, sd = sd_ctrl_hat_boot)
            }
            impute_sample2_right_boot[i_c_boot] = sample2_R_boot
          }
          # Substitute the obs in the main data
          time_ctrl_box_imput_boot[(n_ctrl - n_last_censor_ctrl_boot + 1):n_ctrl] = impute_sample2_right_boot
          cens_stat_ctrl_imput_boot[(n_ctrl - n_last_censor_ctrl_boot + 1):n_ctrl] = 0
        }


        ####### Censor in LEFT tail
        if (n_first_leftcen_ctrl_boot == 0) {
          time_ctrl_box_imput_boot = time_ctrl_box_imput_boot
          cens_stat_ctrl_imput_boot = cens_stat_ctrl_imput_boot
        } else {
          # Impute the censored obs before first event time/right censored
          for (i_c_boot in 1:n_first_leftcen_ctrl_boot) {
            sample2_L_boot = time_ctrl_box_boot[i_c_boot]   # Initial sample
            while (sample2_L_boot >= time_ctrl_box_boot[i_c_boot]) {
              sample2_L_boot = rnorm(1, mean = mu_ctrl_hat_boot, sd = sd_ctrl_hat_boot)
            }
            impute_sample2_left_boot[i_c_boot] = sample2_L_boot
          }

          # Substitute the obs in the main data
          time_ctrl_box_imput_boot[1:n_first_leftcen_ctrl_boot] = impute_sample2_left_boot
          cens_stat_ctrl_imput_boot[1:n_first_leftcen_ctrl_boot] = 0
        }


        # Compute Overlap
        f_trt_boot <- function(x) kde_double_censored(x, time_trt_box_imput_boot, cens_stat_trt_imput_boot)
        f_ctrl_boot <- function(x) kde_double_censored(x, time_ctrl_box_imput_boot, cens_stat_ctrl_imput_boot)

        integrand1 <- function(x) pmin(f_trt_boot(x), f_ctrl_boot(x))
        overlap_boot = integrate(integrand1, -Inf, Inf)$value
        overlap_inv_sample_boot[imput_boot] = 1 + (-overlap_boot)
      }

      overlap_inv_boot[b] = mean(overlap_inv_sample_boot)
    }
  }


  # Calculate p-values
  prop_ovl = sum(overlap_inv_boot > overlap_inv)
  pval_ovl = prop_ovl / length(overlap_inv_boot)


  ####### Result #######
  df_result <- data.frame(
    estimate = sprintf("%.4f", overlap),
    p_value = pval_ovl
  )

  rownames(df_result) <- "OVL"

  # Message text
  message_text <- paste(
    "Hypothesis testing method: Overlap coefficient (OVL)",
    "Significance level: 5% (two-sided)",
    paste0("Number of permutations used for inference: ", boots),
    sep = "\n"
  )

  # Create the list object & Assign a custom class name to this list
  out <- list(message = message_text, result = df_result)
  class(out) <- "survival_comp_roc"
  return(out)
}







######################################################################
############## Our Method with double censoring ###############
#### ------- JOINT ROC Length & OVL -------- #######
######################################################################

joint.roc_ovl_DoubleCenSurvival_test <- function(time1, censor1, time2, censor2, boots, progress = TRUE, plot = FALSE) {

  ##### -------------------------------------------------------------------
  ##### censoring status (0 = event, 1 = right censored, -1 = left censored)

  # double censored data for trt and ctrl
  time_tt = time1
  cens_stat_t = censor1
  n_trt = length(time_tt)

  time_cc = time2
  cens_stat_c = censor2
  n_ctrl = length(time_cc)

  # Log transform
  time_t = log(1 + time_tt)
  time_c = log(1 + time_cc)

  # Sort by time
  idx_t = order(time_t)
  time_trt = time_t[idx_t]
  cens_stat_trt = cens_stat_t[idx_t]

  idx_c = order(time_c)
  time_ctrl = time_c[idx_c]
  cens_stat_ctrl = cens_stat_c[idx_c]

  # Box-Cox (with censoring) to estimate parameters and lambda
  res <- boxcox_censor(time_trt, time_ctrl, cens_stat_trt, cens_stat_ctrl)
  g = res$g
  exitflag = res$exitflag
  gval = res$gval

  lambda = g[5]
  mu_trt_hat  = g[1];  sd_trt_hat  = g[2]
  mu_ctrl_hat = g[3];  sd_ctrl_hat = g[4]

  # Transform data to Box-Cox scale
  time_trt_box  = ((time_trt^lambda)-1) / lambda
  time_ctrl_box = ((time_ctrl^lambda)-1) / lambda


  ###########################################################
  ################# Censored obs imputation #################

  ########### Trt arm
  # Find the number of RIGHT censored obs after last event time OR last right censored
  n_last_censor_trt = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_trt);
  # Find the number of LEFT censored obs before 1st event time OR 1st right censored
  n_first_leftcen_trt = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_trt);

  ########### Ctrl arm
  # Find the number of RIGHT censored obs after last event time OR last right censored
  n_last_censor_ctrl = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_ctrl);
  # Find the number of LEFT censored obs before 1st event time OR 1st right censored
  n_first_leftcen_ctrl = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_ctrl);


  #### Check if both arm's last obs is event

  if (n_last_censor_trt == 0 && n_first_leftcen_trt  == 0 && n_last_censor_ctrl == 0 && n_first_leftcen_ctrl == 0) {

    time_trt_box_imput = time_trt_box
    time_ctrl_box_imput = time_ctrl_box
    cens_stat_trt_imput = cens_stat_trt
    cens_stat_ctrl_imput = cens_stat_ctrl

    # Compute ROC length
    f_trt <- function(x) kde_double_censored(x, time_trt_box_imput, cens_stat_trt_imput)
    f_ctrl <- function(x) kde_double_censored(x, time_ctrl_box_imput, cens_stat_ctrl_imput)
    integrand <- function(x) sqrt(f_trt(x)^2 + f_ctrl(x)^2)

    length_roc = integrate(integrand, -Inf, Inf)$value
    length_roc_diff = length_roc - sqrt(2)

    # Compute Overlap
    integrand1 <- function(x) pmin(f_trt(x), f_ctrl(x))
    overlap = integrate(integrand1, -Inf, Inf)$value
    overlap_inv = 1 + (-overlap)

  } else {

    length_roc_sample = numeric(30)
    overlap_sample = numeric(30)

    #### Otherwise do imputation
    for (imput in 1:30) {

      time_trt_box_imput = time_trt_box
      time_ctrl_box_imput = time_ctrl_box
      cens_stat_trt_imput = cens_stat_trt
      cens_stat_ctrl_imput = cens_stat_ctrl

      impute_sample1_right = numeric(n_last_censor_trt)
      impute_sample1_left = numeric(n_first_leftcen_trt)
      impute_sample2_right = numeric(n_last_censor_ctrl)
      impute_sample2_left = numeric(n_first_leftcen_ctrl)


      ########### Trt arm ##########
      ####### Censor in RIGHT tail
      if (n_last_censor_trt == 0) {
        time_trt_box_imput = time_trt_box_imput
        cens_stat_trt_imput = cens_stat_trt_imput
      } else {
        # Impute the censored obs after last event time/left censored
        for (i_t in 1:n_last_censor_trt) {
          sample1_R = time_trt_box[n_trt - n_last_censor_trt + i_t]   # Initial sample
          while (sample1_R <= time_trt_box[n_trt - n_last_censor_trt + i_t]) {
            sample1_R = rnorm(1, mean = mu_trt_hat, sd = sd_trt_hat)
          }
          impute_sample1_right[i_t] = sample1_R
        }
        # Substitute the obs in the main data
        time_trt_box_imput[(n_trt - n_last_censor_trt + 1):n_trt] = impute_sample1_right
        cens_stat_trt_imput[(n_trt - n_last_censor_trt + 1):n_trt] = 0
      }


      ####### Censor in LEFT tail
      if (n_first_leftcen_trt == 0) {
        time_trt_box_imput = time_trt_box_imput
        cens_stat_trt_imput = cens_stat_trt_imput
      } else {
        # Impute the censored obs before first event time/right censored
        for (i_t in 1:n_first_leftcen_trt) {
          sample1_L = time_trt_box[i_t]   # Initial sample
          while (sample1_L >= time_trt_box[i_t]) {
            sample1_L = rnorm(1, mean = mu_trt_hat, sd = sd_trt_hat)
          }
          impute_sample1_left[i_t] = sample1_L
        }

        # Substitute the obs in the main data
        time_trt_box_imput[1:n_first_leftcen_trt] = impute_sample1_left
        cens_stat_trt_imput[1:n_first_leftcen_trt] = 0
      }


      ########### Ctrl arm ##########
      ####### Censor in RIGHT tail
      if (n_last_censor_ctrl == 0) {
        time_ctrl_box_imput = time_ctrl_box_imput
        cens_stat_ctrl_imput = cens_stat_ctrl_imput
      } else {
        # Impute the censored obs after last event time/left censored
        for (i_c in 1:n_last_censor_ctrl) {
          sample2_R = time_ctrl_box[n_ctrl - n_last_censor_ctrl + i_c]   # Initial sample
          while (sample2_R <= time_ctrl_box[n_ctrl - n_last_censor_ctrl + i_c]) {
            sample2_R = rnorm(1, mean = mu_ctrl_hat, sd = sd_ctrl_hat)
          }
          impute_sample2_right[i_c] = sample2_R
        }

        # Substitute the obs in the main data
        time_ctrl_box_imput[(n_ctrl - n_last_censor_ctrl + 1):n_ctrl] = impute_sample2_right
        cens_stat_ctrl_imput[(n_ctrl - n_last_censor_ctrl + 1):n_ctrl] = 0
      }


      ####### Censor in LEFT tail
      if (n_first_leftcen_ctrl == 0) {
        time_ctrl_box_imput = time_ctrl_box_imput
        cens_stat_ctrl_imput = cens_stat_ctrl_imput
      } else {
        # Impute the censored obs before first event time/right censored
        for (i_c in 1:n_first_leftcen_ctrl) {
          sample2_L = time_ctrl_box[i_c]   # Initial sample
          while (sample2_L >= time_ctrl_box[i_c]) {
            sample2_L = rnorm(1, mean = mu_ctrl_hat, sd = sd_ctrl_hat)
          }
          impute_sample2_left[i_c] = sample2_L
        }

        # Substitute the obs in the main data
        time_ctrl_box_imput[1:n_first_leftcen_ctrl] = impute_sample2_left
        cens_stat_ctrl_imput[1:n_first_leftcen_ctrl] = 0
      }

      # Compute ROC length
      f_trt <- function(x) kde_double_censored(x, time_trt_box_imput, cens_stat_trt_imput)
      f_ctrl <- function(x) kde_double_censored(x, time_ctrl_box_imput, cens_stat_ctrl_imput)
      integrand <- function(x) sqrt(f_trt(x)^2 + f_ctrl(x)^2)

      length_roc_sample[imput] = integrate(integrand, -Inf, Inf)$value

      # Compute Overlap
      integrand1 <- function(x) pmin(f_trt(x), f_ctrl(x))
      overlap_sample[imput] = integrate(integrand1, -Inf, Inf)$value
    }

    length_roc = mean(length_roc_sample)
    length_roc_diff = length_roc - sqrt(2)
    overlap = mean(overlap_sample)
    overlap_inv = 1 - overlap
  }


  ########################################
  #################### Permutation
  ########################################

  time = c(time_t, time_c)    # Merge unsorted values
  cens = c(cens_stat_t, cens_stat_c)
  n = n_trt + n_ctrl

  length_roc_boot_diff = numeric(boots)
  overlap_inv_boot = numeric(boots)

  show_progress <- progress && interactive()
  if (show_progress) {
    cat("Permutation test in progress...\n")
    bar_width <- 50
  }

  for (b in 1:boots) {

    # ----------- Progress bar update ----------

    if (show_progress && (b %% 10 == 0 || b == boots)) {
      prog <- b / boots
      n_stars <- floor(prog * bar_width)
      n_spaces <- bar_width - n_stars

      bar <- paste0(
        "|",
        paste(rep("*", n_stars), collapse = ""),
        paste(rep(" ", n_spaces), collapse = ""),
        "| ",
        sprintf("%3d%%", floor(prog * 100))
      )

      cat("\r", bar, sep = "")
      if (b == boots) cat("\n")  # move to next line at end
      flush.console()
    }

    # ------------------------------------------


    # Generate permutation sample from original sample
    at = sample(1:n, n, replace = FALSE)
    time_t_boot = time[at[1:n_trt]]
    cens_stat_t_boot = cens[at[1:n_trt]]
    time_c_boot = time[at[(n_trt + 1):n]]
    cens_stat_c_boot = cens[at[(n_trt + 1):n]]

    # Sort based on time values
    data1_boot = cbind(time_t_boot, cens_stat_t_boot)
    data_sorted1_boot = data1_boot[order(data1_boot[,1]),]
    time_trt_boot = data_sorted1_boot[, 1]
    cens_stat_trt_boot = data_sorted1_boot[, 2]

    data2_boot = cbind(time_c_boot, cens_stat_c_boot)
    data_sorted2_boot = data2_boot[order(data2_boot[,1]),]
    time_ctrl_boot = data_sorted2_boot[, 1]
    cens_stat_ctrl_boot = data_sorted2_boot[, 2]

    # Perform Boxcox
    res <- boxcox_censor(time_trt_boot,time_ctrl_boot,cens_stat_trt_boot,cens_stat_ctrl_boot)
    g = res$g
    exitflag = res$exitflag
    gval = res$gval

    lambda_boot = g[5]
    mu_trt_hat_boot = g[1];  sd_trt_hat_boot = g[2]
    mu_ctrl_hat_boot = g[3];  sd_ctrl_hat_boot = g[4]

    # Transform the dataset
    time_trt_box_boot = ((time_trt_boot^lambda_boot)-1)/lambda_boot
    time_ctrl_box_boot = ((time_ctrl_boot^lambda_boot)-1)/lambda_boot


    ###########################################################
    ################# Censored obs imputation #################

    ########### Trt arm
    # Find the number of RIGHT censored obs after last event time OR last right censored
    n_last_censor_trt_boot = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_trt_boot)
    # Find the number of LEFT censored obs before 1st event time OR 1st right censored
    n_first_leftcen_trt_boot = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_trt_boot)

    ########### Ctrl arm
    # Find the number of RIGHT censored obs after last event time OR last right censored
    n_last_censor_ctrl_boot = Calc_number_RIGHTcen_after_lasteventORleftcen(cens_stat_ctrl_boot)
    # Find the number of LEFT censored obs before 1st event time OR 1st right censored
    n_first_leftcen_ctrl_boot = Calc_number_LEFTcen_before_1steventORrightcen(cens_stat_ctrl_boot)


    #### Check if both arm's last obs is event
    if (n_last_censor_trt_boot == 0 && n_first_leftcen_trt_boot == 0 && n_last_censor_ctrl_boot == 0 && n_first_leftcen_ctrl_boot == 0) {

      time_trt_box_imput_boot = time_trt_box_boot
      time_ctrl_box_imput_boot = time_ctrl_box_boot
      cens_stat_trt_imput_boot = cens_stat_trt_boot
      cens_stat_ctrl_imput_boot = cens_stat_ctrl_boot

      # Compute ROC length
      f_trt_boot <- function(x) kde_double_censored(x, time_trt_box_imput_boot, cens_stat_trt_imput_boot)
      f_ctrl_boot <- function(x) kde_double_censored(x, time_ctrl_box_imput_boot, cens_stat_ctrl_imput_boot)
      integrand <- function(x) sqrt(f_trt_boot(x)^2 + f_ctrl_boot(x)^2)

      length_roc_boot = integrate(integrand, -Inf, Inf)$value
      length_roc_boot_diff[b] = length_roc_boot - sqrt(2)

      # Compute Overlap
      integrand1 <- function(x) pmin(f_trt_boot(x), f_ctrl_boot(x))
      overlap_boot = integrate(integrand1, -Inf, Inf)$value
      overlap_inv_boot[b] = 1 + (-overlap_boot)

    } else {

      length_roc_sample_boot = numeric(30)
      overlap_inv_sample_boot = numeric(30)

      #### Otherwise do imputation
      for (imput_boot in 1:30) {

        time_trt_box_imput_boot = time_trt_box_boot
        time_ctrl_box_imput_boot = time_ctrl_box_boot
        cens_stat_trt_imput_boot = cens_stat_trt_boot
        cens_stat_ctrl_imput_boot = cens_stat_ctrl_boot

        impute_sample1_right_boot = numeric(n_last_censor_trt_boot)
        impute_sample1_left_boot = numeric(n_first_leftcen_trt_boot)
        impute_sample2_right_boot = numeric(n_last_censor_ctrl_boot)
        impute_sample2_left_boot = numeric(n_first_leftcen_ctrl_boot)


        ########### Trt arm ##########

        ####### Censor in RIGHT tail
        if (n_last_censor_trt_boot == 0) {
          time_trt_box_imput_boot = time_trt_box_imput_boot
          cens_stat_trt_imput_boot = cens_stat_trt_imput_boot
        } else {
          # Impute the censored obs after last event time/left censored
          for (i_t_boot in 1:n_last_censor_trt_boot) {
            sample1_R_boot = time_trt_box_boot[n_trt - n_last_censor_trt_boot + i_t_boot]   # Initial sample
            while (sample1_R_boot <= time_trt_box_boot[n_trt - n_last_censor_trt_boot + i_t_boot]) {
              sample1_R_boot = rnorm(1, mean = mu_trt_hat_boot, sd = sd_trt_hat_boot)
            }
            impute_sample1_right_boot[i_t_boot] = sample1_R_boot
          }
          # Substitute the obs in the main data
          time_trt_box_imput_boot[(n_trt - n_last_censor_trt_boot + 1):n_trt] = impute_sample1_right_boot
          cens_stat_trt_imput_boot[(n_trt - n_last_censor_trt_boot + 1):n_trt] = 0
        }



        ####### Censor in LEFT tail
        if (n_first_leftcen_trt_boot == 0) {
          time_trt_box_imput_boot = time_trt_box_imput_boot
          cens_stat_trt_imput_boot = cens_stat_trt_imput_boot
        } else {
          # Impute the censored obs before first event time/right censored
          for (i_t_boot in 1:n_first_leftcen_trt_boot) {
            sample1_L_boot = time_trt_box_boot[i_t_boot]   # Initial sample
            while (sample1_L_boot >= time_trt_box_boot[i_t_boot]) {
              sample1_L_boot = rnorm(1, mean = mu_trt_hat_boot, sd = sd_trt_hat_boot)
            }
            impute_sample1_left_boot[i_t_boot] = sample1_L_boot
          }

          # Substitute the obs in the main data
          time_trt_box_imput_boot[1:n_first_leftcen_trt_boot] = impute_sample1_left_boot
          cens_stat_trt_imput_boot[1:n_first_leftcen_trt_boot] = 0
        }


        ########### Ctrl arm ##########
        ####### Censor in RIGHT tail
        if (n_last_censor_ctrl_boot == 0) {
          time_ctrl_box_imput_boot = time_ctrl_box_imput_boot
          cens_stat_ctrl_imput_boot = cens_stat_ctrl_imput_boot
        } else {
          # Impute the censored obs after last event time/left censored
          for (i_c_boot in 1:n_last_censor_ctrl_boot) {
            sample2_R_boot = time_ctrl_box_boot[n_ctrl - n_last_censor_ctrl_boot + i_c_boot]   # Initial sample
            while (sample2_R_boot <= time_ctrl_box_boot[n_ctrl - n_last_censor_ctrl_boot + i_c_boot]) {
              sample2_R_boot = rnorm(1, mean = mu_ctrl_hat_boot, sd = sd_ctrl_hat_boot)
            }
            impute_sample2_right_boot[i_c_boot] = sample2_R_boot
          }
          # Substitute the obs in the main data
          time_ctrl_box_imput_boot[(n_ctrl - n_last_censor_ctrl_boot + 1):n_ctrl] = impute_sample2_right_boot
          cens_stat_ctrl_imput_boot[(n_ctrl - n_last_censor_ctrl_boot + 1):n_ctrl] = 0
        }


        ####### Censor in LEFT tail
        if (n_first_leftcen_ctrl_boot == 0) {
          time_ctrl_box_imput_boot = time_ctrl_box_imput_boot
          cens_stat_ctrl_imput_boot = cens_stat_ctrl_imput_boot
        } else {
          # Impute the censored obs before first event time/right censored
          for (i_c_boot in 1:n_first_leftcen_ctrl_boot) {
            sample2_L_boot = time_ctrl_box_boot[i_c_boot]   # Initial sample
            while (sample2_L_boot >= time_ctrl_box_boot[i_c_boot]) {
              sample2_L_boot = rnorm(1, mean = mu_ctrl_hat_boot, sd = sd_ctrl_hat_boot)
            }
            impute_sample2_left_boot[i_c_boot] = sample2_L_boot
          }

          # Substitute the obs in the main data
          time_ctrl_box_imput_boot[1:n_first_leftcen_ctrl_boot] = impute_sample2_left_boot
          cens_stat_ctrl_imput_boot[1:n_first_leftcen_ctrl_boot] = 0
        }


        # Compute ROC length
        f_trt_boot <- function(x) kde_double_censored(x, time_trt_box_imput_boot, cens_stat_trt_imput_boot)
        f_ctrl_boot <- function(x) kde_double_censored(x, time_ctrl_box_imput_boot, cens_stat_ctrl_imput_boot)
        integrand <- function(x) sqrt(f_trt_boot(x)^2 + f_ctrl_boot(x)^2)

        length_roc_sample_boot[imput_boot] = integrate(integrand, -Inf, Inf)$value

        # Compute Overlap
        integrand1 <- function(x) pmin(f_trt_boot(x), f_ctrl_boot(x))
        overlap_boot = integrate(integrand1, -Inf, Inf)$value
        overlap_inv_sample_boot[imput_boot] = 1 + (-overlap_boot)
      }

      length_roc_boot = mean(length_roc_sample_boot)
      length_roc_boot_diff[b] = length_roc_boot - sqrt(2)

      overlap_inv_boot[b] = mean(overlap_inv_sample_boot)
    }
  }


  ############################################################
  ################ Standardized values & Convex Hull #########

  # Perform Standardization of values
  df = rbind(
    cbind(length_roc_boot_diff, overlap_inv_boot),
    c(length_roc_diff, overlap_inv)
  )

  df1 = df[,1]
  df2 = df[,2]

  df11 = sapply(df1, function(x) mean(df1 <= x))
  df22 = sapply(df2, function(x) mean(df2 <= x))

  df_norm = cbind(df11, df22)


  # Separate critical point and Permutation values
  length_roc_diff_norm = df_norm[nrow(df_norm), 1]
  overlap_inv_norm = df_norm[nrow(df_norm), 2]

  df_norm = df_norm[-nrow(df_norm), ]

  length_roc_boot_diff_norm = df_norm[,1]
  overlap_inv_boot_norm = df_norm[,2]


  # Calculate the Euclidean distance for all point in the plot
  euc_distance_boot = numeric(length(overlap_inv_boot_norm))

  for (euc in 1:length(overlap_inv_boot_norm)) {
    euc_distance_boot[euc] = sqrt(
      (overlap_inv_boot_norm[euc] - 0)^2 +
        (length_roc_boot_diff_norm[euc] - 0)^2
    )
  }


  #################################################
  # --------------- Convex Hull Plot --------------
  plot_joint <- NULL

  if (plot) {

    threshold_joint <- quantile(euc_distance_boot, 0.95)

    group_joint <- ifelse(euc_distance_boot <= threshold_joint, "Inside", "Outside")
    df_joint <- data.frame(
      X = overlap_inv_boot_norm,
      Y = length_roc_boot_diff_norm,
      group = factor(group_joint, levels = c("Inside", "Outside"))
    )

    df_in_joint <- df_joint[df_joint$group == "Inside", ]
    hull_indices_joint <- chull(df_in_joint$X, df_in_joint$Y)
    hull_indices_joint <- c(hull_indices_joint, hull_indices_joint[1])
    hull_df_joint <- df_in_joint[hull_indices_joint, ]

    plot_joint <- ggplot(df_joint, aes(x = X, y = Y)) +
      geom_point(data = subset(df_joint, group == "Inside"),
                 color = "#00BFC4", size = 0.6) +
      geom_point(data = subset(df_joint, group == "Outside"),
                 color = "#F8766D", size = 0.6) +
      geom_polygon(data = hull_df_joint,
                   fill = "#00BFC4", alpha = 0.1,
                   color = "#00BFC4", linewidth = 0.6) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      labs(
        x = expression(1 - hat(OVL)),
        y = expression(hat(italic(l))[ROC] - sqrt(2))
      ) +
      theme_classic(base_size = 12) +
      theme(
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      ) +
      geom_point(
        data = data.frame(X = overlap_inv_norm, Y = length_roc_diff_norm),
        aes(x = X, y = Y),
        shape = 8,
        color = "blue",
        size = 1,
        stroke = 1
      )
  }
  # -----------------------------------------------


  # Calculate p-values
  prop_length = sum(length_roc_boot_diff_norm > length_roc_diff_norm)
  pval_length = prop_length / length(length_roc_boot_diff_norm)

  prop_ovl = sum(overlap_inv_boot_norm > overlap_inv_norm)
  pval_ovl = prop_ovl / length(overlap_inv_boot_norm)

  euc_distance_observed = sqrt(
    (overlap_inv_norm - 0)^2 +
      (length_roc_diff_norm - 0)^2
  )

  pval_joint = sum(euc_distance_boot > euc_distance_observed) / length(euc_distance_boot)


  ####### Result #######
  df_result <- data.frame(
    estimate = c(sprintf("%.4f", length_roc), sprintf("%.4f", overlap), "-"),
    p_value = c(pval_length, pval_ovl, pval_joint)
  )

  rownames(df_result) <- c("ROC Length", "OVL", "Joint ROC Length-OVL")

  # Message text
  message_text <- paste(
    "Hypothesis testing methods: ROC length, Overlap coefficient (OVL), and the joint ROC Length-OVL",
    "Significance level: 5% (two-sided)",
    paste0("Number of permutations used for inference: ", boots),
    sep = "\n"
  )

  # Create the list object & Assign a custom class name to this list
  out <- list(message = message_text, result = df_result, plot = plot_joint)
  class(out) <- "survival_comp_roc"
  return(out)
}




############# Output Display Type
#' @method print survival_comp_roc
#' @export
print.survival_comp_roc <- function(x, ...) {
  cat(x$message, "\n\n")
  print(x$result)
  if (!is.null(x$plot)) {
    cat("\n")
    print(x$plot)
  }
  invisible(x)
}

