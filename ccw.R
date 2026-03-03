ccw_sim = function(n, p_e, p_m, alpha, seed = 123){
  toy = data.frame(ID = 1:n)
  
  set.seed(seed)
  toy = toy %>% 
    mutate(exposure_time = rgeom(n, p_e)) %>% 
    mutate(death_time_alpha0 = rgeom(n, p_m)) %>% 
    mutate(death_time = ifelse(
      death_time_alpha0 <=  exposure_time, #does an individual survive until exposure
      death_time_alpha0, #if they don't then their death time is the same
      exposure_time + rgeom(n, p_m+alpha) #using the fact that geometric distribution is memoryless
      #we can generate the time until death using the new probability after exposure
      #note that technically this time can decrease from the original
    ))
  
  
  
  toy = toy %>% 
    mutate(death_0 = F, #all patients start alive
           death_30 = ifelse(death_time < 30, T, F), #check whether the patient has died at every time measurement
           death_60 = ifelse(death_time < 60, T, F),
           death_90 = ifelse(death_time < 90, T, F),
           death_120 = ifelse(death_time < 120, T, F),
           death_150 = ifelse(death_time < 150, T, F)) %>% 
    #people are only in the exposed group if they were prescribed first 90 days and were alive before time of prescription
    mutate(exposed = ifelse(exposure_time < death_time & exposure_time <= 90, T, F)) %>% 
    mutate(exposure_0 = F,
           exposure_30 = ifelse(exposure_time < death_time & exposure_time < 30, T, F),
           exposure_60 = ifelse(exposure_time < death_time & exposure_time < 60, T, F),
           exposure_90 = ifelse(exposure_time < death_time & exposure_time < 90, T, F),
           exposure_120 = exposure_90,
           exposure_150 = exposure_90)
  
  
  
  #create the clones
  toy_clone = toy %>%
    crossing(exposed_clone = c(T, F)) 
  
  
  toy_long = toy_clone %>%
    pivot_longer(
      cols = matches("_(\\d+)$"),
      names_to = c(".value", "time"),
      names_sep = "_"
    ) %>% 
    mutate(time = as.numeric(time))

    #exposed is if they were exposed before 90 days
  #exposure is if they were exposed by the time point
  #exposed_clone is if their clone has been assigned to be exposed
  toy_censor = toy_long %>% 
    mutate(uncensored = ifelse(exposed_clone ==T & exposed == F & time >= 90, F, T)) %>% #the clone is assigned to be exposed but was not by day 90, and it is day 90 or after then it should be censored
    mutate(uncensored = ifelse(exposed_clone == F & exposure == T, F, uncensored)) #the clone is assigned to non-exposed but becomes exposed then it should be censored
  
  #filter out rows after first censored measurement or death 
  toy_filter = toy_censor %>%
    group_by(ID, exposed_clone) %>%
    filter(cumsum(uncensored == F | death == TRUE) <= 1) %>% #remove all rows after first death
    ungroup()
  
  
  
  #need it in wide form to fit the glms for the weights
  toy_wide = toy_filter %>%
    pivot_wider(
      id_cols = c(ID, exposed_clone),
      names_from = time,
      values_from = c(death, exposure, uncensored, exposed),
      names_sep = "_"
    )
  #the denominators of the weights (1 - probability of being censored is same as probability of being uncensored)
  censor_fit60 = glm(uncensored_60 ~ exposure_30 + uncensored_30, data = toy_wide, family = "binomial", na.action = na.exclude) %>% fitted.values()
  censor_fit90 = glm(uncensored_90 ~ exposure_60 + uncensored_60, data = toy_wide, family = "binomial", na.action = na.exclude) %>% fitted.values()

  #numerator of the weights
  censor_60_num = glm(uncensored_60 ~  uncensored_30, data = toy_wide, family = "binomial", na.action = na.exclude) %>% fitted.values()
  censor_90_num = glm(uncensored_90 ~ uncensored_60, data = toy_wide, family = "binomial", na.action = na.exclude) %>% fitted.values()
  
  toy_wide = toy_wide %>% 
    mutate(weights_0 = 1, 
           weights_30 = 1,
           weights_60 = censor_60_num/censor_fit60,
           weights_90 = censor_90_num/censor_fit90 * weights_60,
           weights_120 = weights_90,
           weights_150 = weights_90)
  
  toy_weights = toy_wide %>%
    pivot_longer(
      cols = matches("_(\\d+)$"),
      names_to = c(".value", "time"),
      names_sep = "_"
    ) %>% 
    mutate(time = as.numeric(time)) %>% 
    drop_na()
  
  
  toy_long_notclone = toy %>%
    pivot_longer(
      cols = matches("_(\\d+)$"),
      names_to = c(".value", "time"),
      names_sep = "_"
    ) %>% 
    mutate(time = as.numeric(time))
  
  toy_filter_notclone = toy_censor %>%
    group_by(ID) %>%
    filter(cumsum(death == TRUE) <= 1) %>% #remove all rows after first death
    ungroup()
  
  return(list(data = toy_filter_notclone, dataccw = toy_weights))
  
  
}

#P(T_m <=t | T_E = y)
#probability of death happening before or at t given that exposure occurs at y 
#in the case where exposure can occur after death

P_M_beforeat_t_given_E_at_y = function(t, y, p_m, p_e, alpha){
  if (t<= y){
    return(1 - (1-p_m)^t)
  } else {
    1 - (1-p_m)^y * (1-p_m - alpha)^(t-y)
  }
}


#P(T_m <=t | T_E <= y)
#probability of death happening before or at t given that exposure before or at y 
#in the case where exposure can occur after death


P_M_beforeat_t_given_E_beforeat_y = function(t, y, p_m, p_e, alpha){
  k = 1:y
  # conditional on T_E = k
  P_Tm_given_TE <- ifelse(
    k <= t,
    1 - (1-p_m)^(k-1) * (1-(p_m+alpha))^(t-k+1),  # exposure before t
    1 - (1-p_m)^t                               # exposure after t, death happens before exposure
  )
  # weight by P(T_E = k | T_E <= y)
  probs <- (1-p_e)^(k-1) * p_e / (1 - (1-p_e)^y) * P_Tm_given_TE
  sum(probs)
  
}

#P(T_M<=t|E = 0)
#probability of time to death being less than t given no exposure
P_M_beforeat_t_given_noE = function(t, y, p_m, p_e, alpha){
  return(1 - (1-p_m)^t)
}

logit =  function(p){
  log(p / (1 - p))
}

expit = function(x){
  exp(x)/(1+exp(x))
}

# Extract coefficients
#coefs = coef(model)

#find P(T_M <= t | A) for the pooled logistic models
get_cumulative_prob = function(coefs, exposed = TRUE, times = c(30, 60, 90, 120, 150)){
  h = numeric(length(times))
  for(i in seq_along(times)){
    t = times[i]
    # construct linear predictor for this interval
    lp = coefs[paste0("factor(time)", t)] + coefs["(Intercept)"]
    if(exposed){
      lp = lp + coefs["exposedTRUE"] + coefs[paste0("exposedTRUE:factor(time)", t)] 
    }
    h[i] = exp(lp) / (1 + exp(lp))  # hazard in this interval
  }
  # cumulative probability by the last interval
  1 - prod(1 - h)
}
