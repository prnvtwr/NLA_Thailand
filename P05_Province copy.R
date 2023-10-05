#Province-lvel Models

province_df_list <- split(df,df$Province)



CHIKV_df <- df %>%
  dplyr::select(Year,CHIKV.lag.1,CHIKV.lag.2,rh_mean.lag.1,rh_mean.lag.2,ah_mean.lag.1,ah_mean.lag.2,tp_mean.lag.1,tp_mean.lag.2,
                t2m_mean.lag.1,t2m_mean.lag.2,SO2CMASS_mean.lag.1,SO2CMASS_mean.lag.2, COSC_mean.lag.1,COSC_mean.lag.2,
                PM25_mean.lag.1,PM25_mean.lag.2, logP, CHIKV,Province)
CHIKV_province_list <- split(CHIKV_df,CHIKV_df$Province)

MAL_df <- df %>%
  dplyr::select(Year,MALARIA.lag.1,MALARIA.lag.2,rh_mean.lag.1,rh_mean.lag.2,ah_mean.lag.1,ah_mean.lag.2,tp_mean.lag.1,tp_mean.lag.2,
                t2m_mean.lag.1,t2m_mean.lag.2,SO2CMASS_mean.lag.1,SO2CMASS_mean.lag.2, COSC_mean.lag.1,COSC_mean.lag.2,
                PM25_mean.lag.1,PM25_mean.lag.2, logP, MALARIA,Province)
MAL_province_list <- split(MAL_df, MAL_df$Province)

ENCEP_df <- df %>%
  dplyr::select(Year,ENCEP.lag.1,ENCEP.lag.2,rh_mean.lag.1,rh_mean.lag.2,ah_mean.lag.1,ah_mean.lag.2,tp_mean.lag.1,tp_mean.lag.2,
                t2m_mean.lag.1,t2m_mean.lag.2,SO2CMASS_mean.lag.1,SO2CMASS_mean.lag.2, COSC_mean.lag.1,COSC_mean.lag.2,
                PM25_mean.lag.1,PM25_mean.lag.2, logP, ENCEP,Province)
ENCEP_province_list <- split(ENCEP_df, ENCEP_df$Province)

DENV_df <- df %>%
  dplyr::select(Year,DF.lag.1,DF.lag.2,rh_mean.lag.1,rh_mean.lag.2,ah_mean.lag.1,ah_mean.lag.2,tp_mean.lag.1,tp_mean.lag.2,
                t2m_mean.lag.1,t2m_mean.lag.2,SO2CMASS_mean.lag.1,SO2CMASS_mean.lag.2, COSC_mean.lag.1,COSC_mean.lag.2,
                PM25_mean.lag.1,PM25_mean.lag.2, logP, DF, Province)
DENV_province_list <- split(DENV_df,DENV_df$Province)

#CHIKV Province Level RR
rr_province_CHIKV<- list()
rr_province_CHIKV_var<- list()
numeric_vars <- names(CHIKV_df)[sapply(CHIKV_df,is.numeric)]
for (i in 1:76){
  for (var in numeric_vars[4:17]){
    pred_df <- na.omit(CHIKV_province_list[[i]]) %>%
      mutate(across(!c(var,CHIKV,Province),~mean(.x))) %>%
      mutate(Year = round(Year))
    
    mean_df <-  na.omit(CHIKV_province_list[[i]]) %>%
      mutate(across(!c(CHIKV,Province),~mean(.x))) %>%
      mutate(Year = round(Year))
    mean_pred <- as.vector(exp(predict(c1_fe,mean_df)))
    pred_rest <- predict(c1_fe, pred_df,se.fit = TRUE)
    mean <- exp(pred_rest$fit)
    lo <- exp(pred_rest$fit - (1.96 * pred_rest$se.fit))
    up <- exp(pred_rest$fit + (1.96* pred_rest$se.fit))
    x <- pred_df %>%
      dplyr::select(var)
    out <- cbind(mean/mean_pred,lo/mean_pred,up/mean_pred,x)
    colnames(out) <- c("mean_rr","lo_rr","up_rr", paste0(var))
    rr_province_CHIKV_var[[var]] <- as.data.frame(out) %>%
      filter(up_rr < 'Inf')
    rr_province_CHIKV[[i]] <- rr_province_CHIKV_var
  }
}

#MALARIA Province Level RR
rr_province_MAL<- list()
rr_province_MAL_var<- list()
numeric_vars <- names(MAL_df)[sapply(MAL_df,is.numeric)]
for (i in 1:76){
  for (var in numeric_vars[4:17]){
    pred_df <- na.omit(MAL_province_list[[i]]) %>%
      mutate(across(!c(var,MALARIA,Province),~mean(.x))) %>%
      mutate(Year = round(Year))
    
    mean_df <-  na.omit(MAL_province_list[[i]]) %>%
      mutate(across(!c(MALARIA,Province),~mean(.x))) %>%
      mutate(Year = round(Year))
    mean_pred <- as.vector(exp(predict(c2_fe,mean_df)))
    pred_rest <- predict(c2_fe, pred_df,se.fit = TRUE)
    mean <- exp(pred_rest$fit)
    lo <- exp(pred_rest$fit - (1.96 * pred_rest$se.fit))
    up <- exp(pred_rest$fit + (1.96* pred_rest$se.fit))
    x <- pred_df %>%
      dplyr::select(var)
    out <- cbind(mean/mean_pred,lo/mean_pred,up/mean_pred,x)
    colnames(out) <- c("mean_rr","lo_rr","up_rr", paste0(var))
    rr_province_MAL_var[[var]] <- as.data.frame(out) %>%
      filter(up_rr < 'Inf')
    rr_province_MAL[[i]] <- rr_province_MAL_var
  }
}

#ENCEP Province Level RR
rr_province_ENCEP<- list()
rr_province_ENCEP_var<- list()
numeric_vars <- names(ENCEP_df)[sapply(ENCEP_df,is.numeric)]
for (i in 1:76){
  for (var in numeric_vars[4:17]){
    pred_df <- na.omit(ENCEP_province_list[[i]]) %>%
      mutate(across(!c(var,ENCEP,Province),~mean(.x))) %>%
      mutate(Year = round(Year))
    
    mean_df <-  na.omit(ENCEP_province_list[[i]]) %>%
      mutate(across(!c(ENCEP,Province),~mean(.x))) %>%
      mutate(Year = round(Year))
    mean_pred <- as.vector(exp(predict(c3_fe,mean_df)))
    pred_rest <- predict(c3_fe, pred_df,se.fit = TRUE)
    mean <- exp(pred_rest$fit)
    lo <- exp(pred_rest$fit - (1.96 * pred_rest$se.fit))
    up <- exp(pred_rest$fit + (1.96* pred_rest$se.fit))
    x <- pred_df %>%
      dplyr::select(var)
    out <- cbind(mean/mean_pred,lo/mean_pred,up/mean_pred,x)
    colnames(out) <- c("mean_rr","lo_rr","up_rr", paste0(var))
    rr_province_ENCEP_var[[var]] <- as.data.frame(out) %>%
      filter(up_rr < 'Inf')
    rr_province_ENCEP[[i]] <- rr_province_ENCEP_var
  }
}

#DF Province Level RR
rr_province_DENV<- list()
rr_province_DENV_var<- list()
numeric_vars <- names(DENV_df)[sapply(DENV_df,is.numeric)]
for (i in 1:76){
  for (var in numeric_vars[4:17]){
    pred_df <- na.omit(DENV_province_list[[i]]) %>%
      mutate(across(!c(var,DF,Province),~mean(.x))) %>%
      mutate(Year = round(Year))
    
    mean_df <-  na.omit(DENV_province_list[[i]]) %>%
      mutate(across(!c(DF,Province),~mean(.x))) %>%
      mutate(Year = round(Year))
    mean_pred <- as.vector(exp(predict(c4_fe,mean_df)))
    pred_rest <- predict(c4_fe, pred_df,se.fit = TRUE)
    mean <- exp(pred_rest$fit)
    lo <- exp(pred_rest$fit - (1.96 * pred_rest$se.fit))
    up <- exp(pred_rest$fit + (1.96* pred_rest$se.fit))
    x <- pred_df %>%
      dplyr::select(var)
    out <- cbind(mean/mean_pred,lo/mean_pred,up/mean_pred,x)
    colnames(out) <- c("mean_rr","lo_rr","up_rr", paste0(var))
    rr_province_DENV_var[[var]] <- as.data.frame(out) %>%
      filter(up_rr < 'Inf')
    rr_province_DENV[[i]] <- rr_province_DENV_var
  }
}

save(rr_province_CHIKV, file = paste0(wd,"/CHIKV_RR_province.Rdata"))
save(rr_province_MAL, file = paste0(wd,"/MAL_RR_province.Rdata"))
save(rr_province_ENCEP, file = paste0(wd,"/ENCEP_RR_province.Rdata"))
save(rr_province_DENV, file = paste0(wd,"/DENV_RR_province.Rdata"))


# CHIKV Province Percentages
store_ci_CHIKV <-list()
for (i in 1:76){
  store_ci<- list()
  for(ii in 1:14){
    data <- rr_province_CHIKV[[i]][[ii]] %>%
      arrange(rr_province_CHIKV[[i]][[ii]][,4])%>%
      filter(rr_province_CHIKV[[i]][[ii]][,4] > mean(rr_province_CHIKV[[i]][[ii]][,4]))
    
    data2<- cbind(data[,2],data[,3])
    
    #find rrs with ci_ind > 0  and breaches rr, plot as map
    ci_ind_1 <- apply(data2, 1, function(x) as.numeric((x[1]>= 1 & x[2] >=1)))
    #find rrs with ci_ind > 0 
    ci_ind_2 <- apply(data2, 1, function(x) as.numeric((x[1]<= 1 & x[2] <=1)))
    #plot only mean and colour rr if CI doesnt touch bound
    
    if (sum(ci_ind_1) > sum(ci_ind_2)) {store_ci[[ii]] <- 1}
    if (sum(ci_ind_1) < sum(ci_ind_2)) {store_ci[[ii]] <- -1}
    if (sum(ci_ind_1) == sum(ci_ind_2)) {store_ci[[ii]] <- 0}
    
  }
  store_ci_CHIKV[[i]] <- unlist(store_ci)
}

# MAL Province Percentages
store_ci_MAL <-list()
for (i in 1:76){
  store_ci<- list()
  for(ii in 1:14){
    data <- rr_province_MAL[[i]][[ii]] %>%
      arrange(rr_province_MAL[[i]][[ii]][,4])%>%
      filter(rr_province_MAL[[i]][[ii]][,4] > mean(rr_province_MAL[[i]][[ii]][,4]))
    
    data2<- cbind(data[,2],data[,3])
    
    #find rrs with ci_ind > 0  and breaches rr, plot as map
    ci_ind_1 <- apply(data2, 1, function(x) as.numeric((x[1]>= 1 & x[2] >=1)))
    #find rrs with ci_ind > 0 
    ci_ind_2 <- apply(data2, 1, function(x) as.numeric((x[1]<= 1 & x[2] <=1)))
    #plot only mean and colour rr if CI doesnt touch bound
    
    if (sum(ci_ind_1) > sum(ci_ind_2)) {store_ci[[ii]] <- 1}
    if (sum(ci_ind_1) < sum(ci_ind_2)) {store_ci[[ii]] <- -1}
    if (sum(ci_ind_1) == sum(ci_ind_2)) {store_ci[[ii]] <- 0}
    
  }
  store_ci_MAL[[i]] <- unlist(store_ci)
}

#CHIKV Province Percentages
store_ci_CHIKV <-list()
for (i in 1:76){
  store_ci<- list()
  for(ii in 1:14){
    data <- rr_province_CHIKV[[i]][[ii]] %>%
      arrange(rr_province_CHIKV[[i]][[ii]][,4])%>%
      filter(rr_province_CHIKV[[i]][[ii]][,4] > mean(rr_province_CHIKV[[i]][[ii]][,4]))
    
    data2<- cbind(data[,2],data[,3])
    
    #find rrs with ci_ind > 0  and breaches rr, plot as map
    ci_ind_1 <- apply(data2, 1, function(x) as.numeric((x[1]>= 1 & x[2] >=1)))
    #find rrs with ci_ind > 0 
    ci_ind_2 <- apply(data2, 1, function(x) as.numeric((x[1]<= 1 & x[2] <=1)))
    #plot only mean and colour rr if CI doesnt touch bound
    
    if (sum(ci_ind_1) > sum(ci_ind_2)) {store_ci[[ii]] <- 1}
    if (sum(ci_ind_1) < sum(ci_ind_2)) {store_ci[[ii]] <- -1}
    if (sum(ci_ind_1) == sum(ci_ind_2)) {store_ci[[ii]] <- 0}
    
  }
  store_ci_CHIKV[[i]] <- unlist(store_ci)
}

#ENCEP Province Percentages

store_ci_ENCEP <-list()
for (i in 1:76){
  store_ci<- list()
  for(ii in 1:14){
    data <- rr_province_ENCEP[[i]][[ii]] %>%
      arrange(rr_province_ENCEP[[i]][[ii]][,4])%>%
      filter(rr_province_ENCEP[[i]][[ii]][,4] > mean(rr_province_ENCEP[[i]][[ii]][,4]))
    
    data2<- cbind(data[,2],data[,3])
    
    #find rrs with ci_ind > 0  and breaches rr, plot as map
    ci_ind_1 <- apply(data2, 1, function(x) as.numeric((x[1]>= 1 & x[2] >=1)))
    #find rrs with ci_ind > 0 
    ci_ind_2 <- apply(data2, 1, function(x) as.numeric((x[1]<= 1 & x[2] <=1)))
    #plot only mean and colour rr if CI doesnt touch bound
    
    if (sum(ci_ind_1) > sum(ci_ind_2)) {store_ci[[ii]] <- 1}
    if (sum(ci_ind_1) < sum(ci_ind_2)) {store_ci[[ii]] <- -1}
    if (sum(ci_ind_1) == sum(ci_ind_2)) {store_ci[[ii]] <- 0}
    
  }
  store_ci_ENCEP[[i]] <- unlist(store_ci)
}

# DENV Province Percentages
store_ci_DENV <-list()
for (i in 1:76){
  store_ci<- list()
  for(ii in 1:14){
    data <- rr_province_DENV[[i]][[ii]] %>%
      arrange(rr_province_DENV[[i]][[ii]][,4])%>%
      filter(rr_province_DENV[[i]][[ii]][,4] > mean(rr_province_DENV[[i]][[ii]][,4]))
    
    data2<- cbind(data[,2],data[,3])
    
    #find rrs with ci_ind > 0  and breaches rr, plot as map
    ci_ind_1 <- apply(data2, 1, function(x) as.numeric((x[1]>= 1 & x[2] >=1)))
    #find rrs with ci_ind > 0 
    ci_ind_2 <- apply(data2, 1, function(x) as.numeric((x[1]<= 1 & x[2] <=1)))
    #plot only mean and colour rr if CI doesnt touch bound
    
    if (sum(ci_ind_1) > sum(ci_ind_2)) {store_ci[[ii]] <- 1}
    if (sum(ci_ind_1) < sum(ci_ind_2)) {store_ci[[ii]] <- -1}
    if (sum(ci_ind_1) == sum(ci_ind_2)) {store_ci[[ii]] <- 0}
    
  }
  store_ci_DENV[[i]] <- unlist(store_ci)
}

CHIKV_ci_df <- do.call(rbind, lapply(store_ci_CHIKV, function(x) data.frame(matrix(unlist(x), ncol = 14))))
MAL_ci_df <- do.call(rbind, lapply(store_ci_MAL, function(x) data.frame(matrix(unlist(x), ncol = 14))))
ENCEP_ci_df <- do.call(rbind, lapply(store_ci_ENCEP, function(x) data.frame(matrix(unlist(x), ncol = 14))))
DENV_ci_df <- do.call(rbind, lapply(store_ci_DENV, function(x) data.frame(matrix(unlist(x), ncol = 14))))

#GAM IRR by Province Breakdown

all_rr <- list(rr_province_CHIKV,rr_province_MAL, rr_province_ENCEP, rr_province_DENV)

names_rr <-  c('Relative Humidity\n   [1-Month Lag]',
               'Relative Humidity\n   [2-Month Lag]',
               "Absolute Humidity (g/m³) \n        [1-Month Lag]",
               "Absolute Humidity (g/m³) \n        [2-Month Lag]",
               'Total Precipitation (mm)\n        [1-Month Lag]',
               'Total Precipitation (mm)\n        [2-Month Lag]',
               'Temperature (°C)\n   [1-Month Lag]',
               'Temperature (°C)\n   [2-Month Lag]',
               bquote(atop(SO[2]~'Surface Concentration (mg/m³)', ~'[1-Month Lag]')),
               bquote(atop(SO[2]~'Surface Concentration (mg/m³)', ~'[2-Month Lag]')),
               bquote(atop('CO Surface Concentration (ppb)', ~'[1-Month Lag]')),
               bquote(atop('CO Surface Concentration (ppb)', ~'[2-Month Lag]')),
               bquote(atop(PM[2.5]~'Surface Concentration (µg/m³)', ~'[1-Month Lag]')),
               bquote(atop(PM[2.5]~'Surface Concentration (µg/m³)', ~'[1-Month Lag]')))

panelNames = c('A1:','A2:','B1:','B2:',
               'C1:','C2:','D1:','D2:',
               'E1:','E2:','F1:','F2:',
               'G1:','G2:','H1:', 'H2:')
province_names = unique(df$Province)

pdf(paste0(wd,"/IRR_province_breakdown.pdf"),onefile = TRUE, width = 12, height = 12)
par(las=1,mfrow=c(4,4),cex.lab=1.4,cex.axis=1.2,mar = c(6,6,4,2) + 0.1, mgp = c(5,1,0))
j<-1
for (i in 1:4){
  a<- all_rr[[i]]
  for(ii in 1:76){
    a1<-a[[ii]]
    for(iii in 1:14){
      a2<- a1[[iii]]%>%
        arrange(a1[[iii]][,4])
      
      if (nrow(a2) == 0) {
        next
      } else{
        gamPlotter_rr(m = a2[,1],
                      l = a2[,2],
                      u = a2[,3],
                      x = a2[,4],
                      xlab = names_rr[iii],
                      ylab = paste(DiseaseNames[i]),
                      ci_col = 'steelblue4',
                      ci_col_null = tempcol_ci_null,
                      panelName = paste0(panelNames[[iii]],province_names[[ii]]))
        
      j<-j+1
      }
    }
    par(las=1,mfrow=c(4,4),cex.lab=1.4,cex.axis=1.2,mar = c(6,6,4,2) + 0.1, mgp = c(5,1,0))
  }
}


dev.off()
tiff('test.tiff',units = 'in',width = 12, height = 12, res = 350)
par(las=1,mfrow=c(4,4),cex.lab=1.4,cex.axis=1.2,mar = c(6,6,4,2) + 0.1, mgp = c(5,1,0))

for(iii in 1:14){
  a2<-a1[[iii]] %>%
    arrange(a1[[iii]][,4])
  gamPlotter_rr(m = a2[,1],
                l = a2[,2],
                u = a2[,3],
                x = a2[,4],
                xlab = names_rr[iii],
                ylab = paste(DiseaseNames[k]),
                ci_col = 'steelblue4',
                ci_col_null = tempcol_ci_null,
                panelName = paste0(panelNames[[iii]],province_names[[1]]))
}

de

for (i in 1:4){
  for(ii in 1:76){
    for(iii in 1:14){
      a<- all_rr[[i]][[ii]][[iii]]%>%
        arrange(all_rr[[i]][[ii]][[iii]][,4])}}}