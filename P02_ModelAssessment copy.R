rm(list=ls())
library(tidyverse)
library(stringr)
library(mgcv)
library(MASS)

setwd("/Users/pranav/Documents/JT_Paper/Code")
wd <- "/Users/pranav/Documents/JT_Paper/Code"

#Load dataset and models
load(paste0(wd,'/df.Rdata'))
load(paste0(wd,'/pooled.Rdata'))
load(paste0(wd,'/fixed.Rdata'))

#Model AICs

a <- unlist(lapply(out_full_pooled,AIC))
b <- unlist(lapply(out_full_fixed,AIC))
tabl <- cbind(a,b)
tabl <- tabl[c(1,5,9,13,
               2,6,10,14,
               3,7,11,15,
               4,8,12,16),]
colnames(tabl) <- c('AIC Pooled','AIC FE')
rownames(tabl) <- c(paste(c('Linear','GAM Pollutants','GAM All','GAM All Penalized'),"CHIKV"),
                    paste(c('Linear','GAM Pollutants','GAM All','GAM All Penalized'),"MALARIA"),
                    paste(c('Linear','GAM Pollutants','GAM All','GAM All Penalized'),'ENCEP'),
                    paste(c('Linear','GAM Pollutants','GAM All','GAM All Penalized'),'DF'))

#GAM ALL and GAM All Penalized similar across all diseases
CHIKV_mod <- out_full_fixed[["c1_fe"]]
MAL_mod <- out_full_fixed[['c2_fe']]
ENCEP_mod <- out_full_fixed[['c3_fe']]
DENV_mod <- out_full_fixed [['c4_fe']]

##RRs for all variables

disease<- c("CHIKV","MALARIA","ENCEP","DF")
disease_variables <- list()

for (i in 1:length(disease)){
  disease_variables[[i]] <- paste0(disease[i],".lag.1,",disease[i],".lag.2,")
}
weather_variables <- "rh_mean.lag.1,rh_mean.lag.2,ah_mean.lag.1,ah_mean.lag.2,tp_mean.lag.1,tp_mean.lag.2,t2m_mean.lag.1,t2m_mean.lag.2"
pollutant_variables <- "SO2CMASS_mean.lag.1,SO2CMASS_mean.lag.2, COSC_mean.lag.1,COSC_mean.lag.2,PM25_mean.lag.1,PM25_mean.lag.2"


CHIKV_df <- df %>%
  dplyr::select(Year,CHIKV.lag.1,CHIKV.lag.2,rh_mean.lag.1,rh_mean.lag.2,ah_mean.lag.1,ah_mean.lag.2,tp_mean.lag.1,tp_mean.lag.2,
                t2m_mean.lag.1,t2m_mean.lag.2,SO2CMASS_mean.lag.1,SO2CMASS_mean.lag.2, COSC_mean.lag.1,COSC_mean.lag.2,
                PM25_mean.lag.1,PM25_mean.lag.2, logP, CHIKV,Province)

MAL_df <- df %>%
  dplyr::select(Year,MALARIA.lag.1,MALARIA.lag.2,rh_mean.lag.1,rh_mean.lag.2,ah_mean.lag.1,ah_mean.lag.2,tp_mean.lag.1,tp_mean.lag.2,
                  t2m_mean.lag.1,t2m_mean.lag.2,SO2CMASS_mean.lag.1,SO2CMASS_mean.lag.2, COSC_mean.lag.1,COSC_mean.lag.2,
                  PM25_mean.lag.1,PM25_mean.lag.2, logP, MALARIA,Province)

ENCEP_df <- df %>%
  dplyr::select(Year,ENCEP.lag.1,ENCEP.lag.2,rh_mean.lag.1,rh_mean.lag.2,ah_mean.lag.1,ah_mean.lag.2,tp_mean.lag.1,tp_mean.lag.2,
         t2m_mean.lag.1,t2m_mean.lag.2,SO2CMASS_mean.lag.1,SO2CMASS_mean.lag.2, COSC_mean.lag.1,COSC_mean.lag.2,
         PM25_mean.lag.1,PM25_mean.lag.2, logP, ENCEP,Province)
DENV_df <- df %>%
  dplyr::select(Year,DF.lag.1,DF.lag.2,rh_mean.lag.1,rh_mean.lag.2,ah_mean.lag.1,ah_mean.lag.2,tp_mean.lag.1,tp_mean.lag.2,
                t2m_mean.lag.1,t2m_mean.lag.2,SO2CMASS_mean.lag.1,SO2CMASS_mean.lag.2, COSC_mean.lag.1,COSC_mean.lag.2,
                PM25_mean.lag.1,PM25_mean.lag.2, logP, DF, Province)

#CHIKV RR

store_rr_CHIKV <- list()
numeric_vars <- names(CHIKV_df)[sapply(CHIKV_df,is.numeric)]

for (var in numeric_vars[4:17]) {
  pred_df <- na.omit(CHIKV_df) %>%
    mutate(across(!c(var,CHIKV,Province),~mean(.x))) %>%
    mutate(Year = round(Year))
  
  mean_df <-  na.omit(CHIKV_df) %>%
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
  
  store_rr_CHIKV[[var]] <- as.data.frame(out) %>%
    filter(up_rr < 'Inf')
  
}

#MALARIA RR


store_rr_MAL <- list()
numeric_vars <- names(MAL_df)[sapply(MAL_df,is.numeric)]

for (var in numeric_vars[4:17]) {
  pred_df <- na.omit(MAL_df) %>%
    mutate(across(!c(var,MALARIA,Province),~mean(.x))) %>%
    mutate(Year = round(Year))
  
  mean_df <-  na.omit(MAL_df) %>%
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
  
  store_rr_MAL[[var]] <- as.data.frame(out) %>%
    filter(up_rr < 'Inf')
  
}


#ENCEP RR

store_rr_ENCEP <- list()
numeric_vars <- names(ENCEP_df)[sapply(ENCEP_df,is.numeric)]

for (var in numeric_vars[4:17]) {
  pred_df <- na.omit(ENCEP_df) %>%
    mutate(across(!c(var,ENCEP,Province),~mean(.x))) %>%
    mutate(Year = round(Year))
  
  mean_df <-  na.omit(ENCEP_df) %>%
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
  
  store_rr_ENCEP[[var]] <- as.data.frame(out) %>%
    filter(up_rr < 'Inf')
  
}


#DENV RR

store_rr_DENV <- list()
numeric_vars <- names(DENV_df)[sapply(DENV_df,is.numeric)]

for (var in numeric_vars[4:17]) {
  pred_df <- na.omit(DENV_df) %>%
    mutate(across(!c(var,DF,Province),~mean(.x))) %>%
    mutate(Year = round(Year))
  
  mean_df <-  na.omit(DENV_df) %>%
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
  
  store_rr_DENV[[var]] <- as.data.frame(out) %>%
    filter(up_rr < 'Inf')
  
}

save(store_rr_CHIKV, file = paste0(wd,"/CHIKV_RR.Rdata"))
save(store_rr_MAL, file = paste0(wd,"/MAL_RR.Rdata"))
save(store_rr_ENCEP, file = paste0(wd,"/ENCEP_RR.Rdata"))
save(store_rr_DENV, file = paste0(wd,"/DENV_RR.Rdata"))