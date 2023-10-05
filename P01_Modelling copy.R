rm(list=ls())
library(tidyverse)
library(stringr)
library(mgcv)
library(MASS)

setwd("/Users/pranav/Documents/JT_Paper/Code")
wd <- "/Users/pranav/Documents/JT_Paper/Code"

#Load big dataset (RDS)

store_impute <- do.call(rbind,store_impute)

#Convert normalized case counts to original case counts


colnames(pop_data) <- sub('X', "", colnames(pop_data)) # Remove X from column names

#Population Values
pop_pivot <- pivot_longer(pop_data, cols = !Province, names_to = 'Year', values_to = 'Pop')

#Case counts
case_counts <- read.csv(paste0(wd,"/merged.csv"),header = TRUE) %>%
  dplyr::select(CHIKV,MALARIA,ENCEP,DF,Year,Month,Province)

store_impute_new <- merge(store_impute,pop_pivot) %>%
  dplyr::select(-CHIKV,-MALARIA,-ENCEP,-DF)

store_impute_final <- merge(store_impute_new,case_counts) %>%
  mutate(CHIKV.lag.1 = lag(CHIKV,1),
         CHIKV.lag.2 = lag(CHIKV,2),
         MALARIA.lag.1 = lag(MALARIA,1),
         MALARIA.lag.2 = lag(MALARIA,2),
         ENCEP.lag.1 = lag(ENCEP,1),
         ENCEP.lag.2 = lag(ENCEP,2),
         DF.lag.1 = lag(DF,1),
         DF.lag.2 = lag(DF,2),
         logP = log(Pop))

#Final dataframe for inputting into models
df<- store_impute_final %>%
  dplyr::select(Province,Year,Month,CHIKV,MALARIA,ENCEP,DF,rh_mean,ah_mean, tp_mean, t2m_mean,SO2CMASS_mean,COSC_mean,PM25_mean,Pop, logP) %>%
  mutate(tp_mean = tp_mean * 1000,  # Scale to mm
          SO2CMASS_mean = SO2CMASS_mean * 1000000000, # Scale to x10^-9
          PM25_mean = PM25_mean *1000000000)  %>% #Scale to x10^-9
  mutate(across(CHIKV:PM25_mean,~lag(.x,n=1), .names = "{.col}.lag.1")) %>%
  mutate(across(CHIKV:PM25_mean, ~lag(.x,n=2), .names = '{.col}.lag.2'))

save(df, file = paste0(wd,'/df.Rdata'))

#Modelling

disease <- c("CHIKV","MALARIA","ENCEP","DF")
disease_controls <- list() 
for (i in 1:length(disease)) {
  disease_controls[[i]] <- paste0(disease[i] ,"~", disease[i],".lag.1+", disease[i],".lag.2")} #Set up formula for case counts

weather_controls <- "rh_mean.lag.1 + rh_mean.lag.2 + ah_mean.lag.1 + ah_mean.lag.2 + tp_mean.lag.1 + tp_mean.lag.2 + t2m_mean.lag.1 +
t2m_mean.lag.2"
glm_pollutants <- "SO2CMASS_mean.lag.1 +SO2CMASS_mean.lag.2 + COSC_mean.lag.1+COSC_mean.lag.2+PM25_mean.lag.1+PM25_mean.lag.2"
gam_pollutants <- "s(SO2CMASS_mean.lag.1) +s(SO2CMASS_mean.lag.2) +s(COSC_mean.lag.1)+s(COSC_mean.lag.2)+ s(PM25_mean.lag.1)+
s(PM25_mean.lag.2)"
gam_controls <- "s(rh_mean.lag.1) + s(rh_mean.lag.2) + 
                s(ah_mean.lag.1) + s(ah_mean.lag.2) + 
                s(tp_mean.lag.1) + s(tp_mean.lag.2) + 
                s(t2m_mean.lag.1) + s(t2m_mean.lag.2) +
                s(SO2CMASS_mean.lag.1) +s(SO2CMASS_mean.lag.2) +
                s(COSC_mean.lag.1)+s(COSC_mean.lag.2)+
                s(PM25_mean.lag.1)+s(PM25_mean.lag.2)"

gam_controls_penalized <- "s(rh_mean.lag.1,bs='ts') + s(rh_mean.lag.2,bs='ts') + 
                s(ah_mean.lag.1,bs='ts') + s(ah_mean.lag.2,bs='ts') + 
                s(tp_mean.lag.1,bs='ts') + s(tp_mean.lag.2,bs='ts') + 
                s(t2m_mean.lag.1,bs='ts') + s(t2m_mean.lag.2,bs='ts') +
                s(SO2CMASS_mean.lag.1,bs='ts') +s(SO2CMASS_mean.lag.2,bs='ts') +
                s(COSC_mean.lag.1,bs='ts')+s(COSC_mean.lag.2,bs='ts')+
                s(PM25_mean.lag.1,bs='ts')+s(PM25_mean.lag.2,bs='ts')"

#Linear Models
a1_f <- glm.nb(formula(paste0(disease_controls[[1]], "+",weather_controls ,"+", glm_pollutants, '+ offset(logP)')),data=df,
               na.action = na.omit)
a2_f <- glm(formula(paste0(disease_controls[[2]], "+",weather_controls ,"+", glm_pollutants, '+ offset(logP)')),data=df,
               na.action = na.omit, family = poisson(link = 'log'))
a3_f <- glm.nb(formula(paste0(disease_controls[[3]], "+",weather_controls ,"+", glm_pollutants, '+ offset(logP)')),data=df,
               na.action = na.omit)
a4_f <- glm.nb(formula(paste0(disease_controls[[4]], "+",weather_controls ,"+", glm_pollutants, '+ offset(logP)')),data=df,
               na.action = na.omit)


#gam with spline terms on only pollutants
b1_f <- gam(formula(paste0(disease_controls[[1]], "+",weather_controls ,"+", gam_pollutants, '+ offset(logP)')),data=df,control=list(keepData=T),
            family = nb(),method='REML')
b2_f <- gam(formula(paste0(disease_controls[[2]], "+",weather_controls ,"+", gam_pollutants, '+ offset(logP)')),data=df,control=list(keepData=T),
            family = nb(),method='REML')
b3_f <- gam(formula(paste0(disease_controls[[3]], "+",weather_controls ,"+", gam_pollutants, ' + offset(logP)')),data=df,control=list(keepData=T),
            family = nb(),method='REML')
b4_f <- gam(formula(paste0(disease_controls[[4]], "+",weather_controls ,"+", gam_pollutants, ' + offset(logP)')),data=df,control=list(keepData=T),
            family = nb(),method='REML')

#gam with spline terms on pollutant and meteorological variables
c1_f <- gam(formula(paste0(disease_controls[[1]], "+",gam_controls, '+ offset(logP)')),data=df,control=list(keepData=T),
            family = nb(),method='REML', na.action = na.omit)
c2_f <- gam(formula(paste0(disease_controls[[2]], "+",gam_controls, '+ offset(logP)')),data=df,control=list(keepData=T),
            family = nb(),method='REML', na.action = na.omit)
c3_f <- gam(formula(paste0(disease_controls[[3]], "+",gam_controls, '+ offset(logP)')),data=df,control=list(keepData=T),
            family = nb(),method='REML', na.action = na.omit)
c4_f <- gam(formula(paste0(disease_controls[[4]], "+",gam_controls, '+ offset(logP)')),data=df,control=list(keepData=T),
            family = nb(),method='REML', na.action = na.omit)

#gam with spline terms, all penalized
d1_f <-  gam(formula(paste0(disease_controls[[1]], "+",gam_controls_penalized, '+ offset(logP)')),data=df,control=list(keepData=T),
             family = nb(),method='REML',na.action = na.omit)
d2_f <-  gam(formula(paste0(disease_controls[[2]], "+",gam_controls_penalized, '+ offset(logP)')),data=df,control=list(keepData=T),
             family = nb(),method='REML',na.action = na.omit)
d3_f <-  gam(formula(paste0(disease_controls[[3]], "+",gam_controls_penalized, '+ offset(logP)')),data=df,control=list(keepData=T),
             family = nb(),method='REML',na.action = na.omit)
d4_f <-  gam(formula(paste0(disease_controls[[4]], "+",gam_controls_penalized, '+ offset(logP)')),data=df,control=list(keepData=T),
             family = nb(),method='REML',na.action = na.omit)

#Fixed effects OLS

df$Province <- as.factor(df$Province)

a1_fe <- glm.nb(formula(paste0(disease_controls[[1]], "+",weather_controls ,"+", glm_pollutants, "+","Province",'+ offset(logP)')),
             data=df, na.action = na.omit)
a2_fe <- glm.nb(formula(paste0(disease_controls[[2]], "+",weather_controls ,"+", glm_pollutants, "+","Province",'+ offset(logP)')),
                data=df, na.action = na.omit)
a3_fe <- glm.nb(formula(paste0(disease_controls[[3]], "+",weather_controls ,"+", glm_pollutants, "+","Province",'+ offset(logP)')),
                data=df, na.action = na.omit)
a4_fe <- glm.nb(formula(paste0(disease_controls[[4]], "+",weather_controls ,"+", glm_pollutants, "+","Province",'+ offset(logP)')),
                data=df, na.action = na.omit)

#gam with spline terms on only pollutants
b1_fe <- gam(formula(paste0(disease_controls[[1]], "+",weather_controls ,"+", gam_pollutants,"+","Province", '+ offset(logP)')),
             data=df, family = nb(),control=list(keepData=T),method='REML', na.action = na.omit)
b2_fe <- gam(formula(paste0(disease_controls[[2]], "+",weather_controls ,"+", gam_pollutants,"+","Province", '+ offset(logP)')),
             data=df, family = nb(),control=list(keepData=T),method='REML', na.action = na.omit)
b3_fe <- gam(formula(paste0(disease_controls[[3]], "+",weather_controls ,"+", gam_pollutants,"+","Province", '+ offset(logP)')),
             data=df, family = nb(),control=list(keepData=T),method='REML', na.action = na.omit)
b4_fe <- gam(formula(paste0(disease_controls[[4]], "+",weather_controls ,"+", gam_pollutants,"+","Province", '+ offset(logP)')),
             data=df, family = nb(),control=list(keepData=T),method='REML', na.action = na.omit)

#gam with spline terms on pollutant and meteorological variables
c1_fe <- gam(formula(paste0(disease_controls[[1]], "+", gam_controls,"+","Province", "+ offset(logP)")),data=df,
             control=list(keepData=T),family = nb(),method='REML', na.action = na.omit)
c2_fe <- gam(formula(paste0(disease_controls[[2]], "+", gam_controls,"+","Province", "+ offset(logP)")),data=df,
             control=list(keepData=T),family = nb(),method='REML', na.action = na.omit)
c3_fe <- gam(formula(paste0(disease_controls[[3]], "+", gam_controls,"+","Province", "+ offset(logP)")),data=df,
             control=list(keepData=T),family = nb(),method='REML', na.action = na.omit)
c4_fe <- gam(formula(paste0(disease_controls[[4]], "+", gam_controls,"+","Province", "+ offset(logP)")),data=df,
             control=list(keepData=T),family = nb(),method='REML', na.action = na.omit)

#gam with spline terms, all penalized
d1_fe <-  gam(formula(paste0(disease_controls[[1]], "+", gam_controls_penalized,"+","Province", ' + offset(logP)')),data=df,
              family = nb(),control=list(keepData=T),method='REML', na.action = na.omit)
d2_fe <-  gam(formula(paste0(disease_controls[[2]], "+", gam_controls_penalized,"+","Province", ' + offset(logP)')),data=df,
              family = nb(),control=list(keepData=T),method='REML', na.action = na.omit)
d3_fe <-  gam(formula(paste0(disease_controls[[3]], "+", gam_controls_penalized,"+","Province", ' + offset(logP)')),data=df,
              family = nb(),control=list(keepData=T),method='REML', na.action = na.omit)
d4_fe <-  gam(formula(paste0(disease_controls[[4]], "+", gam_controls_penalized,"+","Province", ' + offset(logP)')),data=df,
              family = nb(),control=list(keepData=T),method='REML', na.action = na.omit)

#GAMM with Random Effects

#GAMM with spline terms on only pollutants
b1_re <- gamm(formula(paste0(disease_controls[[1]], "+",weather_controls ,"+", gam_pollutants, ' + offset(logP)')),
              random=list(Province=~1),data=df,control=list(keepData=T),family = poisson(link = 'log'),method='REML', na.action = na.omit)
b2_re <- gamm(formula(paste0(disease_controls[[2]], "+",weather_controls ,"+", gam_pollutants, ' + offset(logP)')),
              random=list(Province=~1),data=df,control=list(keepData=T),family = poisson(link = 'log'),method='REML', na.action = na.omit)
b3_re <- gamm(formula(paste0(disease_controls[[3]], "+",weather_controls ,"+", gam_pollutants, ' + offset(logP)')),
              random=list(Province=~1),data=df,control=list(keepData=T),family = poisson(link = 'log'),method='REML', na.action = na.omit)
b4_re <- gamm(formula(paste0(disease_controls[[4]], "+",weather_controls ,"+", gam_pollutants, ' + offset(logP)')),
              random=list(Province=~1),data=df,control=list(keepData=T),family = poisson(link = 'log'),method='REML', na.action = na.omit)

#GAMM with spline terms on pollutant and meteorological variables
c1_re <- gamm(formula(paste0(disease_controls[[1]], "+",gam_controls, "+ offset(logP)")),random=list(Province=~1),
              data=df,control=list(keepData=T), family = poisson(link = 'log'),method='REML', na.action = na.omit)
c2_re <- gamm(formula(paste0(disease_controls[[2]], "+",gam_controls, "+ offset(logP)")),random=list(Province=~1),
              data=df,control=list(keepData=T), family = poisson(link = 'log'),method='REML', na.action = na.omit)
c3_re <- gamm(formula(paste0(disease_controls[[3]], "+",gam_controls, "+ offset(logP)")),random=list(Province=~1),
              data=df,control=list(keepData=T), family = poisson(link = 'log'),method='REML', na.action = na.omit)
c4_re <- gamm(formula(paste0(disease_controls[[4]], "+",gam_controls, "+ offset(logP)")),random=list(Province=~1),
              data=df,control=list(keepData=T), family = poisson(link = 'log'),method='REML', na.action = na.omit)

#GAMM with spline terms, all penalized
d1_re <-  gamm(formula(paste0(disease_controls[[1]], "+",gam_controls_penalized, '+ offset(logP)')),random=list(Province=~1),
               data=df,control=list(keepData=T),family = poisson(link = 'log'),method='REML', na.action = na.omit)
d2_re <-  gamm(formula(paste0(disease_controls[[2]], "+",gam_controls_penalized, '+ offset(logP)')),random=list(Province=~1),
               data=df,control=list(keepData=T),family = poisson(link = 'log'),method='REML', na.action = na.omit)
d3_re <-  gamm(formula(paste0(disease_controls[[3]], "+",gam_controls_penalized, '+ offset(logP)')),random=list(Province=~1),
               data=df,control=list(keepData=T),family = poisson(link = 'log'),method='REML', na.action = na.omit)
d4_re <-  gamm(formula(paste0(disease_controls[[4]], "+",gam_controls_penalized, '+ offset(logP)')),random=list(Province=~1),
               data=df,control=list(keepData=T),family = poisson(link = 'log'),method='REML', na.action = na.omit)


out_full_pooled <- list(a1_f=a1_f,
                        a2_f=a2_f,
                        a3_f=a3_f,
                        a4_f=a4_f,
                        b1_f=b1_f,
                        b2_f=b2_f,
                        b3_f=b3_f,
                        b4_f=b4_f,
                        c1_f=c1_f,
                        c2_f=c2_f,
                        c3_f=c3_f,
                        c4_f=c4_f,
                        d1_f=d1_f,
                        d2_f=d2_f,
                        d3_f=d3_f,
                        d4_f=d4_f)

save(out_full_pooled, file = paste0(wd,"/pooled.Rdata"))

out_full_fixed <- list(a1_fe=a1_fe,
                       a2_fe=a2_fe,
                       a3_fe=a3_fe,
                       a4_fe=a4_fe,
                       b1_fe=b1_fe,
                       b2_fe=b2_fe,
                       b3_fe=b3_fe,
                       b4_fe=b4_fe,
                       c1_fe=c1_fe,
                       c2_fe=c2_fe,
                       c3_fe=c3_fe,
                       c4_fe=c4_fe,
                       d1_fe=d1_fe,
                       d2_fe=d2_fe,
                       d3_fe=d3_fe,
                       d4_fe=d4_fe)
save(out_full_fixed,file = paste0(wd,"/fixed.Rdata"))