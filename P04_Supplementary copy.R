rm(list=ls())
library(tidyverse)
library(stringr)
library(mgcv)
library(MASS)
library(ggplot2)
library(patchwork)

setwd("/Users/pranav/Documents/JT_Paper/Code")
wd <- "/Users/pranav/Documents/JT_Paper/Code"


#Load dataset
load(paste0(wd,'/df.Rdata'))


#KDE Plots
#Disease indexes - 4:7
panel_names <- c("A1",'A2','B1','B2')
disease<- c("Chikungunya","Malaria","Japanese Encephalitis","Dengue Fever")
tiff('KDE_plots.tiff',units = 'in',width = 8, height = 8, res = 350)
par(las=1,mfrow=c(2,2),cex.lab=1.2,cex.axis=0.85, mar = c(5,6,4,2) + 0.1)
for (i in 4:7){
  d <- density(df[,i], na.rm = TRUE)
  plot(d, main = disease[i-3], frame = FALSE,col = 'tomato3', ylab = 'Density')
  mtext(text=panel_names[i-3],side=3,adj=0)
  box()
}

dev.off()

#Table 1 (Summary Statistics)
df_incidence <- df %>%
  dplyr::select(CHIKV,MALARIA,ENCEP,DF,Pop)%>%
  mutate(CHIKV_inc = (CHIKV/Pop) * 100000,
         MALARIA_inc = (MALARIA/Pop) * 100000,
         ENCEP_inc = (ENCEP/Pop) * 100000,
         DF_inc = (DF/Pop) * 100000)

# Convert pollutants to microgram

df<-df %>%
  mutate(across(starts_with('SO2CMASS'), ~.x *10)) %>%
  mutate(across(starts_with('SO2CMASS'), ~.x /1000)) #mg/m3

#Distribution of variables
panel_names_dist <- c('A', 'B','C','D','E','F','G')
var_names <- c('Relative Humidity',
               "Absolute Humidity (g/m³)",
               'Total Precipitation (mm)',
               'Temperature (°C)',
               expression("SO"[2]*" Surface Concentration (mg/m³)"),
               'CO Surface Concentration (ppb)',
               expression("PM"[2.5]*' Surface Concentration (µg/m³)'))
tiff('var_dist.tiff',units = 'in',width = 12, height = 8, res = 350)              
par(las=1,mfrow=c(2,4),cex.lab=1.4,cex.axis=1.2, mar = c(6,6,4,2) + 0.1, mgp = c(4,1,0))
for (i in 8:14){
  hist(df[,i], col = 'azure3', border = 'white', breaks = 20, main = '', xlab = var_names[i-7], freq = FALSE)
  lines(density(df[,i]), col = 'tomato3')
  mtext(text = panel_names_dist[i-7], side = 3, adj = 0)
}

dev.off()

#Kolmogorov-Smirnov Test (Test for normal distribution)

ks_results <- list()

for ( i in 8:14){
  print(ks.test(unique(df[,i]), 'pnorm'))
}




