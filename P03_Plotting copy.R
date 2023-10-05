rm(list=ls())
library(tidyverse)
library(stringr)
library(mgcv)
library(MASS)
library(ggplot2)
library(patchwork)

setwd("/Users/pranav/Documents/JT_Paper/Code")
wd <- "/Users/pranav/Documents/JT_Paper/Code"

#RR plot function
tempcol_ci_null <- col2rgb('azure2')
tempcol_ci_null <- rgb(tempcol_ci_null[1],tempcol_ci_null[2],tempcol_ci_null[3],max = 255,alpha=0.5*255)
panelNames = c('A1','A2','A3','A4',
               'B1','B2','B3','B4',
               'C1','C2','C3','C4',
               'D1','D2','D3','D4')
pol_panelNames =  c('A1','A2','A3','A4','A5','A6',
                    'B1','B2','B3','B4','B5','B6',
                    'C1','C2','C3','C4','C5','C6',
                    'D1','D2','D3','D4','D5','D6')
gamPlotter_rr <- function(m,
                          l,
                          u,
                          x,
                          xlab,
                          ylab,
                          ci_col=tempcol_ci,#'azure3',
                          ci_col_null=tempcol_ci_null,#'azure2',
                          panelName){
  df <- cbind(l,u)
  #find CI values which cross 1 
  ci_ind <- apply(df, 1, function(x) as.numeric(!(x[1]<= 1 & x[2] >=1)))
  
  if (max(df[,2])<=20) {
    plot(x=c(min(x,na.rm=T),max(x,na.rm=T)),
         y=c(min(df),max(df)),
         col="white",
         ylab=ylab,
         xlab=xlab)  
  } else{
    plot(x=c(min(x,na.rm=T),max(x,na.rm=T)),
         y=c(0,20),
         col="white",
         ylab=ylab,
         xlab=xlab)  
  }
                                                         #Plotting Area
  
  polygon(x=c(x,max(x),rev(x),min(x)),
          y=c(l,u[length(u)],rev(u),l[1]),col=ci_col_null,border=ci_col_null)   #Entire CI
  #draw regions which do not overlap with CI
  polys <- rle(ci_ind)
  for (i in which(polys$values==1)){
    if (i == 1){
      start_ind = 1
    } else {
      start_ind = cumsum((polys$lengths))[i-1]
      
    }
    end_ind <-  cumsum((polys$lengths))[i]
    
    x_sub <- x[start_ind:end_ind]
    u_sub <- u[start_ind:end_ind]
    l_sub <- l[start_ind:end_ind]
    
    
    polygon(x=c(x_sub,max(x_sub),rev(x_sub),min(x_sub)),
            y=c(l_sub,u_sub[length(u_sub)],rev(u_sub),l_sub[1]),col=ci_col,border=ci_col)
    
  }
  
  tempcol <- col2rgb('orange2')
  tempcol <- rgb(tempcol[1],tempcol[2],tempcol[3],max = 255,alpha=0.5*255)
  abline(v=mean(x,na.rm=T),col='orange2',lwd=2)
  
  tempcol <- col2rgb('orange2')
  tempcol <- rgb(tempcol[1],tempcol[2],tempcol[3],max = 255,alpha=0.5*255)
  
  
  abline(h=1,col="darkgrey",lty=2)
  
  lines(x=x,y=m,lwd=2)
  
  legend(x="topright",
  bty='n',
  legend=c("95% CI","IRR","Mean\n(Exposure)"),
  pch=c(15,NA,NA),
  lty=c(NA,1,1),
  lwd=2,
  col=c('steelblue4',"black","orange2"))  
  
  mtext(text=panelName,side=3,adj=0)
  
}


#Load RR dataframes
load(paste0(wd,"/CHIKV_RR_2.Rdata"))
load(paste0(wd,"/MAL_RR_2.Rdata"))
load(paste0(wd,"/ENCEP_RR_2.Rdata"))
load(paste0(wd,"/DENV_RR_2.Rdata"))

# Prepare for loops
disease_list <- list(store_rr_CHIKV,store_rr_MAL,store_rr_ENCEP,store_rr_DENV)
DiseaseNames = c("Chikugunya","Malaria","Encephalitis","Dengue Fever")
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

#RH and AH plots - Indexes = 1:4
tiff('rhah_irr_3.tiff',units = 'in',width = 12, height = 12, res = 350)
par(las=1,mfrow=c(4,4),cex.lab=1.4,cex.axis=1.2,mar = c(6,6,4,2) + 0.1, mgp = c(4,1,0))

j<-1
for (k in 1:4){
  a <- disease_list[[k]]
  for (i in 1:4){
    a1<-a[[i]]%>%
      arrange(a[[i]][,4])

    gamPlotter_rr(m=a1[,1],
                  l=a1[,2],
                  u=a1[,3],
                  x=a1[,4],
                  xlab = names_rr[i],
                  ylab = paste(DiseaseNames[k]),
                  ci_col='steelblue4',
                  ci_col_null=tempcol_ci_null,
                  panelName=paste0(panelNames[j]))
    j <- j+1
  }
}
legend(x="topright",
       bty='n',
       legend=c("95% CI","IRR","Mean\n(Exposure)"),
       pch=c(15,NA,NA),
       lty=c(NA,1,1),
       lwd=2,
       col=c('steelblue4',"black","orange2"))  
dev.off()

#Precipitation and Temperature Plots - Indexes = 5:8
tiff('tptemp_irr_3.tiff',units = 'in',width = 12, height = 12, res = 350)
par(las=1,mfrow=c(4,4),cex.lab=1.4,cex.axis=1.2,mar = c(6,6,4,2) + 0.1, mgp = c(4,1,0))
j<-1
for(k in 1:4){
  a <- disease_list[[k]]
  for(i in 5:8){
    a1 <- a[[i]] %>%
      arrange(a[[i]][,4])
    
    gamPlotter_rr(m = a1[,1],
                  l = a1[,2],
                  u = a1[,3],
                  x = a1[,4],
                  xlab = names_rr[i],
                  ylab = paste(DiseaseNames[k]),
                  ci_col = 'steelblue4',
                  ci_col_null = tempcol_ci_null,
                  panelName = paste0(panelNames[j]))
    j <- j+1
  }
}
legend(x="topright",
       bty='n',
       legend=c("95% CI","IRR","Mean\n(Exposure)"),
       pch=c(15,NA,NA),
       lty=c(NA,1,1),
       lwd=2,
       col=c('steelblue4',"black","orange2"))  
dev.off()


# Pollutant Plots - Indexes = 9:14
tiff('pol_irr_3.tiff',units = 'in',width = 18, height = 12, res = 350)
par(las=1,mfrow=c(4,6),cex.lab=1.4,cex.axis=1.2,mar = c(6,6,4,2) + 0.1, mgp = c(5,1,0))
j<-1
for (k in 1:4){
  a <- disease_list[[k]]
  for (i in 9:14){
    a1<-a[[i]]%>%
      arrange(a[[i]][,4])
    
    gamPlotter_rr(m=a1[,1],
                  l=a1[,2],
                  u=a1[,3],
                  x=a1[,4],
                  xlab = names_rr[i],
                  ylab = paste(DiseaseNames[k]),
                  ci_col='steelblue4',
                  ci_col_null=tempcol_ci_null,
                  panelName=paste0(pol_panelNames[j]))
    j <- j+1
  }
}
legend(x="topright",
       bty='n',
       legend=c("95% CI","IRR","Mean\n(Exposure)"),
       pch=c(15,NA,NA),
       lty=c(NA,1,1),
       lwd=2,
       col=c('steelblue4',"black","orange2"))  
dev.off()


