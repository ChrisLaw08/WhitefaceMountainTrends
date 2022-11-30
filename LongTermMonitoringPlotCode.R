#Packages needed to run code ####
library(gridExtra)
library(tidyverse)
#library(ggpmisc)
library(latex2exp)
library(ggpubr)
library(gtsummary)
library(deming)
library(Kendall)
library(scales)
library(rstatix)
library(patchwork)
library(IgorR)
library(vroom)
library(tibbletime)
library(lubridate)
library(plotrix)

##To Show available colors
show_col(hue_pal()(6))
colsgg<-(hue_pal()(6))

### Data Clean Up ####
AllCloudData<-read.csv("AllCloudandMetData.csv", header = TRUE)


PaperQuality<-function(...){
  theme_bw()%+replace%
    theme(axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 22, angle = 90),
          axis.text = element_text(size = 22), title = element_text(size = 24, face = "bold"),
          legend.text = element_text(size = 18), legend.title = element_blank(),
          legend.position = "bottom")
}

AllCloudData<-AllCloudData%>%
  dplyr::mutate(TOC = WSOC)%>%
  filter(Year > 1993 & Year<=2021)%>%
  mutate(date = as.POSIXct(date))%>%
  mutate(TOC = case_when(Year > 2017 ~ TOC*(1/.8486), 
                                  TRUE~TOC))%>%
  mutate(Cations = 1e6*10^(-LABPH)+CA+MG+K+Sodium+NH4,
         Anions = NO3+SO4+CL,
         CatAnRat = Cations/Anions)%>%### Create Columns with Cations and Anions ueq/L
  mutate(RPD = 200*(Cations-Anions)/(Cations+Anions),
         Hplus = 1e6*10^(-LABPH))%>% ## RPD Calculations
  mutate(Ratio = Cations/Anions)%>%
  mutate(HCond = (Hplus*349.5)/1000, Regime = ifelse(HCond/SPCOND < 0.35, "Non-Linear", "Linear"),
         PredCond = (Hplus*349.5 + Sodium*50.1+73.5*K+73.5*NH4+59.5*CA+53.5*MG+80*SO4+76.35*CL+71.46*NO3)/1000)%>%
  mutate(Class = ifelse((Cations < 100 & Anions < 100 & abs(RPD) < 100) | (Cations > 100 & abs(RPD) < 25) |(Anions > 100 & abs(RPD)<25), "Valid", "Invalid"))%>%
  mutate(Class = as.factor(Class))%>%
  mutate(Hplus = 10^(-LABPH))%>%
  mutate(SurplusNH4 = NH4-NO3-SO4,
         NH4Surplus = ifelse(NH4 > NO3+SO4, 'Yes', 'No'),
         BUHplus = SO4+NO3-NH4)%>%
  mutate(HCO3Gas = 1e6*(3.4*10^(-2)*410*10^(-6)*10^(-6.36))/Hplus)%>%##Gas phase HCO3
  #mutate(HCO3Fraction = (Hplus*10^(-6.3))/(Hplus^2 + Hplus*10^(-6.36) + 10^(-6.36-10.36)))%>%##Ionizations Fraction of HCO3
  #mutate(HCO3Total = (HCO3Fraction*(MG+CA)/2) + 1e6*HCO3Gas,
  mutate(pHBin = cut(LABPH, 6))%>%
  mutate(IonBalance = Cations-Anions,
         IonBalanceHCO3 = Cations - Anions -HCO3Gas,
         TDpH = -log10(10^(-LABPH)+(CA+MG)/1e6),
         TDHplus = CA+MG+Hplus*1e6,
         LWC = ifelse(is.na(LWC), LWCNew, LWC),
         MajorIonBalance = 1e6*Hplus+NH4+CA+MG-SO4-NO3-CL,
         TDHplusAll = CA+MG+Hplus*1e6+K+Sodium-CL,
         TDpHAll = -log10(TDHplusAll/1e6))



## Figure S2 TOC vs WSOC ####

TOCWSOC<-read.csv("TOCvsWSOC.csv", header = TRUE)

TOCWSOCPlot<-ggplot(TOCWSOC)+
  geom_point(aes(x = TOC, y = WSOC), color = 'forest green', size = 3)+
  geom_smooth(aes(x = TOC, y= WSOC), formula = y~x, method = lm)+
  PaperQuality()+
  stat_regline_equation(aes(x = TOC, y = WSOC,
                            label = paste(..eq.label.., ..rr.label.., 
                                          sep = "~~~~")))+
  labs(x = TeX("\\textbf{TOC ($\\mu$molC L^{-1})}"),
       y = TeX("\\textbf{WSOC ($\\mu$molC L^{-1})}"))
ggsave(TOCWSOCPlot, filename = "FigureS4.png", width =8, height = 8)

## Figure S3 Percent Valid by Month ####

PercentByMonth<-AllCloudData%>%
  select(-TOC)%>% ##Include whole dataset, as TOC start at 2009
  group_by(Year,Class, Month)%>%
  summarise(n = n())%>%
  group_by(Year, Month)%>%
  na.omit()%>%
  mutate(Percent = 100*n/(sum(n)), Percent = ifelse(Month ==9 & Year ==2021,
                                                    NA, Percent))%>%
  mutate(MonthName = case_when(Month == 6 ~ "June",
                               Month == 7 ~ "July",
                               Month == 8 ~ "August",
                               Month == 9 ~ "September"),
         MonthName = factor(MonthName, levels = c("June",
                                                  "July",
                                                  "August",
                                                  'September')))


Monthlist<-c(6,7,8,9)
PercentTrendTable<-lapply(Monthlist, function(x){
  PercentMonthlyMedians<-filter(PercentByMonth, Month == x & Class == "Invalid")%>%
    arrange(Year, Month)
  print(theilsen(data = PercentMonthlyMedians, Percent~Year))
  return(PercentMonthlyMedians)
})

dat_text<-data.frame(
  MonthName = c('June', 'July', 'August', 'September'),
  label = c('y = 2.12x - 4238.0',
            'y = 1.35x - 2684.3', ### Labels are from the regression analysis, 
            'y = 1.00x - 1989.0', ### but I couldn't find a neat way to write them
            'y = 1.244x - 2479.7') ### in each panel without hard coding
)


dat_text<-dat_text%>%
  mutate(MonthName = factor(MonthName, levels = c('June', "July", 
                                                  "August", 'September')))



PercentByMonthPlot<-ggplot(PercentByMonth%>%
                             filter(Month > 5 & Month < 10 & Class == "Invalid")%>%
                             mutate(Month = as.character(Month)))+
  geom_line(aes(x = Year,y = Percent, color = MonthName), size =3)+
  geom_point(aes(x = Year,y = Percent, fill = MonthName), position = 'dodge', color = 'black',
             shape =21, size =3)+
  geom_smooth(aes(x = Year, y = Percent), method = sen, se = FALSE)+
  PaperQuality()+facet_wrap(~MonthName)+
  labs(y = 'Percent Invalid (%)')+
  geom_text(aes(x = 2005, y = 60, label = label),
            data = dat_text)+theme(strip.text = element_text(size = 14, face = 'bold'),
                                   legend.position = "none")


ggsave("PercentValidvsMonth.png",PercentByMonthPlot,
       width = 12, height = 6)


### Figure S4 Unfiltered vs Filtered Ca and Mg ####

CaMgData<-read.csv("Ca_Mg_Rerun_Data.csv")

MgFilterPlot<-ggplot(CaMgData)+
  geom_point(aes(x = IC_Mg/2, y = Unfiltered_IC_Mg/2), size = 3, color = 'black',
             shape = 21, fill = 'blue')+
  geom_smooth(aes(x = IC_Mg/2, y = Unfiltered_IC_Mg/2), 
              color = 'blue', size =1.5, method = lm,
              se = FALSE, formula= y~x+0)+
  PaperQuality()+ggpmisc::stat_poly_eq(aes(x = IC_Mg/2, y = Unfiltered_IC_Mg/2,
                                           label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
                                       formula = y~x+0, size =8)+
  labs(x = TeX('\\textbf{Filtered Mg ($\\mu$mol L$^{-1}$)}'), 
       y= TeX(('\\textbf{Unfiltered Mg ($\\mu$mol L$^{-1}$)}')))+
  ggtitle("Mg")

CaFilterPlot<-ggplot(CaMgData)+
  geom_point(aes(x = IC_Ca/2, y = Unfiltered_IC_Ca/2), size = 3, color = 'black',
             shape = 21, fill = 'red')+
  geom_smooth(aes(x = IC_Ca/2, y = Unfiltered_IC_Ca/2), color = 'red', size =1.5, method = lm,
              se = FALSE, formula = y~x+0)+
  PaperQuality()+ggpmisc::stat_poly_eq(aes(x = IC_Ca/2, y = Unfiltered_IC_Ca/2,
                                           label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
                                       formula = y~x+0, size =8)+
  labs(x = TeX('\\textbf{Filtered Ca ($\\mu$mol L$^{-1}$)}'), 
       y= TeX(('\\textbf{Unfiltered Ca ($\\mu$mol L$^{-1}$)}')))+
  ggtitle("Ca")

CaMgICFilteredPlot<-grid.arrange(CaFilterPlot, MgFilterPlot, ncol = 2)

ggsave(plot =CaMgICFilteredPlot, filename = 'FigureS1.png',width =14, height = 10)

### Figure S5  Total Ca and Mg vs Cationic Ca and Mg ####



CaPlot<-ggplot(CaMgData)+
  geom_point(aes(x = Total_Ca/2, y = IC_Ca/2, fill = pH), size = 3, color = 'black',
             shape = 21)+
  geom_smooth(aes(x = Total_Ca/2, y= IC_Ca/2), color = 'blue', size =1.5, method = lm,
              se = FALSE, formula = y~x+0)+
 # geom_abline(slope =1, intercept = 0, linetype = 'dashed')+
  PaperQuality()+ggpmisc::stat_poly_eq(aes(x = Total_Ca/2, y = IC_Ca/2 ,
                                  label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
                              formula = y~x+0, size =8)+
  labs(x = TeX('\\textbf{Total Ca Measurement ($\\mu$mol L$^{-1}$)}'), 
       y= TeX(('\\textbf{Ion Chromotagraphy Ca Measurement ($\\mu$mol L$^{-1}$)}')),
       fill = 'pH')+
  theme(legend.position = 'right', legend.title = element_text(size = 14, face = 'bold'))+
  scale_fill_viridis_c()+
  ggtitle(TeX('\\textbf{Ca}'))

  
MgPlot<-ggplot(CaMgData)+
  geom_point(aes(x = Total_Mg/2, y = IC_Mg/2, fill = pH), size = 3, color = 'black',
             shape = 21)+
  geom_smooth(aes(x = Total_Mg/2, y= IC_Mg/2), color = 'blue', size =1.5,
              method =lm, formula = y~x+0, se = FALSE)+
  PaperQuality()+ggpmisc::stat_poly_eq(aes(x = Total_Mg/2, y = IC_Mg/2 ,
                                  label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
                              formula = y~x+0, size =8)+

  labs(x = TeX('\\textbf{Total Mg Measurement ($\\mu$mol L$^{-1}$)}'), 
       y= TeX(('\\textbf{Ion Chromotagraphy Mg Measurement ($\\mu$mol L$^{-1}$)}')))+
  theme(legend.position = 'right', legend.title = element_text(size = 14, face = 'bold'))+
  scale_fill_viridis_c()+
  ggtitle(TeX('\\textbf{Mg}'))

CaMgPlot<-grid.arrange(CaPlot, MgPlot, ncol = 2)

ggsave(plot =CaMgPlot, filename = 'FigureS2.png',width =14, height = 10)



### Modeling CaCO3 and MgCO3 Solubility vs pH ####
CaCO3Fun<-function(x){
  2e6*(10^(-x*2)*10^(-8.33))/(10^(-6.36)*10^(-10.3)*410e-6*3.4*10^(-2))
  
}


MgCO3Fun<-function(x){
  2e6*(10^(-x*2)*6.82e-6)/(10^(-6.36)*10^(-10.3)*410e-6*3.4*10^(-2))
  
}

AllCloudData%>%
  filter(Year > 2017)%>%
  summarise(across(where(is.numeric), max, na.rm = TRUE))

### Figure S6 Solubility of CaCO3 and MgCO3 ####
DissolvedCaCO3<-ggplot(AllCloudData)+
  xlim(c(2.5,10))+
  geom_function(fun = CaCO3Fun, aes(color = 'CaCO3'), size = 2)+
  geom_function(fun = MgCO3Fun, aes(color = 'MgCO3'), size = 2)+
  geom_hline(color = 'red', linetype = 'dashed', yintercept = 353.89, size =2)+
  geom_hline(color = 'blue', linetype = 'dashed', yintercept = 2453.59, size =2)+
  annotate(geom = 'text', x = 6, y = 1.5e4, color = 'blue', label = TeX('\\textbf{Maximum Ca$^{2+}$}'), size =10)+
  annotate(geom = 'text', x = 6, y = 1e2, color = 'red', label = TeX('\\textbf{Maximum Mg$^{2+}$}'), size =10)+
  annotate(geom = 'text', x = 8, y = 1e10, color = 'black', label = TeX('\\textbf{Maximum pH}'), size =10)+
  geom_vline(xintercept = 7.075, color = 'black', linetype = 'dashed', size =2)+
  scale_y_log10()+PaperQuality()+labs(x = 'pH', y =TeX('\\textbf{Dissolved Ca$^{2+}$ or Mg$^{2+}$ ($\\mu$eq L$^{-1}$)}'))+
  scale_color_manual(values = c('CaCO3' = 'blue', 'MgCO3' = 'red'))

ggsave(plot = DissolvedCaCO3, filename = 'FigureS3.png', height = 8, width = 12)


### Create Median Trends Separated by Valid and All Data ####

ValidMedians<-AllCloudData%>%
  select(c(Year, LABPH, SPCOND, CA,MG,K,Sodium,
           NH4,NO3,SO4,CL, TOC, WSOC,LWC, Ratio, Class))%>%
  filter(Class == "Valid")%>%
  mutate(SPCOND = log10(SPCOND))%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), 
                   list(median = median, std.error = ~std.error(., na.rm = TRUE)),
                        na.rm = TRUE,
                   .names = "{.col}.{.fn}"))


InvalidMedians<-AllCloudData%>%
  select(c(Year, LABPH, SPCOND, CA,MG,K,Sodium,
           NH4,NO3,SO4,CL, TOC, WSOC,LWC, Ratio, Class))%>%
  filter(Class == "Invalid")%>%
  mutate(SPCOND = log10(SPCOND))%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), 
                   list(median = median, std.error = ~std.error(., na.rm = TRUE)),
                   na.rm = TRUE,
                   .names = "{.col}.{.fn}"))



AllMedians<-AllCloudData%>%
  select(c(Year, LABPH, SPCOND, CA,MG,K,Sodium,
           NH4,NO3,SO4,CL, TOC, WSOC,LWC, Ratio, Class))%>%
  mutate(SPCOND = log10(SPCOND))%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), 
                   list(median = median, std.error = ~std.error(., na.rm = TRUE)),
                   na.rm = TRUE,
                   .names = "{.col}.{.fn}"))

names(ValidMedians)<-paste0(names(ValidMedians), 'Valid')
names(InvalidMedians)<-paste0(names(InvalidMedians), 'Invalid')
AllMedians<-merge(AllMedians, ValidMedians%>%
                    rename(Year = YearValid), by = 'Year')
AllMedians<-merge(AllMedians, InvalidMedians%>%
                    rename(Year = YearInvalid), by = 'Year')
names(ValidMedians)<-str_remove(names(ValidMedians), "Valid")
names(InvalidMedians)<-str_remove(names(InvalidMedians), "Invalid")

## Run the Sen Slope Functions for Valid and Total Data Sets and Create Datatables####

sen <- function(..., weights = NULL) {
  mblm::mblm(..., repeated = FALSE)
} 

TheilSenFunction<-function(df, analytes ,fun){
  df<-df%>%
  #  mutate(Hplus = Hplus*1e6)%>%
    group_by(Year)%>%
    summarise(across(where(is.numeric), fun, na.rm = TRUE))
  df<-as.data.frame(df)
  SenSlope<-list()
  SenInter<-list()
  #Coeff<-list()
  Ken<-list()
  for (i in analytes){
    #print(i)
    # print(theilsen(df[,i]~Year, data = df))
    #SenSlope[[i]]<-theilsen(df[,i]~Year, data = df)
    Slope<-theilsen(df[,i]~Year, data = df)
    SenSlope[[i]]<-Slope$coefficients[[2]]
    SenInter[[i]]<-Slope$coefficients[[1]]
    pvalue<-MannKendall(df[,i])
    Ken[[i]]<-pvalue$sl
    #print(Sen$coefficients[[2]])
    #print("Mann Kendall Result")
    #print(paste("tau:", Ken$tau,"p-value:", Ken$sl, sep = " "))
  }
  SenDF<-tibble(Analyte = unlist(analytes), Slope = unlist(SenSlope), Intercept =unlist(SenInter) , `P-Value` = unlist(Ken))
  SenDF<-SenDF%>%
    mutate(`P-Value` = scientific(`P-Value`, digits = 3),
           Slope = round(Slope, digits = 4),
           Intercept = round(Intercept, digits = 4))
  return(SenDF)
}


ValidSlope<-TheilSenFunction(df = subset(AllCloudData, Class == "Valid"), analytes = list("LABPH", "SPCOND", "SO4", "NO3", "NH4", "TOC", "CA", "MG", "K", 
                                                                                          "Sodium", "CL"), fun = median)
AllSlope<-TheilSenFunction(AllCloudData, analytes = list("LABPH", "SPCOND", "SO4", "NO3", "NH4", "TOC", "CA", "MG", "K", 
                                                         "Sodium", "CL", "Ratio"), fun = median)
ValidSlope$Analyte<-c("pH (units/yr)", "Conductivity (uS/cm yr)", "SO4 (ueq/L yr)", "NO3 (ueq/L yr)", "NH4 (ueq/L yr)", "TOC (umolC/L yr)", "Ca (ueq/L yr)", "Mg (ueq/L yr)", "K (ueq/L yr)", "Na (ueq/L yr)", "Cl (ueq/L yr)")
AllSlope$Analyte<-c("pH (units/yr)", "Conductivity (uS/cm yr)", "SO4 (ueq/L yr)", "NO3 (ueq/L yr)", "NH4 (ueq/L yr)", "TOC (umolC/L yr)", "Ca (ueq/L yr)", "Mg (ueq/L yr)", "K (ueq/L yr)", "Na (ueq/L yr)", "Cl (ueq/L yr)","Ratio")

ValidSlopeNoInter<-ValidSlope%>%
  dplyr::select(-Intercept)%>%
  mutate(`P-Value` = as.numeric(`P-Value`))%>%
  mutate(`P-Value` = ifelse(`P-Value` < 1e-3, "p < 0.001", `P-Value`))
ValidSlopeTable<-ggtexttable(ValidSlopeNoInter, theme = ttheme("light", base_size = 16, colnames.style = colnames_style(face = "italic", size = 16)), rows = NULL)%>%
  tab_add_title(text = paste("a) Theil-Sen Slope and",
                             "Mann Kendall P-Value 
Valid Dataset", sep = "\n"), size = 17, face = "bold",padding = unit(1, "line"),
                just = "left")

AllSlopeNoInter<-AllSlope%>%
  dplyr::select(-Intercept)%>%
  filter(Analyte != "Ratio")%>%
  mutate(`P-Value` = as.numeric(`P-Value`))%>%
  mutate(`P-Value` = ifelse(`P-Value` < 1e-3, "p < 0.001", `P-Value`))
AllSlopeTable<-ggtexttable(AllSlopeNoInter, theme = ttheme("light", base_size = 16, colnames.style = colnames_style(face = "italic", size = 16)), rows = NULL)%>%
  tab_add_title(text = paste("b) Theil-Sen Slope and",
                             "Mann Kendall P-Value
Complete Dataset", sep = "\n"), size = 17, face = "bold",padding = unit(1, "line"),
                just = "left")

## Figure 3 Median Trend Plots ####


ValidpHCondPlot<-ggplot(ValidMedians)+
  geom_line(aes(x = Year, y = LABPH.median),color = 'purple', show.legend = FALSE, size =6)+
  geom_errorbar(aes(x = Year, ymin = LABPH.median - LABPH.std.error, ymax = LABPH.median + LABPH.std.error),
                show.legend = FALSE, color = 'purple')+
  geom_point(aes(x = Year, y = LABPH.median ,fill = 'pH') ,color = 'black', shape = 21, size =6)+
  geom_smooth(aes(x= Year, y= LABPH.median, color = 'pH'), method = sen,
              se = FALSE, show.legend = FALSE, size =2)+
  geom_line(aes(x = Year, y = 3.5+SPCOND.median), color = 'black', show.legend = FALSE, size =6)+
  geom_errorbar(aes(x = Year, ymin = 3.5+SPCOND.median - SPCOND.std.error, ymax = 3.5+SPCOND.median + SPCOND.std.error)
                , show.legend = FALSE, color = 'black')+
  geom_point(aes(x = Year, y = 3.5+SPCOND.median, fill = 'Conductivity'), color = 'black',
             size = 6, shape =21)+
  geom_smooth(aes(x = Year, y = 3.5+SPCOND.median, color = 'Conductivity'), method = sen,
              se = FALSE, show.legend = FALSE, size = 2)+
  scale_x_continuous(breaks = seq(1994, 2022, 4))+
  scale_y_continuous(sec.axis = sec_axis(~(.-3.5),
                                         labels = function(x) comma(round(10^x), accuracy=1),
                                         breaks=seq(0,10,1),
                                         name = TeX("\\textbf{Conductivity ($\\mu$S cm$^{-1}$})")),
                     name = 'pH')+ 
  scale_fill_manual(values = c("pH" = 'purple', 'Conductivity' = "black","SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  scale_color_manual(values = c("pH" = 'purple', 'Conductivity' = "black","SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"), guide = FALSE)+
  annotate(geom = 'text', x = 1998, y = 4.25, size =12, label = 'pH', color = 'purple', 
           fontface = 'bold')+
  annotate(geom = 'text', x = 1998, y = 5.65, size =12, label = 'Conductivity', color = 'black',
           fontface = 'bold')+
  
  annotate(geom = 'text',x = 2010, y = 6, size = 10,
           label = paste(" Conductivity: y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="Conductivity (uS/cm yr)"][[1]],digits = 3),"x +", 
                         signif(ValidSlope$Intercept[ValidSlope$Analyte=="Conductivity (uS/cm yr)"][[1]], digits = 3)), color = "black")+
  annotate(geom = 'text',x = 2010, y = 5.75, size = 10,
           label = paste("pH : y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="pH (units/yr)"][[1]],digits = 3),"x +", 
                         signif(ValidSlope$Intercept[ValidSlope$Analyte=="pH (units/yr)"][[1]], digits = 3)), color = "purple")+
  
  PaperQuality()+theme(axis.text.y.left = element_text(color = 'purple'), axis.title.y.left = element_text(color= 'purple'))+
  labs(title = TeX('\\textbf{b)}'),
       subtitle = TeX("\\textbf{pH and Conductivity Valid Data Only}"))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 34), axis.text.y = element_text(size = 34),
        axis.text.y.right = element_text(size = 34), axis.title.y.right = element_text(size = 34),
        axis.title.x = element_text(size = 34), axis.title.y = element_text(size = 34),
        title = element_text(size = 38))



ValidConcPlots<-ggplot(ValidMedians)+
  geom_errorbar(aes(x = Year, ymin = SO4.median -SO4.std.error, ymax = SO4.median +SO4.std.error), color = 'red')+
  geom_line(aes(x = Year, y = SO4.median), color = "red",size = 6)+
  geom_point(aes(x = Year, y = SO4.median, fill = "SO4"), size = 6, shape = 21)+
  geom_smooth(aes(x = Year, y = SO4.median, color = "SO4"), method = sen, size =2, se = FALSE)+
  geom_errorbar(aes(x = Year, ymin = NO3.median - NO3.std.error, ymax = NO3.median + NO3.std.error), color = 'blue')+
  geom_line(aes(x = Year, y = NO3.median),color = "blue", size = 6)+
  geom_point(aes(x = Year, y = NO3.median, fill = "NO3"), size = 6, shape = 21)+
  geom_smooth(aes(x = Year, y = NO3.median, color = "NO3"), method = sen, size = 2, se = FALSE)+
  geom_errorbar(aes(x = Year, ymin = NH4.median - NH4.std.error, ymax = NH4.median+NH4.std.error), color = 'orange')+
  geom_line(aes(x = Year, y = NH4.median), color = "orange", size = 6)+
  geom_point(aes(x = Year, y = NH4.median, fill = "NH4"), size = 6, shape = 21)+
  geom_smooth(aes(x = Year, y = NH4.median, color = "NH4"), method = sen, size = 2, se = FALSE)+
  geom_line(aes(x = Year, y = WSOC.median/4),color = "green",size = 6)+
  geom_point(aes(x = Year, y = WSOC.median/4, fill = "WSOC"), size = 6, shape =21)+
  geom_errorbar(aes(x = Year, ymin = TOC.median/4 - TOC.std.error/4, ymax = TOC.median/4 + TOC.std.error/4), color = 'forest green')+
  geom_line(aes(x = Year, y = TOC.median/4),color = "forest green",size = 6)+
  geom_point(aes(x = Year, y = TOC.median/4, fill = "TOC"), size = 6, shape =21)+
  geom_smooth(aes(x = Year, y = TOC.median/4, color = "TOC"), method = sen, size = 2, se = FALSE)+
  scale_x_continuous(breaks = seq(1994, 2022, by = 4))+scale_y_continuous(guide = guide_axis(check.overlap = TRUE),
                                                                          sec.axis = sec_axis(trans = ~.x*4, name = TeX("\\textbf{TOC Concentration ($\\mu$molC  L$^{-1}$)}")),
                                                                          name =TeX("\\textbf{Ion Concentration ($\\mu$eq L$^{-1}$)}"))+
  scale_fill_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green", 'WSOC' = "green"))+
  scale_color_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green", 'WSOC' = 'green'), guide = FALSE)+
  PaperQuality()+theme(axis.text.y.right = element_text(color = "forest green"),
                       axis.ticks = element_line(color = "forest green"), axis.title.y.right = element_text(color = "forest green", size = 15))+
  annotate(geom = 'text',x = 2012, y = 200, size = 10,
           label = paste("SO4: y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="SO4 (ueq/L yr)"][[1]],digits = 3),"x +", 
                         signif(ValidSlope$Intercept[ValidSlope$Analyte=="SO4 (ueq/L yr)"][[1]], digits = 3)), color = "red")+
  annotate(geom = 'text',x = 2012, y = 175, size = 10,
           label = paste("NH4: y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="NH4 (ueq/L yr)"][[1]],digits = 3),"x +", 
                         signif(ValidSlope$Intercept[ValidSlope$Analyte=="NH4 (ueq/L yr)"][[1]], digits = 3)), color = "orange")+
  annotate(geom = 'text',x = 2012, y = 150, size = 10,
           label = paste("NO3: y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="NO3 (ueq/L yr)"][[1]],digits = 3),"x +", 
                         signif(ValidSlope$Intercept[ValidSlope$Analyte=="NO3 (ueq/L yr)"][[1]], digits = 3)), color = "blue")+
  annotate(geom = 'text',x = 2012, y = 125, size = 10,
           label = paste("TOC: y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="TOC (umolC/L yr)"][[1]],digits = 3),"x +", 
                         signif(ValidSlope$Intercept[ValidSlope$Analyte=="TOC (umolC/L yr)"][[1]], digits = 3)), color = "forest green")+
  labs(title = TeX('\\textbf{a)}'),
       subtitle = TeX("\\textbf{Major Analytes Valid Data Only}"))+theme(axis.title.y.right = element_text(size =22, face = 'bold'))+
  annotate(geom = 'text', x = 1997, y = 105, size =12, label = TeX('\\textbf{NH$_4^+$}'), color = 'orange', 
           fontface = 'bold')+
  annotate(geom = 'text', x = 1998, y = 175, size =12, label = TeX('\\textbf{SO$_4^{2-}$}'), color = 'red',
           fontface = 'bold')+
  annotate(geom = 'text', x = 2018, y = 90, size =12, label = TeX('\\textbf{TOC}'), color = 'forest green', 
           fontface = 'bold')+
  annotate(geom = 'text', x = 1998, y = 25, size =12, label = TeX('\\textbf{NO$_3^-$}'), color = 'blue',
           fontface = 'bold')+
  annotate(geom = 'text', x = 2020, y = 50, size =12, label = TeX('\\textbf{WSOC}'), color = 'green',
           fontface = 'bold')+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 34), axis.text.y = element_text(size = 34),
        axis.text.y.right = element_text(size = 34), axis.title.y.right = element_text(size = 34),
        axis.title.x = element_text(size = 34), axis.title.y = element_text(size = 34),
        title = element_text(size = 38))


 
ValidPlots<-ggarrange(ValidConcPlots,ValidpHCondPlot, nrow = 1, ncol = 2)


#ggsave(ValidPlots, filename = "Figure3.png",width = 22, height = 11)
#ggsave(ValidSlopeTable, filename = 'Figure3Table.png', height = 8, width = 5)


AllpHCondPlot<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = LABPH.median),color = 'purple', show.legend = FALSE, size =6)+
  geom_errorbar(aes(x = Year, ymin = LABPH.median - LABPH.std.error, ymax = LABPH.median + LABPH.std.error),
                show.legend = FALSE, color = 'purple')+
  geom_point(aes(x = Year, y = LABPH.median ,fill = 'pH') ,color = 'black', shape = 21, size =6)+
  geom_smooth(aes(x= Year, y= LABPH.median, color = 'pH'), method = sen,
              se = FALSE, show.legend = FALSE, size =4)+
  geom_line(aes(x = Year, y = 3.5+SPCOND.median), color = 'black', show.legend = FALSE, size =6)+
  geom_errorbar(aes(x = Year, ymin = 3.5+SPCOND.median - SPCOND.std.error, ymax = 3.5+SPCOND.median + SPCOND.std.error)
                , show.legend = FALSE, color = 'black')+
  geom_point(aes(x = Year, y = 3.5+SPCOND.median, fill = 'Conductivity'), color = 'black',
             size = 6, shape =21)+
  geom_smooth(aes(x = Year, y = 3.5+SPCOND.median, color = 'Conductivity'), method = sen,
              se = FALSE, show.legend = FALSE, size = 4)+
  scale_y_continuous(sec.axis = sec_axis(~(.-3.5),
                                         labels = function(x) comma(round(10^x), accuracy=1),
                                         breaks=seq(0,10,1),
                                         name = TeX("\\textbf{Conductivity ($\\mu$S cm$^{-1}$})")),
                     name = 'pH')+
  scale_x_continuous(breaks = seq(1994, 2022, 4))+
  scale_fill_manual(values = c("pH" = 'purple', 'Conductivity' = "black","SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  scale_color_manual(values = c("pH" = 'purple', 'Conductivity' = "black","SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"), guide = FALSE)+
  annotate(geom = 'text',x = 2010, y = 6, size = 10,
           label = paste("Conductivity: y=",signif(AllSlope$Slope[AllSlope$Analyte=="Conductivity (uS/cm yr)"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="Conductivity (uS/cm yr)"][[1]], digits = 3)), color = "black")+
  annotate(geom = 'text',x = 2010, y = 5.75, size = 10,
           label = paste("pH: y=",signif(AllSlope$Slope[AllSlope$Analyte=="pH (units/yr)"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="pH (units/yr)"][[1]], digits = 3)), color = "purple")+
  
  PaperQuality()+theme(axis.text.y.left = element_text(color = 'purple'), axis.title.y.left = element_text(color= 'purple'))+
  labs(title = TeX('\\textbf{d)}'),
       subtitle = TeX("\\textbf{pH and Conductivity Complete Dataset}"))+
  theme(legend.position = 'none')+
  annotate(geom = 'text', x = 1998, y = 4.25, size =12, label = 'pH', color = 'purple', 
           fontface = 'bold')+
  annotate(geom = 'text', x = 1998, y = 5.65, size =12, label = 'Conductivity', color = 'black',
           fontface = 'bold')+
  theme(axis.text.x = element_text(size = 34), axis.text.y = element_text(size = 34),
        axis.text.y.right = element_text(size = 34), axis.title.y.right = element_text(size = 34),
        axis.title.x = element_text(size = 34), axis.title.y = element_text(size = 34),
        title = element_text(size = 38))




AllConcPlots<-ggplot(AllMedians)+
  geom_errorbar(aes(x = Year, ymin = SO4.median -SO4.std.error, ymax = SO4.median +SO4.std.error), color = 'red')+
  geom_line(aes(x = Year, y = SO4.median), color = "red",size = 6)+
  geom_point(aes(x = Year, y = SO4.median, fill = "SO4"), size = 6, shape = 21)+
  geom_smooth(aes(x = Year, y = SO4.median, color = "SO4"), method = sen, size =2, se = FALSE)+
  geom_errorbar(aes(x = Year, ymin = NO3.median - NO3.std.error, ymax = NO3.median + NO3.std.error), color = 'blue')+
  geom_line(aes(x = Year, y = NO3.median),color = "blue", size = 6)+
  geom_point(aes(x = Year, y = NO3.median, fill = "NO3"), size = 6, shape = 21)+
  geom_smooth(aes(x = Year, y = NO3.median, color = "NO3"), method = sen, size = 2, se = FALSE)+
  geom_errorbar(aes(x = Year, ymin = NH4.median - NH4.std.error, ymax = NH4.median+NH4.std.error), color = 'orange')+
  geom_line(aes(x = Year, y = NH4.median), color = "orange", size = 6)+
  geom_point(aes(x = Year, y = NH4.median, fill = "NH4"), size = 6, shape = 21)+
  geom_line(aes(x = Year, y = WSOC.median/4),color = "green",size = 6)+
  geom_point(aes(x = Year, y = WSOC.median/4, fill = "WSOC"), size = 6, shape =21)+
  geom_smooth(aes(x = Year, y = NH4.median, color = "NH4"), method = sen, size = 2, se = FALSE)+
  geom_errorbar(aes(x = Year, ymin = TOC.median/4 - TOC.std.error/4, ymax = TOC.median/4 + TOC.std.error/4), color = 'forest green')+
  geom_line(aes(x = Year, y = TOC.median/4),color = "forest green",size = 6)+
  geom_point(aes(x = Year, y = TOC.median/4, fill = "TOC"), size = 6, shape =21)+
  geom_smooth(aes(x = Year, y = TOC.median/4, color = "TOC"), method = sen, size = 2, se = FALSE)+
  scale_x_continuous(breaks = seq(1994, 2022, by = 4))+scale_y_continuous(guide = guide_axis(check.overlap = TRUE),
                                                                          sec.axis = sec_axis(trans = ~.x*4, name = TeX("\\textbf{TOC Concentration ($\\mu$molC  L$^{-1}$)}")),
                                                                          name =TeX("\\textbf{Ion Concentration ($\\mu$eq  L$^{-1}$)}"))+
  scale_fill_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green", 'WSOC' = 'green'))+
  scale_color_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green", 'WSOC' = 'green'), guide = FALSE)+
  PaperQuality()+theme(axis.text.y.right = element_text(color = "forest green"),
                       axis.ticks = element_line(color = "forest green"), axis.title.y.right = element_text(color = "forest green", size = 15))+
  annotate(geom = 'text',x = 2012, y = 200, size = 10,
           label = paste("SO4: y=",signif(AllSlope$Slope[AllSlope$Analyte=="SO4 (ueq/L yr)"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="SO4 (ueq/L yr)"][[1]], digits = 3)), color = "red")+
  annotate(geom = 'text',x = 2012, y = 175, size = 10,
           label = paste("NH4: y=",format(signif(AllSlope$Slope[AllSlope$Analyte=="NH4 (ueq/L yr)"][[1]],digits = 3), nsmall = 2),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="NH4 (ueq/L yr)"][[1]], digits = 3)), color = "orange")+
  annotate(geom = 'text',x = 2012, y = 150, size = 10,
           label = paste("NO3: y=",signif(AllSlope$Slope[AllSlope$Analyte=="NO3 (ueq/L yr)"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="NO3 (ueq/L yr)"][[1]], digits = 3)), color = "blue")+
  annotate(geom = 'text',x = 2012, y = 125, size = 10,
           label = paste("TOC: y=",signif(AllSlope$Slope[AllSlope$Analyte=="TOC (umolC/L yr)"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="TOC (umolC/L yr)"][[1]], digits = 3)), color = "forest green")+
  labs(title = TeX('\\textbf{c)}'),
       subtitle = TeX("\\textbf{Major Analytes Complete Dataset}"))+theme(axis.title.y.right = element_text(size =22, face = 'bold'))+
  theme(legend.position = 'none')+
  annotate(geom = 'text', x = 1994, y = 105, size =12, label = TeX('\\textbf{NH$_4^+$}'), color = 'orange', 
           fontface = 'bold')+
  annotate(geom = 'text', x = 1998, y = 175, size =12, label = TeX('\\textbf{SO$_4^{2-}$}'), color = 'red',
           fontface = 'bold')+
  annotate(geom = 'text', x = 2019, y = 135, size =12, label = TeX('\\textbf{TOC}'), color = 'forest green', 
           fontface = 'bold')+
  annotate(geom = 'text', x = 1998, y = 25, size =12, label = TeX('\\textbf{NO$_3^-$}'), color = 'blue',
           fontface = 'bold')+
  annotate(geom = 'text', x = 2020, y = 70, size =12, label = TeX('\\textbf{WSOC}'), color = 'green',
           fontface = 'bold')+
  theme(axis.text.x = element_text(size = 34), axis.text.y = element_text(size = 34),
        axis.text.y.right = element_text(size = 34), axis.title.y.right = element_text(size = 34),
        axis.title.x = element_text(size = 34), axis.title.y = element_text(size = 34),
        title = element_text(size = 38))


AllPlots<-grid.arrange(ValidConcPlots,ValidpHCondPlot, AllConcPlots, AllpHCondPlot, nrow = 2, ncol = 2)
#AllTables<-grid.arrange(ValidSlopeTable, AllSlopeTable, ncol =2)

ggsave(AllPlots, filename = "Figure3.png",width = 34, height = 20)
ggsave(AllTables, filename = 'Figure3Table.png', height = 8, width = 12)




## Figure 4 Percent Valid Plots ####

# Create Precent Valid and Invalid
Percent<-AllCloudData%>%
  select(-TOC)%>% ##Include whole dataset, as TOC start at 2009
  group_by(Year,Class)%>%
  summarise(n = n())%>%
  group_by(Year)%>%
  na.omit()%>%
  summarise(Percent = 100*n/(sum(n)))
Percent<-Percent[seq(1, nrow(Percent), by =2),] ##Filters out data to only look at percent invalid


PercentPlot<-ggplot(Percent)+
  geom_line(aes(Year, y = Percent), color = 'red', size = 4)+
  geom_point(aes(Year, y = Percent), fill = "red", shape = 21, size = 4)+
  PaperQuality()+scale_y_continuous(TeX('\\textbf{Percent Invalid (%)}'))+scale_x_continuous(breaks = seq(1994, 2022, 4))

ggsave(PercentPlot, filename = "Figure4.png", width =6, height = 5)



## Figure S7 Minute vs Hourly Liquid Water Content ####

Minute2016Data<-vroom("2016MinuteData.csv", col_names = TRUE, skip = 1, guess_max = 1e4,
                      delim = ",") ### Read in the data

Minute2016Data[Minute2016Data==-9999.000]<-NA
Minute2016Data<-na.omit(Minute2016Data)

ChrisCloudData2016<-Minute2016Data%>%
  mutate(date = as.POSIXct(paste0(Date, Time), format = "%d/%m/%Y %H:%M:%S"))%>%
  # mutate(LWCChris = `LWC : Value`)%>%
  mutate(LWCChris = ifelse(`CONFIRM : Value` >= 0.4 & `LWC : Value` >= 0.05, `LWC : Value`, NA))%>% 
  as_tbl_time(., index = date)%>%
  collapse_by("12 hourly", start_date =as.POSIXct("2016-05-27 18:00:00"))%>%
  group_by(date)%>%
  summarise(across(where(is.numeric), mean, na.rm =TRUE))%>%
  mutate(date = date+minutes(1), date = as.character(date))%>%
  as.data.frame(.)

ALSC2016DataFiltered<-vroom("WFC12HourMet.csv", col_names = TRUE)%>%
  filter(Year ==2016)%>%
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y %H:%M"))%>%
  as_tbl_time(. , index =date)%>%
  mutate(LWCALSC = LWC)%>%
  mutate(LWCALSC = ifelse(RainCount/60 < .25 & CloudCount/60 >= .25 & LWC >=0.05, LWC, NA))%>%
  collapse_by("12 hourly", start_date = as.POSIXct("2016-05-27 18:00:00"))%>%
  group_by(date)%>%
  summarise(across(where(is.numeric), mean, na.rm=TRUE))%>%
  mutate(date = date+hours(1), date = as.character(date))%>%
  as.data.frame(.)

ALSC2016DataNotFiltered<-vroom("WFC12HourMet.csv", col_names = TRUE)%>%
  filter(Year ==2016)%>%
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y %H:%M"))%>%
  as_tbl_time(. , index =date)%>%
  mutate(LWCALSC = LWC)%>%
  mutate(LWCALSC = ifelse(LWC >= 0.05, LWC, NA))%>%
  collapse_by("12 hourly", start_date = as.POSIXct("2016-05-27 18:00:00"))%>%
  group_by(date)%>%
  summarise(across(where(is.numeric), mean, na.rm=TRUE))%>%
  mutate(date = date+hours(1), date = as.character(date))%>%
  as.data.frame(.)



AllLWCData2016Filtered<-merge(ChrisCloudData2016, ALSC2016DataFiltered)  
AllLWCData2016NotFiltered<-merge(ChrisCloudData2016, ALSC2016DataNotFiltered)  


LWCPlot<-ggplot(AllLWCData2016Filtered)+
  geom_jitter(aes(x = LWCChris, y= LWCALSC, fill = 'Including Deployment and Rain Detection'), shape = 21, 
              color = 'black', size = 3)+
  geom_jitter(aes(x = LWCChris, y = LWCALSC, fill = 'LWC >= 0.05 Only'),
              color = 'black', size =3, data = AllLWCData2016NotFiltered, shape =21)+
  
  geom_smooth(aes(x = LWCChris, y= LWCALSC, color = 'Including Deployment and Rain Detection'), size = 2, method =lm, formula = y~x+0,
              show.legend = FALSE)+
  geom_smooth(aes(x = LWCChris, y = LWCALSC, color = 'LWC >= 0.05 Only'), size =2, data = AllLWCData2016NotFiltered,
              method =lm, formula = y~x+0,
              show.legend = FALSE)+
  stat_regline_equation(aes(x = LWCChris, y= LWCALSC, color = 'Including Deployment and Rain Detection',
                            label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                        formula = y~x+0,
                        size = 4, show.legend = FALSE)+
  
  stat_regline_equation(aes(x = LWCChris, y= LWCALSC, color = 'LWC >= 0.05 Only',
                            label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                        data = AllLWCData2016NotFiltered,
                        label.y.npc = 0.9, formula = y~x+0,
                        size = 4,
                        show.legend = FALSE)+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))+
  PaperQuality()+geom_abline(slope = 1, intercept = 0, linetype = 'dashed',
                             size = 2)+
  labs(x = TeX("\\textbf{LWC from 1-Minute Resolution Data (g m^{-3})}"), y = TeX("\\textbf{LWC from ALSC 1-Hour Resolution Data (g m^{-3})}"))

ggsave(LWCPlot, filename = 'FigureS5.png', height =8, width = 8)

## Figure S8 Comparing Valid and Invalid Medians####

ValidInvalidpH<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = LABPH.medianInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = LABPH.medianInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = LABPH.medianValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = LABPH.medianValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  #scale_fill_manual(values = c("pH" = 'purple', 'Conductivity' = "black","SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  PaperQuality()+theme(axis.text.y.right = element_text(color = "forest green"), axis.ticks = element_line(color = "forest green"), 
                       axis.title.y.left = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))+
  ggtitle("pH")+
  ylab("pH")+scale_x_continuous(breaks = seq(1994, 2022, by = 4))

ValidInvalidCond<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = 10^(SPCOND.medianInvalid)), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = 10^(SPCOND.medianInvalid), fill = "Invalid"), color = 'black', size = 4, shape = 21)+ ### Need to unlog SCPOND to make it fit 
  geom_line(aes(x = Year, y = 10^(SPCOND.medianValid)),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = 10^(SPCOND.medianValid), fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{Conductivity ($\\mu S $ cm^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2022, by = 4))+
  ggtitle("Conductivity")

ValidInvalidSO4<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = SO4.medianInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = SO4.medianInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = SO4.medianValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = SO4.medianValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{SO_4^{-2} ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2022, by = 4))+
  ggtitle(TeX("\\textbf{SO_4^{ 2-}}"))

ValidInvalidNH4<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = NH4.medianInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = NH4.medianInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = NH4.medianValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = NH4.medianValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{NH_4^{ +} ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2022, by = 4))+
  ggtitle(TeX("\\textbf{NH_4^{ +}}"))

ValidInvalidNO3<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = NO3.medianInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = NO3.medianInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = NO3.medianValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = NO3.medianValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{NO_3^{ -} ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2022, by = 4))+
  ggtitle(TeX("\\textbf{NO_3^{ -}}"))

ValidInvalidTOC<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = TOC.medianInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = TOC.medianInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = TOC.medianValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = TOC.medianValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{TOC ($\\mu mol $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(2008, 2022, by = 2), limits = c(2009,2022))+
  ggtitle(TeX("\\textbf{TOC"))

ValidInvalidCA<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = CA.medianInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = CA.medianInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = CA.medianValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = CA.medianValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{Ca^{2+} ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2022, by = 4), limits = c(1994,2022))+
  ggtitle(TeX("\\textbf{Ca^{2+}}"))


ValidInvalidK<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = K.medianInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = K.medianInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = K.medianValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = K.medianValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{K^{+} ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2022, by = 4), limits = c(1994,2022))+
  ggtitle(TeX("\\textbf{K^{+}}"))



ValidInvalidLWC<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = LWC.medianInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = LWC.medianInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = LWC.medianValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = LWC.medianValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{LWC (g m^{-3}$)}"))+scale_x_continuous(breaks = seq(1994, 2022, by = 4), limits = c(1994,2022))+
  ggtitle(TeX("\\textbf{LWC}"))




ValidInvalidAll<-ValidInvalidpH+ValidInvalidCond+ValidInvalidSO4+ValidInvalidNH4+ValidInvalidNO3+ValidInvalidTOC+
  ValidInvalidCA+ValidInvalidK+ValidInvalidLWC&theme(legend.position = "bottom", legend.text = element_text(size = 22))
ValidInvalidAll<-ValidInvalidAll+plot_layout(guides = "collect")
ggsave(filename = "FigureS6.png", ValidInvalidAll,height = 10,width = 20)

## Figure S9  TOC Trends by Month ####

TOCMonthlyMedians<-AllCloudData%>%
  filter(!is.na(TOC) & Month > 5 & Month < 10)%>%
  group_by(Year,Month)%>%
  summarise(TOC = median(TOC,na.rm = TRUE))%>%
  mutate(MonthName = case_when(Month == 6 ~ "June",
                               Month == 7 ~ "July",
                               Month == 8 ~ "August",
                               Month == 9 ~ "September"),
         MonthName = factor(MonthName, levels = c("June",
                                                  "July",
                                                  "August",
                                                  'September')))

Monthlist<-c(6,7,8,9)
TOCTrendTable<-lapply(Monthlist, function(x){
  TOCMonthlyMedians<-filter(TOCMonthlyMedians, Month == x)
  print(TheilSenFunction(TOCMonthlyMedians, analytes = "TOC", fun = median))
  return(TOCMonthlyMedians)
})



dat_text<-data.frame(
  MonthName = c('June', 'July', 'August', 'September'),
  label = c('y = 31.2x - 62581',
            'y = 22.4x - 44585',
            'y = 7.86x - 15511',
            'y = 26.2x - 52548')
)


dat_text<-dat_text%>%
  mutate(MonthName = factor(MonthName, levels = c('June', "July", 
                                                  "August", 'September')))

TOCTrends<-ggplot(TOCMonthlyMedians)+
  geom_line(aes(x = Year, y = TOC, color = as.factor(Month)), 
            size =3)+
  geom_point(aes(x = Year, y=TOC, fill = as.factor(Month)),
             size = 3, shape = 21, color= 'black')+
  geom_smooth(aes(x = Year, y = TOC, color = as.factor(Month)),
              method = sen, se = FALSE, size =2)+
  PaperQuality()+
  facet_wrap(~MonthName)+
  theme(strip.text = element_text(face = 'bold', size =14),
        legend.position = 'none')+
  labs(y = TeX("\\textbf{TOC ($\\mu$molC L$^{-1}$)}"))+
  geom_text(aes(x = 2012, y = 625, label = label),
            data = dat_text)



ggsave('FigureS7.png',TOCTrends, width = 12, height = 7)



## Figure 5 NH4 Neutralization and Cation/Anion Ratio #### 
SurplusNH4<-ggplot(data = subset(AllCloudData))+
  #  stat_summary(aes(x = Year, group = Year, y =LABPH), fun.data = median_IQR,
  # size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = 1e6*10^(-LABPH)), fun.y = "median",
               size = 4 ,color = "red", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = 1e6*10^(-LABPH)),fun.y ="median" ,geom = "point", 
               size = 4, color = "black", fill = "red", shape = 21)+
  stat_summary(aes(x = Year, y = NH4-SO4-NO3), fun.y = "median",
               size = 4 ,color = "blue", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = NH4-SO4-NO3),fun.y ="median" ,geom = "point", 
               size = 4, color = "black", fill = "blue", shape = 21)+
  geom_hline(yintercept = 0, size =3)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Ion Concentrations ($\\mu$eq L$^{-1}$)}"))+scale_x_continuous(breaks = seq(1995,2020,5))+
  ggtitle(TeX('\\textbf{H$^+$ and NH$_4^+$ Neutralization}'))+
  annotate(x = 2010, y = -75, geom = "text", label = TeX("\\textbf{NH$_4^+$ - SO$_4^{-2}$ - NO$_3^-$}"), color = "blue", size = 6)+
  annotate(x = 2010, y = 75, geom = "text", label = TeX("\\textbf{H$^+$}"), color = 'red', size = 6)

CationAnion<-ggplot(data = AllMedians%>%
                      mutate(Ratio.std.error = ifelse(Ratio.std.error > 3 | Ratio.std.error < 0, NA, Ratio.std.error)))+
  geom_errorbar(aes(x = Year, ymin = Ratio.median - Ratio.std.error,
                    ymax = Ratio.median + Ratio.std.error), color = 'gray',
                show.legend = FALSE)+
  geom_line(aes(x = Year, y = Ratio.median), 
               size = 4 ,color = "grey")+
  geom_point(aes(x = Year, y = Ratio.median), 
               size = 4, color = "black", fill = "grey", shape = 21)+

  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Cations/Anions}"),
                     limits = c(0.5,2))+scale_x_continuous(breaks = seq(1995,2020,5))+
  ggtitle('Cation/Anion Ratio')+
  annotate(geom = 'text',x = 1998, y = 1.4, size = 6,
           label = paste("y=",signif(AllSlope$Slope[AllSlope$Analyte=="Ratio"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="Ratio"][[1]], digits = 3)))+
  geom_smooth(aes(x = Year, y = Ratio.median),color = "grey", method = sen, size = 2, se = FALSE)


NH4CationRatio<-grid.arrange(SurplusNH4,CationAnion, layout_matrix = rbind(c(1,NA,1,NA), c(2,NA,2,NA)))

ggsave(filename = "Figure5.png", NH4CationRatio, width = 8, height = 8)

## Figure S10 TOC vs LWC ####

LWCLoadingData<-read.csv("CloudLWCLoading.csv", header = TRUE)%>%
  filter(Year > 2008)%>%
  mutate(Cations = 1e6*10^(-LABPH)+NH4+Sodium+MG+CA+K,
         Anions = SO4+NO3+CL)
LWCLoadingData<-LWCLoadingData%>%
  mutate(date = as.POSIXct(date))%>%
  filter(LWCCalc < 2)%>%
  mutate(LWCCalcBin = cut(LWCCalc, breaks = seq(0.1, 1, .05),
                      include.lowest = TRUE, right = TRUE))%>%
  mutate(LWCBin = cut(LWC, breaks = seq(0.1, 1, .05),
                          include.lowest = TRUE, right = TRUE),
         TIC = Cations+Anions,
         LWCCalcBin = str_replace(LWCCalcBin, pattern = "\\(", "\\["))%>%
  filter(!is.na(LWCCalcBin))%>%
  mutate(InferredpH = -log10(10^(-LABPH)+(CA+MG)/1e6),
         InferredBin = cut(InferredpH, breaks = seq (2, 5, .5)))

  
TOCvsLWCBin<-ggplot(LWCLoadingData)+
  geom_boxplot(aes(x = LWCCalcBin, y = TOC), fill = 'forest green')+
  PaperQuality()+
  scale_y_continuous(name = TeX('\\textbf{TOC ($\\mu molC$ L^{-1})}'))+
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = .5))+
  ggtitle("TOC")+xlab(TeX("\\textbf{LWC (g m^{-3})}"))

TICvsLWCBin<-ggplot(LWCLoadingData)+  
  geom_boxplot(aes(x = LWCCalcBin, y = TIC), fill = "light blue")+
  PaperQuality()+
  scale_y_continuous(name = TeX('\\textbf{TIC ($\\mu eq$L^{-1})}'))+
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = .5))+
  ggtitle("TIC")+xlab(TeX("\\textbf{LWC (g m^{-3})}"))


BinData<-grid.arrange(TICvsLWCBin,TOCvsLWCBin)
ggsave(plot = BinData, filename = "FigureS8.png", width = 12, height = 16)




## Figure 6 Cloud Water Loadings ####


CWLSlope<-TheilSenFunction(LWCLoadingData, analytes = list("LWCCalc","SO4Mass", "NO3Mass", "NH4Mass", "TOCMass", "CAMass", "MGMass", "KMass", 
                                                           "NaMass", "CLMass"), fun = median)
CWLSlope$Analyte<-c("LWC (g/m^3 yr)","SO4 (ug/m^3 yr)","NO3 (ug/m^3 yr)", "NH4 (ug/m^3 yr)", "TOC (ug/m^3 yr)", "Ca (ug/m^3 yr)", 
                    "Mg (ug/m^3 yr)", "K (ug/m^3 yr)", "Na (ug/m^3 yr)", "Cl (ug/m^3 yr)")
CWLNoInter<-CWLSlope%>%
  dplyr::select(-Intercept)%>%
  mutate(`P-Value` = as.numeric(`P-Value`))%>%
  mutate(`P-Value` = ifelse(`P-Value` < 1e-3, "p < 0.001", `P-Value`))%>%
  mutate(Slope = case_when(Analyte == "SO4 (ug/m^3 yr)"~signif(as.numeric(Slope) *(96/2000),digits = 3),
                           Analyte == "NO3 (ug/m^3 yr)"~signif(as.numeric(Slope)*(62/1000), digits = 3),
                           Analyte == "NH4 (ug/m^3 yr)"~signif(as.numeric(Slope)*(18/1000), digits = 3),
                           Analyte == "TOC (ug/m^3 yr)"~signif(as.numeric(Slope)*(12/1000), digits = 3),
                           Analyte == "Ca (ug/m^3 yr)"~signif(as.numeric(Slope)*(40.078/2000), digits = 3),
                           Analyte == "Mg (ug/m^3 yr)"~signif(as.numeric(Slope)*(24.305/2000), digits = 3),
                           Analyte == "K (ug/m^3 yr)"~signif(as.numeric(Slope)*(39.0983 /1000), digits = 3),
                           Analyte == "Na (ug/m^3 yr)"~signif(as.numeric(Slope)*(22.98/1000), digits = 3),
                           Analyte == "Cl (ug/m^3 yr)"~signif(as.numeric(Slope)*(35.453/1000),digits = 3),
                           Analyte == "LWC (g/m^3 yr)"~signif(as.numeric(Slope), digits =3)))
  
CWLTable<-ggtexttable(CWLNoInter, theme = ttheme("light", base_size = 16, colnames.style = colnames_style(face = "italic", size = 16)), rows = NULL)%>%
  tab_add_title(text = paste("Theil-Sen Slope and",
                             "Mann Kendall P-Value
Cloud Water Loadings", sep = "\n"), size = 17, face = "bold",padding = unit(1, "line"),
                just = "left")

CWLMedians<-LWCLoadingData%>%
  group_by(Year)%>% 
  summarise(across(where(is.numeric), 
                                     list(median = median, std.error = ~std.error(., na.rm = TRUE)),
                                     na.rm = TRUE,
                                     .names = "{.col}.{.fn}"))

CWLMassPlots<-ggplot(CWLMedians)+
  geom_errorbar(aes(x = Year, ymin = (SO4Mass.median - SO4Mass.std.error)*(96/2000), ymax = (SO4Mass.median + SO4Mass.std.error)*(96/2000)),
                color = 'red')+
  geom_line(aes(x = Year, y = 96*SO4Mass.median/2000), color = "red",size = 4)+
  geom_point(aes(x = Year, y = 96*SO4Mass.median/2000, fill = "SO4"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = 96*SO4Mass.median/2000, color = "SO4"), method = sen, size =2, se = FALSE)+
  geom_errorbar(aes(x = Year, ymin = (NO3Mass.median - NO3Mass.std.error)*(62/1000), ymax = (NO3Mass.median + NO3Mass.std.error)*(62/1000)),
                color = 'blue')+
  geom_line(aes(x = Year, y = 62*NO3Mass.median/1000),color = "blue", size = 4)+
  geom_point(aes(x = Year, y = 62*NO3Mass.median/1000, fill = "NO3"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = 62*NO3Mass.median/1000, color = "NO3"), method = sen, size = 2, se = FALSE)+
  geom_errorbar(aes(x = Year, ymin = (NH4Mass.median - NH4Mass.std.error)*(18/1000), ymax = (NH4Mass.median + NH4Mass.std.error)*(18/1000)),
                color = 'orange')+
  geom_line(aes(x = Year, y = 18*NH4Mass.median/1000), color = "orange", size = 4)+
  geom_point(aes(x = Year, y = 18*NH4Mass.median/1000, fill = "NH4"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = 18*NH4Mass.median/1000, color = "NH4"), method = sen, size = 2, se = FALSE)+
  geom_errorbar(aes(x = Year, ymin = (TOCMass.median - TOCMass.std.error)*(12/1000), ymax = (TOCMass.median + TOCMass.std.error)*(12/1000)),
                  color = 'forest green')+
  geom_line(aes(x = Year, y = 12*TOCMass.median/1000),color = "forest green",size = 4)+
  geom_point(aes(x = Year, y = 12*TOCMass.median/1000, fill = "TOC"), size = 4, shape =21)+
  geom_smooth(aes(x = Year, y = TOCMass.median*(12/1000), color = 'TOC'),se = FALSE, method = sen, size =2)+### Needed to be adjusted to make there regression line go through the points properly
  scale_x_continuous(breaks = seq(2004, 2022, by = 4))+
  scale_y_continuous(guide = guide_axis(check.overlap = TRUE),
                     sec.axis = sec_axis(trans = ~.x, 
                                         name = TeX("\\textbf{TOC CWL ($\\mu$g C m^{-3})}")),
                                                                          name =TeX("\\textbf{Ion CWL ($\\mu$g m^{-3}$)}"))+
  scale_fill_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  scale_color_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"), guide = FALSE)+
  PaperQuality()+theme(axis.text.y.right = element_text(color = "forest green"),
                       axis.ticks = element_line(color = "forest green"), axis.title.y.right = element_text(color = "forest green", size = 18))+
  annotate(geom = 'text',x = 2012, y = 3.2, size = 6,
           label = paste("SO4: y=",signif(CWLSlope$Slope[CWLSlope$Analyte=="SO4 (ug/m^3 yr)"][[1]]*(96/2000),digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="SO4 (ug/m^3 yr)"][[1]]*(96/2000), digits = 3)), color = "red")+
  annotate(geom = 'text',x = 2012, y = 3.0, size = 6,
           label = paste("NH4: y=",signif(CWLSlope$Slope[CWLSlope$Analyte=="NH4 (ug/m^3 yr)"][[1]]*(18/1000),digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="NH4 (ug/m^3 yr)"][[1]]*(18/1000), digits = 3)), color = "orange")+
  annotate(geom = 'text',x = 2012, y = 2.8, size = 6,
           label = paste("NO3: y=",signif(CWLSlope$Slope[CWLSlope$Analyte=="NO3 (ug/m^3 yr)"][[1]]*(62/1000),digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="NO3 (ug/m^3 yr)"][[1]]*(62/1000), digits = 3)), color = "blue")+
  annotate(geom = 'text',x = 2012, y = 2.6, size = 6,
           label = paste("TOC: y=",signif(CWLSlope$Slope[CWLSlope$Analyte=="TOC (ug/m^3 yr)"][[1]]*(12/1000),digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="TOC (ug/m^3 yr)"][[1]]*(12/1000), digits = 3)), color = "forest green")+
  ggtitle("Major Analytes")+
  annotate(geom = 'text', x = 2007.5, y = 0.25, size =10, label = TeX('\\textbf{NH_4^+}'), color = 'orange', 
           fontface = 'bold')+
  annotate(geom = 'text', x = 2007.5, y = 0.85, size =10, label = TeX('\\textbf{  SO_4^{2-}}'), color = 'red',
           fontface = 'bold')+
  annotate(geom = 'text', x = 2019, y = 3, size =10, label = TeX('\\textbf{TOC}'), color = 'forest green', 
           fontface = 'bold')+
  annotate(geom = 'text', x = 2007.5, y = 0.55, size =10, label = TeX('\\textbf{ NO_3^-}'), color = 'blue',
           fontface = 'bold')+
  theme(axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
        axis.text.y.right = element_text(size = 30), axis.title.y.right = element_text(size = 30),
        axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
        title = element_text(size = 34), legend.position = 'none')

 



ggsave(plot = CWLMassPlots, filename = "Figure6.png",width = 10, height = 8)
ggsave(plot = CWLTable, filename = 'Figure6Table.png', width = 10, height = 6)


## Figure S11 Cation/Anion and Pred/Meas Conductivity vs pH####

CationAnionpH<-ggplot(AllCloudData)+
  geom_point(aes(x = LABPH, y = Cations/Anions, color = Year))+
  geom_hline(yintercept = 1)+
  scale_y_log10(limits  =c(0.01, 10))+
  annotation_logticks(sides = "l")+
  scale_color_gradientn(colors = rainbow(6))+labs(color = "Year", x = "pH",
                                                  y = "Cations/Anions",
                                                  title = "Cation/Anion Ratios vs pH")+
  PaperQuality()+theme(legend.position = "right", legend.title = element_text(size = 12))


PredMeaCondpH<-ggplot(AllCloudData)+
  geom_point(aes(x = LABPH, y = PredCond/SPCOND, color = Year))+
  geom_abline(slope = 0, intercept = 1)+
  geom_abline(slope = 0, intercept = 1.2, linetype = "dashed")+
  geom_abline(slope = 0, intercept = .8, linetype = 'dashed')+
  scale_y_continuous(limits = c(0.1,2))+
  scale_color_gradientn(colors = rainbow(6))+labs(title = "Predicted/Measured Conductivity Ratios vs pH",
                                                  color = "Year", x = "pH",
                                                  y = "Predicted/Measured Conductivty")+
  PaperQuality()+theme(legend.position = "right", legend.title = element_text(size = 12))

RatioPlot<-grid.arrange(CationAnionpH, PredMeaCondpH)

ggsave("FigureS9png", RatioPlot, width = 8, height = 10)



## Figure S12 HCO3 Improvements on Ion Balance ####


HCO3IonBalanceAll<-ggplot(AllCloudData)+
  geom_point(aes(x = Cations, y = Anions, color = "No HCO3"), size = 2)+
  geom_smooth(aes(x= Cations,y = Anions, color = "No HCO3"), method = lm , formula = y~x+0, size =2)+
  geom_point(aes(x = Cations, y = Anions+HCO3Gas, color = "HCO3 Included"), size = 2)+
  geom_smooth(aes(x= Cations,y = Anions+HCO3Gas, color = "HCO3 Included"), method = lm , formula = y~x+0, size = 2)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed" , size = 2)+
  PaperQuality()+scale_x_continuous(name = TeX("\\textbf{Cations ($\\mu eq $ L^{-1}$)}"))+scale_y_continuous(name = TeX("\\textbf{Anions ($\\mu eq $ L^{-1}$)}"))+
  stat_regline_equation(aes(x = Cations, y = Anions, color = "No HCO3", label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), 
                        label.x = 250, label.y = 4750,formula = y~x+0, size = 6)+
  stat_regline_equation(aes(x= Cations, y = Anions, color = "HCO3 Included",label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
                        label.x = 250, label.y = 3750,formula = y~x+0, size = 6)+ggtitle("Ion Balance All Data")

HCO3IonBalanceHighpH<-ggplot(subset(AllCloudData, LABPH > 6))+
  geom_point(aes(x = Cations, y = Anions, color = "No HCO3"), size = 2)+
  geom_smooth(aes(x= Cations,y = Anions, color = "No HCO3"), method = lm , formula = y~x+0, size = 2)+
  geom_point(aes(x = Cations, y = Anions+HCO3Gas, color = "HCO3 Included"), size = 2)+
  geom_smooth(aes(x= Cations,y = Anions+HCO3Gas, color = "HCO3 Included"), method = lm , formula = y~x+0, size = 2)+
  PaperQuality()+scale_x_continuous(name = TeX("\\textbf{Cations ($\\mu eq $ L^{-1}$)}"))+scale_y_continuous(name = TeX("\\textbf{Anions ($\\mu eq $ L^{-1}$)}"))+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed" , size = 2)+
  stat_regline_equation(aes(x = Cations, y = Anions, color = "No HCO3", label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), 
                        label.x = 250, label.y = 3000,formula = y~x+0, size = 6)+
  stat_regline_equation(aes(x= Cations, y = Anions+HCO3Gas, color = "HCO3 Included",label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
                        label.x = 250, label.y = 2250,formula = y~x+0, size = 6)+ggtitle("Ion Balance for pH > 6")



IonBalancePlots<-grid.arrange(HCO3IonBalanceAll, HCO3IonBalanceHighpH)


ggsave(plot = IonBalancePlots, filename = "FigureS10.png", width = 8, height = 8)


## Figure S13 Ion Balance vs TOC ####
#There is a seperate script for calculating the Average LWC Data

TOCBalanceAll<-ggplot(AllCloudData%>%
                        filter(Year > 2008), aes(x = Cations - Anions, y = TOC))+
  geom_point(aes(x= (Cations-Anions), y = TOC), color = 'forest green', size = 3)+
  geom_smooth(aes(x= Cations - Anions, y = TOC, color = "Outliers Included"),method = lm, size = 3)+
  geom_smooth(aes(x= Cations - Anions, y = TOC, color = "Outliers Not Included"),method = lm, size = 3, data = AllCloudData%>%filter(Cations - Anions < 400))+
  ylim(c(0,3500))+ylab(TeX("\\textbf{TOC ($\\mu molC $ L^{-1}$)}"))+
  xlab(TeX("\\textbf{Cations - Anions ($\\mu eq $ L^{-1}$)}"))+
  stat_regline_equation(size =6, label.x = 0, label.y = 3000,
                        aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = "Outliers Not Included"),
                        data = AllCloudData%>%filter(Cations - Anions < 500))+
  stat_regline_equation(size =6, label.x = 0, label.y = 2750,
                          aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = "Outliers Included"),
                          data = AllCloudData)+
  PaperQuality()

ggsave(filename = "FigureS11.png", TOCBalanceAll, width = 8, height = 8)
## Figure S14 Ion Balance by pH ####

IonBalancebypH<-ggplot(subset(AllCloudData, !is.na(pHBin) & SO4+NO3 > NH4))+
  geom_point(aes(x = Cations, y = Anions, color = Year), size = 2)+
  geom_smooth(aes(x= Cations,y = Anions), method = lm , formula = y~x+0, size =2)+
  PaperQuality()+scale_x_continuous(name = TeX("\\textbf{Cations ($\\mu$eq  L$^{-1}$)}"))+
  scale_y_continuous(name = TeX("\\textbf{Anions ($\\mu$eq L$^{-1}$)}"))+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed" , size = 2)+
  facet_wrap(~pHBin)+scale_color_gradientn(colors = rev(rainbow(6)))+
  theme(legend.position = "right", legend.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"))+
  stat_regline_equation(aes(x = Cations, y = Anions, label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), 
                        label.x = 250, label.y = 4500,formula = y~x+0, size = 5)

#ggsave(plot = IonBalancebypH, file = "FigureS5.png", height = 8, width = 12)



ggsave(plot = IonBalancebypH, file = "FigureS13.png", height = 8, width = 12)


## Figure S15 Fraction of TOC that is Ionic ####

TOCData<-filter(AllCloudData, !is.na(TOC))

TOCMedians<-TOCData%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))%>%
  mutate(OA1C = IonBalance/TOC, OA2C = 2*IonBalance/(TOC), OA3C = 3*IonBalance/(TOC),
         OA5C = 5*IonBalance/TOC)

OASlope<-TheilSenFunction(TOCMedians, analytes = c('OA1C', 'OA2C', 'OA3C', 'OA5C'), 
                          fun = median)


OASlopeNoInter<-OASlope%>%
  dplyr::select(-Intercept)%>%
  filter(Analyte != "Ratio")%>%
  mutate(`P-Value` = as.numeric(`P-Value`))%>%
  mutate(`P-Value` = ifelse(`P-Value` < 1e-3, "p < 0.001", `P-Value`))
OASlopeTable<-ggtexttable(OASlopeNoInter, theme = ttheme("light", base_size = 16, colnames.style = colnames_style(face = "italic", size = 16)), rows = NULL)%>%
  tab_add_title(text = paste("b) Theil-Sen Slope and",
                             "Mann Kendall P-Value
Fraction of OA", sep = "\n"), size = 17, face = "bold",padding = unit(1, "line"),
just = "left")


TOCOAPlot<-ggplot(TOCMedians)+
  geom_line(aes(x = Year, y = 100*OA1C, color = '1 C/z'), size =3, show.legend = FALSE)+
  geom_point(aes(x = Year, y = 100*OA1C, fill = '1 C/z'), shape = 21, size =3)+
  geom_smooth(aes(x = Year, y = 100*OA1C, color = '1 C/z'), method = sen,
              se = FALSE)+
  
  geom_line(aes(x = Year, y = 100*OA2C, color = '2 C/z'), size =3, show.legend = FALSE)+
  geom_point(aes(x = Year, y = 100*OA2C, fill = '2 C/z'), shape = 21, size =3)+
  geom_smooth(aes(x = Year, y = 100*OA2C, color = '2 C/z'), method = sen,
              se = FALSE)+
  
  geom_line(aes(x = Year, y = 100*OA3C, color = '3 C/z'), size =3, show.legend = FALSE)+
  geom_point(aes(x = Year, y = 100*OA3C, fill = '3 C/z'), shape = 21, size =3)+
  geom_smooth(aes(x = Year, y = 100*OA3C, color = '3 C/z'), method = sen,
              se = FALSE)+
  
  geom_line(aes(x = Year, y = 100*OA5C, color = '5 C/z'), size =3, show.legend = FALSE)+
  geom_point(aes(x = Year, y = 100*OA5C, fill = '5 C/z'), shape = 21, size =3)+
  geom_smooth(aes(x = Year, y = 100*OA5C, color = '5 C/z'), method = sen,
              se = FALSE)+
  annotate(geom = 'text', x = 2010, y = 7.5, size =12, label = TeX('\\textbf{1 C/z}'), color = 'red', 
           fontface = 'bold')+
  annotate(geom = 'text', x = 2010, y = 12.5, size =12, label = TeX('\\textbf{2 C/z}'), color = 'blue',
           fontface = 'bold')+
  annotate(geom = 'text', x = 2010, y = 20, size =12, label = TeX('\\textbf{3 C/z}'), color = 'orange', 
           fontface = 'bold')+
  annotate(geom = 'text', x = 2010, y = 30, size =12, label = TeX('\\textbf{5 C/z}'), color = 'green',
           fontface = 'bold')+
  PaperQuality()+labs(y = 'Percent of TOC that is Ionic (%)')+
  theme(legend.position = 'none')+
  scale_color_manual(values = c('1 C/z' = 'red', '2 C/z' = 'blue', '3 C/z' = 'orange',
                                '5 C/z' = 'green'))+
  scale_fill_manual(values = c('1 C/z' = 'red', '2 C/z' = 'blue', '3 C/z' = 'orange',
                               '5 C/z' = 'green'))+
  annotate(geom = 'text',x = 2012, y = 45, size = 10,
           label = paste("1 C/z: y=",100*signif(OASlope$Slope[OASlope$Analyte=="OA1C"][[1]],digits = 3),"x +", 
                         100*signif(OASlope$Intercept[OASlope$Analyte=="OA1C"][[1]], digits = 3)), color = "red")+
  annotate(geom = 'text',x = 2012, y = 50, size = 10,
           label = paste("2 C/z: y=",100*signif(OASlope$Slope[OASlope$Analyte=="OA2C"][[1]],digits = 3),"x +", 
                         100*signif(OASlope$Intercept[OASlope$Analyte=="OA2C"][[1]], digits = 3)), color = "blue")+
  annotate(geom = 'text',x = 2012, y = 55, size = 10,
           label = paste("3 C/z: y=",100*signif(OASlope$Slope[OASlope$Analyte=="OA3C"][[1]],digits = 3),"x +", 
                         100*signif(OASlope$Intercept[OASlope$Analyte=="OA3C"][[1]], digits = 3)), color = "orange")+
  annotate(geom = 'text',x = 2012, y = 60, size = 10,
           label = paste("5 C/z: y=",100*signif(OASlope$Slope[OASlope$Analyte=="OA5C"][[1]],digits = 3),"x +", 
                         100*signif(OASlope$Intercept[OASlope$Analyte=="OA5C"][[1]], digits = 3)), color = "green")


ggsave(plot = TOCOAPlot, file = "FigureS12.png", height = 8, width = 12)





## Figure 7a Conductivity vs pH ####

HCondPlot<-ggplot(subset(AllCloudData, !is.na(HCond/SPCOND) & HCond/SPCOND >= 0 & HCond/SPCOND <=1))+
  geom_point(aes(x= LABPH, y = SPCOND, color = HCond/SPCOND),
             size = 3)+
  scale_color_gradientn(colors = rainbow(6), limits = c(0,1))+
  scale_x_continuous(name = TeX("\\textbf{pH}"))+
  geom_abline(slope =-1, intercept = 5.5, linetype = "dotted", size = 2)+
  geom_abline(slope = -1, intercept = 6, linetype ='dashed', size =2)+
  scale_y_log10(limits =c(1,1000),
                name = TeX('\\textbf{Conductivity ($\\mu$S cm$^{-1})$}'))+
  labs(color = TeX("\\textbf{R_{Cond}}"))+PaperQuality()+
  ggtitle('a)')+
  theme(legend.position = "right",
                       legend.title = element_text(size = 20),legend.text = element_text(size = 20),
        axis.title.y = element_text(size = 22), axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22), axis.text.y = element_text(size = 22))
  

ggsave(plot = HCondPlot, filename = "Figure7a.png", 
       height = 6, width = 10)

## Figure 7b Percent in the New Regime #### 

PercentRegime<-AllCloudData%>%
  select(-TOC)%>% ##Include whole dataset, as TOC start at 2009
  group_by(Year,Regime)%>%
  summarise(n = n())%>%
  filter(!is.na(Regime))%>%
  group_by(Year)%>%
  summarise(Percent = 100*n/(sum(n)))%>%
  ungroup()%>%
  add_row(Year = 2004, Percent = 0, .before = 22)



PercentRegime<-PercentRegime[seq(2, nrow(PercentRegime), by =2),] ##Filters out data to only look at percent invalid


PercentPlotRegime<-ggplot(PercentRegime)+
  geom_line(aes(Year, y = Percent), color = 'red', size = 4)+
  geom_point(aes(Year, y = Percent), fill = "red", shape = 21, size = 4)+
  PaperQuality()+scale_y_continuous(TeX('\\textbf{Percent Non-Linear (%)}'))+scale_x_continuous(breaks = seq(1995, 2022, 5))+
  ggtitle("b)")+
  theme(axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size =20))

CondandPerc<-grid.arrange(HCondPlot,PercentPlotRegime, ncol = 2)

ggsave(CondandPerc, filename = "Figure7.png", width =20, height = 6)


## Table 1. Regime Dunn Test ####

RegimeTableData<-AllCloudData%>%
  mutate(`Ion Balance` = Cations - Anions)%>%
  dplyr::select(c(Regime,LABPH:TOC, `Ion Balance`))

names(RegimeTableData)<-paste0(names(RegimeTableData), 
                         c("", "", ' (uS/cm)',' (ueq/L)', ' (ueq/L)',
                           ' (ueq/L)',' (ueq/L)',' (ueq/L)'," (ueq/L)",
                           " (ueq/L)",' (ueq/L)',"", " (g/m^3)",
                           '', '', ' (umolC/L)', ' (ueq/L)'))
                        

DunnTestbyRegime<-RegimeTableData%>%
  #select( dplyr::select(c(Regime,LABPH:TOC, `Ion Balance`)))%>%
 # mutate(`Ion Balance` = Cations - Anions)%>%
  dplyr::select(c(Regime,LABPH:`TOC (umolC/L)`, `Ion Balance (ueq/L)`))%>%
  dplyr::rename(pH = LABPH, `Conductivity (uS/cm) ` = `SPCOND (uS/cm)`, `Na (ueq/L)` = `Sodium (ueq/L)`, 
                `Ca (ueq/L)` = `CA (ueq/L)`,
                `Mg (ueq/L)` = `MG (ueq/L)`, `Cl (ueq/L)` = `CL (ueq/L)`)%>%
  pivot_longer(cols = !Regime, names_to = "Species", values_to = "Conc")%>%
  
  group_by(Species)%>%
  dunn_test(Conc~Regime,p.adjust.method = "none", detailed = TRUE)%>%
  dplyr::select(c("Species", "estimate", "n1", "n2", "p.adj", "p.adj.signif"))



MedianOld<-RegimeTableData%>%
 # mutate(`Ion Balance` = Cations - Anions)%>%
  dplyr::select(c(Regime,LABPH:`TOC (umolC/L)`, `Ion Balance (ueq/L)`))%>%
  dplyr::rename(pH = LABPH, `Conductivity (uS/cm) ` = `SPCOND (uS/cm)`, `Na (ueq/L)` = `Sodium (ueq/L)`, 
                `Ca (ueq/L)` = `CA (ueq/L)`,
                `Mg (ueq/L)` = `MG (ueq/L)`, `Cl (ueq/L)` = `CL (ueq/L)`)%>%
  pivot_longer(cols = !Regime, names_to = "Species", values_to = "Conc")%>%
  filter(Regime == "Linear")%>%
  group_by(Species)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))%>%
  dplyr::rename(`Linear Median` = Conc)


MedianNew<-RegimeTableData%>%
  dplyr::select(c(Regime,LABPH:`TOC (umolC/L)`, `Ion Balance (ueq/L)`))%>%
  dplyr::rename(pH = LABPH, `Conductivity (uS/cm) ` = `SPCOND (uS/cm)`, `Na (ueq/L)` = `Sodium (ueq/L)`, 
                `Ca (ueq/L)` = `CA (ueq/L)`,
                `Mg (ueq/L)` = `MG (ueq/L)`, `Cl (ueq/L)` = `CL (ueq/L)`)%>%
  pivot_longer(cols = !Regime, names_to = "Species", values_to = "Conc")%>%
  filter(Regime == "Non-Linear")%>%
  group_by(Species)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))%>%
  dplyr::rename(`Non-Linear Median` = Conc)



DunnTestbyRegime<-inner_join(DunnTestbyRegime, MedianOld, by = "Species")
DunnTestbyRegime<-inner_join(DunnTestbyRegime, MedianNew, by = "Species")

DunnResults<-DunnTestbyRegime%>%
  dplyr::select(c("Species","n1", "n2", "Non-Linear Median", "Linear Median", "p.adj", "p.adj.signif"))%>%
  arrange(desc(`Non-Linear Median`) )%>%
  dplyr::rename(`Linear n` = n1, `Non-Linear n` = n2)%>%
  mutate(Difference = `Non-Linear Median` - `Linear Median`)%>%
  dplyr::select(c(Species, `Non-Linear Median`,`Linear Median`, Difference, p.adj, `Non-Linear n`, `Linear n`))%>%
  mutate(across(c(`Non-Linear Median`,`Linear Median`, Difference), round, digits = 3),
         p.adj = signif(p.adj, digits = 3))%>%
  mutate(`Non-Linear Median` = format(`Non-Linear Median`, scientific = FALSE),
         `Linear Median` = format(`Linear Median`, scientific = FALSE))%>%
  filter(Species !="Hplus" & Species != "WSOC" & Species != "SampleVolume")%>%
  mutate(`Percent Difference` = 100*(as.numeric(`Non-Linear Median`)-as.numeric(`Linear Median`))/as.numeric(`Non-Linear Median`),
         p.adj = ifelse(p.adj < 0.001, "p < 0.001", p.adj),
         `Percent Difference` = percent(`Percent Difference`/100))%>%
  relocate(`Percent Difference`, .after  = Difference)%>%
  rename(`p-value` = p.adj)
  #mutate(across(where(is.numeric), as.character))%>%
  #mutate(across(`Non-Linear Median`:p, str_remove_all, "0+"))

  
RegimeTable<-ggtexttable(DunnResults, theme = ttheme("light", base_size = 16, colnames.style = colnames_style(face = "italic", size = 16)), rows = NULL)%>%
  tab_add_title(text = paste("Dunn-Test Comparison of Regimes", sep = "\n"), size = 17, face = "bold",padding = unit(1, "line"),
                just = "left")


ggsave(filename = "Table1.png", plot = RegimeTable, width = 14, height = 8)


## Figure S16 Size Resolved Plot of Calcium and Potassium in Toronto ####

CaKAerosolData<-read.pxp('Sample3_TO_090807R_MOUDINEQ_BIN.pxp')
CaKAerosolData<-bind_rows(CaKAerosolData)



CaKAerosolDataNames<-c('Bin', 'DMA', 'TEA', "NH4", "K", "Na", "Mg", "Ca", 
                       "AC", "FM", "Cl", "NO2", "NO3", "SO4", "OX", "DEA")


names(CaKAerosolData)<-CaKAerosolDataNames

CaKAerosolData <- CaKAerosolData[seq(2,nrow(CaKAerosolData),3),] ### Remove triplicate data
CaKAerosolData<-CaKAerosolData%>%
  arrange(Bin)%>%
  mutate(LogDp = log10(Bin/lag(Bin, n = 1)),
         LogDp = ifelse(is.na(LogDp), log10(0.78/0.5), LogDp))%>%
  mutate(TotalMass = 1e-3*((NH4*18.04)+(Na*23)+(K*39.09)+(24.305*Mg)+(DEA*73.14)+(Ca/2*40.08)+(FM*45.03)+(AC*59.05)+(Cl*35.45)+
                            (NO2*46.01)+(NO3*62)+(SO4/2*96.06)+(OX/2*88.03)+
                            (DMA*45.08)+(TEA*101.19)),
         Np = TotalMass/(1.5*(pi/6)*(Bin)^3))%>%
  mutate(MgTot = sum(Mg*LogDp, na.rm = TRUE))%>%
  mutate(MgPer = LogDp*Mg/MgTot)%>%
  mutate(CaTot = sum(Ca*LogDp, na.rm = TRUE))%>%
  mutate(CaPer = LogDp*Ca/CaTot)%>%
  mutate(SO4Tot = sum(SO4*LogDp, na.rm = TRUE))%>%
  mutate(SO4Per = LogDp*SO4/SO4Tot)%>%
  mutate(KTot = sum(K*LogDp, na.rm = TRUE))%>%
  mutate(KPer = LogDp*K/KTot,
         NO3Tot = sum(NO3*LogDp, na.rm = TRUE),
         NO3Per = LogDp*NO3/NO3Tot,
         NpTot = sum(Np),
         NpPer = Np/NpTot)



CaKAerosolMassPlot<-ggplot(CaKAerosolData)+
  geom_line(aes(x= Bin, y = KPer*100, color = "Relative K Mass"), size = 4)+
  geom_point(aes(x= Bin, y = KPer*100, fill = "Relative K Mass"), size = 4, color = "black", shape = 21)+
  geom_line(aes(x= Bin, y = CaPer*100, color = "Relative Ca Mass"), size = 4)+
  geom_point(aes(x= Bin, y = CaPer*100, fill = "Relative Ca Mass"), size = 4, color = "black", shape = 21)+
  geom_line(aes(x = Bin, y = NO3Per*100, color = "Relative NO3 Mass"), size =4)+
  geom_point(aes(x = Bin, y = NO3Per*100, fill = "Relative NO3 Mass"), size =4, shape = 21, color = "black")+
  
  #  geom_line(aes(x= `Sample2$Bin`, y = ClMass/ClTot, color = "Cl"), size = 2)+
  #  geom_point(aes(x= `Sample2$Bin`, y = ClMass/ClTot, color = "Cl"), size = 4)+
  scale_x_log10(name = TeX("\\textbf{Diameter ($\\mu m$)"))+PaperQuality()+
  scale_color_discrete(guide=FALSE)+
  ylab("Percentage of Total Mass (%)")+
  labs(title = 'a)',subtitle = "Percent of Mass of Major Aerosol Species")+
  theme(axis.title.y = element_text(size = 16))

NPPlot<-ggplot(CaKAerosolData)+
  geom_line(aes(x= Bin, y = NpPer*100), color = 'purple',size = 4)+
  geom_point(aes(x= Bin, y = NpPer*100), fill ='purple',size = 4, color = "black", shape = 21)+
  #  geom_line(aes(x= `Sample2$Bin`, y = ClMass/ClTot, color = "Cl"), size = 2)+
  #  geom_point(aes(x= `Sample2$Bin`, y = ClMass/ClTot, color = "Cl"), size = 4)+
  scale_x_log10(name = TeX("\\textbf{Diameter ($\\mu$m)"))+PaperQuality()+
  scale_color_discrete(guide=FALSE)+
  ylab("Percent of Total Number (%)")+
  labs(title= "b)",subtitle = 'Percentage of Total Aerosol Number')+
  theme(axis.title.y = element_text(size = 16))


AllAerosolPlots<-grid.arrange(CaKAerosolMassPlot,NPPlot)

ggsave(filename = "FigureS14.png", plot = AllAerosolPlots, width = 10, height = 12)



## Figure 8. Drop BU and Drop TD pH vs Measured pH ####



TDBUpH<-AllCloudData%>%
  mutate(TDHplus = CA+MG+1e6*10^(-LABPH),
         BUHplus = (SO4+NO3-NH4)*1e-6)%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))%>%
  mutate(TDpH = -log10(TDHplus*1e-6))%>% ## Convert from ueq
  ggplot()+
  geom_line(aes(x= Year, y = TDpH, color = "pH TD",
                ), show.legend = FALSE,size = 4)+
  geom_point(aes(x = Year, y = TDpH, fill = ("pH TD")), size = 4, shape = 21)+
  geom_line(aes(x= Year, y = -log10(Hplus), color = "Measured pH"), 
            size = 4, show.legend = FALSE)+
  geom_point(aes(x = Year, y = -log10(Hplus), fill = "Measured pH"), size = 4,color = 'black', shape = 21)+
  geom_line(aes(x= Year, y = -log10(BUHplus), color = "pH BU"), 
            size = 4, show.legend = FALSE)+
  geom_point(aes(x = Year, y = -log10(BUHplus), fill = "pH BU"), size = 4,color = 'black', shape = 21)+
  PaperQuality()+scale_x_continuous(breaks = seq(1995, 2020, by = 5))+
  ylab('pH')
  
  
ggsave(plot = TDBUpH, filename = "Figure8.png",
       width = 8, height = 8)




##Figure 9 Samples with NH4 > SO4+NO3 ####



PercentNH4Surplus<-AllCloudData%>%
  dplyr::select(-TOC)%>% ##Include whole dataset, as TOC start at 2009
  group_by(Year,NH4Surplus)%>%
  summarise(n = n())%>%
  na.omit()%>%
  ungroup()%>%
  add_row(Year = 2004, NH4Surplus = "Yes", n = 0, .after = 21)%>%
  group_by(Year)%>%
  summarise(Percent = 100*n/(sum(n)))
PercentNH4Surplus<-PercentNH4Surplus[seq(2, nrow(PercentNH4Surplus), by =2),]


FractionalAcidity<-AllCloudData%>%
  mutate(SurplusNH4 = ifelse(SO4+NO3-NH4 > 0, SO4+NO3-NH4, 0))%>%
  mutate(FractionAcidity = (CA+MG+1e6*10^(-LABPH)-SurplusNH4)/(CA+MG+1e6*10^(-LABPH)),
         FractionAcidityZero = ifelse(SurplusNH4 == 0, NA, FractionAcidity))%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))

PercentSurplusMAF<-ggplot(PercentNH4Surplus)+
  geom_line(aes(x = Year, y = Percent/100), color = 'red', show.legend = FALSE ,size =3)+
  geom_point(aes(x = Year, y = Percent/100), fill = 'red',size =3,  show.legend = FALSE,shape =21)+
  geom_line(aes(x = Year, y = 0.5*FractionAcidity, color = 'black'), 
            size =3, data = FractionalAcidity)+
  geom_point(aes(x = Year, y = 0.5*FractionAcidity),fill = 'black',
             size =3, color ='black', shape =21,
             data = FractionalAcidity)+
  geom_line(aes(x = Year, y = 0.5*FractionAcidityZero, color = 'gray'), 
            size =3, data = FractionalAcidity)+
  geom_point(aes(x = Year, y = 0.5*FractionAcidityZero),fill = 'gray',
             size =3, color ='black', shape =21,
             data = FractionalAcidity)+
  geom_hline(yintercept = 0, size =2)+
  scale_x_continuous(breaks = seq(1994, 2022, 4))+
  PaperQuality()+
  scale_y_continuous(sec.axis = sec_axis(~.*2, name = TeX('\\textbf{Missing Acid Fraction}'),
                                         breaks = seq(0,1,0.2)))+
  theme(axis.title.y.left = element_text(color = 'red'),
        axis.text.y.left = element_text(color = 'red'),
        legend.position = 'bottom')+
  scale_color_manual(name = 'Test', labels = c(TeX("MAF (setting H$^+_{BU}$ = 0 w/ NH$_4^+$ > SO$_4^{2-}$+NO$_3^-$)"),
                                               TeX("MAF (excluding samples w/ negative values for H$^+_{BU}$)")), 
                     values = c('black', 'gray'))+
  #scale_fill_manual(name = 'Test', labels = c(TeX("A$^+$"),"B"), values = c('black', 'gray'))+
  labs(y = TeX("\\textbf{Fraction of Samples NH$_4^+$ > SO$_4^{2-}$+NO$_3^-$}"))+
  guides(color = guide_legend(nrow = 2,label.hjust = 0))
  
  
ggsave(filename = 'Figure9.png', plot = PercentSurplusMAF, width = 9, height =7)




PercentNH4SurplusPlot<-ggplot(PercentNH4Surplus)+
  geom_line(aes(Year, y = Percent), color = 'red', size = 4)+
  geom_point(aes(Year, y = Percent), fill = "red", shape = 21, size = 4)+
  PaperQuality()+scale_y_continuous(TeX('\\textbf{Percent NH$_4^+$ > SO$_4^{2-}$ + NO$_3^-$  (%)}'))+scale_x_continuous(breaks = seq(1994, 2022, 4))+
  labs(title = 'a)')



##  Figure S17 Measured H+ vs Predicted H+ ####

LinearCloudTS<-theilsen(Hplus~BUHplus, data = AllCloudData%>%filter(NH4< SO4+NO3 & Regime =='Linear')%>%
                          mutate(Hplus = 1e6*10^(-LABPH)))
NonLinearCloudTS<-theilsen(Hplus~BUHplus, data = AllCloudData%>%filter(NH4< SO4+NO3 & Regime =='Non-Linear')%>%
                             mutate(Hplus = 1e6*10^(-LABPH)))
LinearCloudAdjustTS<-theilsen(TDHplus~BUHplus, data = AllCloudData%>%filter(NH4< SO4+NO3 & Regime =='Linear')%>%
                                mutate(Hplus = 1e6*10^(-LABPH)))
NonLinearCloudAdjustTS<-theilsen(TDHplus~BUHplus, data = AllCloudData%>%filter(NH4< SO4+NO3 & Regime =='Non-Linear')%>%
                                   mutate(Hplus = 1e6*10^(-LABPH)))

MeasuredLinearpHRegime<-ggplot(subset(AllCloudData, Regime == "Linear" & SO4+NO3>NH4))+
  geom_point(aes(x = NO3+SO4-NH4, y = CA+MG+(1e6*10^(-LABPH)), 
                 color = "H+ TD"), size = 3)+
  geom_point(aes(x = NO3+SO4-NH4, y = (1e6*10^(-LABPH)), 
                 color = 'Measured H+'), size = 3)+
  geom_abline(aes(slope = LinearCloudTS$coefficients[[2]], intercept = LinearCloudTS$coefficients[[1]], color = 'Measured H+'),
              size = 2, show.legend = FALSE)+
  geom_abline(aes(slope = LinearCloudAdjustTS$coefficients[[2]], intercept = LinearCloudAdjustTS$coefficients[[1]], color = "H+ TD"),
              size = 2, show.legend = FALSE)+
  geom_abline(slope =1 , intercept = 0, linetype = 'dashed', size = 2)+
  annotate(geom = 'text',x = 900, y = 2900, size = 8,
           label = paste("y=",format(round(LinearCloudTS$coefficients[[2]], digits = 2),nsmall = 2),"x + ",round(LinearCloudTS$coefficients[[1]], digits = 3)),
           color ="#00BFC4")+
  annotate(geom = 'text',x = 900, y = 2600, size = 8,
           label = paste("y=",round(LinearCloudAdjustTS$coefficients[[2]], digits = 2),"x + ",round(LinearCloudAdjustTS$coefficients[[1]], digits = 3)),
           color = "#F8766D")+
  PaperQuality()+
  labs(title = "b)",subtitle    = "Linear Regime",
       x = TeX('\\textbf{(SO$_4^{2-}$+NO$_3^{-}$) - NH$_4^+$ ($\\mu$eq L$^{-1}$)}'),
       y = TeX('\\textbf{H$^+$ ($\\mu$eq L$^{-1}$)}'))

MeasuredNonLinearpHRegime<-ggplot(subset(AllCloudData, Regime == 'Non-Linear' & SO4+NO3>NH4))+
  geom_point(aes(x = NO3+SO4-NH4, y = CA+MG+(1e6*10^(-LABPH)), 
                 color = "H+ TD"), size = 3)+
  geom_point(aes(x = NO3+SO4-NH4, y = (1e6*10^(-LABPH)), 
                 color = 'Measured H+'), size = 3)+
  geom_abline(aes(slope = NonLinearCloudTS$coefficients[[2]], intercept = NonLinearCloudTS$coefficients[[1]], color = 'Measured H+'),
              size = 2, show.legend = FALSE)+
  geom_abline(aes(slope = NonLinearCloudAdjustTS$coefficients[[2]], intercept = NonLinearCloudAdjustTS$coefficients[[1]], color = 'H+ TD'),
              size = 2, show.legend = FALSE)+
  geom_abline(slope =1 , intercept = 0, linetype = 'dashed', size = 2)+
  annotate(geom = 'text',x = 560, y = 1750, size = 8,
           label = paste("y=",format(round(NonLinearCloudTS$coefficients[[2]], digits = 2), nsmall = 2),"x + ",round(NonLinearCloudTS$coefficients[[1]], digits = 1)),
           color = "#00BFC4" )+
  annotate(geom = 'text',x = 575, y = 1550, size = 8,
           label = paste("y=",round(NonLinearCloudAdjustTS$coefficients[[2]], digits = 2),"x + ",round(NonLinearCloudAdjustTS$coefficients[[1]], digits = 1)),
           color ="#F8766D")+
  PaperQuality()+
  labs(title = 'c)',subtitle = "Non-Linear Regime",
       x = TeX('\\textbf{(SO$_4^{2-}$+NO$_3^{-}$) - NH$_4^+$ ($\\mu$eq L$^{-1}$)}'),
       y = TeX('\\textbf{H$^+$ ($\\mu$eq L$^{-1}$)}'))  

EstimatedpHPlots<-grid.arrange(PercentNH4SurplusPlot,MeasuredLinearpHRegime, MeasuredNonLinearpHRegime, layout_matrix = rbind(c(1,1),c(2,3)))


ggsave(filename = "Figure9.png", plot = EstimatedpHPlots,
       width = 12, height = 12)


## Figure 10 Fractional Acidity from Organics ####



## Figure S18 TD pH and Measured pH vs LWC####

LWCSO4PlotInferrred<-ggplot(AllCloudData)+
  geom_point(aes(x = LWC, y = InferredpH, color = SO4), size = 3)+
  scale_color_gradientn(colors = rev(rainbow(6)),  trans= 'log10')+
  PaperQuality()+labs(x  = TeX("\\textbf{LWC (g m$^{-3}$)}"),
                      y = TeX("\\textbf{pH_{TD}}"),
                      color = TeX('\\textbf{SO$_4^{2-}$ ($\\mu$eq L$^{-1}$)}'))+xlim(c(0.05, 1.5))+
  theme(legend.position = 'right', legend.title = element_text(size =14))+
  ylim(c(2,7))


LWCSO4PlotMeasured<-ggplot(AllCloudData)+
  geom_point(aes(x = LWC, y = LABPH, color = SO4), size = 3)+
  scale_color_gradientn(colors = rev(rainbow(6)),  trans= 'log10')+
  PaperQuality()+labs(x  = TeX("\\textbf{LWC (g m$^{-3}$)}"),
                      y = TeX("\\textbf{Measured pH}"),
                      color = TeX('\\textbf{SO$_4^{2-}$ ($\\mu$eq L$^{-1}$)}'))+xlim(c(0.05, 1.5))+
  theme(legend.position = 'right', legend.title = element_text(size =14))+
  ylim(c(2,7))



LWCvsPH<-grid.arrange(LWCSO4PlotInferrred, LWCSO4PlotMeasured)

ggsave("FigureS15.png", plot = LWCvsPH,  width = 8, height = 8)

## Figure 10 TOC vs TDpH

TOCvsMeasuredpH<-ggplot(AllCloudData%>%
                          filter(Year > 2008))+
  geom_point(aes(x = LABPH, y = TOC, color = Year), size = 3)+
  scale_color_gradientn(colors = rev(rainbow(6)))+
  PaperQuality()+labs(y  = TeX("\\textbf{NH$_4^+$ ($\\mu$mol L$^{-}$)}"),
                      x = TeX("\\textbf{pH_{TD}}"),
                      color = TeX('\\textbf{Year}}'))+
  theme(legend.position = 'right', legend.title = element_text())



TOCvsTDpH<-ggplot(AllCloudData%>%
                      filter(Year > 2008))+
  geom_point(aes(x = InferredpH, y = TOC, color = Year), size = 3)+
  scale_color_gradientn(colors = rev(rainbow(6)))+
  PaperQuality()+labs(y  = TeX("\\textbf{NH$_4^+$ ($\\mu$mol L$^{-1}$)}"),
                      x = TeX("\\textbf{pH_{TD}}"),
                      color = TeX('\\textbf{Year}}'))+
  theme(legend.position = 'right', legend.title = element_text())


TOCpHPlot<-grid.arrange(TOCvsMeasuredpH, TOCvsTDpH)

ggsave(TOCpHPlot, file ='Figure10.png', width = 14, height = 8)

## Figure S19 NH4 and NO3 TD pH ####

NH4vsMeasuredpH<-ggplot(AllCloudData)+
  geom_point(aes(x = LABPH, y = NH4, color = Year), size = 3)+
  scale_color_gradientn(colors = rev(rainbow(6)))+
  PaperQuality()+labs(y  = TeX("\\textbf{NH$_4^+$ ($\\mu$mol L$^{-1}$)}"),
                      x = TeX("\\textbf{pH_{Measured}}"),
                      color = TeX('\\textbf{Year}}'))+
  theme(legend.position = 'right', legend.title = element_text())+
  labs(title = "a)",
       subtitle =TeX('\\textbf{NH$_4^+$ vs Measured pH}'))

NO3vsMeasuredpH<-ggplot(AllCloudData)+
  geom_point(aes(x = LABPH, y = NO3, color = Year), size = 3)+
  scale_color_gradientn(colors = rev(rainbow(6)))+
  PaperQuality()+labs(y  = TeX("\\textbf{NO$_3^-$ ($\\mu$mol L$^{-1}$)}"),
                      x = TeX("\\textbf{pH_{Measured}}"),
                      color = TeX('\\textbf{Year}}'))+
  theme(legend.position = 'right', legend.title = element_text())+
  labs(title = "b)",
       subtitle = TeX('\\textbf{NO$_3^-$ vs Measured pH}'))



NH4vsTDpH<-ggplot(AllCloudData)+
  geom_point(aes(x = InferredpH, y = NH4, color = Year), size = 3)+
  scale_color_gradientn(colors = rev(rainbow(6)))+
  PaperQuality()+labs(y  = TeX("\\textbf{NH$_4^+$ ($\\mu$mol L$^{-1}$)}"),
                      x = TeX("\\textbf{pH_{TD}}"),
                      color = TeX('\\textbf{Year}}'))+
  theme(legend.position = 'right', legend.title = element_text())+
  labs(title = "c)",
       subtitle = TeX('\\textbf{NH$_4^+$ vs Inferred Cloud Droplet pH}'))

NO3vsTDpH<-ggplot(AllCloudData)+
  geom_point(aes(x = InferredpH, y = NO3, color = Year), size = 3)+
  scale_color_gradientn(colors = rev(rainbow(6)))+
  PaperQuality()+labs(y  = TeX("\\textbf{NO$_3^-$ ($\\mu$mol L$^{-1}$)}"),
                      x = TeX("\\textbf{pH_{TD}}"),
                      color = TeX('\\textbf{Year}}'))+
  theme(legend.position = 'right', legend.title = element_text())+
  labs(title = "d)",
       subtitle = TeX('\\textbf{NO$_3^-$ vs Inferred Cloud Droplet pH}'))

TDpHConcentrations<-grid.arrange(NH4vsMeasuredpH, NO3vsMeasuredpH,NH4vsTDpH, NO3vsTDpH, ncol = 2)

ggsave(TDpHConcentrations, file ='FigureS16.png', width = 14, height = 8)


## Figure 10 TOC vs Adjusted pH####

ExpFun<-function(x){
  A<-1.7193e5
  Tau<-1/1.3934
  A*exp(-x/Tau)
}

TOCvsBulkpH<-ggplot(AllCloudData%>%
         filter(Year > 2008))+
  geom_point(aes(x = LABPH, y = TOC, color = Year), size =3)+
  geom_function(fun = ExpFun, color = 'black', size =2)+
  scale_color_gradientn(colors = rev(rainbow(6)),
                        breaks = seq(2008, 2022, 2))+
  PaperQuality()+labs(x = TeX("\\textbf{Measured Bulk Cloud Water pH}"),
                      y = TeX("\\textbf{TOC ($\\mu mol L^{-1}$)}"))+
  theme(legend.position = 'right', legend.key.height = unit(1.4, 'cm'))+
  xlim(c(3,7.25))

TOCvsInferredpH<-ggplot(AllCloudData%>%
         filter(Year > 2008))+
  geom_point(aes(x = InferredpH, y = TOC, color = Year), size =3)+
  geom_function(fun = ExpFun, color = 'black', size =2)+
  scale_color_gradientn(colors = rev(rainbow(6)),
                        breaks = seq(2008, 2022, 2))+
  PaperQuality()+labs(x = TeX("\\textbf{Inferred Bulk Cloud Water pH}"),
                      y = TeX("\\textbf{TOC ($\\mu mol L^{-1}$)}"))+
  theme(legend.position = 'right', legend.key.height = unit(1.4, 'cm'))+
  xlim(c(3,7.25))

AllTOCvspH<-grid.arrange(TOCvsInferredpH, TOCvsBulkpH)
  
ggsave(plot = AllTOCvspH, filename = 'Figure11.png', width = 12, height =10)

