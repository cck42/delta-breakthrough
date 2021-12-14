library(tidyverse)
library(gtools)
library(lubridate)
library(readxl)
library(stringr)
library(cowplot)

## Read in all of the Grubaugh lab genomic survillance metadata
viro_Ct_combined <- read_xlsx("/Users/chaneykalinich/Documents/SARS-CoV-2/Vaccine_breakthrough/maintain_spreadsheets/combined_forCt.xlsx")

## Filter lines that don't have diagnostic test platform and Ct for both diagnostic and research lab
combined2 <- viro_Ct_combined %>%
  filter(!is.na(`Test Method.x`)&`Test Method.x`!="NA") %>%
  filter(!is.na(`Yale.N1.FAM.`)) %>%
  filter(!is.na(Value)&`Value`!="NA") %>%
  select(Sample.ID,Yale.N1.FAM.,Yale.N1_GE.mL,`Test Method.x`,Value)
  
#What diagnostic test platforms are represented and how frequently?
counts <- combined2 %>%
  ddply(.(`Test Method.x`),nrow)

##-------------------------------
##Model for YH cobas 6800

#Clean up the data
cobas <- combined2 %>%
  filter(`Test Method.x`=="YH cobas 6800") %>%
  mutate(vals=str_replace(Value,"Detected, ","")) %>%
  mutate(vals=str_replace(vals,"Presumptive Positive, ","")) %>%
  mutate(endat=str_locate(vals,",")) %>%
  mutate(vals2=str_sub(vals,1,4)) %>%
  mutate(vals2=as.numeric((vals2))) %>%
  mutate(Yale_N1=signif(as.numeric(Yale.N1.FAM.),digits=4))

#Visualize the regression
cobasviz <- cobas %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  ggplot(aes(x=vals2,y=Yale_N1)) + geom_abline(intercept=0, slope=1, linetype="dashed", col="grey") + 
    geom_point() + 
    geom_line(stat="smooth", method="lm", fullrange=TRUE) + 
    theme_minimal() + 
    theme(text=element_text(size=14)) + 
    labs(x="Cobas Ct", y="Yale Ct (N1)") +
    #labs(x="",y="")+
    xlim(15,45) +
    ylim(15,45)
plot(cobasviz)
#Run a linear regression
cobasfit <- cobas %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  lm(Yale_N1 ~ vals2,data=.)
#Save the coefficients
cobas_b0 <- cobasfit$coef[[1]]
cobas_b1 <- cobasfit$coef[[2]]

#look at residuals
cobas_resid <- cobas %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  mutate(resid = Yale_N1-(cobas_b0+cobas_b1*vals2)) %>%
  ggplot(aes(x=vals2, y=resid)) + 
  geom_point() + 
  theme_minimal()  + 
  theme(text=element_text(size=14)) + 
  labs(x="Cobas Ct value", y="Yale Ct (N1) residual") +
  #labs(x="",y="")+
  scale_y_continuous(limits=c(-15,15))
plot(cobas_resid)
#Store the SD of the residuals
cobas_sd <- sd(cobas_resid$data$resid)

##-------------------------------
##Model for LMH Genexpert Infinity 48
lmhgenex <- combined2 %>%
  filter(`Test Method.x`=="LMH GENEXPERT INFINITY 48" | `Test Method.x`=="LMH GENEXPERT INFINITY") %>%
  mutate(vals2=as.numeric((Value))) %>%
  mutate(Yale_N1=signif(as.numeric(Yale.N1.FAM.),digits=4))

lmhgenex_viz <- lmhgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  ggplot(aes(x=vals2,y=Yale_N1)) + geom_abline(intercept=0, slope=1, linetype="dashed", col="grey") + 
  geom_point() + 
  geom_line(stat="smooth", method="lm", fullrange=TRUE) + 
  theme_minimal() + 
  theme(text=element_text(size=14)) + 
  labs(x="LMH Genexpert Infinity 48 Ct", y="Yale Ct (N1)")+
  #labs(x="",y="")+
  xlim(15,45) +
  ylim(15,45)
plot(lmhgenex_viz)
lmhgenex_fit <- lmhgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  lm(Yale_N1 ~ vals2,data=.)
lmhgenex_b0 <- lmhgenex_fit$coef[[1]]
lmhgenex_b1 <- lmhgenex_fit$coef[[2]]

#look at residuals
lmhgenex_resid <- lmhgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  mutate(resid = Yale_N1-(lmhgenex_b0+lmhgenex_b1*vals2)) %>%
  ggplot(aes(x=vals2, y=resid)) + 
  geom_point() + 
  theme_minimal()  + 
  theme(text=element_text(size=14)) + 
  labs(x="LMH Genexpert Infinity 48 Ct value", y="Yale Ct (N1) residual")+
  #labs(x="",y="")+
  scale_y_continuous(limits=c(-15,15))
plot(lmhgenex_resid)
lmhgenex_sd <- sd(lmhgenex_resid$data$resid)

##-------------------------------
##Model for SRC Cepheid Genexpert 2
srcgenex <- combined2 %>%
  filter(`Test Method.x`=="SRC CEPHEID GENEXPERT 2" | `Test Method.x`=="SRC CEPHEID GENEXPERT") %>%
  mutate(vals=str_replace(Value,"Presumptive Positive, ","")) %>%
  mutate(vals=str_replace(vals,"Positive, ","")) %>%
  mutate(vals2=str_sub(vals,1,4)) %>%
  mutate(vals2=as.numeric((vals2))) %>%
  mutate(Yale_N1=signif(as.numeric(Yale.N1.FAM.),digits=4))

srcgenex_viz <- srcgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  ggplot(aes(x=vals2,y=Yale_N1)) + geom_abline(intercept=0, slope=1, linetype="dashed", col="grey") + 
  geom_point() + 
  geom_line(stat="smooth", method="lm", fullrange=TRUE) + 
  theme_minimal() + 
  theme(text=element_text(size=14)) + 
  labs(x="SRC Cepheid Genexpert 2 Ct", y="Yale Ct (N1)")+
  #labs(x="",y="")+
  xlim(15,45) +
  ylim(15,45)
plot(srcgenex_viz)
srcgenex_fit <- srcgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  lm(Yale_N1 ~ vals2,data=.)
srcgenex_b0 <- srcgenex_fit$coef[[1]]
srcgenex_b1 <- srcgenex_fit$coef[[2]]

#look at residuals
srcgenex_resid <- srcgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  mutate(resid = Yale_N1-(srcgenex_b0+srcgenex_b1*vals2)) %>%
  ggplot(aes(x=vals2, y=resid)) + 
  geom_point() + 
  theme_minimal()  + 
  theme(text=element_text(size=14)) + 
  labs(x="SRC Cepheid Genexpert 2 Ct value", y="Yale Ct (N1) residual")+
  #labs(x="",y="")+
  scale_y_continuous(limits=c(-15,15))
plot(srcgenex_resid)
srcgenex_sd <- sd(srcgenex_resid$data$resid)

##-------------------------------
##Model for YH Cepheid Genexpert Infinity

yhgenex <- combined2 %>%
  filter(`Test Method.x`=="YH CEPHEID GENEXPERT INFINITY") %>%
  mutate(vals=str_replace(Value,"Detected, ","")) %>%
  mutate(vals=str_replace(vals,"Positive, ","")) %>%
  mutate(vals2=str_sub(vals,1,4)) %>%
  mutate(vals2=as.numeric((vals2))) %>%
  mutate(Yale_N1=signif(as.numeric(Yale.N1.FAM.),digits=4))

yhgenex_viz <- yhgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  ggplot(aes(x=vals2,y=Yale_N1)) + geom_abline(intercept=0, slope=1, linetype="dashed", col="grey") + 
  geom_point() + 
  geom_line(stat="smooth", method="lm", fullrange=TRUE) + 
  theme_minimal() + 
  theme(text=element_text(size=14)) + 
  labs(x="YH Cepheid Genexpert Infinity Ct", y="Yale Ct (N1)")+
  #labs(x="",y="")+
  xlim(15,45) +
  ylim(15,45)
plot(yhgenex_viz)
yhgenex_fit <- yhgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  lm(Yale_N1 ~ vals2,data=.)
yhgenex_b0 <- yhgenex_fit$coef[[1]]
yhgenex_b1 <- yhgenex_fit$coef[[2]]

#look at residuals
yhgenex_resid <- yhgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  mutate(resid = Yale_N1-(yhgenex_b0+yhgenex_b1*vals2)) %>%
  ggplot(aes(x=vals2, y=resid)) + 
  geom_point() + 
  theme_minimal()  + 
  theme(text=element_text(size=14)) + 
  labs(x="YH Cepheid Genexpert Infinity Ct value", y="Yale Ct (N1) residual")+
  #labs(x="",y="")+
  scale_y_continuous(limits=c(-15,15))
plot(yhgenex_resid)
yhgenex_sd <- sd(yhgenex_resid$data$resid)

##-------------------------------
##Model for BH Cepheid Genexpert (Micro)

bhmicrogenex <- combined2 %>%
  filter(`Test Method.x`=="BH CEPHEID GENEXPERT (MICRO)") %>%
  mutate(vals=str_replace(Value,"Presumptive Positive, ","")) %>%
  mutate(vals=str_replace(vals,"Positive, ","")) %>%
  mutate(vals2=str_sub(vals,1,4)) %>%
  mutate(vals2=as.numeric((vals2))) %>%
  mutate(Yale_N1=signif(as.numeric(Yale.N1.FAM.),digits=4))

bhmicrogenex_viz <- bhmicrogenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  ggplot(aes(x=vals2,y=Yale_N1)) + geom_abline(intercept=0, slope=1, linetype="dashed", col="grey") + 
  geom_point() + 
  geom_line(stat="smooth", method="lm", fullrange=TRUE) + 
  theme_minimal() + 
  theme(text=element_text(size=14)) + 
  labs(x="BH CEPHEID GENEXPERT (MICRO) Ct", y="Yale Ct (N1)")+
  #labs(x="",y="")+
  xlim(15,45) +
  ylim(15,45)
plot(bhmicrogenex_viz)
bhmicrogenex_fit <- bhmicrogenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  lm(Yale_N1 ~ vals2,data=.)
bhmicrogenex_b0 <- bhmicrogenex_fit$coef[[1]]
bhmicrogenex_b1 <- bhmicrogenex_fit$coef[[2]]

#look at residuals
bhmicrogenex_resid <- bhmicrogenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  mutate(resid = Yale_N1-(bhmicrogenex_b0+bhmicrogenex_b1*vals2)) %>%
  ggplot(aes(x=vals2, y=resid)) + 
  geom_point() + 
  theme_minimal()  + 
  theme(text=element_text(size=14)) + 
  labs(x="BH CEPHEID GENEXPERT (MICRO) Ct", y="Yale Ct (N1) residual")+
  #labs(x="",y="")+
  scale_y_continuous(limits=c(-15,15))
plot(bhmicrogenex_resid)
bhmicrogenex_sd <- sd(bhmicrogenex_resid$data$resid)

##-------------------------------
##Model for BH Cepheid Genexpert (Chem)

bhchemgenex <- combined2 %>%
  filter(`Test Method.x`=="BH CEPHEID GENEXPERT (CHEM)") %>%
  mutate(vals=str_replace(Value,"Presumptive Positive, ","")) %>%
  mutate(vals=str_replace(vals,"Positive, ","")) %>%
  mutate(vals2=str_sub(vals,1,4)) %>%
  mutate(vals2=as.numeric((vals2))) %>%
  mutate(Yale_N1=signif(as.numeric(Yale.N1.FAM.),digits=4))

bhchemgenex_viz <- bhchemgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  ggplot(aes(x=vals2,y=Yale_N1)) + geom_abline(intercept=0, slope=1, linetype="dashed", col="grey") + 
  geom_point() + 
  geom_line(stat="smooth", method="lm", fullrange=TRUE) + 
  theme_minimal() + 
  theme(text=element_text(size=14)) + 
  labs(x="BH CEPHEID GENEXPERT (CHEM) Ct", y="Yale Ct (N1)")+
  #labs(x="",y="")+
  xlim(15,45) +
  ylim(15,45)
plot(bhchemgenex_viz)
bhchemgenex_fit <- bhchemgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  lm(Yale_N1 ~ vals2,data=.)
bhchemgenex_b0 <- bhchemgenex_fit$coef[[1]]
bhchemgenex_b1 <- bhchemgenex_fit$coef[[2]]

#look at residuals
bhchemgenex_resid <- bhchemgenex %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  mutate(resid = Yale_N1-(bhchemgenex_b0+bhchemgenex_b1*vals2)) %>%
  ggplot(aes(x=vals2, y=resid)) + 
  geom_point() + 
  theme_minimal()  + 
  theme(text=element_text(size=14)) + 
  labs(x="BH CEPHEID GENEXPERT (CHEM) Ct", y="Yale Ct (N1) residual")+
  #labs(x="",y="")+
  scale_y_continuous(limits=c(-15,15))
plot(bhchemgenex_resid)
bhchemgenex_sd <- sd(bhchemgenex_resid$data$resid)


##-------------------------------
##Model for SMC Cepheid Genexpert

smc <- combined2 %>%
  filter(`Test Method.x`=="SMC CEPHEID GENEXPERT") %>%
  mutate(vals=str_replace(Value,"Presumptive Positive, ","")) %>%
  mutate(vals=str_replace(vals,"Positive, ","")) %>%
  mutate(vals2=str_sub(vals,1,4)) %>%
  mutate(vals2=as.numeric((vals2))) %>%
  mutate(Yale_N1=signif(as.numeric(Yale.N1.FAM.),digits=4))

smc_viz <- smc %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  ggplot(aes(x=vals2,y=Yale_N1)) + geom_abline(intercept=0, slope=1, linetype="dashed", col="grey") + 
  geom_point() + 
  geom_line(stat="smooth", method="lm", fullrange=TRUE) + 
  theme_minimal() + 
  theme(text=element_text(size=14)) + 
  labs(x="SMC CEPHEID GENEXPERT Ct", y="Yale Ct (N1)")+
  #labs(x="",y="")+
  xlim(15,45) +
  ylim(15,45)
plot(smc_viz)
smc_fit <- smc %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  lm(Yale_N1 ~ vals2,data=.)
smc_b0 <- smc_fit$coef[[1]]
smc_b1 <- smc_fit$coef[[2]]

#look at residuals
smc_resid <- smc %>%
  filter(!is.na(Yale_N1) & !is.na(vals2)) %>%
  mutate(resid = Yale_N1-(smc_b0+smc_b1*vals2)) %>%
  ggplot(aes(x=vals2, y=resid)) + 
  geom_point() + 
  theme_minimal()  + 
  theme(text=element_text(size=14)) + 
  labs(x="SMC CEPHEID GENEXPERT Ct", y="Yale Ct (N1) residual")+
  #labs(x="",y="")+
  scale_y_continuous(limits=c(-15,15))
plot(smc_resid)
smc_sd <- sd(smc_resid$data$resid)

mydata <- read_xlsx("/Users/chaneykalinich/Downloads/Delta Final 11_3_21_JK_V2.xlsx")

counts2 <- mydata %>%
  ddply(.(`test_method`),nrow)
B0 <- case_when(
  mydata$test_method=="YH cobas 6800"~cobas_b0,
  mydata$test_method=="LMH GENEXPERT INFINITY 48"~lmhgenex_b0,
  mydata$test_method=="SRC CEPHEID GENEXPERT 2"~srcgenex_b0,
  mydata$test_method=="YH CEPHEID GENEXPERT INFINITY"~yhgenex_b0,
  mydata$test_method=="BH CEPHEID GENEXPERT (MICRO)"~bhmicrogenex_b0,
  mydata$test_method=="BH CEPHEID GENEXPERT (CHEM)"~bhchemgenex_b0,
  mydata$test_method=="SMC CEPHEID GENEXPERT"~smc_b0,
  TRUE ~ 0
)
B1 <- case_when(
  mydata$test_method=="YH cobas 6800"~cobas_b1,
  mydata$test_method=="LMH GENEXPERT INFINITY 48"~lmhgenex_b1,
  mydata$test_method=="SRC CEPHEID GENEXPERT 2"~srcgenex_b1,
  mydata$test_method=="YH CEPHEID GENEXPERT INFINITY"~yhgenex_b1,
  mydata$test_method=="BH CEPHEID GENEXPERT (MICRO)"~bhmicrogenex_b1,
  mydata$test_method=="BH CEPHEID GENEXPERT (CHEM)"~bhchemgenex_b1,
  mydata$test_method=="SMC CEPHEID GENEXPERT"~smc_b1,
  TRUE ~ 1
)
mydata2 <- mydata %>%
  bind_cols(B0,B1) %>%
  mutate(Ct_adj=B0+B1*CT_mean) #%>%
mydata3 <- mydata2 %>%
  mutate(viral_copies=(10^((Ct_adj-43.3754126984127)/(-3.76219047619048)))*75*(1000/300)) #%>%

regressionfig <- plot_grid(cobasviz,lmhgenex_viz,bhchemgenex_viz,bhmicrogenex_viz,
                           smc_viz,srcgenex_viz,yhgenex_viz)
regressionfig
residfig <- plot_grid(cobas_resid,lmhgenex_resid,bhchemgenex_resid,bhmicrogenex_resid,
                      smc_resid,srcgenex_resid,yhgenex_resid)
residfig
