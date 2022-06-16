#Script by Dr. Rodrigo Quiroga - rquiroga7 on Github rquiroga777 on Twitter
#Portugal daily mortality data from:
#https://evm.min-saude.pt/
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(readr,MASS,dplyr,scales,data.table,ggplot2,wesanderson,lubridate)

#Declare help functions, notin and moving average
`%notin%` <- Negate(`%in%`)
ma <- function(x, n = 7){stats::filter(x, rep(1 / n, n), sides = 1)}

Sys.setlocale("LC_ALL","Portuguese")
#setwd(dir = "./Portugal_mortality/")
#Read and process data, remove last day due to right censoring
dataPOR <- read_csv("./2022_06_15_Dados_SICO_2022-02-21.csv",guess_max = 250,col_names = TRUE)
#load and process OWID data
dataOWID <- read_csv("./Portugal-owid-covid-data.csv",guess_max = 250,col_names = TRUE)
dataOWID<-dataOWID %>% mutate(year=year(date)) %>% dplyr::select(date,year,new_deaths_smoothed) %>% mutate(Date=as.Date(as.character(paste0("2020-",month(date),"-",day(date))))  )
dataOWID %>% mutate(Date=as.Date(as.character(paste0("2020-",month(date),"-",day(date))))  )
dataOWID$year<-as.factor(as.integer(dataOWID$year))
dataPOR$Data<-as.Date(paste0("2020-",dataPOR$Data),format="%Y-%b-%d")
dataPOR<-dataPOR %>% filter(Data!=as.Date("2020-02-29"))
today<-dataPOR$Data[max(which(!is.na(dataPOR$'2022')))]
dataPOR[which(dataPOR$Data==today),15]<-NA

dataPOR2<-dataPOR
names(dataPOR2)<-c("Date","a2009","a2010","a2011","a2012","a2013","a2014","a2015","a2016","a2017","a2018","a2019","a2020","a2021","a2022")
dataPOR2<-dataPOR2 %>% mutate(prom=(a2009+a2010+a2011+a2012+a2013+a2014+a2015+a2016+a2017+a2018+a2019)/11)
dataPOR3<-dataPOR2 %>% mutate(across(a2009:a2022, ~ .x - prom)) %>% dplyr::select(Date,a2009:a2022)
names(dataPOR3)<-c("Date","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022")
names(dataPOR)<-c("Date","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022")

long <- melt(setDT(dataPOR), id.vars = "Date", variable.name = "year")
long2 <- melt(setDT(dataPOR3), id.vars = "Date", variable.name = "year")
long4 <- left_join(long,dataOWID,by=c("Date","year"))

#long<-long %>% mutate(color=ifelse(as.numeric(year)<=2019,"gray",ifelse(as.numeric(year)==2020,"gold",ifelse(as.numeric(year)==2021,"orange",ifelse(as.numeric(year)==2022,"red",)))))


Sys.setlocale("LC_ALL","English")

#Yearly total all-cause deaths
ggplot(data=long %>% filter(year %notin% c("2009","2010")),aes(x = Date, y = ma(value),group=year))+
 geom_line(aes(color=year),size=1)+
 ylab("Daily mortality")+
 xlab("Date")+
 ggtitle("All cause daily mortality - Portugal")+
 scale_x_date(date_labels = "%b",minor_breaks = NULL,date_breaks = "1 month",expand=c(0,0))+
 scale_y_continuous(expand=c(0,0),breaks=seq(0,700,100))+
 coord_cartesian(ylim=c(200,750))+
 theme_light(base_size=18) +
 theme(,axis.text.x = element_text(angle = 90,hjust = 0.5,vjust=0.35),strip.background =element_rect(fill="gray"), strip.text = element_text(size=18,face="bold",colour = 'black'),legend.position = "none")+
 gghighlight::gghighlight(use_direct_label = FALSE,unhighlighted_colour = alpha("azure3", 0.3)) +
  geom_line(aes(color=year),size=1)+
 labs(caption="Daily mortality graphed as 7 day moving average. Data from: https://evm.min-saude.pt. Graph by Rodrigo Quiroga @rquiroga777 on Twitter.")+
 facet_wrap(~year)
fname<-paste0(today,"_total_mortality_portugal.png");
ggsave(fname, dpi = 600,type="cairo-png",width=15,height=10)

#Yearly total all-cause deaths + reported covid deaths
ggplot(data=long4 %>% filter(year %notin% c("2009","2010")),aes(x = Date, y = ma(value),group=year))+
 geom_line(aes(color=year),size=1)+
 scale_y_continuous(expand=c(0,0),breaks=seq(0,700,100))+
 ylab("Daily mortality")+
 xlab("Date")+
 ggtitle("All cause daily mortality and COVID reported deaths - Portugal")+
 scale_x_date(date_labels = "%b",minor_breaks = NULL,date_breaks = "1 month",expand=c(0,0))+
 coord_cartesian(ylim=c(0,750))+
theme_light(base_size=18) +
 theme(,axis.text.x = element_text(angle = 90,hjust = 0.5,vjust=0.35),strip.background =element_rect(fill="gray"), strip.text = element_text(size=18,face="bold",colour = 'black'),legend.position = "none")+
 gghighlight::gghighlight(use_direct_label = FALSE,unhighlighted_colour = alpha("azure3", 0.3)) +
 geom_line(aes(y=new_deaths_smoothed),size=1,color="black")+
 labs(caption="Daily mortality graphed as 7 day moving average. Covid death data graphed in black lines, obtained from OurWorldInData.\nData from: https://evm.min-saude.pt. Graph by Rodrigo Quiroga @rquiroga777 on Twitter.")+
 facet_wrap(~year)
fname<-paste0(today,"_total_mortality_covid_deaths_portugal.png");
ggsave(fname, dpi = 600,type="cairo-png",width=15,height=10)

#Prepare data for building model
long3<-long
long3$dia<-rep(seq(1,365,1),length(levels(long$year)))
long3$diatotal<-seq(1,nrow(long3),1)
#Subset for training
long3tr<-long3 %>% filter(year %in% seq(2009,2019,1))

#BUILD MODEL WITH A ROBUST REGRESSION WITH YEARLY AND SEASONAL COMPONENTS (SINE+COS+2nd harmonic)
model_sine2_yearly<-rlm(value ~ sin(2*pi*dia/365) + cos(2*pi*dia/365) + sin(4*pi*dia/365) + cos(4*pi*dia/365) + diatotal, data=long3tr, psi = psi.huber)

#Predict expected daily deaths for the whole dataset
predict.model_sine2_yearly<-predict(model_sine2_yearly, newdata=long3, type='response')
todo.df<-data.frame(pred2=predict.model_sine2_yearly,diatotal=long3$diatotal,Date=long3$Date,year=long3$year,deaths=long3$value,prom=long2$value,covid_deaths=long4$new_deaths_smoothed) %>% mutate(excess2=deaths-pred2)

#EXPECTED AND ACTUAL DAILY MORTALITY FOR EACH YEAR
ggplot(data=todo.df %>% filter(year %notin% c("2009","2010")),aes(x = Date, y = pred2,group=year))+
 geom_line(aes(color=year),size=1)+
 #geom_line(data=dataPOR2,x=Date,y=prom,inherit.aes = FALSE)
 ylab("Daily deaths")+
 coord_cartesian(ylim=c(200,500))+
 scale_y_continuous(expand = c(0,0))+
 xlab("Date")+
 ggtitle("Expected (modelled) daily all-cause mortality compared to daily actual deaths- Portugal")+
 scale_x_date(date_labels = "%b",minor_breaks = NULL,date_breaks = "1 month",expand=c(0,0))+
 theme_light(base_size=18) +
 theme(,axis.text.x = element_text(angle = 90,hjust = 0.5,vjust=0.35),strip.background =element_rect(fill="gray"), strip.text = element_text(size=18,face="bold",colour = 'black'),legend.position = "none")+
 gghighlight::gghighlight(use_direct_label = FALSE,n = 1,unhighlighted_colour = alpha("azure3", 0.3)) +
 geom_point(aes(y=deaths,color=year),size=0.8)+
 labs(caption="Daily mortality predicted with robust regression using yearly and seasonal components, trained on data 2009-2019.\nActual daily mortality plotted as dots, expected mortality as lines. Data from: https://evm.min-saude.pt. Graph by Rodrigo Quiroga @rquiroga777 on Twitter.")+
 facet_wrap(~year)
fname<-paste0(today,"_model_trained_portugal.png");
ggsave(fname, dpi = 600,type="cairo-png",width=15,height=10)


#PLOT EXCESS DEATHS
ggplot(data=todo.df %>% filter(year %notin% c("2009","2010")),aes(x = Date, y = ma(excess2),group=year))+
 geom_line(aes(color=ma(excess2)),size=1)+
 scale_color_gradientn(colours = pal, limits = c(0,50),oob=squish) +
 ylab("Daily excess mortality")+
 xlab("Date")+
 ggtitle("Excess mortality - Portugal")+
 scale_x_date(date_labels = "%b",minor_breaks = NULL,date_breaks = "1 month",expand=c(0,0))+
 scale_y_continuous(breaks=seq(-100,300,100),minor_breaks = NULL)+
 theme_light(base_size=18) +
 theme(strip.background =element_rect(fill="gray"), strip.text = element_text(size=18,face="bold",colour = 'black'),legend.position = "none",axis.text.x = element_text(angle = 90,hjust = 0.5,vjust=0.35))+
 gghighlight::gghighlight(use_direct_label = FALSE,unhighlighted_colour = alpha("azure3", 0.3)) +
 facet_wrap(~year)+
 geom_hline(yintercept = 0,colour="black")+
 labs(caption="Moving average (7 days) all-cause excess mortality, with respect to a baseline mortality model that adjusts\nfor yearly and seasonal components trained on 2009-2019 data.Data from: https://evm.min-saude.pt. Graph by Rodrigo Quiroga @rquiroga777 on Twitter.")
fname<-paste0(today,"_mortalidad_portugal_sinepred_ENG.png");
ggsave(fname, dpi = 600,type="cairo-png",width=15,height=10)



#ONLY 2021 + 2022
ggplot(data=todo.df %>% filter(year %in% c("2021","2022")),aes(x = Date, y = ma(excess2),group=year))+
 geom_line(aes(color=ma(excess2)),size=1)+
 scale_color_gradientn(colours = pal, limits = c(0,50),oob=squish) +
 ylab("Daily excess mortality")+
 xlab("Date")+
 ggtitle("Excess mortality - Portugal")+
 scale_x_date(date_labels = "%b",minor_breaks = NULL,date_breaks = "1 month",expand=c(0,0))+
 scale_y_continuous(breaks=seq(-100,300,100),minor_breaks = NULL)+
 theme_light(base_size=18) +
 theme(strip.background =element_rect(fill="gray"), strip.text = element_text(size=22,face="bold",colour = 'black'),legend.position = "none",axis.text.x = element_text(angle = 90,hjust = 0.5,vjust=0.35))+
 gghighlight::gghighlight(use_direct_label = FALSE,unhighlighted_colour = alpha("azure3", 0.3)) +
 facet_wrap(~year)+
 geom_hline(yintercept = 0,colour="black")+
 labs(caption="Moving average (7 days) all-cause excess mortality, with respect to a baseline mortality model that adjusts for yearly and seasonal\ncomponents trained on 2009-2019 data.Data from: https://evm.min-saude.pt. Graph by Rodrigo Quiroga @rquiroga777 on Twitter.")
fname<-paste0(today,"_mortalidad_portugal_sinepred_21_22_ENG.png");
ggsave(fname, dpi = 600,type="cairo-png",width=15,height=10)

#2020, 2021, 2022 + covid deaths
ggplot(data=todo.df %>% filter(year %in% c("2020","2021","2022")),aes(x = Date, y = ma(excess2),group=year))+
 geom_line(aes(color=ma(excess2)),size=1)+
 scale_color_gradientn(colours = pal, limits = c(0,50),oob=squish) +
 ylab("Daily excess mortality")+
 xlab("Date")+
 ggtitle("Excess mortality and reported COVID deaths - Portugal")+
 scale_x_date(date_labels = "%b",minor_breaks = NULL,date_breaks = "1 month",expand=c(0,0))+
 scale_y_continuous(breaks=seq(-100,300,100),minor_breaks = NULL)+
 theme_light(base_size=18) +
 theme(strip.background =element_rect(fill="gray"), strip.text = element_text(size=22,face="bold",colour = 'black'),legend.position = "none",axis.text.x = element_text(angle = 90,hjust = 0.5,vjust=0.35))+
 gghighlight::gghighlight(use_direct_label = FALSE,unhighlighted_colour = alpha("azure3", 0.3)) +
 facet_wrap(~year)+
 #geom_hline(yintercept = 0,colour="black")+
 geom_line(aes(y=covid_deaths),size=1,color="black")+
 labs(caption="Moving average (7d) all-cause excess mortality in color, based on model with yearly and seasonal components\ntrained on 2009-2019 data. Reported covid deaths in black (smoothed deaths from OWID data). Data from: https://evm.min-saude.pt. Graph by Rodrigo Quiroga.")
fname<-paste0(today,"_mortalidad_portugal_sinepred_20_21_22_covid_ENG.png");
ggsave(fname, dpi = 600,type="cairo-png",width=15,height=10)


#PLOT EXCESS %
todo.df<-todo.df %>% mutate(excess2_por=excess2/pred2) %>% mutate(excess2_por_m7=ma(excess2_por))
ggplot(data=todo.df %>% filter(year %notin% c("2009","2010")),aes(x = Date, y = excess2_por_m7,group=year))+
 geom_line(aes(color=ma(excess2)),size=1.1)+
 scale_color_gradientn(colours = pal, limits = c(0,50),oob=squish) +
 ylab("Daily excess mortality")+
 xlab("Date")+
 ggtitle("Excess mortality - Portugal")+
 scale_x_date(date_labels = "%b",minor_breaks = NULL,date_breaks = "1 month",expand=c(0,0))+
 scale_y_continuous(labels=scales::percent_format(accuracy = 5L),breaks=c(-.1,seq(0,1,.10)),minor_breaks = NULL)+
 coord_cartesian(ylim=c(-.1,.9))+
 theme_light(base_size=18) +
 theme(strip.background =element_rect(fill="gray"), strip.text = element_text(size=18,face="bold",colour = 'black'),legend.position = "none",axis.text.x = element_text(angle = 90,hjust = 0.5,vjust=0.35))+
 gghighlight::gghighlight(use_direct_label = FALSE,unhighlighted_colour = alpha("azure3", 0.3)) +
 facet_wrap(~year)+
 geom_hline(yintercept = 0,colour="black")+
 labs(caption="Moving average (7 days) all-cause excess mortality, with respect to a baseline mortality model that adjusts for yearly and seasonal\ncomponents trained on 2009-2019 data.Data from: https://evm.min-saude.pt. Graph by Rodrigo Quiroga @rquiroga777 on Twitter.")
fname<-paste0(today,"_mortalidad_portugal_sinepred_porc_ENG.png");
ggsave(fname, dpi = 600,type="cairo-png",width=15,height=10)

#PLOT EXCESS % 2021-2022
ggplot(data=todo.df %>% filter(year %in% c("2021","2022")),aes(x = Date, y = excess2_por_m7,group=year))+
 geom_line(aes(color=ma(excess2)),size=1.3)+
 scale_color_gradientn(colours = pal, limits = c(0,50),oob=squish) +
 ylab("Daily excess mortality")+
 xlab("Date")+
 ggtitle("All cause daily excess mortality - Portugal: With the current BA.5 wave, Portugal is at the maximum\nexcess mortality levels since the vaccination campaign began")+
 scale_x_date(date_labels = "%b",minor_breaks = NULL,date_breaks = "1 month",expand=c(0,0))+
 scale_y_continuous(labels=scales::percent_format(accuracy = 5L),breaks=c(-.1,seq(0,1,.10)),minor_breaks = NULL)+
 coord_cartesian(ylim=c(-.1,.9))+
 theme_light(base_size=18) +
 theme(strip.background =element_rect(fill="gray"), strip.text = element_text(size=22,face="bold",colour = 'black'),legend.position = "none",axis.text.x = element_text(angle = 90,hjust = 0.5,vjust=0.35))+
 gghighlight::gghighlight(use_direct_label = FALSE,unhighlighted_colour = alpha("azure3", 0.3)) +
 facet_wrap(~year)+
 geom_hline(yintercept = 0,colour="black")+
 labs(caption="Moving average (7 days) all-cause excess mortality, with respect to a baseline mortality model that adjusts for yearly and seasonal\ncomponents trained on 2009-2019 data. Data obtained from: https://evm.min-saude.pt. Graph by Rodrigo Quiroga @rquiroga777 on Twitter.")
fname<-paste0(today,"_mortalidad_portugal_sinepred_porc_21_22_ENG.png");
ggsave(fname, dpi = 600,type="cairo-png",width=15,height=10)
