#
#UNCOMMENT the lines below if you do have the packages already installed
#
install.packages("ggplot2")
install.packages("plyr")
install.packages("splitstackshape")
install.packages("tidyr")

#Necessary Packages to manipulate data and plot values. 
require(plyr)
require(ggplot2)
require(splitstackshape)
require(tidyr)

#Read in  Ct value table
dCt<-read.csv("data/qpcr_ct_values/qpcr_data_consolidated.csv", header=T)

#Split SAMPLE_ID column to create columns for Ploidyulation, treatment, and sample number
dCt<-cSplit(dCt,"Sample", sep= "_", drop=F)

#rename columns appropriately
dCt<-rename(dCt,replace=c("Sample_1"="Ploidy","Sample_2"="Desiccation","Sample_3"="HeatShock","Sample_4"="SampleNum"))

#change NA to 45
dCt[is.na(dCt)] <- 45

#calculate normalized expression of target gene Ct relative to actin Ct using: 2^-(delta Ct)
dCt$HSC70<-2^-(dCt$HSC70-dCt$Actin)
dCt$DNMT1<-2^-(dCt$DNMT1-dCt$Actin)
dCt$MBD2<-2^-(dCt$MBD2-dCt$Actin)
dCt$MeCP2<-2^-(dCt$MeCP2-dCt$Actin)
dCt$HIF1A<-2^-(dCt$HIF1A-dCt$Actin)
dCt$HATHaP2<-2^-(dCt$HATHaP2-dCt$Actin)
dCt$HAT<-2^-(dCt$HAT-dCt$Actin)
dCt$HSP90<-2^-(dCt$HSP90-dCt$Actin)
dCt$SOD<-2^-(dCt$SOD-dCt$Actin)
dCt$ATPsynthetase<-2^-(dCt$ATPsynthetase-dCt$Actin)



#log transform the data to develop normality in data
dCt$HSC70log<-log(dCt$HSC70)
dCt$DNMT1log<-log(dCt$DNMT1)
dCt$MeCP2log<-log(dCt$MeCP2)
dCt$HIF1Alog<-log(dCt$HIF1A)
dCt$HATHaP2log<-log(dCt$HATHaP2)
dCt$HATlog<-log(dCt$HAT)
dCt$HSP90log<-log(dCt$HSP90)
dCt$MBD2log<-log(dCt$MBD2)
dCt$SODlog<-log(dCt$SOD)
dCt$ATPsynthetaselog<-log(dCt$ATPsynthetase)
dCt$COX1log<-log(dCt$COX1)

#############################
### Facet Wrap with ggplot ##
#############################

#First make new data frame with normalized expression of target gene Ct relative to actin Ct using: 2^-(delta Ct)
data <- dCt[,c(14:17, 3:13)]
#Reshape dataframe so that all columns with a transcript name are gathered into one column called "Transcript" and their corresponding delta_Ct values are gathered into one column called "delta_Ct"
STACKED_data <- gather(data,"Transcript", "delta_Ct",5:15)

#plot the data with facets
ggplot(STACKED_data)+geom_boxplot(aes(x=Desiccation, y=delta_Ct,fill=Ploidy))+theme_bw()+
  scale_fill_manual(values=c("#31C4ED", "#62D12B"),
                    labels=c("Diploid","Triploid"))+
  guides(fill=guide_legend(title="Ploidy"))+
  scale_x_discrete(labels=c("Control","Desiccation \n + Elevated Temp.")) +
  labs(x="Treatment", y=expression(paste("Expression (",Delta,"Ct)"))) + 
  facet_wrap(~Transcript, scale = "free")

#Repeat the above steps only for log transformed delta_Ct values
log_data <- dCt[,c(14:17, 18:28)]
STACKED_log_data <- gather(log_data,"Transcript", "log_delta_Ct",5:15)
STACKED_log_data$Transcript <- gsub("log","",STACKED_log_data$Transcript)

ggplot(STACKED_log_data)+geom_boxplot(aes(x=Desiccation, y=log_delta_Ct,fill=Ploidy))+theme_bw()+
  scale_fill_manual(values=c("#31C4ED", "#62D12B"),
                    labels=c("Diploid","Triploid"))+
  guides(fill=guide_legend(title="Ploidy"))+
  scale_x_discrete(labels=c("Control","Desiccation \n + Elevated Temp.")) +
  labs(x="Treatment", y=expression(paste("Expression (log",Delta,"Ct)"))) + 
  facet_wrap(~Transcript, scale = "free")



