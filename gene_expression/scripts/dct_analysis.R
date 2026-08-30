#
# UNCOMMENT the lines below if you do have the packages already installed
#
install.packages("ggplot2")
install.packages("plyr")
install.packages("splitstackshape")

#Necessary Packages to manipulate data and plot values. 
require(plyr)
require(ggplot2)
require(splitstackshape)

#Read in  Ct value table
dCt<-read.csv("data/qpcr_ct_values/qpcr_data_consolidated.csv", header=T)

#Split SAMPLE_ID column to create columns for Ploidyulation, treatment, and sample number
dCt<-cSplit(dCt,"Sample", sep= "_", drop=F)

#rename columns appropriately
dCt<-rename(dCt,replace=c("Sample_1"="Ploidy","Sample_2"="Desiccation","Sample_3"="HeatShock","Sample_4"="SampleNum"))

#calculate normalized expression of target gene Ct relative to actin Ct using: 2^-(delta Ct)
dCt$HSC70<-2^-(dCt$HSC70-dCt$Actin)
dCt$DNMT1<-2^-(dCt$DNMT1-dCt$Actin)
dCt$MBD2<-2^-(dCt$MBD2-dCt$Actin)
dCt$MeCP2<-2^-(dCt$MeCP2-dCt$Actin)
dCt$HIF1A<-2^-(dCt$HIF1A-dCt$Actin)
dCt$HATHaP2<-2^-(dCt$HATHaP2-dCt$Actin)
dCt$HAT<-2^-(dCt$HAT-dCt$Actin)
dCt$HSP90<-2^-(dCt$HSP90-dCt$Actin)


#log transform the data to develop normality in data
dCt$HSC70log<-log(dCt$HSC70)
dCt$DNMT1log<-log(dCt$DNMT1)
dCt$MeCP2log<-log(dCt$MeCP2)
dCt$HIF1Alog<-log(dCt$HIF1A)
dCt$HATHaP2log<-log(dCt$HATHaP2)
dCt$HATlog<-log(dCt$HAT)
dCt$HSP90log<-log(dCt$HSP90)
dCt$MBD2log<-log(dCt$MBD2)

#Run ANOVA's on all log transformed data as well as Tukey's Honestly Significant Difference post hoc test
HSC70<-aov(HSC70log~Ploidy+Desiccation+Ploidy:Desiccation, data=dCt)
HSC70
TukeyHSD(HSC70)
summary(HSC70)

DNMT1<-aov(DNMT1log~Ploidy+Desiccation+Ploidy:Desiccation, data=dCt)
DNMT1
TukeyHSD(DNMT1)
summary(DNMT1)

MeCP2<-aov(MeCP2log~Ploidy+Desiccation+Ploidy:Desiccation, data=dCt)
MeCP2
TukeyHSD(MeCP2)
summary(MeCP2)

HIF1A<-aov(HIF1Alog~Ploidy+Desiccation+Ploidy:Desiccation, data=dCt)
HIF1A
TukeyHSD(HIF1A)
summary(HIF1A)

HATHaP2<-aov(HATHaP2log~Ploidy+Desiccation+Ploidy:Desiccation, data=dCt)
HATHaP2
TukeyHSD(HATHaP2)
summary(HATHaP2)

HAT<-aov(HATlog~Ploidy+Desiccation+Ploidy:Desiccation, data=dCt)
HAT
TukeyHSD(HAT)
summary(HAT)

HSP90<-aov(HSP90log~Ploidy+Desiccation+Ploidy:Desiccation, data=dCt)
HSP90
TukeyHSD(HSP90)
summary(HSP90)



#graph all normalized Ct values to produce boxplots to visualize data

p1 <- ggplot(data=dCt) + 
  geom_boxplot(aes(x=Desiccation, y=HSC70,fill=Ploidy)) + theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Diploid","Triploid"))+
  guides(fill=guide_legend(title="Ploidy"))+
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=20), axis.title.y=element_text(size=20),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.023))+scale_x_discrete(labels=c("Desiccated + Elevated Temp.","Control"))+
  labs(x="Treatment", y=expression(paste("HSC70 Expression (",Delta,"Ct)")))

p1

p1 <- ggplot(data=dCt)+
      geom_boxplot(aes(x=Desiccation, y=DNMT1,fill=Ploidy)) +
      # scale_y_continuous(breaks = seq(0, 0.03, 0.01), limits = c(0, 0.035)) +
      scale_fill_manual(values=c("royalblue1", "orangered1")) +
      # scale_x_discrete(labels=c("Desiccated + Elevated Temp.","Control"))+
      # labs(x="Treatment", y=expression(paste("HSP90 Expression (",Delta,"Ct)"))) +
      theme(line              = element_line(size=1.5),
            rect              = element_rect(size=1.5),
            text              = element_text(size=14,color="black"),
            panel.background  = element_blank(),
            panel.grid.major  = element_blank(), 
            panel.grid.minor  = element_blank(),
            axis.text.x       = element_text(size=16,color="black"),
            axis.text.y       = element_text(size=16,color="black"),
            axis.title.x      = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y      = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.ticks.x      = element_line(color="black"),
            axis.ticks.y      = element_line(color="black"),
            # axis.line         = element_line(color = "black", size = 0.1),
            panel.border      = element_rect(color = "black", fill=NA, size=1.5),
            legend.key        = element_blank(),
            legend.title      = element_blank()) # removes background of legend bullets

p1

ggsave("DNMT1_ploidy.tiff",
       plot   = p1,
       dpi    = 600,
       device = "tiff",
       width  = 5,
       height = 5,
       units  = "in")

ggplot(data=dCt)+geom_boxplot(aes(x=Desiccation, y=MBD2,fill=Ploidy))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Diploid","Triploid"))+
  guides(fill=guide_legend(title="Ploidy"))+
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=20), axis.title.y=element_text(size=20),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.0005))+scale_x_discrete(labels=c("Desiccated + Elevated Temp.","Control"))+
  labs(x="Treatment", y=expression(paste("MBD2 Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Desiccation, y=MBD2,fill=Ploidy))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Diploid","Triploid"))+
  guides(fill=guide_legend(title="Ploidy"))+
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=20), axis.title.y=element_text(size=20),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.0005))+scale_x_discrete(labels=c("Desiccated + Elevated Temp.","Control"))+
  labs(x="Treatment", y=expression(paste("MBD2 Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Desiccation, y=MeCP2,fill=Ploidy))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Diploid","Triploid"))+
  guides(fill=guide_legend(title="Ploidy"))+
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=20), axis.title.y=element_text(size=20),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
<<<<<<< HEAD
  ylim(c(0,0.002))+scale_x_discrete(labels=c("Dessicated + Elevated Temp.","Control"))+
=======
  ylim(c(0,0.00075))+scale_x_discrete(labels=c("Desiccated + Elevated Temp.","Control"))+
>>>>>>> 8b8bf5e62dcbc30895d3a87125404ede2d66328d
  labs(x="Treatment", y=expression(paste("MeCP2 Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Desiccation, y=HIF1A,fill=Ploidy))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Diploid","Triploid"))+
  guides(fill=guide_legend(title="Ploidy"))+
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=20), axis.title.y=element_text(size=20),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.00075))+scale_x_discrete(labels=c("Desiccated + Elevated Temp.","Control"))+
  labs(x="Treatment", y=expression(paste("HIF1A Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Desiccation, y=HATHaP2,fill=Ploidy))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Diploid","Triploid"))+
  guides(fill=guide_legend(title="Ploidy"))+
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=20), axis.title.y=element_text(size=20),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.005))+scale_x_discrete(labels=c("Desiccated + Elevated Temp.","Control"))+
  labs(x="Treatment", y=expression(paste("HATHaP2 Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Desiccation, y=HAT,fill=Ploidy))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Diploid","Triploid"))+
  guides(fill=guide_legend(title="Ploidy"))+
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=20), axis.title.y=element_text(size=20),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.00075))+scale_x_discrete(labels=c("Desiccated + Elevated Temp.","Control"))+
  labs(x="Treatment", y=expression(paste("HAT Expression (",Delta,"Ct)")))

p1 <- ggplot(data=dCt)+
      geom_boxplot(aes(x=Desiccation,y=HSP90,fill=Ploidy))+
      scale_y_continuous(breaks = seq(0, 0.03, 0.01), limits = c(0, 0.035)) +
      scale_fill_manual(values=c("royalblue1", "orangered1")) +
        # scale_x_discrete(labels=c("Desiccated + Elevated Temp.","Control"))+
        # labs(x="Treatment", y=expression(paste("HSP90 Expression (",Delta,"Ct)"))) +
        theme(line              = element_line(size=1.5),
                    rect              = element_rect(size=1.5),
                    text              = element_text(size=14,color="black"),
                    panel.background  = element_blank(),
                    panel.grid.major  = element_blank(), 
                    panel.grid.minor  = element_blank(),
                    axis.text.x       = element_text(size=16,color="black"),
                    axis.text.y       = element_text(size=16,color="black"),
                    axis.title.x      = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                    axis.title.y      = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                    axis.ticks.x      = element_line(color="black"),
                    axis.ticks.y      = element_line(color="black"),
                    # axis.line         = element_line(color = "black", size = 0.1),
                    panel.border      = element_rect(color = "black", fill=NA, size=1.5),
                    legend.key        = element_blank(),
                    legend.title      = element_blank()) # removes background of legend bullets

p1

ggsave("HSP90_ploidy.tiff",
       plot   = p1,
       dpi    = 600,
       device = "tiff",
       width  = 5,
       height = 5,
       units  = "in")
