# Title: GlobalDNAMeth_Polyploids.R
# Author: Matthew George; mattgeorgephd@gmail.com
# Original Author: Shelly Trigg
# Date: 02/01/2021

## clear
rm(list=ls())

## Load R packages
library(readxl)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(broom)

## Grab the WD from the file location and set it
library(rstudioapi)
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path ))
print( getwd() )

#######################################################################################################################################

#Import data
setwd('docs/')
print( getwd() )

data1 <- read_excel("20190206_FirstReadingSTRIGG.xls", sheet = "plate list", col_names = TRUE)
data2 <- read_excel("20190206_2ndReadingSTRIGG.xls", sheet = "plate list", col_names = TRUE)
data3 <- read_excel("20190206_3rdReadingSTRIGG.xls", sheet = "plate list", col_names = TRUE)


Feb2019 <- merge(data1[,c(3,6)],data2[,c(3,6)], by = "Well")
Feb2019 <- merge(Feb2019,data3[,c(3,6)], by = "Well")
colnames(Feb2019) <- c("Well", "Abs1", "Abs2", "Abs3")
Feb2019$mean <- apply(Feb2019[,2:4],1, mean)
Feb2019$sd <- apply(Feb2019[,2:4],1, sd)
max(Feb2019$sd)
#[1] 0.01871462

#SDs are very small between readings so go with data1 readings

#Calculate standard curve

#first read in sample names from plate map list
Feb_plateMap <- read_excel("20190131_Ronit'sSamplesforGlobalDNAMeth.xlsx", sheet = "PlateMap_list", col_names = TRUE)

#make a list of samples to not be considered
rmv_wells <- c("no_sample", "contam", "no_ab")

Feb_data <- merge(Feb_plateMap[which(!(Feb_plateMap$Sample %in% rmv_wells)),], data1[,c(3,6)])
head(Feb_data)

#simplify absorbance column name
colnames(Feb_data)[4] <- "Absorbance"

#Calculate Average absorbance and SD for each sample
Feb_data_avg <- Feb_data %>% group_by(Sample) %>% summarize(avgAbs=mean(Absorbance))
Feb_data_sd <- Feb_data  %>% group_by(Sample) %>% summarize(SD=sd(Absorbance))
Feb_data_avg

Feb_data_avg <- merge(Feb_data_avg,Feb_data_sd, by = "Sample")

Feb_data_avg <- unique(merge(Feb_data[,c("Type","Sample")], Feb_data_avg, by = "Sample"))

#make new data frame with only STDs
curve <- Feb_data_avg[grep("STD", Feb_data_avg$Type),]
print(curve)

#make new column for %meth
curve$perc_meth <- c(0.1,0.2,0.5,1,5,0,10)
#order data by perc_meth
curve <- curve[order(curve$perc_meth),]
print(curve)

#plot curve
ggplot(curve, aes(perc_meth, avgAbs)) + geom_point() 

##last two points are plateauing so remove them
ggplot(curve[which(curve$perc_meth < 5.0),], aes(perc_meth, avgAbs)) + geom_point() + geom_smooth(method = "lm", se = FALSE)

#find regression line
fit <- lm(curve[which(curve$perc_meth < 5.0),"avgAbs"] ~ curve[which(curve$perc_meth < 5.0),"perc_meth"] )
#find the R-squared
rsquared <- summary(fit)$r.squared
rsquared
#[1] 0.9978712

#find the slope (https://www.cyclismo.org/tutorial/R/linearLeastSquares.html)
slope <- fit$coefficients[[2]]
slope

#plot regression line with equation and r-squared (https://www.listendata.com/2016/07/add-linear-regression-equation-and.html)
linear = function(k) {
  z <- list(xx = format(coef(k)[1], digits = 5),
            yy = format(abs(coef(k)[2]), digits = 5),
            r2 = format(summary(k)$r.squared, digits = 5));
  if (coef(k)[1] >= 0)  {
    eq <- substitute(italic(Y) == yy %*% italic(X) + xx*","~~italic(r)^2~"="~r2,z)
  } else {
    eq <- substitute(italic(Y) == yy %*% italic(X) - xx*","~~italic(r)^2~"="~r2,z)  
  }
  as.character(as.expression(eq));              
}

x <- curve[which(curve$perc_meth < 5.0),"perc_meth"]
y <- curve[which(curve$perc_meth < 5.0),"avgAbs"]
fo = y ~ x
linplot <- ggplot(data = curve[which(curve$perc_meth < 5.0),], aes(x = perc_meth, y = avgAbs)) + geom_smooth(method = "lm", se=FALSE, color="black", formula = fo) +  geom_point() + theme_bw()
linplot1 = linplot + annotate("text", x = 0.3, y = 0.7, label = linear(lm(fo)), colour="black", size = 5, parse=TRUE)
linplot1

#calculate perc_meth for samples using equation in EpiGentek MethylFlash kit
#5-mC% = (sample OD - NC OD) / (slope x input sample DNA in ng) * 100%

NC_OD <- Feb_data_avg[grep("5 ul NC", Feb_data_avg$Sample),"avgAbs"]
slopexDNAamt <- slope * 25

Feb_data_avg$avgAbs_NC <- Feb_data_avg$avgAbs - NC_OD
Feb_data_avg$perc_meth <- Feb_data_avg$avgAbs_NC/slopexDNAamt *100

#make new data frame with only sample info
Feb_data_avg_Samples <- Feb_data_avg[which(Feb_data_avg$Type == "Sample"),]

#Add sample description to data frame
for (i in 1:nrow(Feb_data_avg_Samples)){
  if(substr(Feb_data_avg_Samples$Sample[i],1,1)=="D"){
    Feb_data_avg_Samples$ploidy[i] <- "diploid"
  }
  if(substr(Feb_data_avg_Samples$Sample[i],1,1)=="T"){
    Feb_data_avg_Samples$ploidy[i] <- "triploid"
  }
  if(as.numeric(substr(Feb_data_avg_Samples$Sample[i],2,3)) > 10){
    Feb_data_avg_Samples$treatment[i] <- "heat_stress"
  }
  if(as.numeric(substr(Feb_data_avg_Samples$Sample[i],2,3)) < 10){
    Feb_data_avg_Samples$treatment[i] <- "control"
  }
}

Feb_data_avg_Samples$condition <- paste(Feb_data_avg_Samples$ploidy, Feb_data_avg_Samples$treatment, sep = "_")

####################################################################################################################
#plot absorbances x experimental group

head(Feb_data_avg_Samples)

Feb_data_avg_Samples$condition <- factor(Feb_data_avg_Samples$condition, levels=c("diploid_control","triploid_control","diploid_heat_stress","triploid_heat_stress"), ordered=TRUE)


bp1 <- ggplot(Feb_data_avg_Samples, aes(x = condition, y = perc_meth, fill=ploidy)) +  
            geom_boxplot(colour = "black", size = 1.5,outlier.colour="black", outlier.shape=16,
                         outlier.size=2, notch=FALSE) +
            scale_fill_manual(values = c("royalblue1", "orangered1", "darkgreen", "palegreen")) +
            scale_y_continuous(breaks = seq(0, 4, 0.5), limits = c(0, 3)) +
            # scale_x_discrete(labels=c("diploid_control","diploid_head","triploid_control","triploid_heat")) +
            xlab("Sample") + 
            ylab("% 5-mC") +         
            theme(line              = element_line(size=1.5),
                  rect              = element_rect(size=1.5),
                  text              = element_text(size=22,color="black"),
                  panel.background  = element_blank(),
                  panel.grid.major  = element_blank(),
                  axis.ticks        = element_blank(),
                  axis.text         = element_text(size=22,color="black"),
                  panel.grid.minor  = element_blank(),
                  axis.title.x      = element_blank(),
                  axis.line         = element_blank(),
                  panel.border      = element_rect(color = "black", fill=NA, size=2),
                  legend.position   = 'none',
                  legend.key        = element_rect(fill = 'white'))
bp1


## Save figure in figure folder

getwd()
setwd('..')
setwd('figures/')
getwd()

ggsave("[boxplot]5mC_comparison.tiff",
       plot   = bp1,
       dpi    = 1200,
       device = "tiff",
       width  = 7,
       height = 5,
       units  = "in")



head(Feb_data_avg_Samples)

aov_2way <- aov(perc_meth ~ ploidy + treatment + ploidy:treatment, data = Feb_data_avg_Samples)
summary(aov_2way)
aov_2way_model_summary <- glance(aov_2way)
aov_2way_model_summary

#p-value is significant at 0.05 so do a Tukeys test to see which effect is significant
tuk <- TukeyHSD(aov(perc_meth ~ ploidy + treatment + ploidy:treatment, data = Feb_data_avg_Samples))
tuk


library(Johnson)
library(agricolae)
library(nlme)
library(multcomp)

## Transform data
x = Feb_data_avg_Samples$perc_meth
qqnorm(x) # check for normality
qqline(x) # Draw the line
result <- shapiro.test(x) # p-value fail = good, don't need transformation
print(result$p.value)
if(result$p.value<0.05)     {
  x_johnson <- RE.Johnson(x) # transform
  x_transformed = x_johnson$transformed
  qqnorm(x_transformed) # check linearity of tranformed data
  qqline(x_transformed)
  print(shapiro.test(x_transformed))
  x <- x_transformed  
  print("transformed!",quote=FALSE)}
shapiro.test(x)

Feb_data_avg_Samples$perc_meth <- x


# Statistical testings
summary(aov(perc_meth ~ ploidy*treatment, data = Feb_data_avg_Samples))

tx <- with(data, interaction(MN,DMA))
amod <- aov(kPa ~ tx, data=data)
HSD.test(amod, "tx", group=TRUE, console=TRUE)



#based on the Tukey's test results, the diploid heat stress %5-mC methylation is significantly different from 1. the diploid control 2. the triploid control and 3. the triploid heatstress. The triploid heat stress %5-mC is not different than the diploid or triploid controls. So in summary, diploid ctenidia shows a significant increase in global 5-mC methylation in response to heatstress while triploid ctenidia does not. 

