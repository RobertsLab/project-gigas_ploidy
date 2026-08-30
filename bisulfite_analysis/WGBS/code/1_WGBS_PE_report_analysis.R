# Title: WGBS analysis.R
# Author: Matthew George; mattgeorgephd@gmail.com
# Date: 03/30/2021

## clear
rm(list=ls())

## Load R packages
library(readxl)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

## Grab the WD from the file location and set it
library(rstudioapi)
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path ))
print( getwd() )

#######################################################################################################################################

#Import data
data  <- read_excel("data/PE_reports/percent_methylation_summary.xlsx", sheet = "raw", col_names = TRUE)
data2 <- read_excel("data/PE_reports/percent_methylation_summary.xlsx", sheet = "summary", col_names = TRUE)

####################################################################################################################
#plot CpG x ploidy (all have undergone desiccation stress)

bp1 <- ggplot(data) +
      geom_boxplot(aes(x=ploidy, y=CpG, fill=ploidy),color="black") +
      scale_fill_manual(values=c("royalblue1", "orangered1")) +
      # scale_y_continuous(breaks = seq(0, 4, 0.5), limits = c(0, 3)) +
      # scale_x_discrete(labels=c("diploid_control","diploid_head","triploid_control","triploid_heat")) +
      xlab("Treatment") + 
      ylab("%mCpG") +         
      theme(line              = element_line(size=1.5,color="black"),
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
ggsave("figures/[boxplot]mCpG_ploidy.tiff",
       plot   = bp1,
       dpi    = 1200,
       device = "tiff",
       width  = 5,
       height = 5,
       units  = "in")

library(Johnson)
library(agricolae)
library(nlme)
library(multcomp)

## Transform data
x = data$CpG
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

data$CpG <- x


# Statistical testings
summary(aov(CpG ~ ploidy, data = data))

# tx <- with(data, interaction(MN,DMA))
# amod <- aov(kPa ~ tx, data=data)
# HSD.test(amod, "tx", group=TRUE, console=TRUE)

####################################################################################################################
#plot CpG x ploid x heat_shock (45C) following desiccation stress

p2 <- ggplot(data2,aes(x=heat_shock, y=CpG, group=ploidy,color=ploidy)) +
  geom_errorbar(mapping=aes(x=heat_shock,ymin=CpG-se,ymax=CpG+se),width=0.2, size=1,color="black") +
  geom_point(size=5) +
  geom_line(size=1.5) + 
  scale_fill_manual(values=c("royalblue1", "orangered1")) +
  # scale_y_continuous(breaks = seq(0, 4, 0.5), limits = c(0, 3)) +
  # scale_x_discrete(labels=c("diploid_control","diploid_head","triploid_control","triploid_heat")) +
  xlab("Treatment") + 
  ylab("%mCpG") +         
  theme(line              = element_line(size=1.5,color="black"),
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
p2


## Save figure in figure folder
ggsave("figures/[boxplot]mCpG_ploidy_heatshock.tiff",
       plot   = p2,
       dpi    = 1200,
       device = "tiff",
       width  = 5,
       height = 5,
       units  = "in")

library(Johnson)
library(agricolae)
library(nlme)
library(multcomp)

## Transform data
x = data$CpG
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

data$CpG <- x


# Statistical testings
summary(aov(CpG ~ ploidy + heat_shock, data = data))

# tx <- with(data, interaction(MN,DMA))
# amod <- aov(kPa ~ tx, data=data)
# HSD.test(amod, "tx", group=TRUE, console=TRUE)

