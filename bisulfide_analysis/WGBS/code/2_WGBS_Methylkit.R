
# Title: WGBS analysis.R
# Author: Matthew George; mattgeorgephd@gmail.com
# Date: 03/2021

###################################################################################################
## Setup work space

rm(list=ls()) # clear all

## First time? - install these
# install.packages("readr")
# install.packages("tidyverse")
# install.packages("devtools")
# install.packages("vegan")
# install.packages("pheatmap")
# install.packages("gplots")
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install(version = "3.12")
# BiocManager::install(c("GenomicFeatures", "AnnotationDbi","methylKit"))
# browseVignettes("methylKit") #methylKit manual

# Load R packages
library(rstudioapi)
require(tidyverse)
require(devtools)
require(vegan)
require(pheatmap) 
require(gplots)
require(methylKit)
require(readr)

# Grab the WD from the file location and set it
home <- getActiveDocumentContext()$path
setwd(dirname(home)); getwd()
setwd('..')
home <- getwd()
###################################################################################################
## Load previous session data. Return here to save after big unites.

setwd("E:/project-gigas_ploidy/R_sessions")
# save.image("20210528_methylKit_DMLs_cov10.RData") #Save R Data in case R crashes
load("20210528_methylKit_DMLs_cov10.RData") #Load R Data

####################################################################################################
## Find files and generate list

# Put all .cov files into a list for analysis.
setwd("E:/project-gigas_ploidy/cov_files") # set working directory to location of .cov files
analysisFiles <- list("E:/project-gigas_ploidy/cov_files/zr3534_1_R1.fastp-trim.20201202.CpG_report.merged_CpG_evidence.cov",
                      "E:/project-gigas_ploidy/cov_files/zr3534_2_R1.fastp-trim.20201202.CpG_report.merged_CpG_evidence.cov",
                      "E:/project-gigas_ploidy/cov_files/zr3534_3_R1.fastp-trim.20201202.CpG_report.merged_CpG_evidence.cov",
                      "E:/project-gigas_ploidy/cov_files/zr3534_4_R1.fastp-trim.20201202.CpG_report.merged_CpG_evidence.cov",
                      "E:/project-gigas_ploidy/cov_files/zr3534_5_R1.fastp-trim.20201202.CpG_report.merged_CpG_evidence.cov",
                      "E:/project-gigas_ploidy/cov_files/zr3534_6_R1.fastp-trim.20201202.CpG_report.merged_CpG_evidence.cov",
                      "E:/project-gigas_ploidy/cov_files/zr3534_7_R1.fastp-trim.20201202.CpG_report.merged_CpG_evidence.cov",
                      "E:/project-gigas_ploidy/cov_files/zr3534_8_R1.fastp-trim.20201202.CpG_report.merged_CpG_evidence.cov",
                      "E:/project-gigas_ploidy/cov_files/zr3534_9_R1.fastp-trim.20201202.CpG_report.merged_CpG_evidence.cov",
                      "E:/project-gigas_ploidy/cov_files/zr3534_10_R1.fastp-trim.20201202.CpG_report.merged_CpG_evidence.cov")

setwd(dirname(home)); getwd() # reset working directory to home dir

####################################################################################################
## Run Methylkit

# Create data.frame with sample info
sampleMetadata <- data.frame("sampleID" = c("1", "2", "3", "4","5", "6", "7", "8", "9", "10"),
                             "treatment" = c(rep(0, times = 5),rep(1, times = 5))) # Specify which treatment the samples were from. All animals were subjected to desiccation. 0 = diploid, 1 = triploid 
head(sampleMetadata) #Confirm data.frame creation

# Use `methRead` to create a methylation object from the coverage files, and include sample ID and treatment information.
processedFiles <- methylKit::methRead(analysisFiles,
                                      sample.id = list("1", "2", "3", "4","5", "6", "7", "8", "9", "10"),
                                      assembly = "Roslin",
                                      treatment = sampleMetadata$treatment,
                                      pipeline = "bismarkCoverage",
                                      mincov = 2) #Process files. Treatment specified based on ploidy status. Use mincov = 2 to quickly process reads.

# Filter coverage information for minimum 5x coverage and 10x coverage 
# remove PCR duplicates by excluding data in the 99.9th percentile of coverage with hi.perc = 99.9. 
# Normalize coverage between samples to avoid over-sampling reads from one sample during statistical testing

processedFilteredFilesCov1 <- methylKit::filterByCoverage(processedFiles,
                                                           lo.count = 1, lo.perc = NULL,
                                                           high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)

processedFilteredFilesCov5 <- methylKit::filterByCoverage(processedFiles,
                                                          lo.count = 5, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)

processedFilteredFilesCov10 <- methylKit::filterByCoverage(processedFiles,
                                                          lo.count = 10, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)

####################################################################################################
## Create file name lists for plots

# Create data.frame for file names
nFiles       <- 10 #Count number of samples
fileNameCov5  <- data.frame("nameBase"  = rep("../percent-CpG-methylation", times = nFiles),
                            "nameBase2" = rep("../percent-CpG-coverage", times = nFiles),
                            "sample.ID" = c(1:10))
fileNameCov10 <- data.frame("nameBase"  = rep("../percent-CpG-methylation", times = nFiles),
                            "nameBase2" = rep("../percent-CpG-coverage", times = nFiles),
                            "sample.ID" = c(1:10))

fileNameCov5$actualFileName1 <- paste(fileName$nameBase, "-Filtered", "-5xCoverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG methylation plot
fileNameCov5$actualFileName2 <- paste(fileName$nameBase2, "-Filtered", "-5xCoverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG coverage plot

fileNameCov10$actualFileName1 <- paste(fileName$nameBase, "-Filtered", "-10xCoverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG methylation plot
fileNameCov10$actualFileName2 <- paste(fileName$nameBase2, "-Filtered", "-10xCoverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG coverage plot

##########################################################################################################
## Create plots

setwd("E:/project-gigas_ploidy/figures")
wd <- "E:/project-gigas_ploidy/figures"

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileNameCov5$actualFileName1[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getMethylationStats(processedFilteredFilesCov5[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device 
  } #Plot and save %CpG methylation information

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileNameCov5$actualFileName2[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getCoverageStats(processedFilteredFilesCov5[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileNameCov10$actualFileName1[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getMethylationStats(processedFilteredFilesCov10[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device 
} #Plot and save %CpG methylation information

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileNameCov10$actualFileName2[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getCoverageStats(processedFilteredFilesCov10[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information


##########################################################################################################
### Comparative analysis

### Combine all processed files into a single table. Use destrand = FALSE to not destrand. By default only bases with data in all samples will be kept

methylationInformationFilteredCov1 <- methylKit::unite(processedFilteredFilesCov1, 
                                                       destrand = FALSE,
                                                       mc.cores = 2) 

methylationInformationFilteredCov5 <- methylKit::unite(processedFilteredFilesCov5, 
                                                       destrand = FALSE,
                                                       mc.cores = 2) 

methylationInformationFilteredCov10 <- methylKit::unite(processedFilteredFilesCov10, 
                                                       destrand = FALSE,
                                                       mc.cores = 2) 

### Obtain Clustering Information
clusteringInformationFilteredCov5 <- methylKit::clusterSamples(methylationInformationFilteredCov5, dist = "correlation", method = "ward", plot = FALSE) #Save cluster information as a new object
clusteringInformationFilteredCov10 <- methylKit::clusterSamples(methylationInformationFilteredCov10, dist = "correlation", method = "ward", plot = FALSE) #Save cluster information as a new object


### Generate plots in output folder
# jpeg(filename = "Full-Sample-Pearson-Correlation-Plot-FilteredCov5.jpeg", height = 1000, width = 1000) #Save file with designated name
# methylKit::getCorrelation(methylationInformationFilteredCov5, plot = TRUE) #Understand correlation between methylation patterns in different samples
# dev.off()
# 
# jpeg(filename = "Full-Sample-CpG-Methylation-Clustering-FilteredCov5.jpeg", height = 1000, width = 1000) #Save file with designated name
# methylKit::clusterSamples(methylationInformationFilteredCov5, dist = "correlation", method = "ward", plot = TRUE) #Cluster samples based on correlation coefficients
# dev.off()
# 
# jpeg(filename = "Full-Sample-Methylation-PCA-FilteredCov5.jpeg", height = 1000, width = 1000) #Save file with designated name
# methylKit::PCASamples(methylationInformationFilteredCov5) #Run a PCA analysis on percent methylation for all samples
# dev.off() #Turn off plotting device
# 
# jpeg(filename = "Full-Sample-Methylation-Screeplot-FilteredCov5.jpeg", height = 1000, width = 1000) #Save file with designated name
# methylKit::PCASamples(methylationInformationFilteredCov5, screeplot = TRUE) #Run the PCA analysis and plot variances against PC number in a screeplot
# dev.off()


##########################################################################################################
## Differentially methylated loci (DML analysis)

setwd(dirname(home)); getwd() # reset working directory to home dir

# Identify DML

#Code that was used to test calculateDiffMeth parameters
differentialMethylationStatsTrtCov1   <- methylKit::calculateDiffMeth(methylationInformationFilteredCov1)
differentialMethylationStatsTrtCov5   <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5)
differentialMethylationStatsTrtCov10  <- methylKit::calculateDiffMeth(methylationInformationFilteredCov10)

# differentialMethylationStatsTreatment2 <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5, overdispersion = "MN", test = "Chisq")
# head(differentialMethylationStatsTreatment2) #Look at differential methylation statistics

setwd('DML')
write.csv(differentialMethylationStatsTrtCov10, "DML-stats-ploidy-Cov10.csv") #Save table as .csv

diffMethStatsTreatment20 <- methylKit::getMethylDiff(differentialMethylationStatsTrtCov10, difference = 20, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsTreatment20$chr) #Count the number of DML
head(diffMethStatsTreatment20) #Confirm creation

write.csv(diffMethStatsTreatment20, "DML-getMethylDiff-ploidy-Cov10-20.csv") #Save table as .csv

diffMethStatsTreatment50 <- methylKit::getMethylDiff(differentialMethylationStatsTreatment, difference = 50, qvalue = 0.01) #Identify loci that are at least 50% different
length(diffMethStatsTreatment50$chr) #Count the number of DML
head(diffMethStatsTreatment50) #Confirm creation

diffMethStatsTrt50Cov10 <- methylKit::getMethylDiff(differentialMethylationStatsTrtCov10, difference = 50, qvalue = 0.01) #Identify loci that are at least 50% different
length(diffMethStatsTrt50Cov10$chr) #Count the number of DML
head(diffMethStatsTrt50Cov10) #Confirm creation

write.csv(diffMethStatsTrt50Cov10, "DML-getMethylDiff-temp-Cov10-50.csv") #Save table as .csv

##### Hypermethylated DML
diffMethStats50FilteredCov5Hyper <- methylKit::getMethylDiff(differentialMethylationStatsTreatment, difference = 50, qvalue = 0.01, type = "hyper") #Identify hypermethylated loci that are at least 50% different
head(diffMethStats50FilteredCov5Hyper) #Confirm creation
write.csv(diffMethStats50FilteredCov5Hyper, "output/DML-getMethylDiff-temp-Cov5-50-hyper.csv") #Save table as .csv

##### Hypomethylated DML
diffMethStats50FilteredCov5Hypo <- methylKit::getMethylDiff(differentialMethylationStatsTreatment, difference = 50, qvalue = 0.01, type = "hypo") #Identify hypermethylated loci that are at least 50% different
head(diffMethStats50FilteredCov5Hypo) #Confirm creation
write.csv(diffMethStats50FilteredCov5Hypo, "output/DML-getMethylDiff-temp-Cov5-50-hypo.csv") #Save table as .csv

#### DML distribution figure

##### Calculate distribution of DML across chromosomes
DMLchrCounts <-as.data.frame(table(diffMethStatsTreatment50$chr)) #Count the number of DML/chromosome
DMLchrCounts <- DMLchrCounts[-c(11:length(DMLchrCounts$Var1)),]
prefix <- "Chr"
suffix <- seq(1:length(DMLchrCounts$Var1))
my_names <- paste(prefix,suffix,sep="")
DMLchrCounts$chr <- my_names
DMLchrCounts <- DMLchrCounts[,-1] #Remove column with chromosome RefSeq ID
DMLchrCounts <- DMLchrCounts[,c(2,1)] #Reorganize columns
colnames(DMLchrCounts) <- c("chr", "DMLCount") #Rename columns
head(DMLchrCounts) #Confirm formatting
write.csv(DMLchrCounts, "output/DML-per-Chromosome.csv", row.names = FALSE, quote = FALSE) #Save file

##### [*C. gigas* genome information from NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=pacific+oyster)

DMLchrCounts$geneCount <- c(3425, 3526, 3796, 3124, 3915, 4129, 4097, 3617, 1829, 3713) #Create column with number of gene sequences in each chromosome
DMLchrCounts$DMLbyGenes <- (DMLchrCounts$DMLCount / DMLchrCounts$geneCount)*100
head(DMLchrCounts) #Confirm column creation

DMLchrCounts$chrLengthMb <- c(55.79, 73.22, 58.32, 53.13, 73.55, 60.15, 62.11, 58.46, 37.09, 57.54) #Create column with number Mb in each chromosome
DMLchrCounts$DMLbyChrLength <- (DMLchrCounts$DMLCount / DMLchrCounts$chrLengthMb) #Normalize DML counts by chr length in Mb. Axes should indicate that these counts are multiplied by 10e-6
head(DMLchrCounts) #Confirm column creation

DMLchrCounts$chrCpGCounts <- c(1394962, 1237469, 1568195, 1200917, 2008680, 1121671, 1251184, 1633946, 2300335, 740913) #Create column with number of CpGs in each chromosome
DMLchrCounts$DMLbyChrCpG <- (DMLchrCounts$DMLCount / DMLchrCounts$chrCpGCounts) #Normalize DML counts by number of CpGs in each chromosome
range(DMLchrCounts$DMLbyChrCpG) #Look at range of values
head(DMLchrCounts) #Confirm column creation

# 
# ```{bash}
# #Count the number of CGs in each chromosme (excluding MT chromosome)
# grep "NC_047559.1" ../2019-05-13-Generating-Genome-Feature-Tracks/C_gigas-3.0_CG-motif.bed | wc -l
# grep "NC_047560.1" ../2019-05-13-Generating-Genome-Feature-Tracks/C_virginica-3.0_CG-motif.bed | wc -l
# grep "NC_035782.1" ../2019-05-13-Generating-Genome-Feature-Tracks/C_virginica-3.0_CG-motif.bed | wc -l
# grep "NC_035783.1" ../2019-05-13-Generating-Genome-Feature-Tracks/C_virginica-3.0_CG-motif.bed | wc -l
# grep "NC_035784.1" ../2019-05-13-Generating-Genome-Feature-Tracks/C_virginica-3.0_CG-motif.bed | wc -l
# grep "NC_035785.1" ../2019-05-13-Generating-Genome-Feature-Tracks/C_virginica-3.0_CG-motif.bed | wc -l
# grep "NC_035786.1" ../2019-05-13-Generating-Genome-Feature-Tracks/C_virginica-3.0_CG-motif.bed | wc -l
# grep "NC_035787.1" ../2019-05-13-Generating-Genome-Feature-Tracks/C_virginica-3.0_CG-motif.bed | wc -l
# grep "NC_035788.1" ../2019-05-13-Generating-Genome-Feature-Tracks/C_virginica-3.0_CG-motif.bed | wc -l
# grep "NC_035789.1" ../2019-05-13-Generating-Genome-Feature-Tracks/C_virginica-3.0_CG-motif.bed | wc -l
# ```

### Need to update these numbers for c.gigas
DMLchrCounts$chrCpGCounts <- c(1394962, 1237469, 1568195, 1200917, 2008680, 1121671, 1251184, 1633946, 2300335, 740913) #Create column with number of CpGs in each chromosome
DMLchrCounts$DMLbyChrCpG <- (DMLchrCounts$DMLCount / DMLchrCounts$chrCpGCounts) #Normalize DML counts by number of CpGs in each chromosome
range(DMLchrCounts$DMLbyChrCpG) #Look at range of values
head(DMLchrCounts) #Confirm column creation

### Create Figure

pdf("output/2021-05-25-DML-and-Gene-Distribution.pdf", height = 8.5, width = 12) #Save figure as a pdf
par(mar = c(5,7,2,10)) #Change figure boundaries
DMLbarplot <- barplot(as.matrix(t(DMLchrCounts$DMLbyChrCpG)),
                      axes = FALSE, 
                      names.arg = DMLchrCounts$chr,
                      cex.names = 1.5,
                      xlim = c(0.7,11.5),
                      ylim = c(0,60e-05),
                      col = "grey80") #Create a barplot and save as a new object. Use axes = FALSE to remove the y-axis and names.arg to set labels on the x-axis. The object contains x coordinates for bars, so xlim is set at 12 to compensate for maximum value of 11.5
mtext(side = 1, "Chromosome", las=1, line = 3, cex = 1.5) #Add x-axis label
axis(side = 2, line = 1.5, at = seq(0, 60e-05, by = 10e-05), las = 2, col = "grey80", cex.axis = 1.2) #Add y-axis for DML counts
mtext(side = 2, "Number DML per 10,000 CpGs", line = 5, cex = 1.5) #Add y-axis label for DML counts
par(new = TRUE) #Create a new plot
plot(x = DMLbarplot,
     y = DMLchrCounts$geneCount,
     type = "b",
     axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i",
     pch = 16, col = "grey20",
     xlim = c(0,12), ylim = c(0,6500)) #Plot points and lines (type = "b") for gene count by chromosome. Use the coordinates from DMLbarplot (x = DMLbarplot) and set xlim = (0,12) so plots are lined up. Use axes = FALSE to remove both axes. Remove x and y lables (xlab = ""; ylab = ""). Set ylim = (0,6500) to account for max y-values. Use xaxs and yxaxs to remove space between axes.
axis(side = 4, line = 1.5, at = seq(0, 6500, by = 500), las = 2, col = "grey80", cex.axis = 1.2) #Add y-axis for gene sequence counts
mtext(side = 4, "Number of Genes", line = 6, cex = 1.5) #Add y-axis label for gene sequence counts
dev.off() #Turn off plotting device


###### Principal Component Analyses

sample.IDs <- list("1", "2", "3", "4", "5", "6", "7", "8", "9", "10") #Create list of sample IDs
treatmentSpecification <- c(rep(0, times = 5), rep(1, times = 5)) 

plotCustomization <- data.frame(sample = 1:10,
                                treatmentSpecification) #Create dataframe with sample treatment information
head(plotCustomization) #Confirm dataframe creation

allDataPCA <- PCASamples(methylationInformationFilteredCov5, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA) #Look at summary statistics. The first PC explains 15.6% of variation, the second PC explains 12.1% of variation


RColorBrewer::display.brewer.all() #Show all RColorBrewer palettes. I will choose greens.
plotColors <- rev(RColorBrewer::brewer.pal(5, "GnBu")) #Create a color palette for the barplots. Use 5 green-blue shades from RColorBrewer. Reverse theorder so the darkest shade is used first.

pdf("output/2021-05-26-All-Data-PCA.pdf", width = 11, height = 8.5)
par(mar = c(5, 5, 1, 1)) #Specify inner and outer margins
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep(plotColors[2], times = 5), rep(plotColors[4], times = 5)), pch = c(rep(16, times = 5), rep(17, times = 5)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2
#Add multiple white boxes on top of the default black box to manually change the color
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
ordiellipse(allDataPCA, plotCustomization$treatment, show.groups = "1", col = plotColors[4]) #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, plotCustomization$treatment, show.groups = "0", col = plotColors[2]) #Add confidence ellipse around the samples in ambient pCO2
axis(side =  1, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1 (15.6%)", line = 3, cex = 1.5) #Add x-axis label
axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2 (12.1%)", line = 3, cex = 1.5) #Add y-axis label
legend("topleft", 
       pch = c(16, 17), 
       legend = c("diploid", "triploid"), 
       col = c(plotColors[2], plotColors[4]), 
       cex = 1.7, bty = "n") #Add a legend with information about ambient and elevated samples
dev.off()

####### DML

# The first thing I need to do is subset the data so it only includes DML. Then, I can use similar code to what I used above to create PCA plots.

DMLPositions <- rep(0, times = length(diffMethStatsTreatment50$chr)) #Create an empty vector with 598 places to store row numbers
for (i in 1:length(DMLPositions)) {
  DMLPositions[i] <- which(getData(diffMethStatsTreatment50)$start[i] == getData(methylationInformationFilteredCov5)$start)
} #For each DML, save the row number where that DML is found in methylationInformationFilteredCov5
tail(DMLPositions) #Confirm vector was created

DMLMatrix <- methylationInformationFilteredCov5[DMLPositions,] #Subset methylationInformationFilteredCov5Destrand to only include DML and save as a new methylBase object
sum((DMLMatrix$start) == (diffMethStatsTreatment50$start)) == length(diffMethStatsTreatment50$start) #Confirm that start columns are identical. If they are identical, the sum of all TRUE statements should equal the length of the original methylBase object
tail(DMLMatrix) #Confirm methylBase object creation

DMLDataPCA <- PCASamples(DMLMatrix, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(DMLDataPCA) #Look at summary statistics. The first PC explains 47.6% of variation, the second PC explains 9.5% of variation


pdf("output/2021-05-26-DML-Only-PCA.pdf", width = 11, height = 8.5)
par(mar = c(5, 5, 1, 1)) #Specify inner and outer margins
fig.DMLDataPCA <- ordiplot(DMLDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.DMLDataPCA, "sites", col = c(rep(plotColors[2], times = 5), rep(plotColors[4], times = 5)), pch = c(rep(16, times = 5), rep(17, times = 5)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2
#Add multiple white boxes on top of the default black box to manually change the color
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
ordiellipse(DMLDataPCA, plotCustomization$treatment, show.groups = "1", col = plotColors[4]) #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(DMLDataPCA, plotCustomization$treatment, show.groups = "0", col = plotColors[2]) #Add confidence ellipse around the samples in ambient pCO2
axis(side =  1, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1 (61.8%)", line = 3, cex = 1.5) #Add x-axis label
axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2 (7.5%)", line = 3, cex = 1.5) #Add y-axis label
legend("topright", 
       pch = c(16, 17), 
       legend = c("diploid", "triploid"), 
       col = c(plotColors[2], plotColors[4]), 
       cex = 1.7, bty = "n") #Add a legend with information about ambient and elevated samples
dev.off()

#### Multipanel plot

pdf("output/2021-05-26-PCA-Multpanel.pdf", height = 8.5, width = 11) #Save plot
par(mfrow = c(1, 2), oma = c(5, 2, 2, 0), mar = c(0, 3, 0, 5)) #Set up parameters for multipanel plot
#All CpG Loci
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep(plotColors[2], times = 5), rep(plotColors[4], times = 5)), pch = c(rep(16, times = 5), rep(17, times = 5)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2
#Add multiple white boxes on top of the default black box to manually change the color
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
ordiellipse(allDataPCA, plotCustomization$treatment, show.groups = "1", col = plotColors[4]) #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, plotCustomization$treatment, show.groups = "0", col = plotColors[2]) #Add confidence ellipse around the samples in ambient pCO2
axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1 (15.6%)", line = 3, cex = 1.5) #Add x-axis label
axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2 (18.1%)", line = 3, cex = 1.5) #Add y-axis label
mtext(side = 3, line = -5, adj = c(-100,0), text = "    a. All CpG Loci", cex = 1.5)

fig.DMLDataPCA <- ordiplot(DMLDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-20,20), ylim = c(-10,20), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.DMLDataPCA, "sites", col = c(rep(plotColors[2], times = 5), rep(plotColors[4], times = 5)), pch = c(rep(16, times = 5), rep(17, times = 5)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2
#Add multiple white boxes on top of the default black box to manually change the color
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
ordiellipse(DMLDataPCA, plotCustomization$treatment, show.groups = "1", col = plotColors[4]) #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(DMLDataPCA, plotCustomization$treatment, show.groups = "0", col = plotColors[2]) #Add confidence ellipse around the samples in ambient pCO2
axis(side =  1, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1 (61.8%)", line = 3, cex = 1.5) #Add x-axis label
axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2 (7.5%)", line = 3, cex = 1.5) #Add y-axis label
mtext(side = 3, line = -5, adj = c(-100, 0), text = "    b. DML", cex = 1.5) #Add test category
legend(x = 0, y = -22, 
       pch = c(16, 17), 
       legend = c("Control", "Elevated"), 
       col = c(plotColors[2], plotColors[4]),
       y.intersp = 1, x.intersp = 1,
       cex = 1.7, bty = "n") #Add a legend with information about ambient and elevated samples
dev.off()

######################################################################################################
### Heat Map

percMethDML <- percMethylation(DMLMatrix, rowids = TRUE) #Get percent methylation for all samples at DML. Include row IDs (chr, start, end) information
head(percMethDML) #Confirm percent methylation matrix was created

pheatmap(percMethDML, color = rev(plotColors),
         cluster_rows = TRUE, clustering_distance_rows = "euclidean", treeheight_row = 70, show_rownames = FALSE,
         cluster_cols = TRUE, clustering_distance_cols = "euclidean", treeheight_col = 40, show_colnames = FALSE,
         annotation_col = data.frame(pCO2 = factor(rep(c("Ambient","Treatment"), each = 5))),
         annotation_colors = list(pCO2 = c(Ambient = "grey90", Treatment = "grey10")), 
         annotation_legend = FALSE, annotation_names_col = FALSE,
         legend = TRUE) #Create heatmap using pheatmap using percMethDML and plotColors color scheme. Cluster rows and columns using euclidean distances. Adjust the dendogram tree heights and do not show any row or column names. Create a dataframe with treatment information using annotation_col. Use annotation_colors to indicate colors for treatment ("grey90") and ambient ("grey10") samples. Do not include an annotation_legend or name for annotations (annotatino_names_col). Include a legend.

pdf("2021-05-27-DML-Only-Heatmap.pdf", height = 8.5, width = 11)
par(oma = c(0, 1, 0, 0)) #Adjust outer margins
heatmap.2(percMethDML, col = rev(plotColors), scale = "none", margins = c(1,1),
          trace = "none", tracecol = "black",
          labRow = FALSE, labCol = FALSE, 
          ColSideColors = c(rep("grey90", times = 5), rep("grey10", times = 5)),
          key = TRUE, keysize = 1.8, density.info = "density", key.title = "", key.xlab = "% Methylation", key.ylab = "",
          key.par = list(cex.lab = 2.0, cex.axis = 1.5)) #Create heatmap using heatmap.2 from gplots package using percMethDML data. Use plotColors but do not scale data, label rows, or label columns. Use ColSideColors to indicate colors for treatment and ambient samples. Add a legend using key, and adjust keysize. Have key display density data with density.info. Do not add a key title or y-axis label, and label x axis with key.xlab.
mtext("Density", cex = 1.6, las = 3, adj = 0.8, padj = -29) #Manually add y-axis label for key since heatmap.2 doesn't let you change font size
dev.off()


############################################################################################################
### Convert DML list to Bed files

DML <- data.frame(diffMethStatsTreatment50$chr,
                 diffMethStatsTreatment50$start,
                 diffMethStatsTreatment50$end,
                 diffMethStatsTreatment50$meth.diff) # arrange(chr, start) %>% #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML, "2021-05-28-DML-Locations-treatment-50.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file
