#LOAD LIBRARIES
library("readxl")
library("vegan")
library("factoextra")
library("ggplot2")
library("ggfortify")
library("ggrepel")
library("R.utils")

#READ IN DATA
#Load in the “PCA – Normalized Samples and QCs” file as the “varespec” variable and load in the “categories_all_peaks” file as the “categ” variable for making the PCA Plot featuring all detected peaks
#OR
#Load in the “PCA_Amino_Acids_Normalized” file as the “varespec” variable and load in the “categories_amino_acids” file as the “categ” variable for making the PCA Plot featuring just the amino acid data

varespec <- read.csv(file = file.choose(), header = TRUE, na.strings = c("NA")) #this is the name of the dataframe with peak area values. Top row is the peak number label.
categ <- read.csv(file = file.choose(), header = TRUE) #these are labels for the dataframe. First column is the full name ("JPV Day 0 B1_11.cdf area" for example)
#second column is either "Burn" or "Unburn"
#third column is the day ("Day 0", "Day 3", etc.)


#this takes the ZERO values for a given metabolite and assigns the zero values to 10% of the lowest non-missing value (compound minimum)
#this is used for GC-MS data that has been normalized
for (i in 1:ncol(varespec)) {
  minimum_replacement <- min(varespec[,i][which(varespec[,i]>0)], na.rm = TRUE)*0.10
  for (j in 1:nrow(varespec)) {
    if (isZero(varespec[j,i])) { #check to see if the element is zero
      varespec[j,i] <- abs(jitter(minimum_replacement))#replace the zero value with 10% of the minimum of that column and the jitter function adds/subtracts noise from the replaced value
    } 
  }
}



#this is the function that can be used for Pareto Scaling. This is just making the function
paretoscale <- function(data, exclude = T) {
  
  if (exclude == T) {
    # Here we extract numeric data and perform Pareto scaling
    sample_classes <- data[, 1:2]
    x <- data[, 3:dim(data)[2]]
  } else {
    sample_classes <- NULL
    x <- data
  }
  # Here we perform centering
  x.centered <- apply(x, 2, function(x) x - mean(x))
  # Then we perform scaling on the mean-centered matrix
  x.sc <- apply(x.centered, 2, function(x) x/sqrt(sd(x)))
  x.sc <- cbind(sample_classes, x.sc)
  
}


varespec <- paretoscale(varespec, exclude = F) #this apples Pareto Scaling to the dataframe using the Pareto Scaling function that was made above
pca.prcomp <- prcomp(varespec, scale = FALSE, center = FALSE) #the Pareto Scaling function already scales and centers the data. Therefore,scaling and centering the data here is not needed
pca.prcomp

pca.prcomp.data <- data.frame(pca.prcomp$x) #creates principal components as columns, samples as rows
pca.prcomp.data$plotx <- pca.prcomp.data[,1] #assigns first column as plotx
pca.prcomp.data$ploty <- pca.prcomp.data[,2] #assigns second column as ploty

Unburned_or_Burned = as.character(categ$Factor.1) #includes a column with categories
pca.prcomp.data$Unburned_or_Burned <- Unburned_or_Burned #includes a column with categories

Day = as.character(categ$Factor.2)
pca.prcomp.data$Day <- Day


#PCA Plot for Unburned_or_Burned for all peaks (with quality controls)
pca.prcomp.plot <- ggplot(pca.prcomp.data, aes(x=plotx, y=ploty, color = Unburned_or_Burned)) +
  geom_point(size=2.5) +
  scale_color_manual(values=c("red","green4","dodgerblue"))+
  ggtitle("All Detected Peaks (986)") +
  xlab("Principal Component 1 (33.59%)") + ylab("Principal Component 2 \n(9.24%)") +   #YOU WILL HAVE TO MANUALLY PUT IN THESE VALUES FOR EACH NEW DATASET THAT YOU IMPORT 
  labs(color=NULL)+ #this removes the legend title
  coord_fixed(ratio = 1) +
  theme_bw() + theme(aspect.ratio=1) + theme(panel.grid = element_blank()) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0.3))),  #this adds the circle color shading to group the data points
               data = pca.prcomp.data[pca.prcomp.data$Unburned_or_Burned != "versicolor",])+
  scale_shape_discrete(name="", labels=c("Burned","Quality Control","Unburned"))+
  scale_color_manual(name="",labels=c("Burned","Quality Control", "Unburned"),values=c("red","green4","dodgerblue"))+
  theme(axis.text.x=element_text(size=24, colour = "black"))+
  theme(axis.text.y=element_text(size=24, colour = "black"))+
  theme(axis.title.x=element_text(size=32))+
  theme(axis.title.y=element_text(size=32))+
  theme(legend.text=element_text(size=36))+
  theme(plot.title = element_text(hjust =0.5, size=36))

pca.prcomp.plot



#Unburned and Burned without quality controls (amino acid PCA plot)
pca.prcomp.plot <- ggplot(pca.prcomp.data, aes(x=plotx, y=ploty, color = Unburned_or_Burned)) +
  geom_point(size=2.5) +
  scale_color_manual(values=c("red","dodgerblue"))+
  ggtitle("Annotated Amino Acids (6)") +
  xlab("Principal Component 1 (81.19%)") + ylab("Principal Component 2 \n(12.59%)") +   #YOU WILL HAVE TO MANUALLY PUT IN THESE VALUES FOR EACH NEW DATASET THAT YOU IMPORT 
  labs(color=NULL)+ #this removes the legend title
  coord_fixed(ratio = 1) +
  theme_bw() + theme(aspect.ratio=1) + theme(panel.grid = element_blank()) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0.3))),  #this adds the circle color shading to group the data points
               data = pca.prcomp.data[pca.prcomp.data$Unburned_or_Burned != "versicolor",])+
  scale_shape_discrete(name="", labels=c("Burned","Unburned"))+
  scale_color_manual(name="",labels=c("Burned","Unburned"),values=c("red","dodgerblue"))+
  theme(axis.text.x=element_text(size=24, colour = "black"))+
  theme(axis.text.y=element_text(size=24, colour = "black"))+
  theme(axis.title.x=element_text(size=32))+
  theme(axis.title.y=element_text(size=32))+
  theme(legend.text=element_text(size=36))+
  theme(plot.title = element_text(hjust =0.5, size=36))

pca.prcomp.plot




#### PERMANOVA ####

# If you need to install vegan package, use the following:
# install.packages ('vegan')
library (vegan)


#Using Euclidean model:  
varespec_distance <- vegdist(varespec, method = "euclidean") #this generates a distance matrix called "varespec_distance" using the "manhattan method"
varespec_distance
varespec_distance_div <- adonis2(varespec_distance~Factor.1, data=categ, method = "euclidean") #this code then conducts PERMANOVA to determine if there is a significant difference amongst samples according to a certain variable
#Factor.1 refers to burned, unburned, and quality control samples
#"data=categ" refers to the dataframe that contains the categories for my data (e.g., name, burned vs. unburned,sampling date)
#"method" refers to the method being used to create the distance matrix
varespec_distance_div
#Check out the table at the end of this code in the console
#The Factor.1 row refers to burned, unburned, and quality control samples
#The R2 value in the Factor.1 row is the R-squared value which is usually reported
#The Pr(>F) value is the p-value. In this case, the p-value is <=0.001, indicating that there is a significant difference amongst the centroids of the burned, unburned, and quality control groups 
#The next step is to do a post-hoc test to determine which pairs are specifically contributing to those observed significant difference

## Conduct p-value adjustments as needed
## P-value adjustments:
library(stats)
p = c(0.001, 0.001, 0.001)
p.adjusted = p.adjust(p, method = "BH")
p.adjusted

## New method that uses adonis2 (not adonis)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

varespec_distance_div_pairwise <- pairwise.adonis2(varespec_distance~Factor.1, data = categ, sim.method = "euclidean", p.adjust.m = "BH", perm = 999)
#Factor.1 refers to burned, unburned, and quality control samples
#"data=categ" refers to the dataframe that contains the categories for my data (e.g., name, burned vs. unburned,sampling date)
#"sim.method" refers to the method being used to create the distance matrix
varespec_distance_div_pairwise
#In this case, the p-value when comparing Burn vs. Unburn, Burn vs. Quality Control, and Unburn vs. Quality control are all <=0.001 when looking at all peaks
#This result indicates that there is a significant difference between the centroids of burn vs. unburn, burn vs. quality control, and unburn vs. quality control.
