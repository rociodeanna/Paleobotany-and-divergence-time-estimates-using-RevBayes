

#----------------------PACKAGES----------------------#

# Install packages
install.packages("vegan")
install.packages("cluster")
install.packages("ggplot2")

# Call the packages installed
library(vegan)
library(cluster)
library(ggplot2)



#----------------------ANALYSIS----------------------#
  
# Setting working directory
setwd('C:/Users/XXX')

# Read data as CSV
data<-read.csv("Sol_traits.csv", header=T, row.names=1) #

# Make dissimilarity matrix -- this function allows multiple variable types and missing data
# it needed to know what to do with the binary variables -- e.g., sym flags them as symmetric binary
dissim_mat<-daisy(data, metric="gower", type=list(symm=c(6,7,8,9,10,11,12,14,15), ordtype=(13))) # default distance is euclidean

# Perform the NMDS analysis -- you will see each iteration of the NMDS until a solution is reached (i.e., stress was minimized after some number of reconfigurations of the points in 2 dimensions)
mds<-metaMDS(dissim_mat,  trymax = 100, k=3, autotransform = FALSE, wascores = FALSE, noshare = FALSE) # this function automatically detects that the input is a matrix
points<-data.frame(mds$points) # mds reduces the data into two dimensions, MDS1 and MDS2, that's what's in points
rownames(points)<-rownames(data) # this is just to get ready to plot

#This fits vectors (morphological characters) onto an ordination
char = data[(1:15)]
en = envfit(mds, char, permutations = 999, na.rm = TRUE)
en

# The envfit output contains information on the length of the segments for each variable. The segments are scaled to the r2 value, so that the morphological characters with a longer segment are more strongly correlated with the data than those with a shorter segment. 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) # extract the required information from the envfit result
en_coord_cont



#----------------------STRESS----------------------#

stressplot(mds)
print(mds$stress)


# the function NMDS.scree() automatically performs a NMDS for 1-10 dimensions and plots the nr of dimensions vs the stress value
NMDS.scree <- function(x) { # where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1, trymax = 20)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

# Use the function that we just defined to choose the optimal number of dimensions
NMDS.scree(dissim_mat)


#----------------------PLOTTING----------------------#

# plot the results
attach(points)
plot(MDS1, MDS2, xlab = 'MDS1', ylab = 'MDS2', main = 'NMDS Analysis of Extinct and Extant Fruits', cex=0.4)
text(MDS1, MDS2, labels = row.names(points), pos = 4, cex=0.4, offset=0.2)
plot(en)


## plotting better graphics using classifiers for your species
classifier<-read.csv("Sol_classifier.csv", header=T, row.names=1)  # add data with fossil/extant classification
data.scores = as.data.frame(scores(mds$points)) # extract the sample coordinates in the NMDS ordination space and convert to a data frame
data.scores$state = classifier$state # add columns from your original data that contain information that you want (e.g., if it's either fossil or extant)


# plot the NMDS (changes colors, text, etc as you like)
gg = ggplot(data = data.scores, aes(x = MDS1, y = MDS2)) + 
  geom_point(data = data.scores, aes(colour = state), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c("steelblue", "red")) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Sample")

gg

# plot the vectors (morphological characters)
gg = ggplot(data = data.scores, aes(x = MDS1, y = MDS2)) + 
  geom_point(data = data.scores, aes(colour = state), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c("steelblue", "red"))  + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Seed")

gg



