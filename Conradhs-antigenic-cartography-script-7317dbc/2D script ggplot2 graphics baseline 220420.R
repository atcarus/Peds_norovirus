# This script is designed to map the antigenic cartography of neutralization and avidity datasets
# Written by Helen Conrad
# most recently modified 02/07/2022 
#install the package for the antigenic cartography
#install.packages("stringr","devtools", "rainbow", "rgl","reshape2","Rcpp","seqinr","ape")
#devtools::install_github("acorg/Racmacs")
#set the working directory to wherever the data is saved
#setwd(choose.dir())
setwd("C:/Users/conradws/Desktop/General Work/Tech applications/R work/data and results/220211 Norovirus paper/noro data by year and age")
getwd()
#load in the package
library(stringr)
library(Racmacs)
library(reshape2)
library(Rcpp)
library(ape)
library(seqinr)
library(ggplot2)
library(colorspace)

#read in the data into the format accepted by Racmacs and create the map
#df = read.csv(file = "broken out groups BEI data.csv", header = TRUE, sep = ",")
df = read.csv(file = file.choose())
# set the row names to the virus names and "groups" for the row containing the list of sample groups
row.names(df) = df[,c(1)]
df = df[-c(1)]
# select the names of the different groups in the data set and save them as a separate object
groupNames = df[grep("group1",row.names(df)),]
df = df[-grep("group", row.names(df)),]
#save the data as a csv file and then read it back in. This is a little clunky, but Racmacs bugs out if I do it any other way
#write.csv(df, "Cleaned data.csv")
table = df
table[is.na(table)]<-"*"
map = acmap(titer_table = table)

#IMPORTANT: CHECK NEXT INSTRUCTION BEFORE PROCEEDING
#set the dilution step size to whatever was used in the experiments to ensure distances are calculated correctly
dilutionStepsize(map)<- 2

#set the names of the viruses in the acmap object and associate each one with a color
agNames(map) <- row.names(table)
agGroups(map) <- row.names(table)
#set the rng seed to make the results more reproducible 
set.seed(123)
# note that the map is being created for 3 dimensional mds
map = optimizeMap(map = map, number_of_dimensions = 2, number_of_optimizations = 1000, minimum_column_basis = "none")
save.coords(map, "racmacs coordinates.csv")
df1 = read.csv("racmacs coordinates.csv")

# set the names of the antigens and label all the sera as "sera"
# Make a data frame with the antigens and one with the sera
df2 <- df1[str_detect(df1$type, "antigen"), ]
df3 <- df1[str_detect(df1$type, "sera"), ]
df2 = subset(df2, select = -c(1))
df3$type = t(groupNames)
df3 = subset(df3, select = -c(2))
#calculate the antigenic distances
df4 <- df2
rownames(df4) = df4[,1]
df4 = subset(df4, select = -c(1))
AgDistances = melt(as.matrix(dist(df4)))
colnames(AgDistances) = c("Virus 1","Virus 2","Distance")
#write.csv(AgDistances, file.choose())

# set the column names to the same thing so rbind can join them
colnames(df3)<- colnames(df2)
df1<- rbind(df2, df3)
# NOTE: Comment in the line below to remove antigens and only keep sera
#df1=df3
#Keep the antigen names and sera as groups, then cut that column from the data frame
antigenGroups = as.data.frame(unique(df2[,1]))
antigenGroups = antigenGroups[order(antigenGroups[,1]),]
antigenGroups = as.data.frame(antigenGroups)
seraGroups = as.data.frame(unique(df3[,1]))
seraGroups = seraGroups[order(seraGroups[,1]),]
seraGroups = as.data.frame(seraGroups)
df = subset(df1, select = -c(1))
df$Groups = df1$name

# read in the alignment and calculate the distances
#da = read.alignment(file = file.choose(), format = "fasta")
da = read.alignment(file = "Norovirus alignment.fasta", format = "fasta")
da1 = as.matrix(dist.alignment(da, matrix = c("similarity")))
da2 = as.vector(da1)
AgDistanceMatrix = as.matrix(dist(df4))
stristancematrix = as.vector(AgDistanceMatrix)
# scale the distances to make sure it's comparable to the serum data
#da1 = da1*(mean(stristancematrix)/mean(da2))
da1 = da1*7
gencoord = as.data.frame(cmdscale(da1,k = 2))
gencoord$names = row.names(gencoord)
gencoord = gencoord[order(gencoord$names),]

#Add colors to the table of the group names and standardize the column names
antigenGroups$colors = rainbow(nrow(antigenGroups))
seraGroups$colors = rainbow(nrow(seraGroups))
#set the colors based on hex code or color name
antigenGroups[1,2]="springgreen4"
antigenGroups[2,2]="purple"
antigenGroups[3,2]="black"
antigenGroups[4,2]="sienna"
#antigenGroups[5,2]="indianred3"
#antigenGroups[6,2]="#CC0033"
#antigenGroups[7,2]="#339900"
seraGroups[1,2]="red4"
seraGroups[2,2]="red2"
seraGroups[3,2]="violet"
seraGroups[4,2]="thistle"
seraGroups[5,2]="skyblue"
seraGroups[6,2]="royalblue"
colnames(antigenGroups) = c("name", "colors")
colnames(seraGroups)<- colnames(antigenGroups)
Legend<- rbind(seraGroups, antigenGroups)
antigencoord = merge(x = df2, y = antigenGroups, all = FALSE)
colnames(antigencoord) = c("name", "X", "Y", "color")
seracoord = merge(x = df3, y = seraGroups, all = FALSE)
colnames(seracoord) = colnames(antigencoord)
gencoord$colors = antigenGroups$colors
gencoord$type <- "genetic"
antigencoord$type<- "antigen"
seracoord$type<- "sera"

colnames(gencoord) = c("X", "Y", "name", "color", "type")
gencoord = gencoord[c("name", "X", "Y", "color", "type")]
antigencoord = antigencoord[c("name", "X", "Y", "color", "type")]
seracoord = seracoord[c("name", "X", "Y", "color", "type")]
gencoord$name = antigencoord$name
allcoord = rbind(seracoord, gencoord,antigencoord)

allcoord$transparency = allcoord$type
allcoord$transparency[allcoord$transparency == "sera"]<-0.3
allcoord$transparency[allcoord$transparency == "antigen"]<-1
allcoord$transparency[allcoord$transparency == "genetic"]<-1

allcoord$shape = allcoord$type
allcoord$shape[allcoord$shape == "sera"]<-19
allcoord$shape[allcoord$shape == "antigen"]<-17
allcoord$shape[allcoord$shape == "genetic"]<-15

allcoord$size = allcoord$type
allcoord$size[allcoord$size == "sera"]<-1
allcoord$size[allcoord$size == "antigen"]<-2
allcoord$size[allcoord$size == "genetic"]<-2
#Note: run the line below for easy loading if it decides to be impolite and quit working
#allcoord = read.csv(file.choose())

#ggplot2 needs the color names in the form of a named vector or else it won't accept them, trying to submit them iin any other format will cause issues
colorDetails <- setNames(as.character(unique(allcoord$color)),unique(allcoord$name))

ggplot(allcoord, aes())+ 
  geom_point(aes(x = X, y = Y, alpha = transparency, color = name, shape = shape, size = as.numeric(as.character(size))), stroke = 0.7) + 
  scale_colour_manual(labels = unique(allcoord$name), values = colorDetails) + 
  scale_size(range = c(1.5,2.5))+
  scale_alpha_manual(values = c(0.5,1))+
  scale_shape_manual(labels = c("genetic", "antigen", "sera"), values = c(19,17,15))+
  labs(shape = "type", color = "name")+
  theme(panel.background = element_rect(fill = "white", colour = "grey60"), panel.grid.major = element_line(colour = "grey60"), panel.grid.minor = element_line(colour = "grey60")) +
  guides(alpha = "none", size = "none", shape=guide_legend(override.aes=list(size=2), order = 1), color=guide_legend(override.aes=list(size=2)))+
  xlim(-4,4)+ 
  ylim(-4,4)
max(allcoord$X, na.rm = TRUE)
min(allcoord$X, na.rm = TRUE)
max(allcoord$Y, na.rm = TRUE)
min(allcoord$Y, na.rm = TRUE)
write.csv(allcoord, "easy loading coordinates1.csv")

