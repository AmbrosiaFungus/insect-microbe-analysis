library (vegan)
library(ape)
library(ggplot2)


#read in Mapping file with the Details for each Sample
map <- read.csv("/home/robert/Projects/fungal_microbiom_Ai/data/Metadata_Miseq_2017/Metadata2.tsv",sep="\t")

names(map)[1]<-"SampleID"

#read in OTU table which i transposed in Python before

OTU <- read.csv("/home/robert/Projects/bacterial_microbiome_Ai/scripts/forward/Betadiversity/OTU_table_bacteria_filtered.csv", sep=",", header=TRUE)

#OTU table has a weird column name, fix that
names(OTU)[1]<-"SampleID"

#rearrange rows in mapping file to match relative abundance matrix

tab <- map[match(OTU$SampleID, as.character(map$SampleID)),]
summary(OTU$SampleID==tab$SampleID)

#merge the two Dataframes

bacteria <- merge(tab,OTU, by.x="SampleID", by.y = "SampleID")
scaledOTU <- scale(bacteria[1:nrow(bacteria),24:ncol(bacteria)], scale=T)
summary(scaledOTU, display=NULL)
mean(scaledOTU)

#delete the Mock-Community and negatives
bacteria_new <- bacteria[-c(63,64,65,66),]
bacteria_new$Location<-factor(bacteria_new$Location)
bacteria_new$Raffaelea<-factor(bacteria_new$Raffaelea)
#Correspondence analysis
spec.ca <- cca(bacteria_new[1:nrow(bacteria_new),24:ncol(bacteria_new)])
plot(spec.ca)

summary(spec.ca, display=NULL)


#convert OTU table to relative abundances with decostand

normi <- decostand(bacteria_new[1:nrow(bacteria_new),24:ncol(bacteria_new)], method="hellinger")
speci <- wcmdscale(vegdist(normi, method="bray"), eig = T)
plot(speci, cex=0.4)
speci$eig[speci$eig>=0]/sum(speci$eig[speci$eig>=0])
##decostand: pa: Presence/Absence or method "total" for relative abundance
##use hellinger which does a square root transformation before it does the relative abundance , but so it is still proportional

rob <- palette(c("red","blue","green","black", "orange", "yellow")) # set up the colour palette
with(bacteria_new, plot(scores(speci, display='sites'), # identifies coordinates
                  col=rob[bacteria_new$Location], #assign symbols, colours
                 xlab='Dim1 7%', ylab="Dim2 5%", cex=0.9,pch=c(1,2)[bacteria_new$Raffaelea])) # change axis labels
levels(bacteria_new$Raffaelea)
legend('bottomleft', # position the legend
       legend=c('Dampier', 'Dorrigo', 'Mount_Wilson', 'Olney', 'Termeil', 'Wombeyan'), # set legend text
       pch=c(19,19,19,19,19,19), col=rob, # assign symbols and colours
       cex=0.7)
legend('bottomright', # position the legend
       legend=c('no', 'yes'), # set legend text
       pch=c(1,2), col=c('black'), # assign symbols and colours
       cex=0.7)


NMDS = data.frame(Dim1 = scores(speci, choices = c(1), display='sites'), Dim2 = scores(speci, choices = c(2), display='sites'), Gender = bacteria_new$Gender, Location = bacteria_new$Location)

ggplot(NMDS, aes(x=Dim1, y=Dim2, col=Location)) +
  geom_point(aes(shape=bacteria_new$Raffaelea)) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "Bacterial Communties of Austroplatypus incompertus") +
  xlab("Dim1 (9.41%)") +
  ylab("Dim2 (5.36%)") 




#################################################################################################################################################################################################################
###################################################  Look at individual sites ################################################################################

################################################################# Olney ######################################################################################
Olney <- subset(bacteria_new, Location=='Olney')
Olney$Raffaelea<-factor(Olney$Raffaelea)
Olney$Host<-factor(Olney$Host)
Olney$Tree<-factor(Olney$Tree)
Olney$Collected<-factor(Olney$Collected)
Olney$Gender<-factor(Olney$Gender)
Olney$Alive_Dead<-factor(Olney$Alive_Dead)


normi_Olney <- decostand(Olney[1:nrow(Olney),24:ncol(Olney)], method="hellinger")
speci_Olney <- wcmdscale(vegdist(normi_Olney, method="bray"), eig = T)
plot(speci_Olney, cex=0.9)
speci_Olney$eig[speci_Olney$eig>=0]/sum(speci_Olney$eig[speci_Olney$eig>=0])

olney_color <- palette(c("red","blue","green","black", "orange", "yellow")) # set up the colour palette
with(Olney, plot(scores(speci_Olney, display='sites'), # identifies coordinates
                 col=olney_color[Olney$Host], #assign symbols, colours
                 xlab='Dim1 17%', ylab="Dim2 15%", cex=0.9,pch=c(1,2,3,4,5,6,7)[Olney$Tree])) # change axis labels
legend('bottomright', # position the legend
       legend=c('Eucalyptus pilularis', 'Eculayptus piperita'), # set legend text
       pch=c(19,19), col=olney_color, # assign symbols and colours
       cex=0.7)
legend('bottomleft', # position the legend
       legend=c('2','4','12', '13', '14', '16', '21'), # set legend text
       pch=c(1,2,3,4,5,6,7), col=c('black'), # assign symbols and colours
       cex=0.7)
title(main = "Bacterial Community Olney State Forest")



#Variation Partitioning for Olney
ambrosia.env.Olney <- Olney[, c('Gender', 'Tree', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

# The column collected should be a factor
ambrosia.env$Collected <- factor(ambrosia.env$Collected)
# set up full and null models for ordistep
bacteria.cap1.Olney <- capscale(normi_Olney ~ ., data=ambrosia.env.Olney, dist="bray")
bacteria.cap0.Olney <- capscale(normi_Olney ~ 1, data=ambrosia.env.Olney, dist="bray")

summary(ambrosia.env.Olney)

#perform forward and backward selection of explanatory variables 

step.env.Olney <- ordistep(bacteria.cap0.Olney, scope=formula(bacteria.cap1.Olney))

step.env.Olney$anova

OTU.var.Olney <- Olney[1:nrow(Olney),24:ncol(Olney)]

#Do the variation partitioning with these variables

bacteria.var.Olney <- varpart(OTU.var.Olney, ~Gender, ~Tree, data=ambrosia.env.Olney)

plot(bacteria.var.Olney, bg=1:3, Xnames=c('Gender', 'Tree'))

#grep the spatial patterns and transform them into a distance matrix

O.spatial_Ai <- subset(Olney, select=c("Longitude", "Latidude"))
Olney.pcnm <- as.data.frame(scores(pcnm(dist(O.spatial_Ai))))
dim(Olney.pcnm)

ambrosia_Olney.var <- varpart(OTU.var.Olney, Olney.pcnm, ~Tree, data=Olney)

plot(ambrosia_Olney.var, bg=1:3, Xnames=c('space', 'Tree'))

#Redundancy analysis for the site. 

ambrosia_Olney.rda <- rda(OTU.var.Olney ~ Tree + Condition(Olney.pcnm[, 1]) + Condition(Olney.pcnm[, 2]) + Condition(Olney.pcnm[, 3]), data=Olney)
anova(ambrosia_Olney.rda)

################################################## Dampier #########################################################################################

Dampier <- subset(bacteria_new, Location=='Dampier_State_Forest')
Dampier$Raffaelea<-factor(Dampier$Raffaelea)
Dampier$Host<-factor(Dampier$Host)
Dampier$Tree<-factor(Dampier$Tree)
Dampier$Collected<-factor(Dampier$Collected)
Dampier$Gender<-factor(Dampier$Gender)
Dampier$Alive_Dead<-factor(Dampier$Alive_Dead)


normi_Dampier <- decostand(Dampier[1:nrow(Dampier),24:ncol(Dampier)], method="hellinger")
speci_Dampier <- wcmdscale(vegdist(normi_Dampier, method="bray"), eig = T)
plot(speci_Dampier, cex=0.9)
speci_Dampier$eig[speci_Dampier$eig>=0]/sum(speci_Dampier$eig[speci_Dampier$eig>=0])

Dampier_color <- palette(c("red","blue","green","black", "orange", "yellow")) # set up the colour palette
with(Dampier, plot(scores(speci_Dampier, display='sites'), # identifies coordinates
                   col=Dampier_color[Dampier$Host], #assign symbols, colours
                   xlab='Dim1 15%', ylab="Dim2 13%", cex=0.9,pch=c(1,2,3,4,5,6,7)[Dampier$Tree])) # change axis labels
legend('bottomright', # position the legend
       legend=c('Eucalyptus sp. A', 'Eculayptus sp. B'), # set legend text
       pch=c(19,19, 19, 19), col=Dampier_color, # assign symbols and colours
       cex=0.7)
legend('bottomleft', # position the legend
       legend=c('1', '2', '4', '11', '14', '15', '29'), # set legend text
       pch=c(1,2,3,4,5,6,7), col=c('black'), # assign symbols and colours
       cex=0.7)
title(main = "Bacterial Community Dampier State Forest")

#Variation Partitioning for Olney
ambrosia.env.Dampier <- Dampier[, c('Gender', 'Tree', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

# The column collected should be a factor
ambrosia.env.Dampier$Collected <- factor(ambrosia.env.Dampier$Collected)
# set up full and null models for ordistep
bacteria.cap1.Dampier <- capscale(normi_Dampier ~ ., data=ambrosia.env.Dampier, dist="bray")
bacteria.cap0.Dampier <- capscale(normi_Dampier ~ 1, data=ambrosia.env.Dampier, dist="bray")

summary(ambrosia.env.Dampier)

#perform forward and backward selection of explanatory variables 

step.env.Dampier <- ordistep(bacteria.cap0.Dampier, scope=formula(bacteria.cap1.Dampier))

step.env.Dampier$anova

OTU.var.Dampier <- Dampier[1:nrow(Dampier),24:ncol(Dampier)]

#Do the variation partitioning with these variables

bacteria.var.Dampier <- varpart(OTU.var.Dampier, ~Gender, ~Tree, data=ambrosia.env.Dampier)

plot(bacteria.var.Dampier, bg=1:3, Xnames=c('Gender', 'Tree'))

#grep the spatial patterns and transform them into a distance matrix

D.spatial_Ai <- subset(Dampier, select=c("Longitude", "Latidude"))
Dampier.pcnm <- as.data.frame(scores(pcnm(dist(D.spatial_Ai))))
dim(Dampier.pcnm)

ambrosia_Dampier.var <- varpart(OTU.var.Dampier, Dampier.pcnm, ~Tree, data=Dampier)

plot(ambrosia_Dampier.var, bg=1:3, Xnames=c('space', 'Tree'))



############################################################# Dorrigo #######################################################################################################################################
Dorrigo <- subset(bacteria_new, Location=='Dorrigo')
Dorrigo$Raffaelea<-factor(Dorrigo$Raffaelea)
Dorrigo$Host<-factor(Dorrigo$Host)
Dorrigo$Tree<-factor(Dorrigo$Tree)
Dorrigo$Collected<-factor(Dorrigo$Collected)
Dorrigo$Gender<-factor(Dorrigo$Gender)
Dorrigo$Alive_Dead<-factor(Dorrigo$Alive_Dead)


normi_Dorrigo <- decostand(Dorrigo[1:nrow(Dorrigo),24:ncol(Dorrigo)], method="hellinger")
speci_Dorrigo <- wcmdscale(vegdist(normi_Dorrigo, method="bray"), eig = T)
plot(speci_Dorrigo, cex=0.9)
speci_Dorrigo$eig[speci_Dorrigo$eig>=0]/sum(speci_Dorrigo$eig[speci_Dorrigo$eig>=0])

dorrigo_color <- palette(c("red","blue","green","black", "orange", "yellow")) # set up the colour palette
with(Dorrigo, plot(scores(speci_Dorrigo, display='sites'), # identifies coordinates
                   col=dorrigo_color[Dorrigo$Host], #assign symbols, colours
                   xlab='Dim1 18%', ylab="Dim2 14%", cex=0.9,pch=c(1,2,3,4,5,6,7)[Dorrigo$Tree])) # change axis labels
legend('bottomleft', # position the legend
       legend=c('4','6','12'), # set legend text
       pch=c(1,2,3,4,5,6,7), col=c('black'), # assign symbols and colours
       cex=0.7)
title(main = "Bacterial Community Dorrigo National Park")


#Variation Partitioning for Olney
ambrosia.env.Dorrigo <- Dorrigo[, c('Gender', 'Tree', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

# The column collected should be a factor
ambrosia.env.Dorrigo$Collected <- factor(ambrosia.env.Dorrigo$Collected)
# set up full and null models for ordistep
bacteria.cap1.Dorrigo <- capscale(normi_Dorrigo ~ ., data=ambrosia.env.Dorrigo, dist="bray")
bacteria.cap0.Dorrigo <- capscale(normi_Dorrigo ~ 1, data=ambrosia.env.Dorrigo, dist="bray")

summary(ambrosia.env.Dorrigo)

#perform forward and backward selection of explanatory variables 

step.env.Dorrigo <- ordistep(bacteria.cap0.Dorrigo, scope=formula(bacteria.cap1.Dorrigo))

step.env.Dorrigo$anova

OTU.var.Dorrigo <- Dorrigo[1:nrow(Dorrigo),24:ncol(Dorrigo)]

#Do the variation partitioning with these variables

bacteria.var.Dorrigo <- varpart(OTU.var.Dorrigo, ~Gender, ~Tree, data=ambrosia.env.Dorrigo)

plot(bacteria.var.Dorrigo, bg=1:3, Xnames=c('Gender', 'Tree'))

#grep the spatial patterns and transform them into a distance matrix

Dorrigo.spatial_Ai <- subset(Dorrigo, select=c("Longitude", "Latidude"))
Dorrigo.pcnm <- as.data.frame(scores(pcnm(dist(Dorrigo.spatial_Ai, "manhattan"))))
dim(Dorrigo.pcnm)

ambrosia_Dorrigo.var <- varpart(OTU.var.Dorrigo, Dorrigo.pcnm, ~Tree, data=Dorrigo)

plot(ambrosia_Dorrigo.var, bg=1:3, Xnames=c('space', 'Tree'))






######################################################### Termeil ################################################################################################################################################
Termeil <- subset(bacteria_new, Location=='Termeil_State_Forest')
Termeil$Raffaelea<-factor(Termeil$Raffaelea)
Termeil$Host<-factor(Termeil$Host)
Termeil$Tree<-factor(Termeil$Tree)
Termeil$Collected<-factor(Termeil$Collected)
Termeil$Gender<-factor(Termeil$Gender)
Termeil$Alive_Dead<-factor(Termeil$Alive_Dead)


normi_Termeil <- decostand(Termeil[1:nrow(Termeil),24:ncol(Termeil)], method="hellinger")
speci_Termeil <- wcmdscale(vegdist(normi_Termeil, method="bray"), eig = T)
plot(speci_Termeil, cex=0.9)
speci_Termeil$eig[speci_Termeil$eig>=0]/sum(speci_Termeil$eig[speci_Termeil$eig>=0])

termeil_color <- palette(c("red","blue","green","black", "orange", "yellow")) # set up the colour palette
with(Termeil, plot(scores(speci_Termeil, display='sites'), # identifies coordinates
                   col=termeil_color[Termeil$Host], #assign symbols, colours
                   xlab='Dim1 22%', ylab="Dim2 13%", cex=0.9,pch=c(1,2,3,4,5,6,7)[Termeil$Tree])) # change axis labels
levels(Termeil$Tree)

legend('bottomleft', # position the legend
       legend=c('4','5','6','9'), # set legend text
       pch=c(1,2,3,4,5,6,7), col=c('black'), # assign symbols and colours
       cex=0.7)

title(main = "Bacterial Community Termeil State Forest")

###################################################################### Wombeyan Caves ##################################################################################################################################################
Wombeyan <- subset(bacteria_new, Location=='Wombeyan_Caves')
Wombeyan$Raffaelea<-factor(Wombeyan$Raffaelea)
Wombeyan$Host<-factor(Wombeyan$Host)
Wombeyan$Tree<-factor(Wombeyan$Tree)
Wombeyan$Collected<-factor(Wombeyan$Collected)
Wombeyan$Gender<-factor(Wombeyan$Gender)
Wombeyan$Alive_Dead<-factor(Wombeyan$Alive_Dead)


normi_Wombeyan <- decostand(Wombeyan[1:nrow(Wombeyan),24:ncol(Wombeyan)], method="hellinger")
speci_Wombeyan <- wcmdscale(vegdist(normi_Wombeyan, method="bray"), eig = T)
plot(speci_Wombeyan, cex=0.9)
speci_Wombeyan$eig[speci_Wombeyan$eig>=0]/sum(speci_Wombeyan$eig[speci_Wombeyan$eig>=0])

wombeyan_color <- palette(c("red","blue","green","black", "orange", "yellow")) # set up the colour palette
with(Wombeyan, plot(scores(speci_Wombeyan, display='sites'), # identifies coordinates
                    col=wombeyan_color[Wombeyan$Host], #assign symbols, colours
                    xlab='Dim1 21%', ylab="Dim2 11%", cex=0.9,pch=c(1,2,3,4,5,6,7)[Wombeyan$Tree])) # change axis labels
levels(Wombeyan$Tree)
levels(Wombeyan$Host)
legend('topleft', # position the legend
       legend=c('Eucalyptus sp.A', 'Eucalyptus sp.B', 'Eucalyptus sp.C'), # set legend text
       pch=c(19,19), col=olney_color, # assign symbols and colours
       cex=0.7)

legend('bottomleft', # position the legend
       legend=c('2','6','9','15','17'), # set legend text
       pch=c(1,2,3,4,5,6,7), col=c('black'), # assign symbols and colours
       cex=0.7)

title(main = "Bacterial Community Wombeyan Caves National Park")


###################################################################### Mount Wilson ##################################################################################################################################################
MTWilson <- subset(bacteria_new, Location=='Mount_Wilson')
MTWilson$Raffaelea<-factor(MTWilson$Raffaelea)
MTWilson$Host<-factor(MTWilson$Host)
MTWilson$Tree<-factor(MTWilson$Tree)
MTWilson$Collected<-factor(MTWilson$Collected)
MTWilson$Gender<-factor(MTWilson$Gender)
MTWilson$Alive_Dead<-factor(MTWilson$Alive_Dead)


normi_MTWilson <- decostand(MTWilson[1:nrow(MTWilson),24:ncol(MTWilson)], method="hellinger")
speci_MTWilson <- wcmdscale(vegdist(normi_MTWilson, method="bray"), eig = T)
plot(speci_MTWilson, cex=0.9)
speci_MTWilson$eig[speci_MTWilson$eig>=0]/sum(speci_MTWilson$eig[speci_MTWilson$eig>=0])

mtwilson_color <- palette(c("red","blue","green","black", "orange", "yellow")) # set up the colour palette
with(MTWilson, plot(scores(speci_MTWilson, display='sites'), # identifies coordinates
                    col=mtwilson_color[MTWilson$Host], #assign symbols, colours
                    xlab='Dim1 19%', ylab="Dim2 14%", cex=0.9,pch=c(1,2,3,4,5,6,7)[MTWilson$Tree])) # change axis labels
levels(MTWilson$Tree)
levels(MTWilson$Host)
legend('bottomright', # position the legend
       legend=c('Eucalyptus sp.A', 'Eucalyptus sp.B'), # set legend text
       pch=c(19,19), col=olney_color, # assign symbols and colours
       cex=0.7)

legend('bottomleft', # position the legend
       legend=c('1','2','5','6','10'), # set legend text
       pch=c(1,2,3,4,5,6,7), col=c('black'), # assign symbols and colours
       cex=0.7)

title(main = "Bacterial Community  in Mount Wilson")


################################################################################################################################################################################################################

#Shared OTUs between sites
library(VennDiagram)
library(gplots)


Dorrigo <- colnames(OTU[map$Location == "Dorrigo", apply(OTU[map$Location == "Dorrigo",], MARGIN=2, function(x) any(x >0))])
Mt_Wilson <- colnames(OTU[map$Location == "Mount_Wilson", apply(OTU[map$Location == "Mount_Wilson",], MARGIN=2, function(x) any(x >0))])
Olney <- colnames(OTU[map$Location == "Olney", apply(OTU[map$Location == "Olney",], MARGIN=2, function(x) any(x >0))])
Wombeyan <- colnames(OTU[map$Location == "Wombeyan_Caves", apply(OTU[map$Location == "Wombeyan_Caves",], MARGIN=2, function(x) any(x >0))])
Dampier <- colnames(OTU[map$Location == "Dampier_State_Forest", apply(OTU[map$Location == "Dampier_State_Forest",], MARGIN=2, function(x) any(x >0))])
Termeil <- colnames(OTU[map$Location == "Termeil_State_Forest", apply(OTU[map$Location == "Termeil_State_Forest",], MARGIN=2, function(x) any(x >0))])

venn(list(Dampier,Termeil,Wombeyan,Olney,Dorrigo))


#shared between males and females
Male <- colnames(OTU[map$Gender == "Male", apply(OTU[map$Gender == "Male",], MARGIN=2, function(x) any(x >0))])
Female <- colnames(OTU[map$Gender == "Female", apply(OTU[map$Gender == "Female",], MARGIN=2, function(x) any(x >0))])
venn(list(Male,Female))
#shared OTUS between 5 Sites
tmp <- venn(list(Dampier,Termeil,Wombeyan,Olney,Dorrigo))
shared_otus <- attr(tmp, "intersections")

shared_otus
#Get the OTU_names of the different OTUS

index <- match(Female,Male) # solve 1
result1 <- Male[na.omit(index)]

index <- match(Male,Female) # solve 2
result2 <- Female[na.omit(index)]

unique(result1) # or unique(result2) to solve 3


library(venn)
library(tidyverse)
library(stringr)
library(seqinr)


p_th = 0.0;
Dorrigo <- Dorrigo[-1]
Dampier <- Dampier[-1]
Mt_Wilson <- Mt_Wilson[-1]
Termeil <- Termeil[-1]
Olney <- Olney[-1]
Wombeyan <- Wombeyan[-1]
venn = list(Dampier,Termeil,Wombeyan,Olney,Dorrigo,Mt_Wilson)

png("shared_otus_bacteria.png", width = 800, height = 800)

venn.result = venn(venn, snames= c("Dampier", "Termeil", "Wombeyan", "Olney", "Dorrigo", "MT Wilson"), ilabels = TRUE, zcolor = "style", size = 30, cexil = 1.2, cexsn = 1.5)


dev.off()


shared_otus <- attr(venn.result, "intersections")
shared_otus


#get the sequence associated with the indicator OTUs
fas<-read.fasta('/home/robert/Projects/bacterial_microbiome_Ai/scripts/forward/Betadiversity/dna-sequences.fasta')

getSequence(fas['94f21cef0d3a4a2d1e3c8195df4e64f3'], as.string=T) #Order:	Clostridiales
getSequence(fas['9a068d169f0e8f73fe497e704a2b33a7'], as.string=T) #Order:	Sphingomonadales
getSequence(fas['1a39283a8c55f0345152dbc6f2a9a49c'], as.string=T) #Order:	Rhodospirillales
getSequence(fas['5f1b27e4462bc484dd71c1cae42e1826'], as.string=T) #Order:	Burkholderiales #Delftia
getSequence(fas['81dc89dc68af7c7c44d48619f05e70bf'], as.string=T) #Order:	Pseudomonadales
getSequence(fas['344139c7e835ef9511e2d09ccb57867c'], as.string=T) #Order:	Rhodospirillales Acidosoma
getSequence(fas['dc20f8a8af90ba051146e220d06f30ea'], as.string=T) #Order:	Actinomycetales Human Skin Bacteria
getSequence(fas['c7d7102e32d7fd200b282dd544f8be31'], as.string=T) #Order:	Burkholderiales
getSequence(fas['468cab0014f51c350de2580488285bfd'], as.string=T) #Order:	Caulobacterales
getSequence(fas['7d034a80e0a082a542ca50c6d4610059'], as.string=T) #Order:	Burkholderiales
getSequence(fas['873b374057eaf7772ec9a52bbed28af8'], as.string=T) #Quambalaria cyanescens 18S with 16S???????
getSequence(fas['8294fa60adbc14fa938fc39ab1028312'], as.string=T) #Order:	Bacillales
getSequence(fas['2115edb65e8a540a36ebb61172ef37c8'], as.string=T) #Order:	Sphingomonadales
getSequence(fas['0c2932d673807c869c28ef7db7a3fd39'], as.string=T) #Order:	Lactobacillales #Enterococcus
getSequence(fas['5dce03bd2e28eae95acce6cb4daf00f2'], as.string=T) #Order:	Burkholderiales #Cupriavidus
getSequence(fas['be8893edd073b50adb8dbf6406924c5c'], as.string=T) #Order:	Burkholderiales #Delftia
getSequence(fas['31046b33084361c5ae8c75570c067ea1'], as.string=T) #Order:	Burkholderiales #Herbaspirillum sp.
getSequence(fas['9fd72f05128fb87f78b47060670a7083'], as.string=T) #Order:	Enterobacteriales #Escherichia coli
getSequence(fas['39390a7684558dbefc5805272550351b'], as.string=T) #Order:	Burkholderiales #Delftia
getSequence(fas['53454143ac65968ed9605b733ea9ca43'], as.string=T) #Order:	Acidobacteriales #Granulicella
getSequence(fas['bcd7cc25c448519c1320a55b93dcef6d'], as.string=T) #Order:	Rhodospirillales Acidosoma
getSequence(fas['755a96f0a40d76a95dfde6b698dc61ac'], as.string=T) #Order:	Burkholderiales #Ralstonia
getSequence(fas['c52f0e2c597e7ec17af1ca1a6ee142ca'], as.string=T) #Order:	Bacillales #Staphylococcus epidermidis
getSequence(fas['a15c6941507535f32afd3425a017aec7'], as.string=T) #Order:	Actinomycetales #Corynebacterium
getSequence(fas['d61ecd1fa655fabcbe4d4a064e70703f'], as.string=T) #Order:	Clostridiales





























