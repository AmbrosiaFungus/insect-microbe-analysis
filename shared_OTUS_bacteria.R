#Shared OTUs between sites
library(VennDiagram)
library(gplots)

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
