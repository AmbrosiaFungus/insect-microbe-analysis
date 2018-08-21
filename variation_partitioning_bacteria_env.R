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



##########################################################################################################################################################################
######################################################## Variationpartitioning ###########################################################################################

#OTU table merged with the metaData in fungi_new

#normalize the data with decostand

normi <- decostand(bacteria_new[1:nrow(bacteria_new),24:ncol(bacteria_new)], method="hellinger")

#Get the environmental variable I thnk are a key factor

ambrosia.env <- bacteria_new[, c('Location', 'Gender', 'Tree', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

#the tree level has cooccuring patterns between sites so fix that

ambrosia.env$Tree <- factor(with(ambrosia.env, paste(Location, Tree, sep ='_')))
ambrosia.env <- droplevels(ambrosia.env)


# The column collected should be a factor
ambrosia.env$Collected <- factor(ambrosia.env$Collected)
# set up full and null models for ordistep
bacteria.cap1 <- capscale(normi ~ ., data=ambrosia.env, dist="bray")
bacteria.cap0 <- capscale(normi ~ 1, data=ambrosia.env, dist="bray")

summary(ambrosia.env)

#perform forward and backward selection of explanatory variables 

step.env <- ordistep(bacteria.cap0, scope=formula(bacteria.cap1))

step.env$anova

#exclude tree

new.ambrosia.env <- ambrosia.env[, c('Location','Gender', 'Host', 'Raffaelea', 'Collected', 'Alive_Dead')]
new.ambrosia.env <- droplevels(new.ambrosia.env)
summary(new.ambrosia.env)

# set up full and null models for ordistep
bacteria.cap1.new <- capscale(normi ~ ., data=new.ambrosia.env, dist="bray")
bacteria.cap0.new <- capscale(normi ~ 1, data=new.ambrosia.env, dist="bray")


#perform forward and backward selection of explanatory variables 

step.env.new <- ordistep(bacteria.cap0.new, scope=formula(bacteria.cap1.new))

step.env.new$anova


plot(step.env.new)

#Select the OTU table

OTU.var <- bacteria_new[1:nrow(bacteria_new),24:ncol(bacteria_new)]

#Do the variation partitioning with these variables

bacteria.var <- varpart(OTU.var, ~Location, ~Alive_Dead, ~Host, ~Collected, data=new.ambrosia.env)

plot(bacteria.var, bg=1:3, Xnames=c('Location', 'Alive Dead', 'Host', 'Collected'))

title(main = "Variation for Bacteria explained")


#include the distance of trees in the variation partitioning

spatial_Ai <- subset(bacteria_new, select=c("Longitude", "Latidude"))

ambrosia.pcnm <- as.data.frame(scores(pcnm(dist(spatial_Ai))))
dim(ambrosia.pcnm)


ambrosia.var <- varpart(OTU.var, ambrosia.pcnm, ~Location, ~Alive_Dead, ~Host, data=ambrosia.env)

plot(ambrosia.var, bg=1:3, Xnames=c('space', 'Location', 'Alive or Dead', 'Host'))


#test for the significance of phylogeny 

sig.rda.Location <- rda(OTU.var, new.ambrosia.env$Location, cbind(new.ambrosia.env$Alive_Dead, new.ambrosia.env$Host))

anova(sig.rda.Location)

plot(sig.rda.Location)

################################################### Variation partitioning and Spatial Analysis on each site ####################################################################################################

##################################################  Olney State Forest  #################################################################################################

#subset the table for the specific site, in this case Olney

Olney <- subset(bacteria_new, Location=='Olney')

#grep the spatial patterns and transform them into a distance matrix

O.spatial_Ai <- subset(Olney, select=c("Longitude", "Latidude"))
Olney.pcnm <- as.data.frame(scores(pcnm(dist(O.spatial_Ai))))
dim(Olney.pcnm)

#Get only the OTUs

OTU.spatial_olney <- Olney[1:nrow(Olney),24:ncol(Olney)]

#normalize

normi.olney <- decostand(OTU.spatial_olney, method="hellinger")

#Get the variable I am interested in 
ambrosia.env.Olney <- Olney[, c('Gender', 'Tree', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

#the tree level has cooccuring patterns between sites so fix that

ambrosia.env.Olney <- droplevels(ambrosia.env.Olney)

# The column collected should be a factor
ambrosia.env.Olney$Collected <- factor(ambrosia.env.Olney$Collected)
ambrosia.env.Olney$Tree <- factor(ambrosia.env.Olney$Tree)

# set up full and null models for ordistep
bacteria.olney.cap1 <- capscale(normi.olney ~ ., data=ambrosia.env.Olney, dist="bray")
bacteria.olney.cap0 <- capscale(normi.olney ~ 1, data=ambrosia.env.Olney, dist="bray")

summary(ambrosia.env.Olney)

#perform forward and backward selection of explanatory variables 

step.olney.env <- ordistep(bacteria.olney.cap0, scope=formula(bacteria.olney.cap1))

step.olney.env$anova


ambrosia_Olney.var <- varpart(OTU.spatial_olney, Olney.pcnm, ~Tree, data=Olney)

plot(ambrosia_Olney.var, bg=1:3, Xnames=c('space', 'Tree'))

#Redundancy analysis for the site. 

ambrosia_Olney.rda <- rda(OTU.spatial_olney ~ Tree + Condition(Olney.pcnm[, 1]) + Condition(Olney.pcnm[, 2]) + Condition(Olney.pcnm[, 3]), data=Olney)
anova(ambrosia_Olney.rda)


###############################################     Dampier State Forest             ######################################################################

#subset the table for the specific site, in this case Dampier

Dampier <- subset(bacteria_new, Location=='Dampier_State_Forest')

#grep the spatial patterns and transform them into a distance matrix

Dampier.spatial_Ai <- subset(Dampier, select=c("Longitude", "Latidude"))
Dampier.pcnm <- as.data.frame(scores(pcnm(dist(Dampier.spatial_Ai))))
dim(Dampier.pcnm)

#Get only the OTUs

OTU.spatial_dampier <- Dampier[1:nrow(Dampier),24:ncol(Dampier)]

#normalize

normi.dampier <- decostand(OTU.spatial_dampier, method="hellinger")

#Get the variable I am interested in 
ambrosia.env.Dampier <- Dampier[, c('Gender', 'Tree', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

#the tree level has cooccuring patterns between sites so fix that

ambrosia.env.Dampier <- droplevels(ambrosia.env.Dampier)

# The column collected should be a factor
ambrosia.env.Dampier$Collected <- factor(ambrosia.env.Dampier$Collected)
ambrosia.env.Dampier$Tree <- factor(ambrosia.env.Dampier$Tree)

# set up full and null models for ordistep
bacteria.dampier.cap1 <- capscale(normi.dampier ~ ., data=ambrosia.env.Dampier, dist="bray")
bacteria.dampier.cap0 <- capscale(normi.dampier ~ 1, data=ambrosia.env.Dampier, dist="bray")

summary(ambrosia.env.Dampier)

#perform forward and backward selection of explanatory variables 

step.dampier.env <- ordistep(bacteria.dampier.cap0, scope=formula(bacteria.dampier.cap1))

step.dampier.env$anova


ambrosia_dampier.var <- varpart(OTU.spatial_dampier, Dampier.pcnm, ~Tree, data=Dampier)

plot(ambrosia_dampier.var, bg=1:3, Xnames=c('space', 'Tree'))

#Redundancy analysis for the site. 

ambrosia_dampier.rda <- rda(OTU.spatial_dampier ~ Tree + Condition(Dampier.pcnm[, 1]) + Condition(Dampier.pcnm[, 2]) + Condition(Dampier.pcnm[, 3]), data=Dampier)
anova(ambrosia_dampier.rda)


###############################################     Dorrigo National Park            ######################################################################

#subset the table for the specific site, in this case Dorrigo

dorrigo <- subset(bacteria_new, Location=='Dorrigo')

#grep the spatial patterns and transform them into a distance matrix

dorrigo.spatial_Ai <- subset(dorrigo, select=c("Longitude", "Latidude"))
dorrigo.pcnm <- as.data.frame(scores(pcnm(dist(dorrigo.spatial_Ai))))
dim(dorrigo.pcnm)

#Get only the OTUs

OTU.spatial_dorrigo <- dorrigo[1:nrow(dorrigo),24:ncol(dorrigo)]

#normalize

normi.dorrigo <- decostand(OTU.spatial_dorrigo, method="hellinger")

#Get the variable I am interested in 
ambrosia.env.dorrigo <- dorrigo[, c('Gender', 'Tree', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

ambrosia.env.dorrigo <- droplevels(ambrosia.env.dorrigo)

# The column collected should be a factor
ambrosia.env.dorrigo$Collected <- factor(ambrosia.env.dorrigo$Collected)
ambrosia.env.dorrigo$Tree <- factor(ambrosia.env.dorrigo$Tree)

# set up full and null models for ordistep
bacteria.dorrigo.cap1 <- capscale(normi.dorrigo ~ ., data=ambrosia.env.dorrigo, dist="bray")
bacteria.dorrigo.cap0 <- capscale(normi.dorrigo ~ 1, data=ambrosia.env.dorrigo, dist="bray")

summary(ambrosia.env.dorrigo)

#perform forward and backward selection of explanatory variables 

step.dorrigo.env <- ordistep(bacteria.dorrigo.cap0, scope=formula(bacteria.dorrigo.cap1))

step.dorrigo.env$anova


ambrosia_dorrigo.var <- varpart(OTU.spatial_dorrigo, dorrigo.pcnm, ~Gender, data=dorrigo)

plot(ambrosia_dorrigo.var, bg=1:3, Xnames=c('space', 'Gender'))

#Redundancy analysis for the site. 

ambrosia_dorrigo.rda <- rda(OTU.spatial_dorrigo, ambrosia_dorrigo.var$Gender, OTU.spatial_dorrigo, data=dorrigo)
anova(ambrosia_dorrigo.rda)



###############################################     Wombeyan National Park            ######################################################################

#subset the table for the specific site, in this case Wombeyan

wombeyan <- subset(bacteria_new, Location=='Wombeyan_Caves')

#grep the spatial patterns and transform them into a distance matrix

wombeyan.spatial_Ai <- subset(wombeyan, select=c("Longitude", "Latidude"))
wombeyan.pcnm <- as.data.frame(scores(pcnm(dist(wombeyan.spatial_Ai))))
dim(wombeyan.pcnm)

#Get only the OTUs

OTU.spatial_wombeyan <- wombeyan[1:nrow(wombeyan),24:ncol(wombeyan)]

#normalize

normi.wombeyan <- decostand(OTU.spatial_wombeyan, method="hellinger")

#Get the variable I am interested in 
ambrosia.env.wombeyan <- wombeyan[, c('Gender', 'Tree', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

ambrosia.env.wombeyan <- droplevels(ambrosia.env.wombeyan)

# The column collected should be a factor
ambrosia.env.wombeyan$Collected <- factor(ambrosia.env.wombeyan$Collected)
ambrosia.env.wombeyan$Tree <- factor(ambrosia.env.wombeyan$Tree)

# set up full and null models for ordistep
bacteria.wombeyan.cap1 <- capscale(normi.wombeyan ~ ., data=ambrosia.env.wombeyan, dist="bray")
bacteria.wombeyan.cap0 <- capscale(normi.wombeyan ~ 1, data=ambrosia.env.wombeyan, dist="bray")

summary(ambrosia.env.wombeyan)

#perform forward and backward selection of explanatory variables 

step.wombeyan.env <- ordistep(bacteria.wombeyan.cap0, scope=formula(bacteria.wombeyan.cap1))

step.wombeyan.env$anova


ambrosia_wombeyan.var <- varpart(OTU.spatial_wombeyan, wombeyan.pcnm, ~Host, ~Raffaelea, data=wombeyan)

plot(ambrosia_wombeyan.var, bg=1:3, Xnames=c('space', 'Host', 'Raffaelea'))

#Redundancy analysis for the site. 

ambrosia_wombeyan.rda <- rda(OTU.spatial_wombeyan, ambrosia_wombeyan.var$Gender, OTU.spatial_wombeyan, data=wombeyan)
anova(ambrosia_wombeyan.rda)


###############################################     Mount Wilson            ######################################################################

#subset the table for the specific site, in this case Mount Wilson

wilson <- subset(bacteria_new, Location=='Mount_Wilson')

#grep the spatial patterns and transform them into a distance matrix

wilson.spatial_Ai <- subset(wilson, select=c("Longitude", "Latidude"))
wilson.pcnm <- as.data.frame(scores(pcnm(dist(wilson.spatial_Ai))))
dim(wilson.pcnm)

#Get only the OTUs

OTU.spatial_wilson <- wilson[1:nrow(wilson),24:ncol(wilson)]

#normalize

normi.wilson <- decostand(OTU.spatial_wilson, method="hellinger")

#Get the variable I am interested in 
ambrosia.env.wilson <- wilson[, c('Gender', 'Tree', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

ambrosia.env.wilson <- droplevels(ambrosia.env.wilson)

# The column collected should be a factor
ambrosia.env.wilson$Collected <- factor(ambrosia.env.wilson$Collected)
ambrosia.env.wilson$Tree <- factor(ambrosia.env.wilson$Tree)

# set up full and null models for ordistep
bacteria.wilson.cap1 <- capscale(normi.wilson ~ ., data=ambrosia.env.wilson, dist="bray")
bacteria.wilson.cap0 <- capscale(normi.wilson ~ 1, data=ambrosia.env.wilson, dist="bray")

summary(ambrosia.env.wilson)

#perform forward and backward selection of explanatory variables 

step.wilson.env <- ordistep(bacteria.wilson.cap0, scope=formula(bacteria.wilson.cap1))

step.wilson.env$anova


###############################################     Termeil State Forest             ######################################################################

#subset the table for the specific site, in this case Termeil State Forest

termeil <- subset(bacteria_new, Location=='Termeil_State_Forest')

#grep the spatial patterns and transform them into a distance matrix

termeil.spatial_Ai <- subset(termeil, select=c("Longitude", "Latidude"))
termeil.pcnm <- as.data.frame(scores(pcnm(dist(termeil.spatial_Ai))))
dim(termeil.pcnm)

#Get only the OTUs

OTU.spatial_termeil <- termeil[1:nrow(termeil),24:ncol(termeil)]

#normalize

normi.termeil <- decostand(OTU.spatial_termeil, method="hellinger")

#Get the variable I am interested in 
ambrosia.env.termeil <- termeil[, c('Gender', 'Tree', 'Host', 'Raffaelea','Collected', 'Alive_Dead')]

ambrosia.env.termeil <- droplevels(ambrosia.env.termeil)

# The column collected should be a factor
ambrosia.env.termeil$Collected <- factor(ambrosia.env.termeil$Collected)
ambrosia.env.termeil$Tree <- factor(ambrosia.env.termeil$Tree)

# set up full and null models for ordistep
bacteria.termeil.cap1 <- capscale(normi.termeil ~ ., data=ambrosia.env.termeil, dist="bray")
bacteria.termeil.cap0 <- capscale(normi.termeil ~ 1, data=ambrosia.env.termeil, dist="bray")

summary(ambrosia.env.termeil)

#perform forward and backward selection of explanatory variables 

step.termeil.env <- ordistep(bacteria.termeil.cap0, scope=formula(bacteria.termeil.cap1))

step.termeil.env$anova

ambrosia_termeil.var <- varpart(OTU.spatial_termeil, termeil.pcnm, ~Tree, data=termeil)

plot(ambrosia_termeil.var, bg=1:3, Xnames=c('space', 'Tree'))

#Redundancy analysis for the site. 

ambrosia_termeil.rda <- rda(OTU.spatial_wombeyan, ambrosia_wombeyan.var$Tree, OTU.spatial_wombeyan, data=wombeyan)
anova(ambrosia_termeil.rda)




