library(metacoder)
library(vegan)

#read in the taxonomy file created by Qiime2
tax <- read.table("/home/robert/Projects/fungal_microbiom_Ai/scripts/97_otu/biom-taxonomy.tsv", sep='\t', header=TRUE, stringsAsFactors = F)

#read OTU table
otu <- read.csv("/home/robert/Projects/fungal_microbiom_Ai/scripts/97_otu/OTU_table_97.csv", header = T)

#change column name
names(tax)[1]<-"OTUID"


# R transforms some OTUs and put a X infront of them so we have to change that in the taxonomy too
x <- grep('^[0-9]', tax$OTUID)
tax[x, 'OTUID'] <- paste("X", tax[x, 'OTUID'], sep='')

all(names(otu[, -1]) %in% tax$OTUID)

#Have to transpose the otu-table

trans_otu <- t(otu)

colnames(trans_otu) = trans_otu[1, ] # the first row will be the header

trans_otu = trans_otu[-1, ]          # removing the first row.

#extract the taxon from tax and add to OTU table
new <- as.data.frame(tax[,c("Taxon")])

k <- cbind(trans_otu,new)

names(k)

#change a column name

colnames(k)[159] <- "lineage"


#remove Unassigned

k<-k[!(k$lineage=="Unassigned"),]

k <- k[,c(159,1:ncol(k))]

#remove last column
k <- k[-160]   # same as above



obj <- parse_tax_data(k, class_cols = "lineage", class_sep = ";",class_key = c(tax_rank = "info", tax_name = "taxon_name"),class_regex = "^(.+)__(.+)$")

print(obj)

obj$data$tax_data <- zero_low_counts(obj, "tax_data", min_count = 5)


no_reads <- rowSums(obj$data$tax_data[, hmp_samples$sample_id]) == 0
sum(no_reads)
















heat_tree(obj)


