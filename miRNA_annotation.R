# Clear Global Environment
rm(list=ls())

# Set working directory
setwd(".../miRNA_Pathvisio_Lisa/")

#import data from file
file <- file.path(getwd(),"TNFA_acute/TNFA_acute_Lisa.csv")
dat <- read.delim(file, header = TRUE, sep = ",", dec = ".")


# Load needed libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", "miRBaseConverter")

library(dplyr)
library(miRBaseConverter)

# Split columns
# Split name column into firstname and last name
library(tidyr)
dat2 <- dat %>% separate(X, c("MiRNA", "Second", "Other"), extra = "drop", sep = "/")

# Download accessions
data(miRNATest)
Accessions = miRNATest$Accession

# Check the version 
miRNANames = dat2$MiRNA
version=checkMiRNAVersion(miRNANames, verbose = TRUE)

# Convert mature miRNA names to precursor 
result1=miRNA_MatureToPrecursor(miRNANames)
datPrecursor <- cbind(dat2, precursor=result1$Precursor)


# In case of duplicates, remove the least significant ones
dat.ord <- datPrecursor[order(datPrecursor$p01),]
dat.ord2 <- dat.ord %>% distinct(precursor, .keep_all = TRUE)


## Add accessions of precursor 
miRNANames2 = dat.ord2$precursor
result3 = miRNA_NameToAccession(miRNANames2,version = "v22")
dat.ann <- cbind(dat.ord2, accession=result3[,2])

write.table(dat.ann, "DE_list_acute_annotatedmiRBase.txt", sep = "\t", dec = ".", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
