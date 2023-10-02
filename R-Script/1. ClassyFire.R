library(data.table)
library(dplyr)
library(ChemmineR)

#Set Variables
setwd("G:\\Mi unidad\\CICESE\\Tesis\\2. Analysis\\2. Results\\2. FBMN\\0. Enrichment\\1. GNPS\\1. PreResults\\ClassyFire")
WD <- getwd()

smiles <- read.csv(paste0(WD, "//0. Source//SMILES.csv"))

#Make TSV to analyze in Classyfire
write.table(smiles, file = ".//0. Source//SMILES_ClassyFire.tsv", 
            row.names=FALSE, sep="\t",col.names=TRUE) 

#When the classyfire analysis finishes, download the sdf file
#https://www.bioconductor.org/packages/devel/bioc/vignettes/ChemmineR/inst/doc/ChemmineR.html#Export_of_Compounds

#Read ClassyFire result in sdf format
sdfset <- read.SDFset((paste0(WD, "//0. Source//11227002.sdf")))

blockmatrix <- datablock2ma(datablocklist=datablock(sdfset))
metID <- sdfid(sdfset)
blockmatrix <- cbind(metID,blockmatrix)
blockmatrix <- as.data.frame(blockmatrix[-1,])

#Insert ID to dataframe
Superclass <- sub(" __ > <Superclass>","",blockmatrix$Superclass)
blockmatrix$Superclass <- Superclass

Class <- sub(" __ > <Class>","",blockmatrix$Class)
blockmatrix$Class <- Class

Subclass <- sub(" __ > <Subclass>","",blockmatrix$Subclass)
blockmatrix$Subclass <- Subclass

Classyfire_result <- blockmatrix %>% 
select(metID, Superclass, Class, Subclass)
rownames(Classyfire_result) <- NULL	
Classyfire_result <- merge(Classyfire_result,smiles, by ="metID")

write.csv(Classyfire_result, file = "ClassyFire.csv", 
          row.names = FALSE, col.names = TRUE)
