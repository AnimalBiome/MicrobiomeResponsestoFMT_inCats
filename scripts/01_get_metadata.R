################################################################################
#
#               Fecal microbiome responses to FMT in Cats
#                      
#       Rojas et al 2022. Microbiome responses to fecal microbiota 
#             transplantation in cats with chronic digestive issues.
#
#                     Code Created By: Connie A.Rojas
#                     Created On: 22 Sept 2021
#                     Last updated: 6 December 2022
#
################################################################################

## CODE FOR: 
#       A) setting sample metadata factors; 
#       B) saving ASV abundance table after removing ASVs classified as Unknown / Chloroplast

source(file="scripts/00_configure.R"); #set up environment


################################################################################
#             1.  Load sample metadata, ASV table, and ASV taxonomy             
################################################################################

#load metadata
meta=read.csv("data/00_cats_metadata.csv", stringsAsFactors = F);

#load ASV table output by DADA2 R tutorial 
#https://benjjneb.github.io/dada2/tutorial.html
load("data/00_seqtab_dada2.Rdata");

#load ASV taxonomy output by DADA2 R tutorial
#silva v. 138
tax=read.csv("data/00_ASV_taxonomy_silva138.csv", header=T,row.names=1); 


################################################################################
#             2.  Tidy ASV table           
################################################################################

#transpose ASV table so that ASVs are rows
asvdf=as.data.frame(seqtab.nochim);
asvdf=as.data.frame(t(asvdf));

#remove duplicate samples from ASV table
remo=c("9QVBFPd","TLDR72d","KS28KNd","TND32Wd");
asvdf=asvdf[,!colnames(asvdf) %in% remo];

#remove ASVs classified as Chloroplast / Unknown Bacteria
asvdf=asvdf[rownames(asvdf) %in% rownames(tax),];

#remove ASVs with 0 counts across samples
asvdf=asvdf[rowSums(asvdf) > 0,];


################################################################################
#             3.  Set metadata factors                
################################################################################

meta$sex<-factor(meta$sex, levels=c("Female","Male"));
meta$IBD=factor(meta$IBD, levels=c("Yes","No","Maybe"));
meta$antibiotics=factor(meta$antibiotics, levels=c("Yes","No"));

meta$symptom=factor(meta$symptom, levels=c("Diarrhea","Constipation",
                                                   "VomDiarr","Other"));

meta$diet_combo=factor(meta$diet_combo, levels=c("Dry","Wet","Raw","Wet.Dry",
                                                         "Wet.Raw","Other"));

meta$response2=factor(meta$response2, levels=c("Responder","Non-Responder"));

meta=meta[order(meta$sampleID),];


################################################################################
#             4. save formatted metadata file, and ASV table
################################################################################

save(meta, file="data/01_cats_meta.Rdata");
save(asvdf, file="data/01_ASV_table.Rdata");