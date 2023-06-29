################################################################################
#
#               Fecal microbiome responses to FMT in Cats
#                      
#       Rojas et al 2023. Microbiome responses to fecal microbiota 
#             transplantation in cats with chronic digestive issues.
#
#                     Code Created By: Connie A.Rojas
#                     Created On: 22 Sept 2021
#                     Last updated: 28 April 2023
#
################################################################################

## CODE FOR: 
#       A) calculating microbiome alpha-diversity
#       B) running liner models on microbiome alpha-diversity
#       C) plotting microbiome alpha-diversity (boxplots)

source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#             1. Load filtered ASV abundance table and 
#                       sample metadata                 
################################################################################

load("data/01_ASV_table.Rdata");
load("data/01_cats_meta.Rdata");


################################################################################
#             2. Create phylogenetic tree of ASV sequences for
#                 calculating Faith's phylogenetic diversity
################################################################################

asvf=t(asvdf)
seqs <- getSequences(as.matrix(asvf));
names(seqs) <- seqs;
alignment <- AlignSeqs(DNAStringSet(seqs),
                       anchor=NA,verbose=FALSE);
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA");
dm <- dist.ml(phangAlign);
treeNJ <- NJ(dm);
fit = pml(treeNJ, data=phangAlign);
fitGTR <- update(fit, k=4, inv=0.2);
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE,
                    optGamma=TRUE, rearrangement = "stochastic",
                    control = pml.control(trace = 0));

# save the output as an .Rdata file
save(fitGTR, file="data/02_ASV_phylotree.Rdata");
load("data/02_ASV_phylotree.Rdata");


################################################################################
#             3. Calculate microbiome alpha-diversity 
################################################################################

# save sampleIDs because we will need these to restore them later
sampleID=colnames(asvdf);
noXs=gsub("X","",sampleID);
snames=data.frame(sampleID, noXs);

# only keep samples from FMT recipients 
fmt=meta$sampleID[meta$group=="recipient"];
asvdf=asvdf[,colnames(asvdf) %in% fmt];

# remove ASVs that do not appear in dataset
asvdf=asvdf[rowSums(asvdf)>=1,];

# make a phyloseq object and calculate richness
ps<- phyloseq(otu_table(asvdf, taxa_are_rows=TRUE));
alpha=estimate_richness(ps,split = TRUE, 
                        measures = c("Observed","Chao1","Shannon"));

# restore sampleIDs
rownames(alpha) <- gsub("X","",row.names(alpha))
rownames(alpha)=snames$sampleID[
  match(rownames(alpha),snames$noXs)];

# calculate Pielou's evenness and append to alpha df
pev=evenness(ps, 'pielou')
alpha=merge(alpha, pev, by="row.names", all=TRUE); 
rownames(alpha)=alpha$Row.names; alpha$Row.names=NULL;

# calculate Faith's phylogenetic diversity and append to alpha df
faith=pd(t(asvdf), phy_tree(fitGTR$tree), include.root=F);
alpha=merge(alpha, faith, by="row.names", all=TRUE); 
colnames(alpha)[1]="sampleID";


################################################################################
#             4. Run linear models on microbiome alpha-diversity
################################################################################

# append sample metadata to alphadiv table
am=inner_join(alpha, meta, by="sampleID");

# LMM microbiome alpha-div ~ host predictors
bdf=am[am$type=="before",];
adf=am[am$type=="after",];

mymetrics=c("Observed", "pielou","PD");

for(i in 1:3)
{
  print(paste("LMM preFMT ~ host predictors for:", mymetrics[i]));
  bmod=lmer(bdf[,mymetrics[i]]~
              bdf$response2+
              bdf$antibiotics+
              bdf$symptom+
              bdf$dry_food+
              (1|round(bdf$age_yrs))+
              (1|bdf$sex),
            na.action="na.exclude");
  print(Anova(bmod));
};

for(i in 1:3)
{
  print(paste("LMM postFMT ~ host predictors for:", mymetrics[i]));
  amod=lmer(adf[,mymetrics[i]]~
              adf$response2+
              adf$antibiotics+
              adf$symptom+
              adf$dry_food+
              (1|round(adf$age_yrs))+
              (1|adf$sex),
            na.action="na.exclude");
  print(Anova(amod));
};

# post-hoc comparisons
mod1=lmer(pielou~response2+antibiotics+symptom+dry_food+
            (1|round(bdf$age_yrs)), data=bdf);
summary(glht(mod1, linfct = mcp(dry_food = "Tukey")))


################################################################################
#             5. boxplots of microbiome alpha-diversity
################################################################################

# preFMT Pielou's evenness ~ Clinical signs
bdf$symptom=factor(bdf$symptom, levels=c("Diarrhea","VomDiarr",
                                         "VomConstip","Constipation"));
box1=ggplot(data=bdf, 
          mapping=aes(x=symptom,y=pielou))+
  geom_boxplot(lwd=0.7)+
  theme_bw()+ 
  labs(x = "Clinical signs",
       y = "Pielou's Evenness",  
       title="Pre-FMT")+
  theme(axis.title = element_text(size = 13, face="bold"),
        plot.title=element_text(size = 13, face="bold"),
        axis.text = element_text(size = 11.4),
        legend.position="none",
        panel.background = element_rect(colour = "black", 
                                        linewidth=2)); plot(box1);

# preFMT Pielou's evenness ~ dry kibble consumption
bdf$dry_food=factor(bdf$dry_food,levels=c("Yes","No"));
box2=ggplot(data=bdf, 
          mapping=aes(x=dry_food,y=pielou))+
  geom_boxplot(lwd=0.7)+
  theme_bw()+ 
  labs(x = "Dry Kibble",
       y = "Pielou's Evenness",  
       title="Pre-FMT")+
  theme(axis.title = element_text(size = 13, face="bold"), 
        plot.title=element_text(size = 13, face="bold"),
        axis.text = element_text(size = 13),
        legend.position="none",
        panel.background = element_rect(colour = "black", 
                                        linewidth=2)); plot(box2);

# save plots
ggsave(filename="02_alphadiv_clinicalSigns.pdf",
       device="pdf",path="./figures",
       plot=box1,
       width=4.5, 
       height=4,
       units="in",
       dpi=500);

ggsave(filename="02_alphadiv_diet.pdf",
       device="pdf",path="./figures",
       plot=box2,
       width=4.5, 
       height=4,
       units="in",
       dpi=500);


