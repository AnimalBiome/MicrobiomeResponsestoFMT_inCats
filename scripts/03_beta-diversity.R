################################################################################
#
#               Fecal microbiome responses to FMT in Cats
#                      
#       Rojas et al 2023. Microbiome responses to fecal microbiota 
#             transplantation in cats with chronic digestive issues.
#
#                     Code Created By: Connie A.Rojas
#                     Created On: 02 Oct 2021
#                     Last updated: 28 April 2023
#
################################################################################

## CODE FOR: 
#       A) testing microbiome beta-diversity ~ host predictors (PERMANOVA)
#       B) testing microbiome dispersion ~ host predictors (PERMDISP)
#       C) PCoA ordinations of microbiome variation ~ host predictors


source(file="scripts/00_configure.R"); #set up environment


################################################################################
#             1. Load ASV table, ASV phylogenetic tree, and sample metadata                 
################################################################################  

load("data/01_cats_meta.Rdata");
load("data/01_ASV_table.Rdata");

# retain samples from FMT recipients 
meta2=meta[meta$group=="recipient",];
before=meta2[meta2$type=="before",];
after=meta2[meta2$type=="after",];

# subset to before samples, after samples, or both [PICK ONE]
asv=asvdf[,colnames(asvdf) %in% before$sampleID];  
asv=asvdf[,colnames(asvdf) %in% after$sampleID]; 


# remove singleton/doubleton ASVs
asv=asv[rowSums(asv)>2,];
asv=asv[,order(colnames(asv))];


################################################################################
#             2. Calculate Bray-Curtis and Weighted Unifrac distances 
#                           and run PERMANOVAS
################################################################################

# get distances
bray<-apply(asv, 2, function(i) (i/sum(i))*100);
ps1<- phyloseq(otu_table(bray, taxa_are_rows=TRUE));
bray.dist=phyloseq::distance(ps1, method="bray");

ps2=phyloseq(otu_table(asv,taxa_are_rows=TRUE));
ps2=microbiome::transform(ps2,"clr",target="OTU");
clr.dist=phyloseq::distance(ps2, method="euclidean");

# microbiome beta-diversity ~ host predictors
mydist=list(bray.dist, clr.dist);
names=c("Bray-Curtis","Aitchison");
met=c("bray","euclidian"); 

for(i in 1:2)
{
  print(paste("PERMANOVA ~ host predictors, N=46",names[i]));
  print(adonis2(mydist[[i]]~
                  response2+
                  antibiotics+  
                  symptom+
                  dry_food+
                  round(age_yrs)+
                  sex,
                data=after,   #before or after
                method = met[i],
                by="margin",
                permutations = 999));
};
# library("pairwiseAdonis");
for(i in 1:3)
{
  print(pairwise.adonis(mydist[[i]],
                        as.factor(after$symptom), #before or after
                        sim.method = met[i],
                        p.adjust.m = "fdr"));
};




################################################################################
#             3. run PERMDISP tests of microbiome dispersion
################################################################################

# clinical signs
for(i in 1:2)
{
  print(paste("PERMDISP test, N=67",names[i]));
  bdisper<-with(after,    # before or after
                betadisper(mydist[[i]],symptom));
  print(anova(bdisper));
  #print(TukeyHSD(bdisper));
};

# dry kibble consumption
for(i in 1:2)
{
  print(paste("PERMDISP test, N=67",names[i]));
  bdisper<-with(after,    # before or after
                betadisper(mydist[[i]],dry_food));
  print(anova(bdisper));
  #print(TukeyHSD(bdisper));
};



################################################################################
#             4. generate PC coordinates for PCoA
################################################################################

# calculate coordinates for PCoA
pcoa_dec=cmdscale(clr.dist, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "sampleID");
pcoa_met=merge(pcoa,after,by="sampleID");

# calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

# set order offactors
pcoa_met$symptom=factor(pcoa_met$symptom, levels=c("Diarrhea","VomDiarr",
                                               "VomConstip","Constipation"));
pcoa_met$dry_food=factor(pcoa_met$dry_food,levels=c("Yes","No"));


################################################################################
#             5. plot PCoAs!
################################################################################

# PCoA color-coded by clinical signs
mycol=c("#1f78b4","#ffd92f","palegreen","#f781bf"); 
pcoa1=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=symptom),  
             size = 2.7,
             pch=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="PostFMT - Clinical Signs")+   
  scale_fill_manual(values=mycol)+  
  stat_ellipse(geom = "polygon", alpha =0.12, aes(fill = symptom))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", linewidth =2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa1);

# PCoA color-coded by diet category
mycol=c("#1b9e77", "#d95f02", "darkorchid2","steelblue1","#e6ab02","#b3e2cd");
pcoa2=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=dry_food),  
             size = 2.7,
             pch=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="PostFMT - Dry Kibble")+   
  scale_fill_manual(values=mycol)+  
  stat_ellipse(geom = "polygon", alpha =0.12, aes(fill = dry_food))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", linewidth =2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa2);


################################################################################
#             6. save plots
################################################################################

mps=arrangeGrob(pcoa1,pcoa2,nrow=1);
ggsave(filename="03_PCoAs_hostpredictors_aitchison.pdf",
       device="pdf",path="./figures",
       plot=mps,
       width=11,
       height=4,
       units="in",
       dpi=500);

