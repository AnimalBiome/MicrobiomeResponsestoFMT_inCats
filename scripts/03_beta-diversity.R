################################################################################
#
#               Fecal microbiome responses to FMT in Cats
#                      
#       Rojas et al 2022. Microbiome responses to fecal microbiota 
#             transplantation in cats with chronic digestive issues.
#
#                     Code Created By: Connie A.Rojas
#                     Created On: 02 Oct 2021
#                     Last updated: 4 November 2022
#
################################################################################

## CODE FOR: 
#       A) testing microbiome beta-diversity ~ host predictors (PERMANOVA)
#       B) testing microbiome dispersion ~ host predictors (PERMDISP)


source(file="scripts/00_configure.R"); #set up environment


################################################################################
#             1. Load ASV table, ASV phylogenetic tree, and sample metadata                 
################################################################################  

load("data/01_cats_meta.Rdata");
load("data/01_ASV_table.Rdata");
load("data/02_ASV_phylotree.Rdata");

# retain samples from FMT recipients 
meta2=meta[meta$group=="recipient",];
meta2=meta2[complete.cases(meta2$antibiotics),];
before=meta2[meta2$type=="before",];
after=meta2[meta2$type=="after",];

# subset to before samples, after samples, or both [PICK ONE]
asv=asvdf[,colnames(asvdf) %in% before$sampleID];  
asv=asvdf[,colnames(asvdf) %in% after$sampleID]; 
asv=asvdf[,colnames(asvdf) %in% meta2$sampleID]; 

# remove singleton/doubleton ASVs
asv=asv[rowSums(asvdf)>2,];
asv=asv[,order(colnames(asv))];


################################################################################
#             2. Calculate Bray-Curtis and Weighted Unifrac distances 
#                           and run PERMANOVAS
################################################################################

# get distances
bray<-apply(asv, 2, function(i) (i/sum(i))*100);

ps1<- phyloseq(otu_table(bray, taxa_are_rows=TRUE));
bray.dist=phyloseq::distance(ps1, method="bray");

ps2<- phyloseq(otu_table(asv, taxa_are_rows=TRUE),
               phy_tree(fitGTR$tree));
set.seed(711);
phy_tree(ps2) <- root(phy_tree(ps2), sample(taxa_names(ps2), 1), resolve.root = TRUE); 
wun.dist=phyloseq::distance(ps2, method="wunifrac");

# microbiome beta-diversity ~ host predictors
mydist=list(bray.dist,wun.dist);
names=c("Bray-Curtis","WUNIFRAC");
met=c("bray","euclidian"); 

for(i in 1:2)
{
  print(paste("PERMANOVA ~ host predictors, N=67",names[i]));
  print(adonis2(mydist[[i]]~     
                  IBD+
                  symptom+
                  antibiotics+
                  response2+
                  diet_combo+
                  age_yrs+
                  sex,
                data=after,   # before or after 
                method = met[i],
                by="terms",
                permutations = 999));
};

# preFMT vs. postFMT
meta2$name[(meta2$sampleID=="2HH288")|
             (meta2$sampleID=="NHHXZN")]="Tabbytha2";
for(i in 1:2)
{
  print(paste("PERMANOVA preFMT vs postFMT, N=67, using:", names[i]));
  print(adonis2(mydist[[i]]~
                  type+
                  name,
                data=meta2,
                method = met[i],
                by="terms",
                permutations = 999));
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

# symtoms
for(i in 1:2)
{
  print(paste("PERMDISP test, N=67",names[i]));
  bdisper<-with(after,    # before or after
                betadisper(mydist[[i]],diet_combo));
  print(anova(bdisper));
  #print(TukeyHSD(bdisper));
};

# IBD
for(i in 1:2)
{
  print(paste("PERMDISP test, N=67",names[i]));
  bdisper<-with(after,   # before or after
                betadisper(mydist[[i]],IBD));
  print(anova(bdisper));
  #print(TukeyHSD(bdisper));
};


################################################################################
#             4. generate PC coordinates for PCoA
################################################################################

# calculate coordinates for PCoA
pcoa_dec=cmdscale(bray.dist, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "sampleID");
pcoa_met=merge(pcoa,after,by="sampleID");

# calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);


################################################################################
#             5. plot PCoAs!
################################################################################

# PCoA color-coded by IBD
mycol=c("#fa9fb5","#000000","lightgoldenrod") #for IBD
pcoa1=ggplot(pcoa_met, aes(Axis1, Axis2))+
geom_point(mapping=aes(fill=IBD),  
           size = 2.5,
           pch=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="PostFMT - IBD")+   
  scale_fill_manual(values=mycol)+  
  stat_ellipse(geom = "polygon", alpha =0.12, aes(fill = IBD))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa1);

# PCoA color-coded by clinical symptoms
mycol=c("#1f78b4","#ffd92f","palegreen","#f781bf"); 
pcoa2=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=symptom),  
             size = 2.5,
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
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa2);

# PCoA color-coded by diet category
mycol=c("#1b9e77", "#d95f02", "darkorchid2","steelblue1","#e6ab02","#b3e2cd");
pcoa3=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=diet_combo),  
             size = 2.5,
             pch=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="PostFMT - Diet")+   
  scale_fill_manual(values=mycol)+  
  stat_ellipse(geom = "polygon", alpha =0.12, aes(fill = diet_combo))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa3);

# barplot showing % variance explained by each host predictor
cats=c("Diet","ClinicalSigns","IBD","Antibiotics","ResptoTmt","Age","Sex");
cats=c(cats,cats); stype=c(rep("PreFMT",7), rep("PostFMT",7));
vals1=c(11,10,0.1,0.1,0.1,0.1,0.1,11,12,5,0.1,0.1,0.1,0.1);
mdf=data.frame(cats, vals1,stype);

mdf$cats=factor(mdf$cats, levels=c("Diet","ClinicalSigns","IBD",
                                   "Antibiotics","ResptoTmt","Age","Sex"));
mdf$stype=factor(mdf$stype, levels=c("PreFMT","PostFMT"));

barcol=c("#bebada", "#7fcdbb");
b1=ggplot(mdf, aes(fill=stype, y=vals1, x=cats))+
  geom_bar(position="dodge",stat="identity",
           color="black")+
  labs(y="% Variance Explained",
       x="",
       fill="")+  
  scale_fill_manual(values=barcol)+ 
  theme_bw()+
  theme(legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        axis.text.x=element_text(size=13, angle=45, vjust = 0.5),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(b1);


################################################################################
#             6. save plots
################################################################################

mps=arrangeGrob(b1,pcoa1,pcoa2,pcoa3,nrow=2);
ggsave(filename="03_PCoAs_hostpredictors.pdf",
       device="pdf",path="./figures",
       plot=mps,
       width=12.5,
       height=9,
       units="in",
       dpi=500);