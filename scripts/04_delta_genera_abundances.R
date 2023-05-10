################################################################################
#
#               Fecal microbiome responses to FMT in Cats
#                      
#       Rojas et al 2023. Microbiome responses to fecal microbiota 
#             transplantation in cats with chronic digestive issues.
#
#                     Code Created By: Connie A.Rojas
#                     Created On: 20 Jan 2022
#                     Last updated: 8 May 2023
#
################################################################################

## CODE FOR: 
#       A) testing which bacterial genera correlate with host predictors
#             i) taxa that are part of the core healthy cat microbiome
#            ii) taxa that are potentially pathogenic to cats
#       B) plotting the relative abundances of these genera over time

source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#             1.  Load filtered ASV table, ASV taxonomy, sample metadata, 
#                 table of core taxa, and table of pathogenic taxa
################################################################################

load("data/01_ASV_table.Rdata");
tax=read.csv("data/00_ASV_taxonomy_silva138.csv", header=T,row.names=1); 
load("data/01_cats_meta.Rdata");
core=read.csv("data/00_core_taxa.csv", header=T);
pato=read.csv("data/00_pathogenic_taxa.csv",header=T);

# attach taxonomy to the ASV table 
tax$Genus[(tax$Genus=="Ruminococcus_2")|
            (tax$Genus=="Ruminococcus_1")]="Ruminococcus"
#asv_tax=merge(asv,tax,by="row.names"); for new cats 48
asv_tax=merge(asvdf,tax,by="row.names"); 
colnames(asv_tax)[1]="ASV";

# make a vector of FMT sample IDs
meta2=meta[meta$group=="recipient",]; 
samples=meta2$sampleID;
rownames(meta2)=meta2$sampleID; meta2$sampleID=NULL;


################################################################################
#             2. Calculate the changes in abundance for the 21 core 
#                         bacterial genera (postFMT - preFMT)
################################################################################

# select bacterial taxonomic rank 
abun=asv_tax[,which(names(asv_tax) 
                    %in% c(samples, "Genus"))];
colnames(abun)[ncol(abun)]="taxa";

# calculate taxa relative abundances 
abun=aggregate(.~taxa, abun, sum);  
abun[,-1] <- lapply(abun[,-1], function(x) (x/sum(x))*100);

# subset df to only taxa that form part of core in healthy cats
rownames(abun)=abun$taxa; abun$taxa=NULL;
abun=as.data.frame(t(abun));
abun2=abun[,colnames(abun) %in% core$Genus];

# confirm that all core taxa have a prevalence > 10%
jac=(abun2>0)*1;jac=data.frame(t(jac));
jac=jac[rowSums(jac)>5,]; 
#jac=jac[rowSums(jac)>13,];
nrow(jac)==ncol(abun2);

# append sample metadata
abunmeta=merge(abun2, meta2, by="row.names");
colnames(abunmeta)[1]="sampleID";   

# wrangle df so columns are  CatName | Taxa | Before | After 
abunmeta2=abunmeta[,c(2:22,24,25)]; 
abunmeta2=melt(abunmeta2, id=c("name","type")); 
befores=abunmeta2[abunmeta2$type=="before",];
afters=abunmeta2[abunmeta2$type=="after",];
colnames(befores)[4]="before"; befores$type=NULL;
colnames(afters)[4]="after"; afters$type=NULL;
ab=merge(befores,afters,by=c("name","variable"))
ab$tchange=ab$after-ab$before;

# reattach metadata
ab=ab[,c(1,2,5)];
ab= spread(ab, variable, tchange);
ab2= merge(ab, meta2[meta2$type=="after",], by="name");


################################################################################
#             3. Calculate the changes in abundance for the pathogenic 
#                         bacterial genera (postFMT - preFMT)
################################################################################

# select bacterial taxonomic rank 
abun=asv_tax[,which(names(asv_tax) 
                    %in% c(samples, "Genus"))];
colnames(abun)[ncol(abun)]="taxa";

# calculate taxa relative abundances 
abun=aggregate(.~taxa, abun, sum);  
abun[,-1] <- lapply(abun[,-1], function(x) (x/sum(x))*100);

# subset df to only include potentially pathogenic taxa of cats 
rownames(abun)=abun$taxa; abun$taxa=NULL;
abun=as.data.frame(t(abun));
abun=abun[,colnames(abun) %in% pato$Genus];

# only retain taxa with prevalence > 10%
jac=(abun>0)*1;jac=data.frame(t(jac));
jac=jac[rowSums(jac)>5,]; 
abun3=abun[,colnames(abun)%in%rownames(jac)];

# add sample metadata
abunmeta=merge(abun3, meta2, by="row.names");
colnames(abunmeta)[1]="sampleID";

# wrangle df so columns are  CatName | Taxa | Before | After 
abunmeta2=abunmeta[,c(2:7,9,10)]; 
abunmeta2=melt(abunmeta2, id=c("name","type"));
befores=abunmeta2[abunmeta2$type=="before",];
afters=abunmeta2[abunmeta2$type=="after",];
colnames(befores)[4]="before"; befores$type=NULL;
colnames(afters)[4]="after"; afters$type=NULL;
ab=merge(befores,afters,by=c("name","variable"))
ab$tchange=ab$after - ab$before;

# reattach metadata
ab=ab[,c(1,2,5)];
ab= spread(ab, variable, tchange);
ab3= merge(ab, meta2[meta2$type=="after",], by="name");


################################################################################
#             4. Construct linear models to evaluate whether 
#                 taxa abundances vary with host factors 
################################################################################

# use for loop to test each core bacterial genera [21 in total]

for(i in 2:22)
{
  print(paste("Linear mixed model for:", colnames(ab2)[i]));
  bmod=glm(ab2[,i]~
             factor(ab2$response2)+
             factor(ab2$antibiotics)+
             factor(ab2$symptom)+
             factor(ab2$dry_food),
            data=ab2, na.action=na.omit);
  #print(summary(bmod));
  print("Are host predictors statistically significant?");
  print(Anova(bmod));
};

# quick break for means
aggregate(ab2$Blautia, list(ab2$symptom), FUN=mean);
aggregate(ab2$Parabacteroides, list(ab2$response2), FUN=mean);
aggregate(ab2$Peptococcus, list(ab2$symptom), FUN=mean);
aggregate(ab2$Ruminococcus, list(ab2$response2), FUN=mean);

# use for loop to test each pathogenic bacteria genera [6 for 48 cats]

for(i in 2:7)
{
  print(paste("Linear mixed model for:", colnames(ab3)[i]));
  bmod=glm(ab3[,i]~
             ab3$response2+
             ab3$antibiotics+
             ab3$symptom+
             ab3$dry_food,
            data=ab3, na.action=na.omit);
  #print(summary(bmod));
  print("Are host predictors statistically significant?");
  print(Anova(bmod));
};

# quick break for means
aggregate(ab3$Desulfovibrio, list(ab3$symptom), FUN=mean);
aggregate(ab3$`Escherichia-Shigella`, list(ab3$antibiotics), FUN=mean);
aggregate(ab3$Desulfovibrio, list(ab3$response2), FUN=mean);

# for posthoc testing
ab3$symptom=factor(ab3$symptom);
ab3$response2=factor(ab3$response2);
ab3$antibiotics=factor(ab3$antibiotics);
ab3$dry_food=factor(ab3$dry_food);
mod1=glm(ab3$`Escherichia-Shigella`~response2+antibiotics+symptom+dry_food,
         data=ab3);
summary(glht(mod1, linfct = mcp(antibiotics = "Tukey")));


################################################################################
#             5. Plot the abundances of the bacterial genera that 
#              yielded significant p-values according to the models above
################################################################################

# Plot five genera color-coded by clinical signs
CG=abun2[,colnames(abun2) %in% c("Clostridium ss1","Collinsella", "Negativibacillus",
                                 "Subdoligranulum")];
patho=data.frame(abun3[,1:2]);patho$Clostridioides=NULL; 

CG=merge(CG, patho,by="row.names");  colnames(CG)[1]="sampleID"
CG<-reshape2::melt(CG, id.vars="sampleID",value.name = "abun");
CG=merge(CG,meta,by="sampleID");
CG$type=factor(CG$type, levels=c("before","after"));
CG$symptom=factor(CG$symptom,levels=c("Diarrhea","VomDiarr","VomConstip","Constipation"))

mycol=c("#1f78b4","#ffd92f","palegreen","#f781bf"); 

p1=ggplot(data=CG, 
          mapping=aes(x=type,y=abun))+
  facet_wrap(~variable, scales="free_y",ncol=3)+
  geom_line(aes(group = name,
                color=symptom))+
  scale_color_manual(values=mycol)+
  scale_x_discrete(labels=c("before" = "preFMT", 
                            "after" = "postFMT"),
                   expand = c(0.07, 0.07))+
  theme_bw()+ 
  labs(x = "",
       y="Relative Abun. (%)",
       color="Clinical Signs")+
  theme(legend.title=element_text(size=11,face="bold"),
        text = element_text(size=12),
        legend.position="top",
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size=12),
        legend.margin=margin(t=-10),
        legend.box.spacing = unit(0, "pt"),
        plot.margin = unit(c(0.1,0.5,0.1,0.5), "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"),
        strip.text = element_text(size =11, face="bold"))+
  guides(color = guide_legend(override.aes = list(size = 3))); plot(p1);

# Plot five genera color-coded by to Treatment
CG=abun2[,colnames(abun2) %in% c("sampleID","Butyricicoccus","Megamonas","Peptococcus",
                                 "Ruminococcus")];
patho=data.frame(abun3[,5:6]);patho$Sarcina=NULL; 

CG=merge(CG, patho,by="row.names");  colnames(CG)[1]="sampleID"
CG<-reshape2::melt(CG, id.vars="sampleID",value.name = "abun");
CG=merge(CG,meta,by="sampleID");
CG$type=factor(CG$type, levels=c("before","after"));
CG$response2=factor(CG$response2, levels=c("Responder","Non-Responder"));

mycol=c("#88419d","#FDBF6F");
p2=ggplot(data=CG, 
          mapping=aes(x=type,y=abun))+
  facet_wrap(~variable, scales="free_y",ncol=3)+
  geom_line(aes(group = name,
                color=response2))+
  scale_color_manual(values=mycol)+
  scale_x_discrete(labels=c("before" = "preFMT", 
                            "after" = "postFMT"),
                   expand = c(0.07, 0.07))+
  theme_bw()+ 
  labs(x = "",
       y="Relative Abun. (%)",
       color="Response to Tmt")+
  theme(legend.title=element_text(size=11,face="bold"),
        text = element_text(size=12),
        legend.position="top",
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size=12),
        legend.margin=margin(t=-14),
        legend.box.spacing = unit(0, "pt"),
        plot.margin = unit(c(0.1,0.5,0.1,0.5), "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"),
        strip.text = element_text(size =11, face="bold"))+
  guides(color = guide_legend(override.aes = list(size = 3))); plot(p2);


# Plot three genera color-coded by Abx
CG=abun2[,colnames(abun2) %in% c("sampleID","Negativibacillus","Megamonas")];
CG$Megamonas=NULL;
patho=data.frame(abun3[,2:3]);

CG=merge(CG, patho,by="row.names");  colnames(CG)[1]="sampleID"
CG<-reshape2::melt(CG, id.vars="sampleID",value.name = "abun");
CG=merge(CG,meta,by="sampleID");
CG$type=factor(CG$type, levels=c("before","after"));
CG$antibiotics=factor(CG$antibiotics, levels=c("Yes","No"));

mycol=c("#41b6c4","#f768a1");
p3=ggplot(data=CG, 
          mapping=aes(x=type,y=abun))+
  facet_wrap(~variable, scales="free_y",nrow=1)+
  geom_line(aes(group = name,
                color=antibiotics))+
  scale_color_manual(values=mycol)+
  scale_x_discrete(labels=c("before" = "preFMT", 
                            "after" = "postFMT"),
                   expand = c(0.07, 0.07))+
  theme_bw()+ 
  labs(x = "",
       y="Relative Abun. (%)",
       color="Prior Abx Use")+
  theme(legend.title=element_text(size=11,face="bold"),
        text = element_text(size=12),
        legend.position="top",
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size=12),
        legend.margin=margin(t=-14),
        legend.box.spacing = unit(0, "pt"),
        plot.margin = unit(c(0.1,0.5,0.1,0.5), "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"),
        strip.text = element_text(size =11, face="bold"))+
  guides(color = guide_legend(override.aes = list(size = 3))); plot(p3);


# Plot three genera color-coded by Dry Kibble consumption
CG=abun2[,colnames(abun2) %in% c("sampleID","Peptoclostridium","Subdoligranulum")];
CG$sampleID=rownames(CG);
CG<-reshape2::melt(CG, id.vars="sampleID",value.name = "abun");
CG=merge(CG,meta,by="sampleID");
CG$type=factor(CG$type, levels=c("before","after"));
CG$dry_food=factor(CG$dry_food, levels=c("Yes","No"));

mycol=c("#1b9e77", "#d95f02");
p4=ggplot(data=CG, 
          mapping=aes(x=type,y=abun))+
  facet_wrap(~variable, scales="free_y",nrow=1)+
  geom_line(aes(group = name,
                color=dry_food))+
  scale_color_manual(values=mycol)+
  scale_x_discrete(labels=c("before" = "preFMT", 
                            "after" = "postFMT"),
                   expand = c(0.07, 0.07))+
  theme_bw()+ 
  labs(x = "",
       y="Relative Abun. (%)",
       color="Dry Kibble")+
  theme(legend.title=element_text(size=11,face="bold"),
        text = element_text(size=12),
        legend.position="top",
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size=12),
        legend.margin=margin(t=-14),
        legend.box.spacing = unit(0, "pt"),
        plot.margin = unit(c(0.1,0.5,0.1,0.5), "cm"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"),
        strip.text = element_text(size =11, face="bold"))+
  guides(color = guide_legend(override.aes = list(size = 3))); plot(p4);


################################################################################
#             6. save plots!
################################################################################

mypts=arrangeGrob(p1,p2,p3,p4, nrow=4,
                  heights = c(1/3, 1/3, 1/4.6, 1/4.6),
                  layout_matrix=rbind(c(1,1,1,1),c(2,2,2,2),
                                      c(3,3,3,3),c(4,4,4,NA)));

ggsave(filename="04_RelAbunTaxa_Plots.pdf",
       device="pdf",path="./figures",
       plot=mypts,
       width=7,
       height=12,
       units="in",
       dpi=500);

