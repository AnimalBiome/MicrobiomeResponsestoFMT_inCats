################################################################################
#
#               Fecal microbiome responses to FMT in Cats
#                      
#       Rojas et al 2022. Microbiome responses to fecal microbiota 
#             transplantation in cats with chronic digestive issues.
#
#                     Code Created By: Connie A.Rojas
#                     Created On: 20 October 2022
#                     Last updated: 8 May 2023
#
################################################################################

## CODE FOR: 
#       A) determining which ASVs "engrafted" in the FMT recipient
#       B) making plots showing taxonomy of engrafted ASVs

source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#             1. ASV table, ASV taxonomy, and sample metadata                
################################################################################  
load("data/01_cats_meta.Rdata");
load("data/01_ASV_table.Rdata");
tax=read.csv("data/00_ASV_taxonomy_silva138.csv", header=T,row.names=1); 

# convert ASV table to presence/absence
adf=merge(asvdf,tax[,7:8],by="row.names");
rownames(adf)=adf$asv_species;
adf$Row.names=NULL; adf$Species=NULL; adf$asv_species=NULL;
jac=(adf>0)*1;

# subset data frame to keep samples from stool donors and FMT recipients
meta=meta[meta$group!="healthy",];
meta$donor=as.character(meta$donor);


################################################################################
#           2. see which ASVs engrafted in each individual cat and
#                     calculate ASV engrafment rate 
################################################################################
# ASV engraftment rate = 
#         number of ASVs shared between postFMT samples and their stool donors 
#         (excluding taxa shared between preFMT samples and donors) DIVIDED BY
#         ----------------------------------------------------------------
#         number of ASVs found in stool donor's sample(s)
#         (excluding taxa shared between preFMT samples and donors) 

mycats=unique(meta$name[meta$group=="recipient"]);
taxlab=c();
for(i in 1:length(mycats))
{
  #shared between preFMT and donor aka "pre-shared"
  dname=meta$donor[meta$name==mycats[i]];
  dID=meta$sampleID[(meta$donor==dname)&(meta$group=="donor")];
  ipre=meta$sampleID[(meta$name==mycats[i])&(meta$type=="before")];
  ipost=meta$sampleID[(meta$name==mycats[i])&(meta$type=="after")];
  
  graf1=data.frame(jac[,colnames(jac) %in% c(dID,ipre)],check.names=F);
  colnames(graf1)=meta$group[match(colnames(graf1),meta$sampleID)];
  graf1=graf1[graf1$recipient>0,];
  graf1=graf1[rowSums(graf1)>=2,]; 
  
  # shared between postFMT and donor (subtract pre-shared)
  graf2=data.frame(jac[,colnames(jac)%in% c(dID,ipost)],check.names=F);
  colnames(graf2)=meta$group[match(colnames(graf2),meta$sampleID)];
  graf2=graf2[graf2$recipient>0,];
  graf2=graf2[rowSums(graf2)>=2,]; 
  graf2=graf2[!rownames(graf2)%in% rownames(graf1),];
  taxlab=c(taxlab, rownames(graf2));  #saves taxonomic labels of ASVs
  
  # donor bacteria that could engraft (subtract preshared)
  graf3=data.frame(jac[,colnames(jac)%in% c(dID,ipost)],check.names=F);
  colnames(graf3)=meta$group[match(colnames(graf3),meta$sampleID)];
  graf3$recipient=NULL;
  graf3$rcounts=rowSums(graf3);graf3=graf3[graf3$rcounts>0,];
  graf3=graf3[!rownames(graf3)%in% rownames(graf1),];
  
  # ASV engraftment rate
  print(mycats[i]);
  #print(nrow(graf3));
  #print(nrow(graf2))
  print((nrow(graf2)/nrow(graf3))*100);
};

# I typed up the values shown in the console and saved file as "05_engraftment_rate"


################################################################################
#         3. Linear model of ASV engraftment rate ~ host predictors                
################################################################################

erate=read.csv("data/05_engraftment_rate.csv",header=T);
edf=merge(erate,meta[meta$type=="before",],by="name");

# engrafment rate ~ host predictors
mod1=lmer(edf$engraf.rate~edf$response2+ edf$antibiotics+
            edf$symptom+edf$dry_food+
            (1|round(edf$age_yrs))+
            (1|edf$sex),
          na.action="na.exclude");
Anova(mod1);

# engrafment rate ~ donor name
kruskal.test(edf$engraf.rate~edf$donor);


################################################################################
#         4. Plots of ASV engraftment rates                
################################################################################

# bar plot of ASV engrafment rates for the 68 FMT recipients
erate=erate[erate$name %in% edf$name,];
ep1=ggplot(erate,aes(reorder(name,engraf.rate), y=engraf.rate))+
  geom_bar(position="dodge",stat="identity",linewidth=1.2, color="#525252",fill="#bdbdbd")+ 
  coord_cartesian(ylim=c(0,40))+
  labs(y="ASV Engraftment Rate (%)",
       x="FMT Recipients")+
  theme_classic() + 
  theme(axis.title = element_text(size = 14, face="bold"), 
        axis.text = element_text(size = 16),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none"); plot(ep1);

# boxplot of ASV engrafment rate ~ donor
mycol=c("#02818a",'#7570b3',"#CAB2D6",
          '#e78ac3','#fc9272',"black",'#FDBF6F',"darkgoldenrod");
ep2=ggplot(data=edf, 
          mapping=aes(x=donor,y=engraf.rate,fill=donor))+
  geom_boxplot(lwd=0.7)+
  coord_cartesian(ylim=c(0,30))+
  scale_x_discrete(labels=c("D1*","D2","D3","D4*","D5","D6","D7","D8"))+
  theme_classic()+ 
  scale_fill_manual(values=mycol)+
  labs(x = "Donor ID",
       y = "ASV Engraftment Rate in Recipients")+
  theme(axis.title = element_text(size = 14, face="bold"), 
        axis.text = element_text(size = 16),
        legend.position="none"); plot(ep2);
aggregate(edf$engraf.rate, list(edf$donor), FUN=mean);


################################################################################
#         5. Plot of the taxonomy of engrafted ASVs               
################################################################################

# combine the taxonomic labels of ASVs that engrafted across FMT recipients
# and collapse by Genus
mydf=data.frame(table(taxlab)); 
mydfb=as.data.frame(str_split_fixed(mydf$taxlab, "\\|",2));
mydf$taxo=mydfb$V2;
mydf$taxlab=NULL;
genus=aggregate(.~taxo, mydf, sum);  
genus$abun=(genus$Freq/sum(genus$Freq))*100;

# plot
genus$dummy="cats";
genus=genus[genus$Freq>15,];

b1=ggplot(genus, aes(y=abun, fill=dummy,reorder(taxo,abun)))+
  geom_bar(position="dodge",stat="identity",
           color="black")+
  coord_flip() +
  labs(y="Percentage of Engrafted ASVs",
       x="",
       fill="")+  
  scale_fill_manual(values="#fa9fb5")+  
  theme_bw()+
  theme(legend.position="none",
        legend.text=element_text(size=13),
        legend.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=14,face="bold"),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=14,face="bold"),
        axis.title.y=element_text(size=13, face="bold"));plot(b1);

## save plots
combp=arrangeGrob(ep1,ep2,b1, nrow=2,
                  layout_matrix=rbind(c(1,1,2,2),c(3,3,3,3)));

ggsave(filename="figures/05_ASV_engraftment_plots.pdf",
       device="pdf",path=".",
       plot=combp,
       width=10,
       height=10,
       units="in",
       dpi=500);

