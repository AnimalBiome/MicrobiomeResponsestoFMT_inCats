################################################################################
#
#               Fecal microbiome responses to FMT in Cats
#                      
#       Rojas et al 2022. Microbiome responses to fecal microbiota 
#             transplantation in cats with chronic digestive issues.
#
#                     Code Created By: Connie A.Rojas
#                     Created On: 15 Feb 2022
#                     Last updated: 4 November 2022
#
################################################################################

## CODE FOR: 
#       A) calculating how similar FMT samples are to the healthy reference set
#       B) correlating those microbiome similarity values ~ host predictors
#       C) plotting those microbiome similarity values ~ host predictors

source(file="scripts/00_configure.R"); #set up environment


################################################################################
#             1. Load ASV table, ASV phylogenetic tree, and sample metadata                 
################################################################################  

load("data/01_cats_meta.Rdata");
load("data/01_ASV_table.Rdata");
load("data/02_ASV_phylotree.Rdata");

# retain samples from FMT recipients and healthy animals
meta2=meta[meta$group!="donor",];
asv=asvdf[,colnames(asvdf) %in% meta2$sampleID];  

# remove singleton/doubleton ASVs
asv=asv[rowSums(asvdf)>2,];
asv=asv[,order(colnames(asv))];


################################################################################
#             2. Generate Bray-Curtis and Weighted Unifrac distances 
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

# clr (center log ratio) for Aitchison distance ("compositions")
clra=data.frame(clr(asv));
ps3<- phyloseq(otu_table(clra, taxa_are_rows=TRUE));
clr.dist=phyloseq::distance(ps3, method="euclidean");


################################################################################
#             3. how similar are the fecal microbiomes of FMT recipients to
#                   those of healthy animals after treatment?
################################################################################

#  melt dissimilarity distance and append sample metadata
df1=reshape2::melt(as.matrix(clr.dist));   # bray.dist OR wun.dist
colnames(df1)[1]="sampleID";
df1=inner_join(df1, meta2[,c(1,3,4,10,16,20,23:25,30,31)], by="sampleID"); 
colnames(df1)=c("sam1","sampleID","val","type1","name1","bod1","antibiotics1",
                "symptom1","wet1","raw1","dry1","age1","response1");
df1=inner_join(df1, meta2[,c(1,3,4,10,16,20,23:25,30,31)], by="sampleID");

# retain FMTrecipient -- age-matched healthy animal comparisons
# the two cats must be less than 1 year apart in their ages
df1$dissim=df1$val;
df2=df1[df1$type=="healthy",];
df2=df2[df2$sam1!=df2$sampleID,];
df2=df2[df2$type1!="healthy",];
df2$agediff=abs(df2$age1-df2$age_yrs);
df2=df2[df2$agediff<=1.1,]; #1.05 before

# new diet controlled
df2=df2[df2$dry1==df2$dry_food,];
#df2=df2[df2$raw1==df2$raw_food,];

# new variable that saves the names of the two cats being compared
df2$comp=paste(df2$name1,"-",df2$name);

# calculate postFMT similarity - preFMT similarity (aka DELTA similarity)
dfA=df2[df2$type1=="after",];
dfB=df2[df2$type1=="before",];
mdf=merge(dfA,dfB[,c(24,26)],by="comp");
#mdf$sdelta=mdf$simi.x-mdf$simi.y; 
mdf$sdelta=(-1)*(mdf$dissim.x-mdf$dissim.y);

# run linear model to test DELTA similarity ~ host predictors
full.model<- glm(sdelta ~ factor(response1)+factor(antibiotics1)+
                   factor(symptom1)+factor(dry1),
                 data = mdf,na.action = "na.omit");

Anova(full.model);

# get mean similarity for a group
aggregate(mdf$sdelta, list(mdf$dry1), FUN=mean);


################################################################################
#             4. Plots of delta similarity ~ host predictors
#                            Means and SE
###############################################################################

# delta similarity ~ response to treatment
b=aggregate(sdelta~response1,mdf,mean); colnames(b)[2]="mean_abun";   
c=aggregate(sdelta~response1, mdf, sd); colnames(c)[2]="sd_abun";     
d=merge(b,c,by="response1");                                          

d$size=table(mdf$response1);                                          
d$err=d$sd_abun/sqrt(d$size);
d$upper=d$mean_abun+ d$err;
d$lower=d$mean_abun-d$err;

p1=ggplot(d,aes(x=response1, y=mean_abun,ymax=0.01))+    
  geom_point(size=2, color="#fc9272")+
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    colour="#fc9272", width = 0)) +
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "#737373", size=0.5)+
  theme_classic()+
  scale_x_discrete(labels=c("Responder" = "Resp.", 
                            "Non-Responder" = "Non-Resp."))+
  labs(y="Mean Δ Similarity to Healthy",
       x="")+                                                                         
  theme(legend.title=element_blank(),
        text = element_text(size=12),
        legend.position="none",
        plot.title = element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold")); plot(p1); 

  
# delta similarity ~ clinical signs
b=aggregate(sdelta~symptom1,mdf,mean); colnames(b)[2]="mean_abun";   
c=aggregate(sdelta~symptom1, mdf, sd); colnames(c)[2]="sd_abun";     
d=merge(b,c,by="symptom1");                                          

d$size=table(mdf$symptom1);                                          
d$err=d$sd_abun/sqrt(d$size);
d$upper=d$mean_abun+ d$err;
d$lower=d$mean_abun-d$err;

d$symptom1=factor(d$symptom1,levels=c("Diarrhea","VomDiarr",
                                      "VomConstip","Constipation"));

p2=ggplot(d,aes(x=symptom1, y=mean_abun, ymin=-0.06))+    
  geom_point(size=2, color="#fc9272")+
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    colour="#fc9272", width = 0)) +
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "#737373", size=0.5)+
  theme_classic()+
  labs(y="Mean Δ Similarity to Healthy",
       x="")+                                                                         
  theme(legend.title=element_blank(),
        text = element_text(size=12),
        legend.position="none",
        plot.title = element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold")); plot(p2); 

# delta similarity ~ antibiotic use
b=aggregate(sdelta~antibiotics1,mdf,mean); colnames(b)[2]="mean_abun";   
c=aggregate(sdelta~antibiotics1, mdf, sd); colnames(c)[2]="sd_abun";     
d=merge(b,c,by="antibiotics1");                                          

d$size=table(mdf$antibiotics1);                                          
d$err=d$sd_abun/sqrt(d$size);
d$upper=d$mean_abun+ d$err;
d$lower=d$mean_abun-d$err;

d$antibiotics1=factor(d$antibiotics1,levels=c("Yes","No"));

p4=ggplot(d,aes(x=antibiotics1, y=mean_abun, ymin=-0.03))+    
  geom_point(size=2, color="#fc9272")+
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    colour="#fc9272", width = 0)) +
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "#737373", size=0.5)+
  theme_classic()+
  labs(y="Mean Δ Similarity to Healthy",
       x="")+                                                                         
  theme(legend.title=element_blank(),
        text = element_text(size=12),
        legend.position="none",
        plot.title = element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold")); plot(p4); 

# delta similarity ~ diet category
b=aggregate(sdelta~dry1,mdf,mean); colnames(b)[2]="mean_abun";   
c=aggregate(sdelta~dry1, mdf, sd); colnames(c)[2]="sd_abun";     
d=merge(b,c,by="dry1");                                          

d$size=table(mdf$dry1);                                          
d$err=d$sd_abun/sqrt(d$size);
d$upper=d$mean_abun+ d$err;
d$lower=d$mean_abun-d$err;

d$dry1=factor(d$dry1,levels=c("Yes","No"));

p5=ggplot(d,aes(x=dry1, y=mean_abun, ymax=0.02))+    
  geom_point(size=2, color="#fc9272")+
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    colour="#fc9272", width = 0)) +
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "#737373", size=0.5)+
  theme_classic()+
  labs(y="Mean Δ Similarity to Healthy",
       x="")+                                                                         
  theme(legend.title=element_blank(),
        text = element_text(size=12),
        legend.position="none",
        plot.title = element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold")); plot(p5); 


################################################################################
#             5. save plots
###############################################################################

myps=grid.arrange(p1,p2,p4,p5,
                    layout_matrix = rbind(c(1,1,2,2,2),
                                          c(4,4,5,5,NA)));

ggsave(filename="06_deltaSimilarity_hostfactors_aitchison.pdf",
       device=cairo_pdf,path="./figures",
       plot=myps,
       width=9,
       height=7,
       units="in",
       dpi=500);