################################################################################
#
#               Fecal microbiome responses to FMT in Cats
#                      
#       Rojas et al 2023. Microbiome responses to fecal microbiota 
#             transplantation in cats with chronic digestive issues.
#
#                     Code Created By: Connie A.Rojas
#                     Created On: 20 April 2023
#                     Last updated: 28 April 2023
#
################################################################################

## CODE FOR: plotting microbiome composition for 8 FMT participants

source(file="scripts/00_configure.R"); #set up environment


################################################################################
#             1. Load ASV table, ASV taxonomy, and sample metadata                 
################################################################################  

load("data/01_ASV_table.Rdata");
tax=read.csv("data/00_ASV_taxonomy_silva138.csv", header=T,row.names=1); 
load("data/01_cats_meta.Rdata");

# retain samples from FMT recipients
meta2=meta[meta$group=="recipient",];
asv=asvdf[,colnames(asvdf) %in% meta2$sampleID]; 
asv=asv[rowSums(asv) > 0,];

# merge ASV table with ASV taxonomy
asv2=merge(asv,tax,by="row.names");
asv2$Row.names=NULL;


################################################################################
#             2. Make plots of bacterial genus abundances                 
###############################################################################

# delete unwanted columns
dgen=asv2[,c(1:92,98)];
colnames(dgen)[ncol(dgen)]="taxa";
dgen=aggregate(.~taxa,dgen, sum);
dgen=dgen[rowSums(dgen[-1])>0,];

# calculate proportions
dgen[,-1] <- lapply(dgen[,-1], function(x) (x/sum(x))*100);
print(colSums(dgen[-1]));
dgen$AVG=rowMeans(dgen[,-1]); 

# keep taxa >1.6% relative abundance across samples
dgen=dgen[dgen$AVG>1.655,]; 
dgen$AVG=NULL;

# denote the rest of phyla as "Other"
newrow=c(NA, 100-colSums(dgen[2:ncol(dgen)])); 
dgen=rbind(dgen, newrow); 
dgen$taxa=as.character(dgen$taxa);
dgen[nrow(dgen),1]="Other";

# melt data frame for ggplot
genbar<-reshape2::melt(dgen, id.vars="taxa",value.name = "abun");
colnames(genbar)[2]="sampleID";
genbar=merge(genbar, meta2, by="sampleID");

# set up color palette
fam_col1=c("#771155", "#AA4488", "#CC99BB", "#114477","#77AADD", 
           "#117777", "#88CCAA",
           "#777711",  
           "#DDDD77","#DDAA77","#bdbdbd","#774411", "#AA7744", "#AA4455", "grey39", "black");

fam_col2=c("#771155", "#AA4488", "#CC99BB", "#114477","#77AADD", 
           "#117777", "#88CCAA","#774411",
           "#777711",  
           "#DDDD77","#DDAA77","#bdbdbd", "#AA7744", "#AA4455", "grey39", "black");


# subset df to pick eight random cats
want=c("Chelsea","Leonidas","Goldy","Ender",
       "Morrigan","Spanky",
       "Poe","Cinders", "Archie",
       "Frankie","Misha","Tenga");
genbar=genbar[genbar$name %in% want,];

# anonimize cat names
genbar$name=factor(genbar$name,levels=c("Chelsea","Leonidas","Goldy","Ender",
                                        "Morrigan","Spanky",
                                        "Poe","Cinders", "Archie",
                                        "Frankie","Misha","Tenga"));
supp.labs <- c("Cat1","Cat2","Cat3","Cat4","Cat5","Cat6","Cat7","Cat8",
               "Cat9","Cat10","Cat11","Cat12");
names(supp.labs) <- levels(genbar$name);


# plot
barp3=ggplot(data=genbar, 
             aes(type,y=abun, fill=taxa))+
  geom_bar(stat = "identity")+
  facet_wrap(~name,scales="free",ncol=4,
             labeller=labeller(name=supp.labs))+
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial Genus",
       title="")+
  scale_fill_manual(values=fam_col2)+
  theme(legend.position="right", 
        legend.text = element_text(size=14),
        legend.title = element_text(size=13, face="bold"),
        legend.key.size = unit(1.2,"line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size=13),
        axis.text.x=element_text(size=12, angle=85, vjust=0.5),
        strip.text = element_text(size =13,face="bold"),
        axis.title.x = element_text(size=13, face="bold")); plot(barp3);

ggsave(filename="07_microbiome_composition_last12cats.pdf",
       device="pdf",path="./figures",
       plot=barp3,
       width=10,
       height=9.5,
       units="in",
       dpi=500);