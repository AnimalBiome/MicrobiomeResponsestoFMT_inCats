################################################################################
#
#               Fecal microbiome responses to FMT in Cats
#                      
#       Rojas et al 2022. Microbiome responses to fecal microbiota 
#             transplantation in cats with chronic digestive issues.
#
#                     Code Created By: Connie A.Rojas
#                     Created On: 22 Sept 2021
#                     Last updated: 4 December 2022
#
################################################################################

## CODE FOR: 
#       configuring R workspace and printing R version and package versions
#       for reader


################################################################################
#             1.  Configure the workspace for subsequent R project scripts                 
################################################################################

# set conditions for R session
rm(list=ls());
options(scipen=999);
options(stringsAsFactors = FALSE) ;

# load necessary packages
library(pacman);
pacman::p_load("car","MASS","dplyr","tidyr","reshape2","vegan","ggplot2",
               "picante","indicspecies","lme4","lmtest","multcomp","grid",
               "phyloseq","Biostrings","QsRutils","phangorn","ape","pheatmap",
               "stringr","dada2","DECIPHER","gridExtra","decontam", "sjPlot",
               "lubridate","purrr","stringi","MuMIn","chron","philentropy",
               "microbiome","cowplot","Maaslin2");


################################################################################
#             2. Communicate the R version and package versions to reader                 
################################################################################

print("This code was developed with R version 3.6.2");

print("The packages used and their versions were: pheatmap_1.0.12 | 
phangorn_2.5.5|QsRutils_0.1.4| Biostrings_2.54.0| phyloseq_1.30.0| 
multcomp_1.4-15| lmtest_0.9-38| lme4_1.1-26| indicspecies_1.7.9| picante_1.8.2| 
nlme_3.1-151| ape_5.4-1| ggplot2_3.3.3| vegan_2.5-7| permute_0.9-5| 
reshape2_1.4.4| tidyr_1.1.2| dplyr_1.0.3| MASS_7.3-53| car_3.0-10| pacman_0.5.1| 
stringr_1.4.0| DECIPHER_2.14.0| dada2_1.14.1|gridExtra_2.3| decontam_1.6.0 | 
sjPlot_2.8.7| lubridate_1.7.9.2| purrr_0.3.4|stringi_1.5.3| MuMIn_1.43.17 | 
chron_2.3-56 |philentropy_0.4.0 | microbiome_1.15.2 | cowplot_1.1.1 |
      Maaslin2_1.7.3");

