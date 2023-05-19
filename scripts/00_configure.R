################################################################################
#
#               Fecal microbiome responses to FMT in Cats
#                      
#       Rojas et al 2023. Microbiome responses to fecal microbiota 
#             transplantation in cats with chronic digestive issues.
#
#                     Code Created By: Connie A.Rojas
#                     Created On: 22 Sept 2021
#                     Last updated: 19 May 2023
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
pacman::p_load("car","MASS","dplyr","tidyr","vegan","ggplot2",
               "picante","lme4","lmtest","multcomp", "reshape2",
               "phyloseq","phangorn","ape","stringr","dada2",
               "gridExtra","purrr","stringi",
               "compositions","pairwiseAdonis");


################################################################################
#             2. Communicate the R version and package versions to reader                 
################################################################################

print("This code was developed with R version 4.3.0");

print("The packages used and their versions were: car_3.1-2| MASS_7.3-60| 
      dplyr_1.1.2| tidyr_1.3.0| vegan_2.6-4| ggplot2_3.4.2| picante_1.8.2|
      lme4_1.1-33| lmtest_0.9-40 | multcomp_1.4-23| reshape2_1.4.4| 
      phyloseq_1.44.0| phangorn_2.11.1| ape_5.7-1 | stringr_1.5.0| 
      dada2_1.28.0| gridExtra_2.3| purrr_1.0.1| stringi_1.7.12| 
      compositions_2.0-6 | pairwiseAdonis_0.4.1");

