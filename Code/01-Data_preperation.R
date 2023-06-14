#######################################
# Purpose: Loading the SOEP data      #
# Author: Henrik-Alexander Schubert   #
# Date: 14.06.2023                    #
# E-Mail: schubert@demogr.mpg.de      #
# Pre-Requisites: full repository     #
#######################################

# Load the functions
source("Functions/packages.R")
source("Functions/functions.R")
source("Functions/graphics.R")

### Load the long individual tracking file ------------------------------

# Load the tracking file
#tracking <- read_dta(file = "SOEP_V36/Stata/ppathl.dta")

# remove the dataset
#remove(tracking)

### Clean individual Identification Data ------------------------------

# Load the data
id <- read_dta(file = "SOEP_V36/Stata/ppfad.dta")

# Residence during reunification by sex
tab(id$sex, id$loc1989)

# Create a Treatment Variable
id <- id |> 
  mutate(treatment = case_when(loc1989 == 1 & corigin ==1 ~ 1,
                               loc1989 == 2 | loc1989 ==3 ~ 0)) 

# Select the variables
id <- id |> select(birthregion, loc1989, treatment, sex, gebmonat, gebjahr, pid, persnr,  ends_with("hhnr"))

# reshape dataset into a long version
id.long <- id %>% pivot_longer(cols = 10:45,
                    names_to = "wave",
                    values_to = "household")

# create a clean number of waves
id.long$wave <- id.long$wave %>% str_remove(pattern = "hhnr")

#create a survey year variable based on the wave information
id.long <- id.long %>% group_by(persnr) %>% mutate(syear = 1984:2019) %>% ungroup()

### Household Questionnary -------------------------------------------
# Create a Variable indicating residence in East Germany

# Load household data
hh <- read_dta(file = "SOEP_V36/Stata/hgen.dta")

# Preperations: state-codes for Eastern states
selection <- c(1, 4, 8, 13, 14, 16)

# Finally creating the variable
hh <- hh |> 
  mutate(East = case_when(hgnuts1 %in% selection ~ 1,
                          hgnuts1 >= 0 & hgnuts1 %notin% selection  ~ 0))

# Remove other variables
hh <- hh |> select(hgnuts1, East, syear, hghinc, cid)


### Merge Household and Individual Data ############
df <- inner_join(hh, id.long, by = c("syear", "cid" = "household"))

# drop the old data files
remove(hh, id, id.long)

# recode treatment variable
df <- df %>% mutate(treatment2 = treatment, 
              treatment = case_when(treatment2 ==0 ~ 0,
                                    treatment2 == 1 & East == 0 ~ 1)) %>% select(persnr, treatment, syear)



### individual Data --------------------------------------

# Load the data
pers <- read.delim("SOEP_V36/Stata/pl_reduced.txt")


### Relationship Data  ############################

# Load relationship calendar
rel <- read_dta("SOEP_V36/Stata/biocouplm.dta")

# Drop observations which cannot be connected to a couple
rel <-subset(rel, subset = coupid > 0 | spelltyp >= 0)


### Relationship

# relationship retrospective
rel.retr <- read_dta("SOEP_V36/Stata/biocouply.dta")






