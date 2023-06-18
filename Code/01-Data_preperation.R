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

# Select variables
id <- subset(id, select = c(persnr, pid, birthregion, loc1989, austritt))

# Clean the birthregion
id$birthregion <- ifelse(id$birthregion %in% 11:16, "East",
                         ifelse(id$birthregion %in% 1:10, "West", NA_character_))

# Impute birth region if missing using location in 1989
id$birthregion <- ifelse(is.na(id$birthregion) & id$loc1989 == 2, "West",
                         ifelse(is.na(id$birthregion) & id$loc1989 == 1, "East", id$birthregion))

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
#pers <- read.delim("SOEP_V36/Stata/pl_reduced.txt")


#


### Biographical data -----------------------------------

# Load the data
bio <- read_dta("SOEP_V36/Stata/biol.dta")

# 
### Activity calendar ----------------------------------

# Load the data
act <- read_dta("SOEP_V36/Stata/artkalen.dta")


#Join the data
act <- left_join(act, id, by = c("pid", "persnr"))


#
act2 <- act |> 
  filter(spellnr == 1) |> 
  mutate(Total = n()) |> 
  group_by(spelltyp) |> 
  summarise(prop = n() / unique(Total),
            data = "full")


# 
act3 <-act |> 
  filter(is.na(birthregion) & spellnr == 1) |> 
  mutate(Total = n()) |> 
  group_by(spelltyp) |> 
  summarise(prop = n() / unique(Total),
            data = "missing")


#
bind_rows(act2, act3) |> 
  mutate(activity = case_when(spelltyp == 1 ~ "Voll erwerbstätig",
                            spelltyp == 2 ~ "kurzarbeit",
                            spelltyp == 3 ~ "teilzeit",
                            spelltyp == 4 ~ "betr. Ausbildung",
                            spelltyp == 5 ~ "arbeitslos",
                            spelltyp == 6 ~ "Rent",
                            spelltyp == 7 ~ "elternzeit",
                            spelltyp == 8 ~ "studium",
                            spelltyp == 9 ~ "Zivildienst",
                            spelltyp == 10 ~ "Hausfrau(-mann)",
                            spelltyp == 11 ~ "nebenberufl. Tätigkeit",
                            spelltyp == 12 ~ "sonstiges",
                            spelltyp == 13 ~ "in Lehre",
                            spelltyp == 14 ~ "Fortbildung",
                            spelltyp == 15 ~ "Minijob")) |> 
  ggplot(aes(x = activity, y = prop, fill = data)) +
  geom_col(position = position_dodge()) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(MPIDRblue, MPIDRred)) +
  coord_flip() +
  ylab("") + xlab("")

# Save the figure
ggsave(last_plot(), filename = "Figures/selection_birthreg.pdf")



### Relationship Data  ############################

# Load relationship calendar
rel <- read_dta("SOEP_V36/Stata/biocouplm.dta")

# Drop observations which cannot be connected to a couple
rel <-subset(rel, subset = coupid > 0 | spelltyp >= 0)


### Relationship

# relationship retrospective
rel.retr <- read_dta("SOEP_V36/Stata/biocouply.dta")






