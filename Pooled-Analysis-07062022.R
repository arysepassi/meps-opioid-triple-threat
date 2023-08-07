#################################################################################
# Title:          Pooled Analysis                                               #
# Project:        MEPS Opioid Triple Threat Study                               #
# Programmer:     Aryana Sepassi                                                #
# Date:           06 July 2022                                                  #
# Updated:        NA                                                            #
# Updated by:     NA                                                            #
# Notes:          NA                                                            #
#################################################################################

################################################################################
# Load libraries                                                               #
################################################################################
devtools::install_github("tzoltak/prediction")
devtools::install_github("tzoltak/margins")

library(dplyr)            #For mutate function/data manipulation (required)
library(tidyr)            #For miscellaneous data manipulation functions 
library(maditr)           #For dcast function (to reshape data from long to wide)
library(foreign)
library(ggplot2)
library(viridis)
library(survey)
library(magrittr)
library(gee)
library(MASS)
library(ggsurvey)
library(prediction)
library(margins)
library(data.table)
library(logiBin)
library(sandwich)
library(sjstats)
library(lmtest)
library(emmeans)


options(scipen =999)      #Removes scientific notation from R output (easier to interpret)
options(digits =2)

setwd("C:/Users/aryse/OneDrive/Desktop/R Scripts/MEPS Opioid Study")

################################################################################
# Load in 2009-2019 pooled data file                                           #
################################################################################

load(file="pooled_meps.Rdata")

################################################################################
# Edits/Additions to Data Prior to Descriptive Statistics                      #
################################################################################

#Drop and recode the double panel/triple panel/single panel/gaba opioid variables: -------------------------------------------

#Code for opioid ONLY: The person ONLY reported an opioid and NO OTHER STUDY DRUGS: 
pool <- pool %>% 
  mutate(opioid_only=ifelse((narc_total>0 & nonnarc_total==0 & smr_total==0 & gaba_total==0 & bzd_total==0),1,0))

table(pool$single_tot) #25,348 patients with only an opioid reported 

#Code for BZD ONLY: The person ONLY reported a BZD and NO OTHER STUDY DRUGS:
pool <- pool %>%
  mutate(bzd_only=ifelse((bzd_total>0 & narc_total==0 & nonnarc_total==0 & smr_total==0 & gaba_total==0),1,0))

table(pool$bzd_only) #5,408 with only a BZD reported 

#Code for SMR ONLY: The person ONLY reported an SMR and NO OTHER STUDY DRUGS: 
pool <- pool %>%
  mutate(smr_only=ifelse((smr_total>0 & narc_total==0 & nonnarc_total==0 & bzd_total==0 & gaba_total==0),1,0))

table(pool$smr_only) #2,219 with only an SMR reported 

#Code for NONARC ONLY: The person ONLY reported a nonnarcotic analgesic and NO OTHER STUDY DRUGS: 
pool <- pool %>%
  mutate(nonnarc_only=ifelse((nonnarc_total>0 & narc_total==0 & bzd_total==0 & gaba_total==0 & smr_total==0),1,0))


table(pool$nonnarc_only) #22,658 with only a NONNARC reported 

#Code for GABA only: The person ONLY reported gabapentin with NO OTHER STUDY DRUGS: 
pool <- pool %>%
  mutate(gaba_only=ifelse((gaba_total>0 & narc_total==0 & nonnarc_total==0 & bzd_total==0 & smr_total==0),1,0))

table(pool$gaba_only) #2,504 with only a gabapentinoid reported 

#To flag for those with an opioid and benzodiazepine, they must have reported BOTH IN THE SAME PANEL with NO nonnarcotics, 
#SMRs, or GABA 

pool <- pool %>%
  mutate(opioid_bzd=ifelse( (((narc_panel1>0 & bzd_panel1>0) | (narc_panel2>0 & bzd_panel2>0) | (narc_panel3>0 & bzd_panel3>0) | (narc_panel4>0 & bzd_panel4>0) | (narc_panel5>0 & bzd_panel5>0)) & 
                              nonnarc_total==0 & smr_total==0 & gaba_total==0),1,0))


table(pool$opioid_bzd) #1820 

#To flag for those with an opioid and benzodiazepine AND SMR, they must have reported BOTH IN THE SAME PANEL with NO nonnarcotics, 
#SMRs, or GABA 

pool <- pool %>%
  mutate(opioid_bzd_smr=ifelse( (((narc_panel1>0 & bzd_panel1>0 & smr_panel1>0) | (narc_panel2>0 & bzd_panel2>0 & smr_panel2>0) | (narc_panel3>0 & bzd_panel3>0 & smr_panel3>0) | (narc_panel4>0 & bzd_panel4>0 & smr_panel4>0) | (narc_panel5>0 & bzd_panel5>0 & smr_panel5>0)) & 
                               nonnarc_total==0 & gaba_total==0),1,0))

table(pool$opioid_bzd_smr) #442
table(pool$year, pool$opioid_bzd_smr) 

#To flag for those with an opioid and a gabapentin: 
pool <- pool %>%
  mutate(opioid_gaba=ifelse( (((narc_panel1>0 & gaba_panel1>0) | (narc_panel2>0 & gaba_panel2>0) | (narc_panel3>0 & gaba_panel3>0) | (narc_panel4>0 & gaba_panel4>0) | (narc_panel5>0 & gaba_panel5>0)) & 
                               nonnarc_total==0 & smr_total==0 & bzd_total==0),1,0))


table(pool$opioid_gaba) #1254
table(pool$year, pool$opioid_gaba, useNA='always') 



#Create a variable flag for if the person did not have ANY of the study drugs (control group)
pool$bzd_total[is.na(pool$bzd_total)] <- 0
pool$narc_total[is.na(pool$narc_total)] <- 0
pool$nonnarc_total[is.na(pool$nonnarc_total)] <- 0
pool$smr_total[is.na(pool$smr_total)] <- 0
pool$gaba_total[is.na(pool$gaba_total)] <- 0

pool <- pool %>%
  mutate(none=ifelse((bzd_total==0 & narc_total==0 & nonnarc_total==0 & smr_total==0 & gaba_total==0),1,0))

table(pool$none, useNA='always')
#300,558 out of 376,738 unweighted participants did not have any of the study drugs 

#Change year to a factor 
table(pool$year)

#Create one variable that categorizes the different drug groups into one: ---------------------------------------------------
pool <- pool %>%
  mutate(drug_group=ifelse(opioid_only>0,"Opioid Only",
                           ifelse(bzd_only>0,"Benzodiazepine Only",
                                  ifelse(smr_only>0,"Skeletal Muscle Relaxant Only",
                                         ifelse(nonnarc_only>0, "Non-Narcotic Analgesic Only",
                                                ifelse(gaba_only>0, "Gabapentin Only",
                                                       ifelse(opioid_bzd>0,"Opioid+Benzodiazepine Only",
                                                              ifelse(opioid_bzd_smr>0,"Opioid + BZD + SMR Only",
                                                                     ifelse(opioid_gaba>0, "Opioid+Gabapentin Only",NA)))))))))


pool$drug_group <- as.factor(pool$drug_group)

table(pool$drug_group)

#Change Cost Data to Reflect Adjustments for Inflation: ---------------------------------------------------------------------

#MEPS recommends using PHCE or PCE-Health Total Indices for Pooling Total Expenditures: https://meps.ahrq.gov/about_meps/Price_Index.shtml 
#CPI-M for Pooling Out of Pocket Expenditures 
#Trends in Expenditures overall with use of PCE 
#Trends on out of pocket expenditures using CPI 

#MEPS recommends use of CPI for trends in pooling total expenditures and PCHE for pooling expenditure by type of service 

#Pooling expenditures requires a price index that is specific to health care services. 

### Part 1: Adjust for Total Expenditure using PCE (table 2 in link)

#The PCE for 2009 = 94.062. The PCE for 2021 = 115.530. 2021 price = 2009 price x (2009 PCE/1950 PCE)
#https://www.meps.ahrq.gov/about_meps/Price_Index.shtml
#PCE2009 = 94.062
#PCE2010 = 95.747
#PCE2011 = 98.170
#PCE2012 = 100
#PCE2013 = 101.354 
#PCE2014 = 102.887
#PCE2015 = 103.116
#PCE2016 = 104.146
#PCE2017 = 106.051
#PCE2018 = 108.318
#PCE2019 = 109.922
#PCE2021 = 115.530 

#PCE2019/PCE2021

PCE <- data.frame(year=c(2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019),
                  pce_multiplier=c(0.8141781,0.8287631,0.849736,0.865576,0.8772959,0.8905652,0.8925474,0.9014628,0.917952,0.9375747,0.9514585))


test <- merge(pool, PCE, by='year', all.x=TRUE)


test <- test %>%
  mutate(totexp_adj=pce_multiplier*totexp)

summary(pool$totexp)
summary(pool$totexp_adj)

pool <- test

pool <- pool %>%
  mutate(rxexp_adj=pce_multiplier*rxexp)

pool <- pool %>%
  mutate(optexp_adj=pce_multiplier*optexp)

pool <- pool %>%
  mutate(ertexp_adj=pce_multiplier*ertexp)

pool <- pool %>%
  mutate(iptexp_adj=pce_multiplier*iptexp)


### Calculate Charlson Comorbidity Index Score from vcom variables ---------------------------------------------------------
pool$vcom1[is.na(pool$vcom1)] <- 0
pool$vcom2[is.na(pool$vcom2)] <- 0
pool$vcom3[is.na(pool$vcom3)] <- 0
pool$vcom4[is.na(pool$vcom4)] <- 0
pool$vcom5[is.na(pool$vcom5)] <- 0
pool$vcom6[is.na(pool$vcom6)] <- 0
pool$vcom7[is.na(pool$vcom7)] <- 0
pool$vcom8[is.na(pool$vcom8)] <- 0
pool$vcom9[is.na(pool$vcom9)] <- 0
pool$vcom10[is.na(pool$vcom10)] <- 0
pool$vcom11[is.na(pool$vcom11)] <- 0
pool$vcom12[is.na(pool$vcom12)] <- 0
pool$vcom13[is.na(pool$vcom13)] <- 0
pool$vcom14[is.na(pool$vcom14)] <- 0
pool$vcom15[is.na(pool$vcom15)] <- 0
pool$vcom16[is.na(pool$vcom16)] <- 0

#CCI Scoring Algorithm: 
# Age: 50-59 years = +1, 60-69 years = +2, 70-79 years = +3, 80 years or more = +4
# MI: +1
#CHF: +1
# PVD: +1
# CVA or TIA: +1 
# Dementia: +1 
# COPD: +1 
#Connective tissue disease: +1
#PUD: +1
#Liver Dis: +1 mild, +3 moderate
#DM: +1 uncomplicated, +2 end organ damange, 
#hemiplegia: +2 
#Moderate-severe CKD: +2
#Tumors: Localized +2, metastatic +6
#Leukemia: +2 
#Lymphoma: +2
#AIDS: +6

pool <- pool %>%
  mutate(age_cci=ifelse(age<50,0,
                        ifelse(49<age & age<60,1,
                               ifelse(59<age & age<70,2,
                                      ifelse(69<age & age<80,3,
                                             ifelse(age>79,4, NA))))))

pool <- pool %>%
  mutate(mi=ifelse(vcom1>0, 1,0))

pool <- pool %>%
  mutate(chf=ifelse(vcom2>0,1,0))

pool <- pool %>%
  mutate(pvd=ifelse(vcom3>0,1,0))

pool <- pool %>%
  mutate(cva=ifelse(vcom4>0,1,0))

pool <- pool %>%
  mutate(dementia=ifelse(vcom5>0,1,0))

pool <- pool %>%
  mutate(copd=ifelse(vcom6>0,1,0))

pool <- pool %>%
  mutate(connect=ifelse(vcom7>0,1,0))

pool <- pool %>%
  mutate(pud=ifelse(vcom8>0,1,0))

pool <- pool %>%
  mutate(mildliverdis=ifelse(vcom9>0,1,0))

pool <- pool %>%
  mutate(diab=ifelse(vcom10>0,1,0))

pool <- pool %>%
  mutate(hemi=ifelse(vcom11>0,2,0))

pool <- pool %>%
  mutate(ckd=ifelse(vcom12>0,2,0))

pool <- pool %>%
  mutate(localcancer=ifelse(vcom13>0,2,0))

pool <- pool %>%
  mutate(modliverdis=ifelse(vcom14>0,3,0))

pool <- pool %>%
  mutate(metastatic=ifelse(vcom15>0, 6,0))

pool <- pool %>%
  mutate(aids=ifelse(vcom16>0,6,0))


#Calculate final CCI score: 
pool <- pool %>%
  mutate(cci=age_cci+mi+chf+pvd+cva+dementia+copd+connect+pud+mildliverdis+diab+hemi+ckd+localcancer+modliverdis+metastatic+aids)

summary(pool$cci)

### Create Age Category  ---------------------------------------------------------

pool <- pool %>%
  mutate(age_cat=ifelse(17<age&age<25,1,
                        ifelse(24<age&age<45,2,
                               ifelse(44<age&age<65,3,
                                      ifelse(age>64,4,NA)))))

table(pool$age_cat,useNA='always')

################################################################################
# Generate weighted  unmatched Table 1 Demographics                            #
################################################################################

#Apply exclusion criteria: Not flagged for any of the study drugs and age < 18 
colnames(pool)
table(pool$drug_group, useNA='always')

pool <- pool %>%
  mutate(studyflag=ifelse(is.na(drug_group) | age<18,0,1))

table(pool$studyflag) #48,216 out of 328,522 participants total are included in the study 

#Save pool dataset for counts and future analysis (PICK UP FROM HERE)-------------------------
save(pool, file="pool.Rdata")

load(file="pool.Rdata")

#Age, continuous, age/category, gender, race, ethnicity, marital status, education, region, poverty status, insurance coverage, 

#Define the survey design 

mepsdsn = svydesign(
  id=~varpsu,
  strata=~varstr,
  weights=~poolwt,
  data=pool,
  nest=TRUE
)




####NOTE: Survey functions CANNOT BE USED ON A SUBSET OF THE DATA. The subset group must be specified within the design function. 

table(pool$studyflag)

svymean(~age, design=subset(mepsdsn, studyflag))

tbl <- svytable(~age_cat, design=subset(mepsdsn, studyflag))
prop.table(tbl)*100

tbl<-svytable(~male, design=subset(mepsdsn, studyflag))
tbl
prop.table(tbl)*100

tbl<-svytable(~racev1x, design=subset(mepsdsn, studyflag))
tbl
prop.table(tbl)*100

table(test$racev1x, useNA='always')

tbl<-svytable(~hispanx, design=subset(mepsdsn, studyflag))
tbl
prop.table(tbl)*100

tbl<-svytable(~marry, design=subset(mepsdsn, studyflag))
tbl
prop.table(tbl)*100

tbl<-svytable(~neweducode, design=subset(mepsdsn, studyflag))
tbl
prop.table(tbl)*100

tbl<-svytable(~region, design=subset(mepsdsn, studyflag))
tbl
prop.table(tbl)*100

tbl<-svytable(~povcat, design=subset(mepsdsn, studyflag))
tbl
prop.table(tbl)*100

tbl<-svytable(~inscov, design=subset(mepsdsn, studyflag))
tbl
prop.table(tbl)*100

svymean(~cci, design=subset(mepsdsn, studyflag))

svymean(~totexp, design=subset(mepsdsn, studyflag))
svyby(~totexp, ~drug_group, design=subset(mepsdsn, studyflag), svymean)

svymean(~rxexp, design=subset(mepsdsn, studyflag))
svyby(~rxexp, ~drug_group, design=subset(mepsdsn, studyflag), svymean)

svymean(~optexp, design=subset(mepsdsn, studyflag))
svyby(~optexp, ~drug_group, design=subset(mepsdsn, studyflag), svymean)

svymean(~ertexp, design=subset(mepsdsn, studyflag))
svyby(~ertexp, ~drug_group, design=subset(mepsdsn, studyflag), svymean)

svymean(~iptexp, design=subset(mepsdsn, studyflag))
svyby(~iptexp, ~drug_group, design=subset(mepsdsn, studyflag), svymean)

svymean(~rxtot, design=subset(mepsdsn, studyflag))
svyby(~rxtot, ~drug_group, design=subset(mepsdsn, studyflag), svymean)

svymean(~obtotv, design=subset(mepsdsn, studyflag))
svyby(~obtotv, ~drug_group, design=subset(mepsdsn, studyflag), svymean)

svymean(~ertot, design=subset(mepsdsn, studyflag))
svyby(~ertot, ~drug_group, design=subset(mepsdsn, studyflag), svymean)

svymean(~ipdis, design=subset(mepsdsn, studyflag))
svyby(~ipdis, ~drug_group, design=subset(mepsdsn, studyflag), svymean)

svymean(~ipngtd, design=subset(mepsdsn, studyflag))
svyby(~ipngtd, ~drug_group, design=subset(mepsdsn, studyflag), svymean)


#Generate a second study flag only to account for the five study groups: opioid only, opioid+BZD, opioid+BZD+SMR, opioid+gabapentin, non-narcotic analgesic 
table(pool$studyflag, pool$drug_group)

pool <- pool %>%
  mutate(flag2 = ifelse(studyflag==1&(drug_group=="Opioid Only"| drug_group=="Opioid+Benzodiazepine Only"| drug_group=="Opioid + BZD + SMR Only"| drug_group=="Opioid+Gabapentin Only"|drug_group=="Non-Narcotic Analgesic Only"),1,0))
table(pool$flag2, pool$drug_group)

svytable(~drug_group, design=subset(mepsdsn, flag2))

svyby(~age, ~drug_group, design=subset(mepsdsn, flag2), svymean)

svytable(~age_cat+drug_group, design=subset(mepsdsn, flag2))

svytable(~male+drug_group, design=subset(mepsdsn, flag2))

svytable(~male+drug_group, design=subset(mepsdsn, flag2))

svytable(~racev1x+drug_group, design=subset(mepsdsn, flag2))

svytable(~hispanx+drug_group, design=subset(mepsdsn, flag2))

svytable(~marry+drug_group, design=subset(mepsdsn, flag2))

svytable(~neweducode+drug_group, design=subset(mepsdsn, flag2))

svytable(~region+drug_group, design=subset(mepsdsn, flag2))

svytable(~povcat+drug_group, design=subset(mepsdsn, flag2))

svytable(~inscov+drug_group, design=subset(mepsdsn, flag2))

svyby(~cci, ~drug_group, design=subset(mepsdsn, flag2), svymean)

#Perform unweighted comparative statistical analysis 
#Continuous values = ANOVA. Categorical = chi-squared test 

svychisq(~age_cat+drug_group, design=subset(mepsdsn, flag2), statistic="adjWald")
svychisq(~male+drug_group, design=subset(mepsdsn, flag2), statistic="adjWald")
svychisq(~racev1x+drug_group, design=subset(mepsdsn, flag2), statistic="adjWald")
svychisq(~hispanx+drug_group, design=subset(mepsdsn, flag2), statistic="adjWald")
svychisq(~marry+drug_group, design=subset(mepsdsn, flag2), statistic="adjWald")
svychisq(~neweducode+drug_group, design=subset(mepsdsn, flag2), statistic="adjWald")
svychisq(~region+drug_group, design=subset(mepsdsn, flag2), statistic="adjWald")
svychisq(~povcat+drug_group, design=subset(mepsdsn, flag2), statistic="adjWald")
svychisq(~inscov+drug_group, design=subset(mepsdsn, flag2), statistic="adjWald")


#Save pool as final unmatched pooled dataset ---------
save(pool, file="pool_unmatch.Rdata")

library(dplyr)            #For mutate function/data manipulation (required)
library(tidyr)            #For miscellaneous data manipulation functions 
library(maditr)           #For dcast function (to reshape data from long to wide)
library(foreign)
library(ggplot2)
library(viridis)
library(survey)
library(magrittr)
library(gee)
library(MASS)
library(ggsurvey)
library(margins)          #To estimate marginal effects for model objects 
library(prediction)
library(svyVGAM)
library(emmeans)
library(sjstats)          #To perform survey weighted negative binomial regression

options(scipen =999)      #Removes scientific notation from R output (easier to interpret)
options(digits =2)

setwd("C:/Users/aryse/OneDrive/Desktop/R Scripts/MEPS Opioid Study")


load(file="pool_unmatch.Rdata")

#Create a flag to indicate only the study drugs to be included in the analysis (treatGrp) 
pool <- pool %>%
  mutate(treatGrp = ifelse(drug_group=="Non-Narcotic Analgesic Only", "Non-Narcotic",
                           ifelse(drug_group=="Opioid + BZD + SMR Only", "OpioidBZDSMR",
                                  ifelse(drug_group=="Opioid Only","Opioid",
                                         ifelse(drug_group=="Opioid+Benzodiazepine Only","OpioidBZD",
                                                ifelse(drug_group=="Opioid+Gabapentin Only", "OpioidGABA",NA))))))
table(pool$treatGrp)
#Check to see that treatGrp truly reflects the drug_group values of interest: 
table(pool$treatGrp, pool$drug_group) #treatGrp only includes non-narcotic analgesics, opioid+BZD+SMR, opioid only, opioid+BZD, and op+GABA 

#Create a flag to generate inclusion for the study group - if the person is also over age 18 and is NOT NA for drug_group, then flag = 1: 

pool <- pool %>%
  mutate(flag2 = ifelse(is.na(treatGrp)|age<18,0,1))

table(pool$flag2)
table(pool$flag2, pool$treatGrp) #Correct ------ 


##Recode the negative marital status category to a positive level (negative integers are not allowed) 

table(pool$marry)

pool <- pool %>%
  mutate(marry2 = ifelse(marry==1, 1,
                         ifelse(marry==2,2,
                                ifelse(marry==3,3,
                                       ifelse(marry==4,4,
                                              ifelse(marry==5,5,
                                                     ifelse(marry==-1,6, NA)))))))

#Generate a variable that indicates whether person had a zero or non-zero expenditure (for the two-part model)

pool <- pool %>%
  mutate(zero=ifelse(totexp_adj==0,0,
                     ifelse(totexp_adj>0 & totexp_adj<9999999999999,1,NA)))

table(pool$zero)


#Region: 
pool <- pool %>%
  mutate(region2 = ifelse(region=="1 Northeast","Northeast",
                          ifelse(region=="2 Midwest","Midwest",
                                 ifelse(region=="3 South","South",
                                        ifelse(region=="4 West","West",NA)))))
table(pool$region2, useNA='always') 


#Re-Code the education variable (neweducode was incorrectly coded): 
table(pool$hideg)
table(pool$eduyrdg)

pool <- pool %>%
  mutate(pool,
         education=case_when(
           hideg== -15 ~ 1,
           hideg== -9 ~ 1, 
           hideg== -8 ~ 1,
           hideg== -7 ~ 1, 
           hideg== -1 ~ 1, 
           eduyrdg==-9 ~1,
           eduyrdg==-8 ~1,
           eduyrdg==-7 ~1,
           eduyrdg==-1 ~1,
           eduyrdg== 1 ~1,
           eduyrdg==10 ~1,
           hideg== 1   ~2,
           eduyrdg== 2 ~2,
           hideg== 2   ~3,
           hideg== 3   ~3,
           eduyrdg== 3 ~3,
           eduyrdg== 4 ~3,
           hideg== 7   ~4,
           eduyrdg== 5 ~4,
           eduyrdg== 6 ~4,
           eduyrdg== 7 ~4,
           hideg== 4   ~5,
           eduyrdg== 8 ~5,
           hideg== 5   ~6,
           hideg== 6   ~6,
           eduyrdg== 9 ~6
         ))

pool$education <- factor(pool$education,
                         levels=c(1,2,3,4,5,6),
                         labels=c("1 Not Ascertained or <18 Y/O",
                                  "2 No Diploma",
                                  "3 GED or HS",
                                  "4 Assoc, Other, or <4 Yrs College",
                                  "5 BS",
                                  "6 MS, Doctorate, Profesional"))

table(pool$education, useNA='always')

table(pool$treatGrp)

18379+1820+442+1254 

#21895 opioid users 
#354843 opioid non-users 

#save(pool, file="C:/Users/aryse/OneDrive/Desktop/pool.Rdata")
#write.dta(pool, file="C:/Users/aryse/OneDrive/Desktop/pool.dta")

#Generate the survey design object
colnames(pool)

table(pool$treatGrp, useNA='always')


mepsdsn = svydesign(
  id=~varpsu,
  strata=~varstr,
  weights=~pooledwt, 
  data=pool, 
  nest=TRUE
)

study_dsn <- subset(mepsdsn, flag2==1) #Flag 2 denotes a flag for the study drugs 

###### Results, Table 1: Generate Demographic Contingency Tables: -----------------------------------------------------------------------
#Note: using the svyby command and svytable commands below generate the same contingency tables: 
colnames(pool)

#Total counts for each drug category: 
svytable(~treatGrp, design=study_dsn)
15916017+1802395+437755+14899352+1186319


#Age Category by drug group: 
#For entire cohort:
svytable(~age_cat, design=study_dsn)

svytable(~age_cat+treatGrp, design=study_dsn)
round(prop.table(svytable(~age_cat+treatGrp, design=study_dsn),margin=2)*100,digits=2)
svychisq(~treatGrp+age_cat, design=study_dsn, statistic="adjWald")

#Gender by drug group: 
svytable(~male, design=study_dsn)

svytable(~male+treatGrp, design=study_dsn)
round(prop.table(svytable(~male+treatGrp, design=study_dsn),margin=2)*100,digits=2)
svychisq(~treatGrp+male, design=study_dsn, statistic="adjWald")

#Race by drug group: 
svytable(~racev1x, design=study_dsn)

svytable(~racev1x+treatGrp, design=study_dsn, exclude=NULL)
round(prop.table(svytable(~racev1x+treatGrp, design=study_dsn),margin=2)*100,digits=2)
svychisq(~treatGrp+racev1x, design=study_dsn, statistic="adjWald")

#Ethnicity by drug group: 
svytable(~hispanx, design=study_dsn)
svytable(~hispanx+treatGrp, design=study_dsn)
round(prop.table(svytable(~hispanx+treatGrp, design=study_dsn),margin=2)*100,digits=2)
svychisq(~hispanx+treatGrp, design=study_dsn, statistic="adjWald")

#Marital status by drug group: 
svytable(~marry2, design=study_dsn)
svytable(~marry2+treatGrp, design=study_dsn)
round(prop.table(svytable(~marry2+treatGrp, design=study_dsn),margin=2)*100,digits=2)
svychisq(~marry2+treatGrp, design=study_dsn, statistic="adjWald")

#Education by drug group: 
svytable(~education, design=study_dsn)
svytable(~education+treatGrp, design=study_dsn)
round(prop.table(svytable(~education+treatGrp, design=study_dsn),margin=2)*100,digits=2)
svychisq(~education+treatGrp, design=study_dsn, statistic="adjWald")

#Region by drug group: 
svytable(~region2, design=study_dsn)
svytable(~region2+treatGrp, design=study_dsn)round(prop.table(svytable(~education+treatGrp, design=study_dsn),margin=2)*100,digits=2)
round(prop.table(svytable(~region2+treatGrp, design=study_dsn),margin=2)*100,digits=2)
svychisq(~region2+treatGrp, design=study_dsn, statistic="adjWald")

#Income category by drug group: 
svytable(~povcat, design=study_dsn)
svytable(~povcat+treatGrp, design=study_dsn)
round(prop.table(svytable(~povcat+treatGrp, design=study_dsn),margin=2)*100,digits=2)
svychisq(~povcat+treatGrp, design=study_dsn, statistic="adjWald")

#Insurance coverage by drug group: 
svytable(~inscov, design=study_dsn)
svytable(~inscov+treatGrp, design=study_dsn)
round(prop.table(svytable(~inscov+treatGrp, design=study_dsn),margin=2)*100,digits=2)
svychisq(~inscov+treatGrp, design=study_dsn, statistic="adjWald")

#Charlson Comorbidity Index by drug group: 
svyby(~cci, ~treatGrp, design=study_dsn, svymean)
svyranktest(~cci+treatGrp, design=study_dsn, test='KruskalWallis')

####Results, Table 2: Outcomes: Unadjusted Total Expenditure --------------------------------------------------------------------------------------
#This code generates results that counts mean total expenditure and expenditure of inpatient stays, outpatient stays, office based
#stays, ER stays, and prescription medications across the five  treatment groups: 

svyby(~totexp_adj, ~treatGrp, design=study_dsn, svymean)
svyby(~iptexp_adj, ~treatGrp, design=study_dsn, svymean)
svyby(~optexp_adj, ~treatGrp, design=study_dsn, svymean)
svyby(~obvexp_adj, ~treatGrp, design=study_dsn, svymean)
svyby(~ertexp_adj, ~treatGrp, design=study_dsn, svymean)
svyby(~rxexp_adj, ~treatGrp, design=study_dsn, svymean)

svyby(~totexp_adj, ~treatGrp, design=study_dsn, svymean)
summary(svyglm(totexp_adj~treatGrp, design=study_dsn))






#Expenditure is not normally distributed; use the Kruskal-Wallis test as a replacement for ANOVA: 
svyranktest(~totexp_adj+treatGrp, design=study_dsn, test='KruskalWallis') 
svyranktest(~iptexp_adj+treatGrp, design=study_dsn, test='KruskalWallis')
svyranktest(~optexp_adj+treatGrp, design=study_dsn, test='KruskalWallis')
svyranktest(~obvexp_adj+treatGrp, design=study_dsn, test='KruskalWallis')
svyranktest(~ertexp_adj+treatGrp, design=study_dsn, test='KruskalWallis')
svyranktest(~rxexp_adj+treatGrp, design=study_dsn, test='KruskalWallis')

####Results, Table 2: Outcomes: Unadjusted Total Resource Utilization --------------------------------------------------------------------------------------
#This code generates results that count mean resource utilization across different resource categories and the 
#five treatment groups 
svyby(~ipdis, ~treatGrp, design=study_dsn, svymean)
svyby(~optotv, ~treatGrp, design=study_dsn, svymean)
svyby(~obtotv, ~treatGrp, design=study_dsn, svymean)
svyby(~ertot, ~treatGrp, design=study_dsn, svymean)
svyby(~rxtot, ~treatGrp, design=study_dsn, svymean)

#Resource Utilization is not normally distributed; use the Kruskal-Wallis test as a replacement for ANOVA: 
svyranktest(~ipdis+treatGrp, design=study_dsn, test='KruskalWallis') 
svyranktest(~optotv+treatGrp, design=study_dsn, test='KruskalWallis')
svyranktest(~obtotv+treatGrp, design=study_dsn, test='KruskalWallis')
svyranktest(~ertot+treatGrp, design=study_dsn, test='KruskalWallis')
svyranktest(~rxtot+treatGrp, design=study_dsn, test='KruskalWallis')

####Results, Table 4: Average Marginal Effects of Expenditure Outcomes--------------------------------------------------------------------

                                                #####################
################################################ # Total Expenditure #####################################################################
                                                #####################

#Recall that only positive integers can be used with a gamma distribution. 
#Therefore, first, evaluate how many participants in each treatment group had a total expenditure of $0 for the year: 
#First, create a variable for total expenditure; if totexp_adj > 0, then 1. Else, 0: 
pool <- pool %>% 
  mutate(totexp_adj_tpm = ifelse(totexp_adj==0,0,
                                 ifelse(totexp_adj>0, 1, NA)))

table(pool$totexp_adj_tpm, pool$treatGrp)
#Only one participant from the opioid group had a total expenditure of $0. We will exclude this person from the analysis 
#Subset to exclude this observation for the analysis 
mepsdsn = svydesign(
  id=~varpsu,
  strata=~varstr,
  weights=~pooledwt,
  data=pool,
  nest=TRUE
)
study_dsn <- subset(mepsdsn, flag2==1) 
study_dsn1 <- subset(study_dsn, totexp_adj_tpm==1)


#Generalized Linear Model using Gamma distribution with a log link fit  
m1 <- svyglm(totexp_adj~treatGrp+factor(age_cat)+factor(year.x)+factor(male)+factor(hispanx)+
               factor(marry2)+factor(education)+factor(region2)+factor(povcat)+factor(inscov)+cci,
             design=study_dsn1,
             family=Gamma (link=log))
summary(m1)

#Exponentiate the coefficients to produce non-log results: 
x <- exp(coef(m1))
format(round(x, 2), nsmall=2)

y <- exp(confint(m1))
format(round(y, 2), nsmall=2)

#Produce regression results with robust standard errors 
cov.m1 <- vcov(m1, type="HC0")
std.err <- sqrt(diag(cov.m1))
r.est <- cbind(Estimate= coef(m1), "Robust SE" = std.err,
               "Pr(>|z|)" = 2 * pnorm(abs(coef(m1)/std.err), lower.tail=FALSE),
               LL = coef(m1) - 1.96 * std.err,
               UL = coef(m1) + 1.96 * std.err)

r.est


#Average Marginal Effects: 
margins(m1, design=study_dsn1) %>% summary()

                                                ########################
################################################# Inpatient Expenditure######################################################################
                                                ########################

#Recall that only positive integers can be used with a gamma distribution. 
#Therefore, first, evaluate how many participants in each treatment group had a total expenditure of $0 for the year: 
pool <- pool %>% 
  mutate(iptexp_adj_tpm = ifelse(iptexp_adj==0,0,
                                 ifelse(iptexp_adj>0, 1, NA)))

table(pool$iptexp_adj_tpm, pool$treatGrp)
mepsdsn = svydesign(
  id=~varpsu,
  strata=~varstr,
  weights=~poolwt,
  data=pool,
  nest=TRUE
)
study_dsn <- subset(mepsdsn, flag2==1) 
#There are a lot of zero values for inpatient expenditure; therefore, a two part model would be ideal: 
####### The two part models will be generated in stata: 
#Export pool dataset as a stata dataset: 
write.dta(pool, "C:/Users/aryse/OneDrive/Desktop/mepsopioid.dta")

#STATA was used to generate output for AMEs and Regression Coefficients for Inpatient Expenditures, Outpatient Expenditures, 
#Office-Based Expenditures, and ED Expenditures 

                                                ########################
#################################################   Rx Expenditure     ######################################################################
                                                ########################


pool <- pool %>% 
  mutate(rxexp_adj_tpm = ifelse(rxexp_adj==0,0,
                                 ifelse(rxexp_adj>0, 1, NA)))
table(pool$rxexp_adj_tpm, pool$treatGrp)
#Only Non-Narcotic and Opioid Groups have people in zero expenditure (total = 31 participants)

#Exclude from analysis 
mepsdsn = svydesign(
  id=~varpsu,
  strata=~varstr,
  weights=~poolwt,
  data=pool,
  nest=TRUE
)
study_dsn <- subset(mepsdsn, flag2==1) 
study_dsn1 <- subset(study_dsn, rxexp_adj_tpm==1)


#Generalized Linear Model using Gamma distribution with a log link fit  
m2 <- svyglm(rxexp_adj~treatGrp+factor(age_cat)+factor(year.x)+factor(male)+factor(hispanx)+
               factor(marry2)+factor(education)+factor(region2)+factor(povcat)+factor(inscov)+cci,
             design=study_dsn1,
             family=Gamma (link=log))

summary(m2)


#Exponentiate the coefficients to produce non-log results: 
x <- exp(coef(m1))
format(round(x, 2), nsmall=2)

y <- exp(confint(m1))
format(round(y, 2), nsmall=2)

#Produce regression results with robust standard errors 
cov.m1 <- vcov(m1, type="HC0")
std.err <- sqrt(diag(cov.m1))
r.est <- cbind(Estimate= coef(m1), "Robust SE" = std.err,
               "Pr(>|z|)" = 2 * pnorm(abs(coef(m1)/std.err), lower.tail=FALSE),
               LL = coef(m1) - 1.96 * std.err,
               UL = coef(m1) + 1.96 * std.err)

r.est

#Average Marginal Effects: 
margins(m2, design=study_dsn1) %>% summary()

#Part 2: Resource Utilization Results-------------------------------------------------
#Resource use variables are count variables with overdispersion; a series of negative binomial models were fit in this case. 

#This code checks for overdispersion (conditional means and variances) 
with(pool, tapply(rxtot, drug_group, function(x){
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))

#Conditional means and variances are not equal, indicating that over-dispersion is present (SD for certain groups is higher)
#Therefore, we may fit a negative binomial distribution. As an alternative, quasi-poisson to account for changes in variance may also work.
#Model 3: Number of inpatient stay discharges (ipdis)
m3 <- svyglm.nb(ipdis~factor(treatGrp)+factor(age_cat)+factor(year.x)+factor(male)+factor(hispanx)+factor(marry2)+
                  factor(education)+factor(region2)+factor(povcat)+factor(inscov)+factor(racev1x),
                design=study_dsn)

m3

#Model 4: Outpatient Visits (optotv)
m4 <- svyglm.nb(optotv~factor(treatGrp)+factor(age_cat)+factor(year.x)+factor(male)+factor(hispanx)+factor(marry2)+factor(education)+
                  factor(region2)+factor(povcat)+factor(inscov)+factor(racev1x)+cci,
                design=study_dsn)

m4

#Model 5: Office-Based Visits (obtotv) using svyglm.nb 
m5 <- svyglm.nb(obtotv~factor(treatGrp)+factor(age_cat)+factor(year.x)+factor(male)+factor(hispanx)+factor(marry2)+
                  factor(education)+factor(region2)+factor(povcat)+factor(inscov),
                design=study_dsn)

m5

#Model 6: ER Visits (ertot)
m6 <- svyglm.nb(ertot~factor(treatGrp)+factor(age_cat)+factor(year.x)+factor(male)+factor(hispanx)+factor(marry2)+
                  factor(education)+factor(region2)+factor(povcat)+factor(inscov)+factor(racev1x)+cci,
                design=study_dsn)

m6


#Model 7: Prescription Drug Use 
m7 <- svyglm.nb(rxtot~factor(treatGrp)+factor(age_cat)+factor(year.x)+factor(male)+factor(hispanx)+factor(marry2)+
                  factor(education)+factor(region2)+factor(povcat)+factor(inscov)+cci,
                design=study_dsn)

m7


#### Part 3: Changes in outcome categories year over year compared to non-narcotic analgesic only group: 
#These are all linear regression models adjusted for covariates using an interaction term to evaluate changes 
#over time (year x drug group)

#Model 1: Average total expenditure per year change 
m1 <- svyglm(totexp_adj~factor(year.x)+age+factor(racev1x)+
               factor(hispanx)+factor(marry2)+factor(region2)+factor(povcat)+
               factor(inscov)+factor(education)+factor(male)+cci+treatGrp*year.x,
             design=study_dsn)

summary(m1)
confint(m1)

#Model 2: Average inpatient expenditure per year change 
m2 <- svyglm(iptexp_adj~factor(year.x)+age+factor(racev1x)+
               factor(hispanx)+factor(marry2)+factor(region2)+factor(povcat)+
               factor(inscov)+factor(education)+factor(male)+cci+treatGrp*year.x,
             design=study_dsn)

summary(m2)
confint(m2)

#Model 3: Average outpatient expenditure per year change 
m3 <- svyglm(optexp_adj~factor(year.x)+age+factor(racev1x)+
               factor(hispanx)+factor(marry2)+factor(region2)+factor(povcat)+
               factor(inscov)+factor(education)+factor(male)+cci+treatGrp*year.x,
             design=study_dsn)

summary(m3)
confint(m3)

#Model 4: Average office-based expenditure per year change 
m4 <- svyglm(obvexp_adj~factor(year.x)+age+factor(racev1x)+
               factor(hispanx)+factor(marry2)+factor(region2)+factor(povcat)+
               factor(inscov)+factor(education)+factor(male)+cci+treatGrp*year.x,
             design=study_dsn)

summary(m4)
confint(m4)

#Model 5: Average ER expenditure per year change 
m5 <- svyglm(ertexp_adj~factor(year.x)+age+factor(racev1x)+
               factor(hispanx)+factor(marry2)+factor(region2)+factor(povcat)+
               factor(inscov)+factor(education)+factor(male)+cci+treatGrp*year.x,
             design=study_dsn)

summary(m5)
confint(m5)



###########################################################################################################################################################################         
### Aim 2 Results 
################################################################################################################################################
#Figure 1a: Average of the Total Medical Expenditure per Person per Year##################################################

#Step 1: Generate a table of the average total medical expenditure per person per drug type using the survey package: 
df <- svyby(~totexp_adj, ~treatGrp+year.x, design=study_dsn, svymean)
df
#Drop standard error column
df <- subset(df, select=-c(se))
      
#Step 2: Plot in ggplot2: 
fig1a <- ggplot(data=df,
                aes(x=year.x, y=totexp_adj, group=treatGrp, color=treatGrp)) + 
                geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
                geom_point() +
                scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
               # theme_bw() +
                theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),panel.background=element_blank()) +
                theme(axis.line.x=element_line(color='black',size=0.5),
                      axis.line.y=element_line(color='black',size=0.5))+
  theme(axis.text.x=element_text(size=12, color='black'), axis.text.y=element_text(size=16, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=12)) + 
  theme(axis.title.y=element_text(color='black', size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  xlab('Year') +
  ylab('Average Total Medical Expenditure Per Person') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=10)) +
  theme(legend.text=element_text(size=10)) +
  scale_y_continuous(labels=scales::dollar_format(), limits=c(0,40000))+
  theme(legend.key=element_rect(fill="white"))+ #removes gray background in legend icons 
  annotate("text",x=2019,y=40000,label="A", size=8)

  

fig1a

ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure1a.tiff", width=8, height=8, device='tiff', dpi=700)

#Figure 1b: Average of the Total Inpatient Expenditure per Person per Year##################################################

#Step 1: Generate a table of the average total medical expenditure per person per drug type using the survey package: 
df <- svyby(~iptexp_adj, ~treatGrp+year.x, design=study_dsn, svymean)
df
#Drop standard error column
df <- subset(df, select=-c(se))

#Step 2: Plot in ggplot2: 
fig1b <- ggplot(data=df,
                aes(x=year.x, y=iptexp_adj, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Average Total Inpatient Expenditure Per Person') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=10)) +
  scale_y_continuous(labels=scales::dollar_format(), limits=c(0,15000)) +
  theme(legend.key=element_rect(fill="white"))+ #removes gray background in legend icons 
  annotate("text",x=2019,y=15000,label="B", size=8)


fig1b
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure1b.tiff", width=8, height=8, device='tiff', dpi=700)

#Figure 1c: Average of the Total Outpatient Expenditure per Person per Year##################################################
#Step 1: Generate a table of the average total medical expenditure per person per drug type using the survey package: 
df <- svyby(~optexp_adj, ~treatGrp+year.x, design=study_dsn, svymean)
df

#Step 2: Plot in ggplot2: 
fig1c <- ggplot(data=df,
                aes(x=year.x, y=optexp_adj, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Average Total Outpatient Expenditure Per Person') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8))+
  scale_y_continuous(labels=scales::dollar_format(), limits=c(0,4000)) +
  theme(legend.key=element_rect(fill="white"))+ #removes gray background in legend icons 
  annotate("text",x=2019,y=4000,label="C", size=8)


fig1c
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure1c.tiff", width=8, height=8, device='tiff', dpi=700)

#Figure 1d: Average of the Total Office-Based Expenditure per Person per Year##################################################

#Step 1: Generate a table of the average total medical expenditure per person per drug type using the survey package: 
df <- svyby(~obvexp_adj, ~treatGrp+year.x, design=study_dsn, svymean)
df

#Step 2: Plot in ggplot2: 
fig1d <- ggplot(data=df,
                aes(x=year.x, y=obvexp_adj, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Average Total Office-Based Expenditure Per Person') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8))+
  scale_y_continuous(labels=scales::dollar_format(), limits=c(0,6000)) +
  theme(legend.key=element_rect(fill="white"))+ #removes gray background in legend icons 
  annotate("text",x=2019,y=6000,label="D", size=8)


fig1d
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure1d.tiff", width=8, height=8, device='tiff', dpi=700)

remove(fig1d)
#Figure 1e: Average of the Total ED Expenditure per Person per Year##################################################

#Step 1: Generate a table of the average total medical expenditure per person per drug type using the survey package: 
df <- svyby(~ertexp_adj, ~treatGrp+year.x, design=study_dsn, svymean)
df

#Step 2: Plot in ggplot2: 
fig1e <- ggplot(data=df,
                aes(x=year.x, y=ertexp_adj, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Average Total ED Expenditure Per Person') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8))+
  scale_y_continuous(labels=scales::dollar_format(), limits=c(0,1300)) +
  theme(legend.key=element_rect(fill="white"))+ #removes gray background in legend icons 
  annotate("text",x=2019,y=1300,label="E", size=8)

fig1e
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure1e.tiff", width=8, height=8, device='tiff', dpi=700)
remove(fig1e)

#Figure 1f: Average of the Total Prescription Drug Expenditure per Person per Year##################################################
#Step 1: Generate a table of the average total medical expenditure per person per drug type using the survey package: 
df <- svyby(~rxexp_adj, ~treatGrp+year.x, design=study_dsn, svymean)
df

#Step 2: Plot in ggplot2: 
fig1f <- ggplot(data=df,
                aes(x=year.x, y=rxexp_adj, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Average Total Prescription Drug Expenditure Per Person') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8))+
  scale_y_continuous(labels=scales::dollar_format(), limits=c(0,20000)) +
  theme(legend.key=element_rect(fill="white"))+ #removes gray background in legend icons 
  annotate("text",x=2019,y=20000,label="F", size=8)

fig1f

ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure1f.tiff", width=8, height=8, device='tiff', dpi=700)

remove(fig1f)

###### FIGURE 2 ##############################################################################################
#Figure 2a: Average Inpatient Discharges per Person per Year##################################################

#Step 1: Generate a table of the average inpatient discharges per person per drug type using the survey package: 
df <- svyby(~ipdis, ~treatGrp+year.x, design=study_dsn, svymean)
df

fig2a <- ggplot(data=df,
                aes(x=year.x, y=ipdis, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Average Inpatient Discharges Per Person') +
  ylim(0,1) +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8))+
  theme(legend.key=element_rect(fill="white"))+ #removes gray background in legend icons 
  annotate("text",x=2019,y=1,label="A", size=8)  

fig2a

ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure2a.tiff", width=8, height=8, device='tiff', dpi=700)
remove(fig2a)

#Figure 2b: Average Outpatient Visits per Person per Year##################################################
#Step 1: Generate a table of the average inpatient discharges per person per drug type using the survey package: 
df <- svyby(~optotv, ~treatGrp+year.x, design=study_dsn, svymean)
df

fig2b <- ggplot(data=df,
                aes(x=year.x, y=optotv, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Average Outpatient Visits Per Person') +
  ylim(0,6) +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8))+
  theme(legend.text=element_text(size=8))+
  theme(legend.key=element_rect(fill="white"))+ #removes gray background in legend icons 
  annotate("text",x=2019,y=6,label="B", size=8)  


fig2b

ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure2b.tiff", width=8, height=8, device='tiff', dpi=700)

remove(fig2b)

#Figure 2c: Average Office-Based Visits per Person per Year##################################################
#Step 1: Generate a table of the average inpatient discharges per person per drug type using the survey package: 
df <- svyby(~obtotv, ~treatGrp+year.x, design=study_dsn, svymean)
df

fig2c <- ggplot(data=df,
                aes(x=year.x, y=obtotv, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Average Office-Based Visits Per Person') +
  ylim(5,30) +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.key=element_rect(fill="white"))+ #removes gray background in legend icons 
  annotate("text",x=2019,y=30,label="C", size=8)  


fig2c

ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure2c.tiff", width=8, height=8, device='tiff', dpi=700)

remove(fig2c)

#Figure 2d: Average ER Visits per Person per Year##################################################
#Step 1: Generate a table of the average inpatient discharges per person per drug type using the survey package: 
df <- svyby(~ertot, ~treatGrp+year.x, design=study_dsn, svymean)
df

fig2d <- ggplot(data=df,
                aes(x=year.x, y=ertot, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Average Emergency Department Visits Per Person') +
  ylim(0,1.5) +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) + 
  theme(legend.text=element_text(size=8)) +
  theme(legend.key=element_rect(fill="white"))+ #removes gray background in legend icons 
  annotate("text",x=2019,y=1.5,label="D", size=8)  


fig2d

ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure2d.tiff", width=8, height=8, device='tiff', dpi=700)

remove(fig2d)

#Figure 2e: Average Number of Prescription Medications per Person per Year##################################################
#Step 1: Generate a table of the average inpatient discharges per person per drug type using the survey package: 
df <- svyby(~rxtot, ~treatGrp+year.x, design=study_dsn, svymean)
df

fig2e <- ggplot(data=df,
                aes(x=year.x, y=rxtot, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Average Number of Prescription Medications Per Person') +
  ylim(10,70) +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.key=element_rect(fill="white"))+ #removes gray background in legend icons 
  annotate("text",x=2019,y=70,label="E", size=8)  


fig2e

ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure2e.tiff", width=8, height=8, device='tiff', dpi=700)

remove(fig2e)

###### FIGURE 3 ##############################################################################################
#Figure 3: Generate a line graph of the number of unique participants per drug category per year: 
df <- svyby(~treatGrp, ~year.x, design=study_dsn, svytotal)
#Use gather function to change from wide to long format: 
library(tidyr)
library(scales)

df_long <- gather(df, treatGrp, n, `treatGrpNon-Narcotic`:treatGrpOpioidGABA, factor_key=TRUE)
df_long <- subset(df_long, select=-c(`se.treatGrpNon-Narcotic`, se.treatGrpOpioid, se.treatGrpOpioidBZD, se.treatGrpOpioidBZDSMR,se.treatGrpOpioidGABA))
df_long

fig3 <- ggplot(data=df_long,
                aes(x=year.x, y=n, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(axis.line.x=element_line(color='black',size=0.5),
        axis.line.y=element_line(color='black',size=0.5))+
  theme(axis.text.x=element_text(size=12, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=12)) + 
  theme(axis.title.y=element_text(color='black', size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  xlab('Year') +
  ylab('Total Number of Participants') +
#  ylim(0,2500000) +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  scale_y_continuous(label=comma, limits=c(0,2000000))+
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=10)) + #1314-1316 adjust font sizes for legend 
  theme(legend.text=element_text(size=10)) +
  theme(legend.title=element_text(size=10)) +
  theme(legend.key=element_rect(fill="white")) #removes gray background in legend icons 

fig3

ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figure1.tiff", width=8, height=5, device='tiff', dpi=700)


###### Supplementary Material Tables ##############################################################################################
#Table S4: Number of unique participants per drug group per year: 
svytable(~year.x+treatGrp, design=study_dsn)

#Table S5: Average Expenditure per person per drug group per year: 
#Total expenditure: 
df_ipt <- svyby(~totexp_adj, ~year.x+treatGrp, design=study_dsn, svymean)
df_ipt

#Inpatient expenditure: 
df <- svyby(~iptexp_adj, ~year.x+treatGrp, design=study_dsn, svymean)
df

#Outpatient expenditure: 
df <- svyby(~optexp_adj, ~year.x+treatGrp, design=study_dsn, svymean)
df

#Office-Based expenditure: 
df <- svyby(~obvexp_adj, ~year.x+treatGrp, design=study_dsn, svymean)
df

#ED expenditure: 
df <- svyby(~ertexp_adj, ~year.x+treatGrp, design=study_dsn, svymean)
df

#Drug expenditure: 
df <- svyby(~rxexp_adj, ~year.x+treatGrp, design=study_dsn, svymean)
df

#####Table S5: Average resource utilization per person per drug group per year, 2009-2019
#Inpatient Discharges: 
df <- svyby(~ipdis, ~year.x+treatGrp, design=study_dsn, svymean)
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2
colnames(pool)

#Outpatient Visits: 
df <- svyby(~optotv, ~year.x+treatGrp, design=study_dsn, svymean)
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#Office-Based Visits: 
df <- svyby(~obtotv, ~year.x+treatGrp, design=study_dsn, svymean)
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#ED Visits: 
df <- svyby(~ertot, ~year.x+treatGrp, design=study_dsn, svymean)
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#Prescription Medicines: 
df <- svyby(~rxtot, ~year.x+treatGrp, design=study_dsn, svymean)
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2


#Supplemental Figure 3: Total number of resources utilized per drug group per year: 

#Figure S3A: Total Inpatient Discharges per Drug Group per Year 
df <- svyby(~ipdis, ~treatGrp+year.x, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df

figS3A <- ggplot(data=df,
                aes(x=year.x, y=ipdis, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Total Number of Inpatient Discharges') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  scale_y_continuous(labels=comma, limits=c(0,600000)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=10)) + 
  annotate("text",x=2019,y=600000,label="A", size=8)

figS3A
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figureS3A.tiff", width=8, height=8, device='tiff', dpi=700)
remove(figS3A)

#Figure S3B: Total Outpatient Visits per Drug Group per Year 
df <- svyby(~optotv, ~treatGrp+year.x, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df

figS3B <- ggplot(data=df,
                 aes(x=year.x, y=optotv, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Total Number of Outpatient Visits') +
  scale_y_continuous(label=comma, limits=c(0,2500000))+
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=10))+
  annotate("text",x=2019,y=2500000,label="B", size=8)


figS3B
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figureS3B.tiff", width=8, height=8, device='tiff', dpi=700)
remove(figS3B)

#Figure S3C: Total Office-Based Visits per Drug Group per Year 
df <- svyby(~obtotv, ~treatGrp+year.x, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df

figS3C <- ggplot(data=df,
                 aes(x=year.x, y=obtotv, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Total Number of Office-Based Visits') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  scale_y_continuous(label=comma,limits=c(0,20000000)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8))+
  annotate("text",x=2019,y=20000000,label="C", size=8)

figS3C
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figureS3C.tiff", width=8, height=8, device='tiff', dpi=700)
remove(figS3C)

#Figure S3D: Total Emergency Department Visits per Drug Group per Year 
df <- svyby(~ertot, ~treatGrp+year.x, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df

figS3D <- ggplot(data=df,
                 aes(x=year.x, y=ertot, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Total Number of Emergency Department Visits') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  scale_y_continuous(label=comma,limits=c(0,1000000))+
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=10)) +
  annotate("text",x=2019,y=1000000,label="D", size=8)

figS3D
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figureS3D.tiff", width=8, height=8, device='tiff', dpi=700)
remove(figS3D)

#Figure S3E: Total Prescription Medicines per Drug Group per Year 
df <- svyby(~rxtot, ~treatGrp+year.x, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df

figS3E <- ggplot(data=df,
                 aes(x=year.x, y=rxtot, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Total Number of Prescription Medicines') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  scale_y_continuous(label=comma, limits=c(0,35000000))+
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=10)) +
  annotate("text",x=2019,y=35000000,label="E", size=8)


figS3E
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figureS3E.tiff", width=8, height=8, device='tiff', dpi=700)
remove(figS3E)

#Supplemental Figure 4: Sum total expenditure per drug group per year: 

#Figure S4A: Total medical expenditure: 
df <- svyby(~totexp_adj, ~treatGrp+year.x, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df

figS4A <- ggplot(data=df,
                 aes(x=year.x, y=totexp_adj, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Total Medical Expenditure') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8))+
  scale_y_continuous(labels=scales::dollar_format(),limits=c(0,20000000000))+
  annotate("text",x=2019,y=20000000000,label="A", size=8)

figS4A

ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figureS4A.tiff", width=8, height=8, device='tiff', dpi=700)

remove(figS4A)

#Figure S4B: Total inpatient expenditure: 
df <- svyby(~iptexp_adj, ~treatGrp+year.x, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df

figS4B <- ggplot(data=df,
                 aes(x=year.x, y=iptexp_adj, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Total Inpatient Expenditure') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8))+
  scale_y_continuous(labels=scales::dollar_format(), limits=c(0,7000000000)) +
  annotate("text",x=2019,y=7000000000,label="B", size=8)

figS4B
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figureS4B.tiff", width=8, height=8, device='tiff', dpi=700)
remove(figS4B)

#Figure S4C: Total outpatient visit expenditure:
df <- svyby(~optexp_adj, ~treatGrp+year.x, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df

figS4C <- ggplot(data=df,
                 aes(x=year.x, y=optexp_adj, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Total Outpatient Visit Expenditure') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=))+
  scale_y_continuous(labels=scales::dollar_format(), limits=c(0,3500000000))+
  annotate("text",x=2019,y=3500000000,label="C", size=8)

figS4C
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figureS4C.tiff", width=8, height=8, device='tiff', dpi=700)
remove(figS4C)

#Figure S4D: Total office-based visit expenditure:
df <- svyby(~obvexp_adj, ~treatGrp+year.x, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df

figS4D <- ggplot(data=df,
                 aes(x=year.x, y=obvexp_adj, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Total Office-Based Visit Expenditure') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8))+
  scale_y_continuous(labels=scales::dollar_format(), limits=c(0,4500000000))+
  annotate("text",x=2019,y=4500000000,label="D", size=8)

figS4D
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figureS4D.tiff", width=8, height=8, device='tiff', dpi=700)
remove(figS4D)

#Figure S4E: Total emergency department expenditure:
df <- svyby(~ertexp_adj, ~treatGrp+year.x, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df

figS4E <- ggplot(data=df,
                 aes(x=year.x, y=ertexp_adj, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Total Emergency Department Expenditure') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8))+
  scale_y_continuous(labels=scales::dollar_format(), limits=c(0,1000000000))+
  annotate("text",x=2019,y=1000000000,label="E", size=8)

figS4E
ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figureS4E.tiff", width=8, height=8, device='tiff', dpi=700)
remove(figS4E)

#Figure S4F: Total prescription medicine expenditure:
df <- svyby(~rxexp_adj, ~treatGrp+year.x, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df

figS4F <- ggplot(data=df,
                 aes(x=year.x, y=rxexp_adj, group=treatGrp, color=treatGrp)) + 
  geom_line(lwd=0.8) +                                                  #LWD changes line thickness. 
  geom_point() +
  scale_color_manual("Drug Group",values=c("coral3","darkgoldenrod2","darkturquoise","darkolivegreen3","deepskyblue4"),labels=c("Non-Narcotic","Opioid","Opioid+BZD","Opioid+BZD+SMR","Opioid+GABA")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(size=10, color='black'), axis.text.y=element_text(size=15, color='black')) + 
  theme(axis.title.x=element_text(color="black", size=10)) + 
  theme(axis.title.y=element_text(color='black', size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  xlab('Year') +
  ylab('Total Prescription Medicine Expenditure') +
  scale_x_continuous(breaks=c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=8))+
  scale_y_continuous(labels=scales::dollar_format(), limits=c(0,4000000000))+
  annotate("text",x=2019,y=4000000000,label="F", size=8)


figS4F

ggsave(filename="C:/Users/aryse/OneDrive/Desktop/figureS4F.tiff", width=8, height=8, device='tiff', dpi=700)
remove(figS4F)
remove(df)

#Supplemental Table 7: Total Resource Utilization per Group per Year: 
#Inpatient Discharges 
df <- svyby(~ipdis, ~year.x+treatGrp, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#Outpatient Visits 
df <- svyby(~optotv, ~year.x+treatGrp, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#Office-Based Visits 
df <- svyby(~obtotv, ~year.x+treatGrp, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#Emergency Department Visits 
df <- svyby(~ertot, ~year.x+treatGrp, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#Prescription Medications 
df <- svyby(~rxtot, ~year.x+treatGrp, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2


#Supplemental Table 8: Total Expenditure per Group per Year: 
#Total Medical
df <- svyby(~totexp_adj, ~year.x+treatGrp, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2




#Inpatient
df <- svyby(~iptexp_adj, ~year.x+treatGrp, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#Outpatient
df <- svyby(~optexp_adj, ~year.x+treatGrp, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#Office-based visits 
df <- svyby(~obvexp_adj, ~year.x+treatGrp, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#Emergency Department Visits 
df <- svyby(~ertexp_adj, ~year.x+treatGrp, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#Prescription Medications  
df <- svyby(~rxexp_adj, ~year.x+treatGrp, design=study_dsn, svytotal)
df <- subset(df, select=-c(se))
df2 <- reshape(df, idvar="treatGrp", timevar="year.x", direction="wide")
df2

#library(openxlsx)
#write.xlsx(df2, 'table.xlsx')



