install.packages("forestplot")

library (survival)
library(survminer)
library (forestplot)
library(ggplot2)
packageVersion("survminer") 
load("sarcopenia.RData")

#Read clinical data and CT measurements
CTresults <- read.csv("CT Results.csv")
HH_ov_clinical <- read.csv("HH_ov_clinical.csv")
same_index <- Reduce (intersect, list (CTresults$Patient_ID,HH_ov_clinical$PatientID))
HH_ov_clinical <- HH_ov_clinical [HH_ov_clinical$PatientID %in% same_index,]
HH_ov_clinical_results <- merge (HH_ov_clinical, CTresults, by.x = "PatientID", by.y = "Patient_ID", all.x = T)

#SMA cut-off value for sarcopenia
shapiro.test(HH_ov_clinical_results$lite_ai_area_cm2)
SMA_cutoff <- as.numeric (summary (HH_ov_clinical_results$lite_ai_area_cm2)[2])

#Sarcopenia patients identified by SMA
sarcopenia_SMA <- HH_ov_clinical_results [which (HH_ov_clinical_results$lite_ai_area_cm2 < SMA_cutoff), ]
nonsarcopenia_SMA <- HH_ov_clinical_results [which (HH_ov_clinical_results$lite_ai_area_cm2 >= SMA_cutoff), ]
nrow(sarcopenia_SMA)
nrow (nonsarcopenia_SMA)


#SMD cut-off value for sarcopenia 
shapiro.test(HH_ov_clinical_results$lite_ai_density_HU)
SMD_mean <- mean(HH_ov_clinical_results$lite_ai_density_HU)
SMD_sd <- sd (HH_ov_clinical_results$lite_ai_density_HU)
SMD_cutoff <- SMD_mean - SMD_sd

#sarcopenia patients identified by SMD
sarcopenia_SMD <- HH_ov_clinical_results[which(HH_ov_clinical_results$lite_ai_density_HU < SMD_cutoff),]
nonsarcopenia_SMD <- HH_ov_clinical_results[which(HH_ov_clinical_results$lite_ai_density_HU >= SMD_cutoff),]

#Clinical data distribution
#Age 
mean (HH_ov_clinical_results$Age.at.initial.diagnosis)
sd (HH_ov_clinical_results$Age.at.initial.diagnosis)
mean(nonsarcopenia_SMA$Age.at.initial.diagnosis)
sd (nonsarcopenia_SMA$Age.at.initial.diagnosis)
mean (sarcopenia_SMA$Age.at.initial.diagnosis)
sd (sarcopenia_SMA$Age.at.initial.diagnosis)
mean (nonsarcopenia_SMD$Age.at.initial.diagnosis)
sd (nonsarcopenia_SMD$Age.at.initial.diagnosis)
mean (sarcopenia_SMD$Age.at.initial.diagnosis)
sd (sarcopenia_SMD$Age.at.initial.diagnosis)

t.test (nonsarcopenia_SMA$Age.at.initial.diagnosis, sarcopenia_SMA$Age.at.initial.diagnosis)
t.test(nonsarcopenia_SMD$Age.at.initial.diagnosis, sarcopenia_SMD$Age.at.initial.diagnosis)

#FIGO stage
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Stage == 1),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Stage == 2),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Stage == 3),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Stage == 4),])

nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Stage == 1),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Stage == 2),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Stage == 3),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Stage == 4),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Stage == 1),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Stage == 2),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Stage == 3),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Stage == 4),])

#Chi-square test: FIGO stage SMA-sarcopenia: early stage vs. late stage
table.sarc.sma.stage <- matrix (c(14,7,130,42),nrow = 2, ncol = 2)
table.sarc.sma.stage

chisq.test(table.sarc.sma.stage)

nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Stage == 1),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Stage == 2),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Stage == 3),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Stage == 4),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Stage == 1),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Stage == 2),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Stage == 3),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Stage == 4),])

table.sarc.smd.stage <- matrix (c(19,2,139,33),nrow = 2, ncol = 2)
table.sarc.smd.stage
fisher.test(table.sarc.smd.stage)

#Tumour grade
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Grade_di == 0),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Grade_di == 1),])

nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Grade_di == 0),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Grade_di == 1),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Grade_di == 0),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Grade_di == 1),])

table.sma.tum.gd <- matrix(c(143,48,4,1), nrow = 2, ncol = 2)
table.sma.tum.gd
fisher.test(table.sma.tum.gd)

nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Grade_di == 0),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Grade_di == 1),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Grade_di == 0),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Grade_di == 1),])

table.smd.tum.gd <- matrix(c(157,34,4,1), nrow = 2, ncol = 2)
table.smd.tum.gd
fisher.test(table.smd.tum.gd)

#Initial treatment
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Surgery_type == 0),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Surgery_type == 1),])

nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Surgery_type == 0),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Surgery_type == 1),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Surgery_type == 0),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Surgery_type == 1),])
table.sma.treatment <-  matrix(c(123,36,22,13), nrow = 2, ncol = 2)
table.sma.treatment
chisq.test(table.sma.treatment)

nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Surgery_type == 0),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Surgery_type == 1),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Surgery_type == 0),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Surgery_type == 1),])
table.smd.treatment <-  matrix(c(134,25,26,9), nrow = 2, ncol = 2)
table.smd.treatment
chisq.test(table.smd.treatment)

#Surgery outcome
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Surgery.outcome == 0),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Surgery.outcome == 1),])

nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Surgery.outcome == 0),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Surgery.outcome == 1),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Surgery.outcome == 0),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Surgery.outcome == 1),])

table.sma.surg.outcm <-  matrix(c(108,36,29,12), nrow = 2, ncol = 2)
table.sma.surg.outcm
chisq.test(table.sma.surg.outcm)

nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Surgery.outcome == 0),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Surgery.outcome == 1),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Surgery.outcome == 0),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Surgery.outcome == 1),])

table.smd.surg.outcm <-  matrix(c(118,26,33,8), nrow = 2, ncol = 2)
table.smd.surg.outcm
chisq.test(table.smd.surg.outcm)

#Chemotherapy outcome
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Primary.chemotherapy.outcome == "complete.response"),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Primary.chemotherapy.outcome == "partial.response"),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Primary.chemotherapy.outcome == "stable.disease"),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Primary.chemotherapy.outcome == "progressive.disease"),])

nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Primary.chemotherapy.outcome == "complete.response"),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Primary.chemotherapy.outcome == "partial.response"),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Primary.chemotherapy.outcome == "stable.disease"),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Primary.chemotherapy.outcome == "progressive.disease"),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Primary.chemotherapy.outcome == "complete.response"),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Primary.chemotherapy.outcome == "partial.response"),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Primary.chemotherapy.outcome == "stable.disease"),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Primary.chemotherapy.outcome == "progressive.disease"),])

table.sma.chemo <- matrix (c(54,19,21,8,5,4,5,0), nrow = 2, ncol = 4)
table.sma.chemo
fisher.test(table.sma.chemo)

nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Primary.chemotherapy.outcome == "complete.response"),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Primary.chemotherapy.outcome == "partial.response"),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Primary.chemotherapy.outcome == "stable.disease"),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Primary.chemotherapy.outcome == "progressive.disease"),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Primary.chemotherapy.outcome == "complete.response"),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Primary.chemotherapy.outcome == "partial.response"),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Primary.chemotherapy.outcome == "stable.disease"),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Primary.chemotherapy.outcome == "progressive.disease"),])

table.smd.chemo <- matrix (c(58,15,24,5,6,3,4,1), nrow = 2, ncol = 4)
table.smd.chemo
fisher.test(table.smd.chemo)

#CA125
summary (HH_ov_clinical_results$Ca125)
summary (nonsarcopenia_SMA$Ca125)
shapiro.test(nonsarcopenia_SMA$Ca125)
summary (sarcopenia_SMA$Ca125)
shapiro.test(sarcopenia_SMA$Ca125)
wilcox.test(nonsarcopenia_SMA$Ca125, sarcopenia_SMA$Ca125)

summary (nonsarcopenia_SMD$Ca125)
shapiro.test(nonsarcopenia_SMD$Ca125)
summary (sarcopenia_SMD$Ca125)
shapiro.test(sarcopenia_SMD$Ca125)
wilcox.test(nonsarcopenia_SMD$Ca125, sarcopenia_SMD$Ca125)

#Tumour content
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Cellularity == 0),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$Cellularity == 1),])

nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Cellularity == 0),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$Cellularity == 1),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Cellularity == 0),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$Cellularity == 1),])

table.sma.celluarity <- matrix (c(80,27,45,14), nrow = 2, ncol = 2)
table.sma.celluarity
chisq.test(table.sma.celluarity)

nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Cellularity == 0),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$Cellularity == 1),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Cellularity == 0),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$Cellularity == 1),])

table.smd.celluarity <- matrix (c(90,17,51,8), nrow = 2, ncol = 2)
table.smd.celluarity
chisq.test(table.smd.celluarity)

#SMA
summary (HH_ov_clinical_results$lite_ai_area_cm2)
summary (nonsarcopenia_SMA$lite_ai_area_cm2)
summary (sarcopenia_SMA$lite_ai_area_cm2)
wilcox.test(nonsarcopenia_SMA$lite_ai_area_cm2, sarcopenia_SMA$lite_ai_area_cm2)

summary (nonsarcopenia_SMD$lite_ai_area_cm2)
summary (sarcopenia_SMD$lite_ai_area_cm2)
wilcox.test(nonsarcopenia_SMD$lite_ai_area_cm2, sarcopenia_SMD$lite_ai_area_cm2)

#SMD
summary (HH_ov_clinical_results$lite_ai_density_HU)
sd (HH_ov_clinical_results$lite_ai_density_HU)

summary (nonsarcopenia_SMA$lite_ai_density_HU)
sd (nonsarcopenia_SMA$lite_ai_density_HU)
summary (sarcopenia_SMA$lite_ai_density_HU)
sd (sarcopenia_SMA$lite_ai_density_HU)
t.test(nonsarcopenia_SMA$lite_ai_density_HU, sarcopenia_SMA$lite_ai_density_HU)

summary (nonsarcopenia_SMD$lite_ai_density_HU)
sd (nonsarcopenia_SMD$lite_ai_density_HU)
summary (sarcopenia_SMD$lite_ai_density_HU)
sd (sarcopenia_SMD$lite_ai_density_HU)
t.test(nonsarcopenia_SMD$lite_ai_density_HU, sarcopenia_SMD$lite_ai_density_HU)

#OS events
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$OS.event == 0),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$OS.event == 1),])

nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$OS.event == 0),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$OS.event == 1),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$OS.event == 0),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$OS.event == 1),])
table.sma.os <- matrix(c(59,15,87,34), nrow = 2, ncol = 2)
table.sma.os
chisq.test(table.sma.os)

nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$OS.event == 0),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$OS.event == 1),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$OS.event == 0),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$OS.event == 1),])
table.smd.os <- matrix(c(60,14,100,21), nrow = 2, ncol = 2)
table.smd.os
chisq.test(table.smd.os)

#Survival time
shapiro.test(HH_ov_clinical_results[which(HH_ov_clinical_results$OS.event == 1),]$Overall.survival..days.)
summary (HH_ov_clinical_results[which(HH_ov_clinical_results$OS.event == 1),]$Overall.survival..days.)
summary (nonsarcopenia_SMA[which(nonsarcopenia_SMA$OS.event == 1),]$Overall.survival..days.)
summary (sarcopenia_SMA[which(sarcopenia_SMA$OS.event == 1),]$Overall.survival..days.)
wilcox.test(nonsarcopenia_SMA[which(nonsarcopenia_SMA$OS.event == 1),]$Overall.survival..days., sarcopenia_SMA[which(sarcopenia_SMA$OS.event == 1),]$Overall.survival..days.)

summary (nonsarcopenia_SMD[which(nonsarcopenia_SMD$OS.event == 1),]$Overall.survival..days.)
summary (sarcopenia_SMD[which(sarcopenia_SMD$OS.event == 1),]$Overall.survival..days.)
wilcox.test(nonsarcopenia_SMD[which(nonsarcopenia_SMD$OS.event == 1),]$Overall.survival..days., sarcopenia_SMD[which(sarcopenia_SMD$OS.event == 1),]$Overall.survival..days.)

#PFS event
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$PFS.event == 0),])
nrow (HH_ov_clinical_results[which(HH_ov_clinical_results$PFS.event == 1),])

nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$PFS.event == 0),])
nrow (nonsarcopenia_SMA[which(nonsarcopenia_SMA$PFS.event == 1),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$PFS.event == 0),])
nrow (sarcopenia_SMA[which(sarcopenia_SMA$PFS.event == 1),])
table.sma.pfs <- matrix(c(68,26,65,22), nrow = 2, ncol = 2)
table.sma.pfs
chisq.test(table.sma.pfs)

nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$PFS.event == 0),])
nrow (nonsarcopenia_SMD[which(nonsarcopenia_SMD$PFS.event == 1),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$PFS.event == 0),])
nrow (sarcopenia_SMD[which(sarcopenia_SMD$PFS.event == 1),])
table.smd.pfs <- matrix(c(74,20,73,14), nrow = 2, ncol = 2)
table.smd.pfs
chisq.test(table.smd.pfs)

#PFS time 
shapiro.test(HH_ov_clinical_results[which(HH_ov_clinical_results$PFS.event == 1),]$Progression.free.survival..days.)
summary (HH_ov_clinical_results[which(HH_ov_clinical_results$PFS.event == 1),]$Progression.free.survival..days.)
summary (nonsarcopenia_SMA[which(nonsarcopenia_SMA$PFS.event == 1),]$Progression.free.survival..days.)
summary (sarcopenia_SMA[which(sarcopenia_SMA$PFS.event == 1),]$Progression.free.survival..days.)
wilcox.test(nonsarcopenia_SMA[which(nonsarcopenia_SMA$PFS.event == 1),]$Progression.free.survival..days., sarcopenia_SMA[which(sarcopenia_SMA$PFS.event == 1),]$Progression.free.survival..days.)

summary (nonsarcopenia_SMD[which(nonsarcopenia_SMD$PFS.event == 1),]$Progression.free.survival..days.)
summary (sarcopenia_SMD[which(sarcopenia_SMD$PFS.event == 1),]$Progression.free.survival..days.)
wilcox.test(nonsarcopenia_SMD[which(nonsarcopenia_SMD$PFS.event == 1),]$Progression.free.survival..days., sarcopenia_SMD[which(sarcopenia_SMD$PFS.event == 1),]$Progression.free.survival..days.)

##Survival analysis
#Univariate cox regression
#OS and PFS object
os.object <- Surv(HH_ov_clinical_results$Overall.survival..days., HH_ov_clinical_results$OS.event)
pfs.object <- Surv(HH_ov_clinical_results$Progression.free.survival..days., HH_ov_clinical_results$PFS.event)

#Age
cox.age.os <- coxph (os.object~Age.at.initial.diagnosis, data = HH_ov_clinical_results)
cox.age.pfs <- coxph (pfs.object~Age.at.initial.diagnosis, data = HH_ov_clinical_results)
summary(cox.age.os)
summary (cox.age.pfs)

#FIGO stage
cox.figo.os <- coxph (os.object~Stage, data = HH_ov_clinical_results)
cox.figo.pfs <- coxph (pfs.object~Stage, data = HH_ov_clinical_results)
summary (cox.figo.os)
summary (cox.figo.pfs)

#Tumour grade
cox.grade.os <- coxph (os.object~Grade_di, data = HH_ov_clinical_results)
cox.grade.pfs <- coxph (pfs.object~Grade_di, data = HH_ov_clinical_results)
summary (cox.grade.os)
summary (cox.grade.pfs)

#Surgery outcome
cox.surg.os <- coxph (os.object~Surgery.outcome, data = HH_ov_clinical_results)
cox.surg.pfs <- coxph (pfs.object~Surgery.outcome, data = HH_ov_clinical_results)
summary (cox.surg.os)
summary (cox.surg.pfs)

#Chemotherapy outcome
#Divide to complete chemotherapy and incomplete chemotherapy
HH_ov_clinical_results$chemo.di <- NA
for (i in 1:nrow (HH_ov_clinical_results)){
  if (!is.na(HH_ov_clinical_results$Primary.chemotherapy.outcome[i])){
  if (HH_ov_clinical_results$Primary.chemotherapy.outcome[i] == "complete.response"){
    HH_ov_clinical_results$chemo.di[i] <- 0
  }
  else if (HH_ov_clinical_results$Primary.chemotherapy.outcome[i] == "partial.response" ||
           HH_ov_clinical_results$Primary.chemotherapy.outcome[i] == "stable.disease" ||
           HH_ov_clinical_results$Primary.chemotherapy.outcome[i] == "progressive.disease"){
    HH_ov_clinical_results$chemo.di[i] <- 1
  }
  }
}

cox.chemo.os <- coxph (os.object~chemo.di, data = HH_ov_clinical_results)
cox.chemo.pfs <- coxph (pfs.object~chemo.di, data = HH_ov_clinical_results)
summary (cox.chemo.os)
summary (cox.chemo.pfs)

#Tumour cellularity
cox.cellularity.os <- coxph(os.object~Cellularity, data = HH_ov_clinical_results)
cox.cellularity.pfs <- coxph(pfs.object~Cellularity, data = HH_ov_clinical_results)
summary (cox.cellularity.os)
summary (cox.cellularity.pfs)

#Initial treatment
cox.ini.os <- coxph (os.object~Surgery_type, data = HH_ov_clinical_results)
cox.ini.pfs <- coxph (pfs.object~Surgery_type, data = HH_ov_clinical_results)
summary (cox.ini.os)
summary (cox.ini.pfs)

#Skeletal muscle area
cox.sma.os <- coxph (os.object~lite_ai_area_cm2, data = HH_ov_clinical_results)
cox.sma.pfs <- coxph (pfs.object~lite_ai_area_cm2, data = HH_ov_clinical_results)
summary (cox.sma.os)
summary (cox.sma.pfs)

#Skeletal muscle density
cox.smd.os <- coxph (os.object~lite_ai_density_HU, data = HH_ov_clinical_results)
cox.smd.pfs <- coxph (pfs.object~lite_ai_density_HU, data = HH_ov_clinical_results)
summary (cox.smd.os)
summary (cox.smd.pfs)

#SMA-sarcopenia
HH_ov_clinical_results$sarcopenia_SMA <- NA
for (i in 1:nrow(HH_ov_clinical_results)){
  if (HH_ov_clinical_results$lite_ai_area_cm2[i] < SMA_cutoff){
    HH_ov_clinical_results$sarcopenia_SMA[i] <- 1
  }else {
    HH_ov_clinical_results$sarcopenia_SMA[i] <- 0
  }
}

cox.smasarc.os <- coxph (os.object~sarcopenia_SMA, data = HH_ov_clinical_results)
cox.smasarc.pfs <- coxph (pfs.object~sarcopenia_SMA, data = HH_ov_clinical_results)
summary (cox.smasarc.os)
summary (cox.smasarc.pfs)

#SMD-sarcopenia
HH_ov_clinical_results$sarcopenia_SMD <- NA
for (i in 1:nrow(HH_ov_clinical_results)){
  if (HH_ov_clinical_results$lite_ai_density_HU[i] < SMD_cutoff){
    HH_ov_clinical_results$sarcopenia_SMD[i] <- 1
  }else {
    HH_ov_clinical_results$sarcopenia_SMD[i] <- 0
  }
}

cox.smdsarc.os <- coxph (os.object~sarcopenia_SMD, data = HH_ov_clinical_results)
cox.smdsarc.pfs <- coxph (pfs.object~sarcopenia_SMD, data = HH_ov_clinical_results)
summary (cox.smdsarc.os)
summary (cox.smdsarc.pfs)

#Multivariate analysis
#OS:FIGO stage+surgery outcome+chemotherapy outcome+SMA
cox.multi.os <- coxph (os.object~Stage+Surgery.outcome+chemo.di+lite_ai_area_cm2, data = HH_ov_clinical_results)
summary (cox.multi.os)

#PFS: FIGO stage+surgery outcome+chemotherapy outcome+SMA+SMD-sarcopenia
cox.multi.pfs <- coxph (pfs.object~Stage+Surgery.outcome+chemo.di+lite_ai_area_cm2+sarcopenia_SMD, data = HH_ov_clinical_results)
summary (cox.multi.pfs)

#Forest plot of HR in Multivariate Cox Regression analysis
cox_reg_os <- read.csv("cox_reg.csv", header = F)
png("Cox forest plot.png", res = 300, width = 2500, height = 2000)
forestplot(labeltext = as.matrix (cox_reg_os[,1:3]),
           mean = cox_reg_os$V4,
           lower = cox_reg_os$V5,
           upper = cox_reg_os$V6,
           is.summary = c(T,T,F,F,F,F,T,F,F,F,F,F),
           hrzl_lines = list ("2" = gpar (col = "black", lwd = 2), "13" = gpar(col = "black", lwd = "2")),
           xlab = "Hazard Ratio",
           zero = 1,
           align=c("l"),
           lwd.zero = 2,
           clip = c(0.8,3),
           xticks = c (1,1.5,2.5,3,3.5),
           txt_gp = fpTxtGp(label = gpar(fontfamily = "serif"),ticks = gpar(cex = 1), xlab = gpar(cex = 1), cex = 1),
           boxsize = 0.3,
           graph.pos = 4,
           fn.ci_norm = fpDrawDiamondCI,
           col = fpColors(line = "black", box = "red", zero = "cadetblue",),
           lty.ci = 7,
           lwd.ci = 2,
           colgap = unit(0.05,"npc"),
           graphwidth = unit (0.2,"npc"),
           family="sans",font=2
          )
dev.off()

which (HH_ov_clinical_results$sarcopenia_SMD == 1)
#KM plot
fit.os.sma <- survfit(os.object~sarcopenia_SMA, data = HH_ov_clinical_results)
#SMA-sarcopenia OS KM plot
fit.os.sma.plot <- ggsurvplot(fit.os.sma, risk.table = TRUE, pval = TRUE,surv.median.line = "hv",
                              xlab = "Survival time (days)",legend = c(0.8,0.9), conf.int = TRUE,
                              legend.title = "", 
                              legend.labs = c("Non-sarcopenia", "Sarcopenia"),
                              break.x.by = 1000,
                              palette = "aaas",
                              data = HH_ov_clinical_results)
png(filename = "SMA_sarcopenia_os.png",res = 300, width = 2000, height = 2000)
fit.os.sma.plot
dev.off()
survdiff (os.object~sarcopenia_SMA, data = HH_ov_clinical_results)

#SMD-sarcopenia OS KM plot
fit.os.smd <- survfit(os.object~sarcopenia_SMD, data = HH_ov_clinical_results)
survdiff (os.object~sarcopenia_SMD, data = HH_ov_clinical_results)
fit.os.smd.plot <- ggsurvplot(fit.os.smd, risk.table = TRUE, pval = TRUE,surv.median.line = "hv",
                              xlab = "Survival time (days)",legend = c(0.8,0.9),conf.int = TRUE,
                              legend.title = "", 
                              legend.labs = c("Non-sarcopenia", "Sarcopenia"),
                              break.x.by = 1000,
                              palette = "aaas")
png(filename = "SMD_sarcopenia_os.png",res = 300, width = 2000, height = 2000)
fit.os.smd.plot
dev.off()

#SMA-sarcopenia PFS KM plot
fit.pfs.sma <- survfit(pfs.object~sarcopenia_SMA, data = HH_ov_clinical_results)
fit.pfs.sma.plot <- ggsurvplot(fit.pfs.sma, risk.table = TRUE, pval = TRUE,surv.median.line = "hv",
                              xlab = "Progression-free time (days)",legend = c(0.8,0.9),conf.int = TRUE,
                              legend.title = "", 
                              legend.labs = c("Non-sarcopenia", "Sarcopenia"),
                              break.x.by = 1000,
                              palette = "aaas")

png(filename = "SMA_sarcopenia_pfs.png",res = 300, width = 2000, height = 2000)
fit.pfs.sma.plot
dev.off()
survdiff (pfs.object~sarcopenia_SMA, data = HH_ov_clinical_results)

#SMD-sarcopenia PFS KM plot
fit.pfs.smd <- survfit(pfs.object~sarcopenia_SMD, data = HH_ov_clinical_results)
fit.pfs.smd.plot <- ggsurvplot(fit.pfs.smd, risk.table = TRUE, pval = TRUE,surv.median.line = "hv",
                               xlab = "Progression-free time (days)",legend = c(0.8,0.9),conf.int = TRUE,
                               legend.title = "", 
                               legend.labs = c("Non-sarcopenia", "Sarcopenia"),
                               break.x.by = 1000,
                               palette = "aaas")
png(filename = "SMD_sarcopenia_pfs.png",res = 300, width = 2000, height = 2000)
fit.pfs.smd.plot
dev.off()
survdiff (pfs.object~sarcopenia_SMD, data = HH_ov_clinical_results)

save.image("sarcopenia.RData")
