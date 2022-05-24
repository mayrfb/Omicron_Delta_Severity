machine<-"P:/ORD_Shaikh_202011058D/"
waning_ve3 <- "projects/research/severity/"

# main patient file, curated via dsp1_make_code.R script
library(lubridate)
library(dplyr)
library(data.table)

set.seed(12345)

dsp11 <- readRDS(paste0(machine, waning_ve3, "dsp11_jan2022_v2.RDS"))
dsp11$inftype <- ifelse(dsp11$pos_test_date>=18901 & dsp11$pos_test_date <=18973,"delta",
                           ifelse(dsp11$pos_test_date>=18986 & dsp11$pos_test_date <=19010, "omicron", "other"))

dsp12 <- filter(dsp11, inftype %in% c("delta", "omicron"))

#fix covariate coding
dsp12$gender<-as.factor(dsp12$Gender)
dsp12$cci_cat<-dsp12$CCI_2yrs
dsp12$cci_cat<-as.factor(ifelse(dsp12$cci_cat>2,3,dsp12$cci_cat))
dsp12$race_cat<-ifelse(dsp12$Race %in% c("WHITE", "WHITE NOT OF HISP ORIG"),"white",
                            ifelse(dsp12$Race=="BLACK OR AFRICAN AMERICAN","black","other"))
dsp12$race_cat[is.na(dsp12$race_cat)]<-"other"
dsp12$race_cat<-as.factor(dsp12$race_cat)
dsp12$age_cat<-as.factor(ifelse(dsp12$Age<40,0,
                                     ifelse(dsp12$Age<50,1,
                                            ifelse(dsp12$Age<60,2,
                                                   ifelse(dsp12$Age<70,3,
                                                          ifelse(dsp12$Age<80,4,
                                                                 ifelse(dsp12$Age<90,5,6)))))))
dsp12$madi_cat<-as.factor(cut(dsp12$madi,breaks=c(0,seq(40,130,10),200),labels=1:11,include.lowest=TRUE))

dsp12$Series2Sta3n <- as.factor(dsp12$Series2Sta3n)
dsp12$calweekSeries2Date <- lubridate::week(dsp12$Series2Date)
dsp12$calweekSeries3Date <- lubridate::week(dsp12$Series3Date)

#identifying vaccination 
dsp13 <- filter(dsp12, Series1Vaccine=="Pfizer" & Series2Vaccine=="Pfizer" & Series1Date >="2020-12-11" | (Series1Vaccine=="Moderna"& Series2Vaccine =="Moderna" & Series1Date>="2020-12-18") | is.na(Series1Vaccine) | is.na(Series2Vaccine) | is.na(Series3Vaccine))
dsp13$booster <- ifelse(dsp13$Series3Date >= "2021-09-22",1,0)

dsp13$Series1Date2 <- as.numeric(as_date(dsp13$Series1Date))
dsp13$Series2Date2 <- as.numeric(as_date(dsp13$Series2Date))
dsp13$Series3Date2 <- as.numeric(as_date(dsp13$Series3Date))

# replacing missing dates with 19357 (2022-12-31)

dsp13$Series1Date2[is.na(dsp13$Series1Date2)] <- 19357
dsp13$Series2Date2[is.na(dsp13$Series2Date2)] <- 19357
dsp13$Series3Date2[is.na(dsp13$Series3Date2)] <- 19357

dsp13$timevacc2toinfec <- dsp13$pos_test_date-dsp13$Series2Date2
dsp13$timevacc3toinfec <- dsp13$pos_test_date-dsp13$Series3Date2

dsp13$vaccines <- ifelse(dsp13$Series3Date2<=dsp13$pos_test_date & dsp13$Series2Date2 <=dsp13$pos_test_date & dsp13$Series1Date2 <= dsp13$pos_test_date, 3,
                         ifelse(dsp13$Series2Date2 <=dsp13$pos_test_date & dsp13$Series1Date2<=dsp13$pos_test_date & dsp13$Series3Date2 > dsp13$pos_test_date,2,
                                ifelse(dsp13$Series1Date2 <= dsp13$pos_test_date & dsp13$Series2Date2>dsp13$pos_test_date & dsp13$Series3Date2 > dsp13$pos_test_date,1,
                                       ifelse(dsp13$Series3Date2 > dsp13$pos_test_date & dsp13$Series2Date2 > dsp13$pos_test_date & dsp13$Series1Date2>dsp13$pos_test_date, 0, NA))))

dsp14 <- filter(dsp13, Series3Date2 >=18892) #booster start 9/22/2021
                         
# replace NA for calweekSeries2Date to 0 for matching
dsp14$calweekSeries2Date[is.na(dsp14$calweekSeries2Date)] <- 0
dsp14$calweekSeries3Date[is.na(dsp14$calweekSeries3Date)] <- 0

library(MatchIt)
library(officer)
library(gtsummary)
library(kableExtra)
library(flextable)

# read in comorbidity data
charlson <- readRDS(paste0(machine,waning_ve3, "charlson_dsp11.RDS"))
elix <- readRDS(paste0(machine,waning_ve3, "elixh_dsp11.RDS"))
colnames(charlson)[1] <- "PatientICN"
colnames(elix)[1] <- "PatientICN"

dsp14 <- merge(dsp14,elix, by="PatientICN", all.x=T)
dsp14 <- merge(dsp14,charlson, by="PatientICN", all.x=T)
dsp14$cci_cat2 <- as.factor(ifelse(dsp14$score.x <=2, dsp14$score.x, 3))
dsp14$cci_cat2[is.na(dsp14$cci_cat2)] <- 0

#create matches for omicron and delta
m.out <- matchit(inftype ~ age_cat +  race_cat+ gender+ cci_cat2 + calweekSeries2Date + Sta3n + vaccines + madi_cat, 
                 data = dsp14[dsp14$multvis==1,], 
                 method = "cem",
                 k2k = TRUE,
                 k2k.method=NULL)

odmdata <- match.data(m.out)
saveRDS(odmdata, paste0(machine,waning_ve3, "odmdata_revisions_cohort_jan22.RDS"))

#odmdata <- readRDS(paste0(machine, waning_ve3, "odmdata_revisions_cohort_jan22.RDS"))

#pic <- readRDS(paste0(machine,waning_ve2,"preindexconditions_jan2022.RDS"))
#odmdata2 <- merge(odmdata, pic[, c(1,18,22,23,30,36,39,44,60,61,95)], by="PatientICN", all.x=T, all.y=F)

# read in from disk
#odmdata2 <- readRDS(paste0(machine, "odmcohort_jan2022.RDS"))
odmdata_mrn <- as.list(odmdata$PatientICN)
odm_old <- readRDS(paste0(machine, waning_ve3, "odmcohort_jan2022.RDS"))
odm_old <- odm_old[odm_old$PatientICN %in% odmdata2_mrn,]

odmdata %>%
  dplyr::select(Age, age_cat, inftype, gender, race_cat, CCI_2yrs, chf.x, hypunc, copd, diabunc, rf, lymph, canc, metacanc.x, rheumd.x, alcohol, drug, depre, score.x) %>%
  tbl_summary(
    by=inftype,
    sort = all_categorical() ~ "frequency",
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = c(all_continuous2() ~2,all_categorical()~1),
    label = list(age_cat ~ "Age category", gender ~ "Sex", race_cat ~ "Race"),
    missing ="no",
    type = list(age_cat ~ "categorical", gender ~ "categorical", race_cat ~ "categorical")) %>%
  as_flex_table() %>%
  autofit() %>%
  flextable::save_as_docx(odmdata2, path = paste0(machine, waning_ve3,"table1.docx"))

length(odmdata2$vaccines[odmdata2$vaccines==2 & odmdata2$timevacc2toinfec >=90 & odmdata2$inftype=="omicron"])

# timing of infection related to vaccination
table(odmdata$inftype[odmdata$vaccines==2 & between(odmdata$pos_test_date-odmdata$Series2Date2, -1000, 14)])
table(odmdata$inftype[odmdata$vaccines==3 & (odmdata$pos_test_date > odmdata$Series2Date2+14 & odmdata$pos_test_date <odmdata$Series3Date2+14)])
table(odmdata$inftype[odmdata$vaccines==3 & (odmdata$pos_test_date > odmdata$Series3Date2+14)])

### REDO
# vaccination status at time of infection
table(odmdata$vaccines[odmdata$vaccines==0], odmdata$inftype[odmdata$vaccines==0])
table(odmdata$vaccines[odmdata$vaccines>0], odmdata$inftype[odmdata$vaccines>0])
table(odmdata$vaccines[odmdata$vaccines==1], odmdata$inftype[odmdata$vaccines==1])
table(odmdata$inftype[odmdata$vaccines==2 & odmdata$timevacc2toinfec <84], odmdata$vaccines[odmdata$vaccines==2 & odmdata$timevacc2toinfec <84])
table(odmdata$inftype[odmdata$vaccines==2 & odmdata$timevacc2toinfec >=84], odmdata$vaccines[odmdata$vaccines==2 & odmdata$timevacc2toinfec >=84])
table(odmdata$inftype[odmdata$vaccines==3 & odmdata$timevacc3toinfec >= 14], odmdata$vaccines[odmdata$vaccines==3 & odmdata$timevacc3toinfec >= 14])


# Infectionstatus related to vaccination


odmdata$infecprevacc2 <- ifelse(odmdata$vaccines %in% c(1,2) & odmdata$pos_test_date < odmdata$Series2Date2+14,1,0)
odmdata$infecprevacc3 <- ifelse(odmdata$vaccines %in% c(2,3) & odmdata$pos_test_date >= odmdata$Series2Date2+14 & odmdata$pos_test_date < odmdata$Series3Date2+14,1,0)
odmdata$infecpostvacc3 <-ifelse(odmdata$vaccines %in% c(3) & odmdata$pos_test_date >= odmdata$Series3Date2+14,1,0)


table(odmdata$inftype[odmdata$vaccines %in% c(1,2)], odmdata$infecprevacc2[odmdata$vaccines %in% c(1,2)])
table(odmdata$inftype[odmdata$vaccines %in% c(2,3)], odmdata$infecprevacc3[odmdata$vaccines %in% c(2,3)])
table(odmdata$inftype[odmdata$vaccines %in% c(3)], odmdata$infecpostvacc3[odmdata$vaccines %in% c(3)])

# defining disease outcomes
odmdata$infec <- ifelse(odmdata$pos_test_date<19357 & odmdata$covid_hosp_date >=19357 & odmdata$death_time >odmdata$pos_test_date+28,1,0)
odmdata$hosp <- ifelse(odmdata$covid_hosp_date <19357 & odmdata$covid_icu_date>=19357 & odmdata$death_time >odmdata$pos_test_date+28,1,0)
odmdata$crit <- ifelse(odmdata$pos_test_date < 19357 & (odmdata$covid_icu_date < 19357 | odmdata$death_time <19357 & odmdata$death_time <=odmdata$pos_test_date+28),1,0)

### END REDO

######################################################################################################################
# vaccination status at time of infection
table(odmdata$inftype[odmdata$vaccines==2 & odmdata$timevacc2toinfec >=84], odmdata$vaccines[odmdata$vaccines==2 & odmdata$timevacc2toinfec >=84])
table(odmdata$inftype[odmdata$vaccines==3 & odmdata$timevacc3toinfec >= 14], odmdata$vaccines[odmdata$vaccines==3 & odmdata$timevacc3toinfec >= 14])


# Infectionstatus related to vaccination

odmdata$infecprevacc2 <- ifelse(odmdata$pos_test_date < odmdata$Series2Date2+14,1,0)
odmdata$infecprevacc3 <- ifelse(odmdata$pos_test_date >= odmdata$Series2Date2+14 & odmdata$pos_test_date < odmdata$Series3Date2+14,1,0)
odmdata$infecpostvacc3 <- ifelse(odmdat2$pos_test_date >= odmdata$Series3Date2+14,1,0)


table(odmdata$inftype[odmdata$vaccines %in% c(1,2)], odmdata$infecprevacc2[odmdata$vaccines %in% c(1,2)])
table(odmdata$inftype[odmdata$vaccines %in% c(1,2,3)], odmdata$infecprevacc3[odmdata$vaccines %in% c(1,2,3)])
table(odmdata$inftype[odmdata$vaccines %in% c(1,2,3)], odmdata$infecpostvacc3[odmdata$vaccines %in% c(1,2,3)])
#########################################################################################################################

saveRDS(odmdata, paste0(machine, waning_ve3, "odmcohort_jan2022_revisions_ver2.RDS"))
odmdata <- readRDS(paste0(machine,waning_ve3, "odmcohort_jan2022_revisions_ver2.RDS"))

# table 2
table(odmdata$inftype, odmdata$infec)
table(odmdata$inftype, odmdata$hosp)
table(odmdata$inftype, odmdata$crit)

# sensitivity analysis by omicron:delta daily admission ratio (definition in line 668-669)
#ratio <3.25
table(odmdata$inftype[odmdata$Sta3n %in% hosp_low_ratio], odmdata$infec[odmdata$Sta3n %in% hosp_low_ratio])
table(odmdata$inftype[odmdata$Sta3n %in% hosp_low_ratio], odmdata$hosp[odmdata$Sta3n %in% hosp_low_ratio])
table(odmdata$inftype[odmdata$Sta3n %in% hosp_low_ratio], odmdata$crit[odmdata$Sta3n %in% hosp_low_ratio])

# ratio >=3.25
table(odmdata$inftype[odmdata$Sta3n %in% hosp_high_ratio], odmdata$infec[odmdata$Sta3n %in% hosp_high_ratio])
table(odmdata$inftype[odmdata$Sta3n %in% hosp_high_ratio], odmdata$hosp[odmdata$Sta3n %in% hosp_high_ratio])
table(odmdata$inftype[odmdata$Sta3n %in% hosp_high_ratio], odmdata$crit[odmdata$Sta3n %in% hosp_high_ratio])

# table 3
table(odmdata$infec[odmdata$infecprevacc3==1], odmdata$inftype[odmdata$infecprevacc3==1])
table(odmdata$hosp[odmdata$infecprevacc3==1], odmdata$inftype[odmdata$infecprevacc3==1])
table(odmdata$crit[odmdata$infecprevacc3==1], odmdata$inftype[odmdata$infecprevacc3==1])

table(odmdata$infec[odmdata$infecpostvacc3==1], odmdata$inftype[odmdata$infecpostvacc3==1])
table(odmdata$hosp[odmdata$infecpostvacc3==1], odmdata$inftype[odmdata$infecpostvacc3==1])
table(odmdata$crit[odmdata$infecpostvacc3==1], odmdata$inftype[odmdata$infecpostvacc3==1])

odmdata$hosp_crit <- as.factor(ifelse(odmdata$hosp==1 | odmdata$crit==1, 1,0))


odmdata$out_cat <- ifelse(odmdata$hosp==1,2,
                           ifelse(odmdata$crit==1,3,
                                  ifelse(odmdata$infec==1,1,99)))


# clr models
odmdata$inftype2 <- as.factor(ifelse(odmdata$inftype=="delta", 0,1))
odmdata$vac_cat <- as.factor(ifelse(odmdata$vaccines==0, 0,
                           ifelse(odmdata$vaccines==1, 1,
                                  ifelse(odmdata$vaccines==2 & odmdata$timevacc2toinfec < 84, 2,
                                         ifelse(odmdata$vaccines==2 & odmdata$timevacc2toinfec >=84, 3,
                                                ifelse(odmdata$vaccines==3 & odmdata$timevacc3toinfec >14, 4, 99))))))


odmdata$age_cat2<- as.factor(ifelse(odmdata$Age<40,0,
                                ifelse(odmdata$Age<50,1,
                                       ifelse(odmdata$Age<60,2,
                                              ifelse(odmdata$Age<70,3,
                                                     ifelse(odmdata$Age<80,4,
                                                            ifelse(odmdata$Age<90,5,6)))))))

odmdata$sex<- as.factor(ifelse(odmdata$Gender=="F", 0,1))

odmdata$race_cat2 <- as.factor(ifelse(odmdata$race_cat=="white", 0,
                             ifelse(odmdata$race_cat=="black", 1,2)))
library(gtools)
odmdata$madi_cat2 <- quantcut(odmdata$madi, q=4)
odmdata$madi_cat3 <- quantcut(odmdata$madi, q=3)
lr_odmdata <- odmdata[, c("inftype2", "vac_cat", "age_cat2", "sex", "race_cat2", "cci_cat2", "madi_cat3", "Sta3n", "hosp", "crit", "infec", "hosp_crit")]

# label variables
library(survival)
library(lme4)
library(stats)
library(MASS)
library(sjPlot)
library(expss)
library(sjlabelled)
library(ggplot2)
library(cowplot)

odmdata2lr <- odmdata[odmdata$vac_cat %in% c(0,2,3,4),]

odmdata2lr$Variant <- factor(odmdata2lr$inftype2, labels = c("Delta variant", "Omicron variant"))
odmdata2lr$Vaccination_status <- factor(odmdata2lr$vac_cat, labels = c("Unvaccinated","2nd vaccine dose < 3 months", "2nd vaccine dose > 3 months", "3rd vaccine dose >14 days"))
odmdata2lr$Age_group <- factor(odmdata2lr$age_cat2, labels = c("Age 18-39", "Age 40-49", "Age 50-59", "Age 60-69", "Age 70-79", "Age 80-89", "Age 90+"))
odmdata2lr$Sex <- factor(odmdata2lr$sex, labels = c("Female, Male"))
odmdata2lr$Racial_group <- factor(odmdata2lr$race_cat2,labels = c("White", "Black", "Other"))
odmdata2lr$CCI_group <- factor(odmdata2lr$cci_cat2, labels = c("No comorbidities", "1 comorbidity", "2 comorbidites", "3+ comorbidities"))
odmdata2lr$Madi_group <- factor(odmdata2lr$madi_cat2, labels = c("30-85", "86-93", "94-99", "100-172"))
odmdata2lr$VAMC <- factor(odmdata2lr$Sta3n)

set_label(odmdata2lr$hosp) <- "Moderate disease"
set_label(odmdata2lr$crit) <- "Severe/Critical disease"



model <- glm(hosp ~ Variant + Vaccination_status + Age_group + Sex + Racial_group + CCI_group+Madi_group+Sta3n,
             family = binomial(link="logit"),
             data = odmdata2lr)

exp(cbind(coef(model), confint(model)))
summary(model)

model2a <- glm(crit ~ Variant + Vaccination_status + Age_group + Sex + Racial_group + CCI_group+Madi_group+Sta3n,
               family = binomial(link="logit"),
               data = odmdata2lr)

exp(cbind(coef(model2a), confint(model2a)))
summary(model2a)


model3 <- glm(hosp_crit ~ Variant + Vaccination_status + Age_group + Sex + Racial_group+ CCI_group+Madi_group + Sta3n,
              family = binomial(link="logit"),
              data = odmdata2lr)

exp(cbind(coef(model3), confint(model3)))
summary(model3)

model3_unadj <- glm(hosp_crit ~ Variant,
              family = binomial(link="logit"),
              data = odmdata2lr)


exp(cbind(coef(model3_unadj), confint(model3_unadj)))
summary(model3_unadj)



# plot models for table 4

sjPlot::set_theme(theme="scatterw")

p1 <- plot_model(model,
           #group.terms=c(1,2,2,2,3,3,3,3,3,3,4,5,5,6,6,6),
           show.values = F,
           value.offset = .4,
           value.size = 3,
           dot.size = 3,
           line.size = 1.0,
           vline.color="grey",
           prefix.labels = "label")

p11 <- p1+theme_sjplot() + coord_flip()

p2 <- plot_model(model2a,
           
           show.values = F,
           value.offset = .4,
           value.size = 3,
           dot.size = 3,
           line.size = 1.0,
           vline.color = "grey",
           axis.lim = c(0.1, 30))

p22 <- p2+theme_sjplot()+coord_flip()

plot_grid(p11, p22, labels="AUTO")

odmdata2 <- odmdata
odmdata2 = apply_labels(odmdata2,
                        inftype2 = c("Delta variant" = 0,
                                     "Omicron variant" = 1),
                        vac_cat = c("No vaccine" =0,
                                    "One vaccine dose" =1,
                                    "2nd vaccine dose <3 Months"=2,
                                    "2nd vaccine dose >3 Months"=3,
                                    "3rd vaccine dose > 14 days"=4,
                                    "Other" = 99),
                        age_cat2= c("Age 18-39" =0,
                                    "Age 40-49" =1,
                                    "Age 50-59" =2,
                                    "Age 60-69" =3,
                                    "Age 70-79" =4,
                                    "Age 80-89" =5,
                                    "Age 90+" = 6),
                                    
                        sex = c("female"=0,
                                "male"= 1),
                        race_cat2 = c("White"=0,
                                      "Black" =1,
                                      "Other" =2),
                        cci_cat2 = c("No comorbidities" =0,
                                    "1  comorbidity" = 1,
                                    "2  comorbidities" = 2,
                                    "3+ comorbidities" =3))


# supplement tables

dsp14[dsp14$multvis==1,] %>%
  dplyr::select(Age, age_cat, inftype, gender, race_cat, CCI_2yrs, cci_cat2, chf.x, hypunc, copd, diabunc, rf, lymph, canc, metacanc.x, rheumd.x, alcohol, drug, depre, score.x) %>%
  tbl_summary(
    by=inftype,
    sort = all_categorical() ~ "frequency",
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous2() ~1, all_categorical()~1),
    label = list(Age ~ "Age", gender ~ "Sex", race_cat ~ "Race"),
    missing ="no",
    type = list(Age ~ "continuous", gender ~ "categorical", race_cat ~ "categorical")) %>%
  as_flex_table()


# timing of infection related to vaccination
table(dsp14$inftype[dsp14$vaccines==2 & between(dsp14$pos_test_date-dsp14$Series2Date2, -1000, 14)])
table(dsp14$inftype[dsp14$vaccines==3 & (dsp14$pos_test_date > dsp14$Series2Date2+14 & dsp14$pos_test_date <dsp14$Series3Date2+14)])
table(dsp14$inftype[dsp14$vaccines==3 & (dsp14$pos_test_date > dsp14$Series3Date2+14)])



# vaccination status at time of infection
dsp15 <- dsp14[dsp14$multvis==1,]
table(dsp15$vaccines[dsp15$vaccines==0], dsp15$inftype[dsp15$vaccines==0])
table(dsp15$vaccines[dsp15$vaccines==1], dsp15$inftype[dsp15$vaccines==1])
table(dsp15$inftype[dsp15$vaccines==2 & dsp15$timevacc2toinfec <84], dsp15$vaccines[dsp15$vaccines==2 & dsp15$timevacc2toinfec <84])
table(dsp15$inftype[dsp15$vaccines==2 & dsp15$timevacc2toinfec >=84], dsp15$vaccines[dsp15$vaccines==2 & dsp15$timevacc2toinfec >=84])
table(dsp15$inftype[dsp15$vaccines==3 & dsp15$timevacc3toinfec >= 14], dsp15$vaccines[dsp15$vaccines==3 & dsp15$timevacc3toinfec >= 14])


# Infectionstatus related to vaccination


dsp15$infecprevacc2 <- ifelse(dsp15$vaccines %in% c(1,2) & dsp15$pos_test_date < dsp15$Series2Date2+14,1,0)
dsp15$infecprevacc3 <- ifelse(dsp15$vaccines %in% c(2,3) & dsp15$pos_test_date >= dsp15$Series2Date2+14 & dsp15$pos_test_date < dsp15$Series3Date2+14,1,0)
dsp15$infecpostvacc3 <-ifelse(dsp15$vaccines %in% c(3) & dsp15$pos_test_date >= dsp15$Series3Date2+14,1,0)


table(dsp15$inftype[dsp15$vaccines %in% c(1,2)], dsp15$infecprevacc2[dsp15$vaccines %in% c(1,2)])
table(dsp15$inftype[dsp15$vaccines %in% c(2,3)], dsp15$infecprevacc3[dsp15$vaccines %in% c(2,3)])
table(dsp15$inftype[dsp15$vaccines %in% c(3)], dsp15$infecpostvacc3[dsp15$vaccines %in% c(3)])


dsp15$infec <- ifelse(dsp15$pos_test_date<19357 & dsp15$covid_hosp_date >=19357 & dsp15$death_time >dsp15$pos_test_date+28,1,0)
dsp15$hosp <- ifelse(dsp15$covid_hosp_date <19357 & dsp15$covid_icu_date>=19357 & dsp15$death_time >dsp15$pos_test_date+28,1,0)
dsp15$crit <- ifelse(dsp15$pos_test_date < 19357 & (dsp15$covid_icu_date < 19357 | dsp15$death_time <19357 & dsp15$death_time <=dsp15$pos_test_date+28),1,0)

table(dsp15$inftype, dsp15$infec)
table(dsp15$inftype, dsp15$hosp)
table(dsp15$inftype, dsp15$crit)

table(dsp15$infec[dsp15$infecprevacc3==1], dsp15$inftype[dsp15$infecprevacc3==1])
table(dsp15$hosp[dsp15$infecprevacc3==1], dsp15$inftype[dsp15$infecprevacc3==1])
table(dsp15$crit[dsp15$infecprevacc3==1], dsp15$inftype[dsp15$infecprevacc3==1])

table(dsp15$infec[dsp15$infecpostvacc3==1], dsp15$inftype[dsp15$infecpostvacc3==1])
table(dsp15$hosp[dsp15$infecpostvacc3==1], dsp15$inftype[dsp15$infecpostvacc3==1])
table(dsp15$crit[dsp15$infecpostvacc3==1], dsp15$inftype[dsp15$infecpostvacc3==1])

dsp15$hosp_crit <- as.factor(ifelse(dsp15$hosp==1 | dsp15$crit==1, 1,0))

# LR Models for matched cohort stratified by variant predominant time period

library(survival)
library(lme4)
library(stats)
library(MASS)


model4 <- glm(hosp ~  Vaccination_status + Age_group + Sex + Racial_group + CCI_group + Madi_group + Sta3n,
             family = binomial(link="logit"),
             data = odmdata2lr[odmdata2lr$inftype=="delta",])

exp(cbind(coef(model4), confint(model4)))
summary(model4)

model5 <- glm(crit ~  Vaccination_status + Age_group + Sex + Racial_group + CCI_group + Madi_group + Sta3n,
              family = binomial(link="logit"),
              data = odmdata2lr[odmdata2lr$inftype=="delta",])

exp(cbind(coef(model5), confint(model5)))
summary(model5)

model6 <- glm(hosp_crit ~  Vaccination_status + Age_group + Sex + Racial_group + CCI_group + Madi_group + Sta3n,
              family = binomial(link="logit"),
              data = odmdata2lr[odmdata2lr$inftype=="delta",])


exp(cbind(coef(model6), confint(model6)))
summary(model6)

model6_unadj <- glm(hosp_crit ~  Vaccination_status,
              family = binomial(link="logit"),
              data = odmdata2lr[odmdata2lr$inftype=="delta",])


exp(cbind(coef(model6_unadj), confint(model6_unadj)))
summary(model6_unadj)


# supplementary table 4

severity <- read.csv(paste0(machine, waning_ve3, "postindexseverity_jan22.csv"))
colnames(severity)[1] <- "PatientICN"

odmdata_sev <- merge(odmdata, severity, by="PatientICN", all.x=T, all.y=F)
odmdata_sev2 <- filter(odmdata_sev, hosp_crit==1)

chisq.test(table(odmdata_sev2$Vasopressor30d, odmdata_sev2$inftype))
chisq.test(table(odmdata_sev2$Dialysis30d, odmdata_sev2$inftype))
chisq.test(table(odmdata_sev2$Ventilator30d, odmdata_sev2$inftype))
chisq.test(table(odmdata_sev2$OxygenLowFlow30d, odmdata_sev2$inftype))
chisq.test(table(odmdata_sev2$OxygenHighFlow30d, odmdata_sev2$inftype))

# supplementary table 5

#alternative variant time period definition

#delta 1/10/2021 to 12/4/2021 (18901 to 18965)
#omicron 1/2/2022 to 1/18/2022 (18994 to 19010)


#dsp11 <- readRDS(paste0(machine, waning_ve3, "dsp11_jan2022_v2.RDS"))
set.seed(12345)

dsp14_alt <- dsp14
dsp14_alt$inftype3 <- ifelse(dsp14_alt$pos_test_date>=18901 & dsp14_alt$pos_test_date <=18965,"delta",
                         ifelse(dsp14_alt$pos_test_date>=18994 & dsp14_alt$pos_test_date <=19010, "omicron", "other"))
dsp14_alt <- filter(dsp14_alt, inftype3 %in% c("delta", "omicron"))

#create matches for omicron and delta
m.out_alt <- matchit(inftype3 ~ age_cat +  race_cat+ gender+ cci_cat2 + madi_cat+ calweekSeries2Date + Sta3n + vaccines, 
                 data = dsp14_alt, 
                 method = "cem",
                 k2k = TRUE,
                 k2k.method=NULL)

odmdata_alt <- match.data(m.out_alt)

odmdata_alt$infec <- ifelse(odmdata_alt$pos_test_date<19357 & odmdata_alt$covid_hosp_date >=19357 & odmdata_alt$death_time >odmdata_alt$pos_test_date+28,1,0)
odmdata_alt$hosp <- ifelse(odmdata_alt$covid_hosp_date <19357 & odmdata_alt$covid_icu_date>=19357 & odmdata_alt$death_time >odmdata_alt$pos_test_date+28,1,0)
odmdata_alt$crit <- ifelse(odmdata_alt$pos_test_date < 19357 & (odmdata_alt$covid_icu_date < 19357 | odmdata_alt$death_time <19357 & odmdata_alt$death_time <=odmdata_alt$pos_test_date+28),1,0)
odmdata_alt$hosp_crit <- as.factor(ifelse(odmdata_alt$hosp==1 | odmdata_alt$crit==1, 1,0))

table(odmdata_alt$inftype3, odmdata_alt$infec)

odmdata_alt$out_cat <- ifelse(odmdata_alt$hosp==1,2,
                          ifelse(odmdata_alt$crit==1,3,
                                 ifelse(odmdata_alt$infec==1,1,99)))

# aggregate admissions per day per sta3n - redo using dsp14 dataset, not matched odmdata2 dataset

agg_omicron <- do.call(data.frame,
               aggregate(covid_hosp_date ~ Sta3n, dsp14[dsp14$inftype=="omicron"& dsp14$covid_hosp_date<19357,],
                         function(x) c(Total_Admits = length(x), Total_Admits_Sta3n = sum(x), 
                                       Sta3n_prop = sum(x)/length(x))))

agg_delta <- do.call(data.frame,
                       aggregate(covid_hosp_date ~ Sta3n, dsp14[dsp14$inftype=="delta" & dsp14$covid_hosp_date<19357,],
                                 function(x) c(Total_Admits = length(x), Total_Admits_Sta3n = sum(x), 
                                               Sta3n_prop = sum(x)/length(x))))
agg_omicron <- agg_omicron[, c(1:2)]
agg_delta <- agg_delta[, c(1:2)]

colnames(agg_omicron)[2] <- "omicron_admits"
colnames(agg_delta)[2] <- "delta_admits"

agg_admits <- merge(agg_omicron, agg_delta, by="Sta3n")
agg_admits$avg_omicron_per_day <- agg_admits$omicron_admits/24
agg_admits$avg_delta_per_day <- agg_admits$delta_admits/72
agg_admits$ratio <- agg_admits$avg_omicron_per_day/agg_admits$avg_delta_per_day

# acute care beds per Sta3n
library(stringr)
beds <- fread(paste0(machine, waning_ve3, "cum_beds_station.csv"))
facility_list <- as.list(unique(odmdata2$Sta3n))
beds <- beds[beds$Sta3n %in% facility_list]
wards <- fread(paste0(machine, waning_ve3, "wards_vamc.csv"))

#beds <- beds[beds$RoomBedDescription %in% c("BED", "COVID", "CRITICAL", "INTERMEDIATE", "PRIVATE", "ISOLATION", "IM CARE",
#                                            "SURGICAL INTERMEDIATE", "SEMI-PRIVATE", "SEMI-PRIV", "MED/SURG", "STEPDOWN",
#                                            "MED/SURG/NEURO", "ORTHO/INF", "ICU/MICU", "PATIENT CARE", "TEMPORARY MEDICAL BEDS",
#                                            "MEDICINE", "MED/SURG", "CCU", "MICU", "PROGRESSIVE CARE UNIT", "2 BED ROOM",
#                                            "SINGLE ROOM", "ICU", "NEGATIVE PRESSURE", "UNIT", "INTERMEDIATE",
#                                            "SICU", "MEDICINE", "OBSERVATION", "PULM", "ED-BED", "DOUBLE ROOM",
#                                            "INDIVIDUAL", "3 WEST", "WARD", "MEDICAL", "TELEMETRY", "IMCU", "OXYGEN", "BEDROOM",
#                                            "ACUTE MED", "OVERFLOW BED", "GEN CARE", "MEDICINE", "OX", "STEP DOWN", "HOSP BED", "HOSP. BED",
#                                            "HOSPITAL BED", "REGULAR", "BED", "PACU OVERFLOW", "SINGLE ROOM")]

wards <- wards[wards$Specialty %in% c("CARDIAC-STEPDOWN UNIT", "CARDIAC INTENSIVE CARE UNIT", "CARDIAC STEPDOWN UNIT", "CARDIAC SURGERY",
                                      "EAR, NOSE, THROAT (ENT)", "ED-OBSERVATION", "ENDOCRINOLOGY", "GASTROENTEROLOGY", "GEM ACUTE MEDICINE",
                                      "GEM INTERMEDIATE CARE", "GEN MEDICINE (ACUTE)", "GENERAL SURGERY", "GENERAL (ACUTE MEDICINE)", 
                                      "GERONTOLOGY", "HEMATOLOGY/ONCOLOGY", "INTERMEDIATE MEDICINE", "MEDICAL ICU", "MEDICAL OBSERVATION",
                                      "MEDICAL STEPDOWN", "METABOLIC", "NEUROLOGY", "NEUROLOGY OBSERVATION", "NEUROSURGERY", "OPHTHALMOLOGY",
                                      "ORAL SURGERY", "ORTHOPEDIC", "PERIPHERAL VASCULAR", "PLASTIC SURGERY", "PODIATRY", "PROCTOLOGY", 
                                      "PULMONARY, NON-TB", "PULMONARY, TUBERCULOSIS", "RESPITE CARE (MEDICINE)", "SURGICAL ICU", "SURGICAL OBSERVATION",
                                      "SURGICAL STEPDOWN", "TELEMETRY", "THORACIC SURGERY", "TRANSPLANTATION", "UROLOGY", "VASCULAR")]

beds <- beds[beds$WardLocationSID %in% wards$WardLocationSID,]
beds$period <- ifelse(beds$CensusDate <="2021-12-12", "delta", "omicron")
beds <- beds[beds$CensusDate <="2021-12-12" | (beds$CensusDate >="2021-12-25" & beds$CensusDate <="2022-01-17"),]

agg_beds <- beds %>%
  group_by(Sta3n, CensusDate) %>%
  arrange(Sta3n, CensusDate) %>%
  mutate(total_beds_operated=sum(BedsOperating)) %>%
  mutate(total_beds_authorized = sum(BedsAuthorized)) %>%
  mutate(total_transfers_in = sum(CumTransfersIn)) %>%
  mutate(total_transfers_out = sum(CumTransfersOut)) %>%
  ungroup()

agg_beds2 <- agg_beds %>%
  group_by(Sta3n, period) %>%
  slice_min(total_beds_operated) %>%
  distinct(., Sta3n, total_beds_operated, total_beds_authorized, period, .keep_all = F )

saveRDS(agg_beds2, paste0(machine,waning_ve3,"agg_beds.RDS"))
saveRDS(agg_admits, paste0(machine, waning_ve3, "agg_admits.RDS"))

agg_admits2 <- merge(agg_admits, agg_beds2[, c(1,3)], by="Sta3n", all.x=T )
agg_admits2 <- agg_admits2 %>%
  distinct(., Sta3n, omicron_admits, delta_admits, avg_delta_per_day, avg_omicron_per_day, ratio, total_beds_operated, .keep_all = T)
agg_admits2 <- agg_admits2[complete.cases(agg_admits2),]
agg_admits2$total_omicron_to_totalbeds <- agg_admits2$omicron_admits/agg_admits2$total_beds_operated
agg_admits2$total_delta_to_totalbeds <- agg_admits2$delta_admits/agg_admits2$total_beds_operated
saveRDS(agg_admits2, paste0(machine, waning_ve3, "agg_admits2.RDS"))

hosp_low_ratio <- as.list(agg_admits2$Sta3n[agg_admits2$ratio < 3.2500])
hosp_high_ratio <- as.list(agg_admits2$Sta3n[agg_admits2$ratio >= 3.2500])

library(ggplot2)

ggplot(agg_beds2[agg_beds2$total_beds_operated>50,], aes(x=total_beds_operated, y=Sta3n))+
  geom_segment(aes(yend = Sta3n), xend = 0, colour = "grey50") +
  geom_point(size = 3, aes(colour = period)) +
  scale_colour_brewer(palette = "Set1", limits = c("omicron", "delta")) +
  theme_bw() + 
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = c(1,0.55),
    legend.justification = c(1,0.5)
  )


# monoclonal antibody treatment
paxlovid <- fread(paste0(machine, waning_ve3, "paxlovid.csv"))
paxlovid <- paxlovid[, c(2,6:8)]
colnames(paxlovid)[1] <- "PatientICN"
colnames(paxlovid)[2] <- "MAB"
colnames(paxlovid)[3] <- "Date"
colnames(paxlovid)[4] <- "Sta3n"

odmdata2_mrn <- as.list(odmdata2$PatientICN)

paxlovid2 <- paxlovid[paxlovid$PatientICN %in% odmdata2_mrn,]

mabs <- readRDS(paste0(machine, waning_ve3, "covid_mabs.RDS"))
mabs2 <- mabs[mabs$patienticn %in% odmdata2_mrn,]
mabs_mrn <- as.list(unique(mabs$patienticn))
pax_mrn <- as.list(unique(paxlovid$PatientICN))


# re-run matching with excluding patients treated with monoclonal antibodies
`%!in%` <- Negate(`%in%`)
dsp14_2 <- dsp14[dsp14$PatientICN %!in% c(mabs_mrn, pax_mrn),]
mab_treated <- dsp14[dsp14$PatientICN %in% c(mabs_mrn, pax_mrn),]

m.out2 <- matchit(inftype ~ age_cat +  race_cat+ gender+ cci_cat + calweekSeries2Date + Sta3n + vaccines + madi_cat, 
                 data = dsp14_2[dsp14_2$multvis==1,], 
                 method = "cem",
                 k2k = TRUE,
                 k2k.method=NULL)

odmdata_nomab <- match.data(m.out2)
saveRDS(odmdata_nomab, paste0(machine,waning_ve3, "odmdata_revisions_nomab.RDS"))

# read in nomab file from drive
odmdata_nomab <- readRDS(paste0(machine,waning_ve3, "odmdata_revisions_nomab.RDS"))

odmdata_nomab$infec <- ifelse(odmdata_nomab$pos_test_date<19357 & odmdata_nomab$covid_hosp_date >=19357 & odmdata_nomab$death_time >odmdata_nomab$pos_test_date+28,1,0)
odmdata_nomab$hosp <- ifelse(odmdata_nomab$covid_hosp_date <19357 & odmdata_nomab$covid_icu_date>=19357 & odmdata_nomab$death_time >odmdata_nomab$pos_test_date+28,1,0)
odmdata_nomab$crit <- ifelse(odmdata_nomab$pos_test_date < 19357 & (odmdata_nomab$covid_icu_date < 19357 |odmdata_nomab$death_time <19357 & odmdata_nomab$death_time <=odmdata_nomab$pos_test_date+28),1,0)

odmdata_nomab$hosp_crit <- as.factor(ifelse(odmdata_nomab$hosp==1 | odmdata_nomab$crit==1, 1,0))


odmdata_nomab$out_cat <- ifelse(odmdata_nomab$hosp==1,2,
                          ifelse(odmdata_nomab$crit==1,3,
                                 ifelse(odmdata_nomab$infec==1,1,99)))

table(odmdata_nomab$inftype, odmdata_nomab$out_cat)
chisq.test(table(odmdata_nomab$inftype, odmdata_nomab$out_cat))

# read in comorbidity data
charlson <- readRDS(paste0(machine,waning_ve3, "charlson_dsp11.RDS"))
elix <- readRDS(paste0(machine,waning_ve3, "elixh_dsp11.RDS"))
colnames(charlson)[1] <- "PatientICN"
colnames(elix)[1] <- "PatientICN"

odmdata_nomab2 <- merge(odmdata_nomab,elix, by="PatientICN", all.x=T)
odmdata_nomab2 <- merge(odmdata_nomab2,charlson, by="PatientICN", all.x=T)
saveRDS(odmdata_nomab2, paste0(machine,waning_ve3, "odmdata_revisions_cohort_jan22_nomab.RDS"))

mab_treated2 <- merge(mab_treated, elix, by="PatientICN", all.x=T)
mab_treated2 <- merge(mab_treated2, charlson, by="PatientICN", all.x=T)
saveRDS(mab_treated2, paste0(machine, waning_ve3, "mab_treated_cohort_.RDS"))

odmdata_nomab2 <- readRDS(paste0(machine, waning_ve3, "odmdata_revisions_cohort_jan22_nomab.RDS"))

#pic <- readRDS(paste0(machine,waning_ve2,"preindexconditions_jan2022.RDS"))
#odmdata2 <- merge(odmdata, pic[, c(1,18,22,23,30,36,39,44,60,61,95)], by="PatientICN", all.x=T, all.y=F)

# characteristics of patients treated with mabs
library(dplyr)
library(gtsummary)

mab_treated2 %>%
  dplyr::select(Age, age_cat, inftype, gender, race_cat, score.x) %>% #chf.x, hypunc, copd, diab, rf, canc, rheumd.y, score.x) %>%
  tbl_summary(
    by=inftype,
    sort = all_categorical() ~ "frequency",
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous2() ~2,
    label = list(age_cat ~ "Age category", gender ~ "Sex", race_cat ~ "Race"),
    missing ="no",
    type = list(age_cat ~ "categorical", gender ~ "categorical", race_cat ~ "categorical")) %>%
  as_flex_table() %>%
  autofit() %>%
  flextable::save_as_docx(mab_treated2, path = paste0(machine, waning_ve3,"table1_mab.docx"))

length(odmdata2$vaccines[odmdata2$vaccines==2 & odmdata2$timevacc2toinfec >=90 & odmdata2$inftype=="omicron"])

# NOMABS timing of infection related to vaccination
table(odmdata_nomab2$inftype[odmdata_nomab2$vaccines==2 & between(odmdata_nomab2$pos_test_date-odmdata_nomab2$Series2Date2, -1000, 14)])
table(odmdata_nomab2$inftype[odmdata_nomab2$vaccines==3 & (odmdata_nomab2$pos_test_date > odmdata_nomab2$Series2Date2+14 & odmdata_nomab2$pos_test_date <odmdata_nomab2$Series3Date2+14)])
table(odmdata_nomab2$inftype[odmdata_nomab2$vaccines==3 & (odmdata_nomab2$pos_test_date > odmdata_nomab2$Series3Date2+14)])

### NO MABS
# vaccination status at time of infection
table(odmdata_nomab2$vaccines[odmdata_nomab2$vaccines==0], odmdata_nomab2$inftype[odmdata_nomab2$vaccines==0])
table(odmdata_nomab2$vaccines[odmdata_nomab2$vaccines>0], odmdata_nomab2$inftype[odmdata_nomab2$vaccines>0])
table(odmdata_nomab2$vaccines[odmdata_nomab2$vaccines==1], odmdata_nomab2$inftype[odmdata_nomab2$vaccines==1])
table(odmdata_nomab2$inftype[odmdata_nomab2$vaccines==2 & odmdata_nomab2$timevacc2toinfec <84], odmdata_nomab2$vaccines[odmdata_nomab2$vaccines==2 & odmdata_nomab2$timevacc2toinfec <84])
table(odmdata_nomab2$inftype[odmdata_nomab2$vaccines==2 & odmdata_nomab2$timevacc2toinfec >=84], odmdata_nomab2$vaccines[odmdata_nomab2$vaccines==2 & odmdata_nomab2$timevacc2toinfec >=84])
table(odmdata_nomab2$inftype[odmdata_nomab2$vaccines==3 & odmdata_nomab2$timevacc3toinfec >= 14], odmdata_nomab2$vaccines[odmdata_nomab2$vaccines==3 & odmdata_nomab2$timevacc3toinfec >= 14])


# Infectionstatus related to vaccination

odmdata_nomab2$infecprevacc2 <- ifelse(odmdata_nomab2$vaccines %in% c(1,2) & odmdata_nomab2$pos_test_date <= odmdata_nomab2$Series2Date2+14,1,0)
odmdata_nomab2$infecprevacc3 <- ifelse(odmdata_nomab2$vaccines %in% c(2,3) & odmdata_nomab2$pos_test_date > odmdata_nomab2$Series2Date2+14 & odmdata_nomab2$pos_test_date <= odmdata_nomab2$Series3Date2+14,1,0)
odmdata_nomab2$infecpostvacc3 <-ifelse(odmdata_nomab2$vaccines %in% c(3) & odmdata_nomab2$pos_test_date >= odmdata_nomab2$Series3Date2+14,1,0)


table(odmdata_nomab2$inftype[odmdata_nomab2$vaccines %in% c(1,2)], odmdata_nomab2$infecprevacc2[odmdata_nomab2$vaccines %in% c(1,2)])
table(odmdata_nomab2$inftype[odmdata_nomab2$vaccines %in% c(2,3)], odmdata_nomab2$infecprevacc3[odmdata_nomab2$vaccines %in% c(2,3)])
table(odmdata_nomab2$inftype[odmdata_nomab2$vaccines %in% c(3)], odmdata_nomab2$infecpostvacc3[odmdata_nomab2$vaccines %in% c(3)])

# defining disease outcomes
odmdata_nomab2$infec <- ifelse(odmdata_nomab2$pos_test_date<19357 & odmdata_nomab2$covid_hosp_date >=19357 & odmdata_nomab2$death_time >odmdata_nomab2$pos_test_date+28,1,0)
odmdata_nomab2$hosp <- ifelse(odmdata_nomab2$covid_hosp_date <19357 & odmdata_nomab2$covid_icu_date>=19357 & odmdata_nomab2$death_time >odmdata_nomab2$pos_test_date+28,1,0)
odmdata_nomab2$crit <- ifelse(odmdata_nomab2$pos_test_date < 19357 & (odmdata_nomab2$covid_icu_date < 19357 | odmdata_nomab2$death_time <19357 & odmdata_nomab2$death_time <=odmdata_nomab2$pos_test_date+28),1,0)


saveRDS(odmdata_nomab2, paste0(machine, waning_ve3, "odmcohort_jan2022_revisions_nomab.RDS"))

# table 2
table(odmdata_nomab2$inftype, odmdata_nomab2$infec)
table(odmdata_nomab2$inftype, odmdata_nomab2$hosp)
table(odmdata_nomab2$inftype, odmdata_nomab2$crit)
