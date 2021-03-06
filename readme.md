# Relationship between blood markers of inflammation and PD: an observational and MR study
- Authors: Mel Jensen, Ben Jacobs, Ruth Dobson, Sara Bandres-Ciga, Cornelis Blauwendraat, Anette Schrag, The International Parkinsons Disease Genomics Consortium, Alastair J Noyce
- Code: Ben J
- Preprint: [here](https://www.medrxiv.org/content/10.1101/2020.09.13.20189530v3)
- Updated: Ben J, 10/12/20

## Overview
- Convert and extract UKB phenotype data
- Observational study in UKB
- Mendelian randomisation

## Convert and extract UKB phenotype data

First we extract the column names
````unix
# navigate to the working directory
cd /data/Wolfson-UKBB-Dobson/ukb_pheno_0204/

# get headers
head -n1 ukb_pheno.tsv > colnames
````
Then we rename the columns to make life easier
````R
library(readr)
library(dplyr)
setwd("/data/Wolfson-UKBB-Dobson/ukb_pheno_0204/")

bd = read_tsv("colnames")
print("read in pheno data as df")

# get the column names from bd (in format f.etc)
colnames = colnames(bd)
print("got colnames")
# get data dictionary
data_dic = readr::read_csv("/data/Wolfson-UKBB-Dobson/ukb_phenotype_data/Data_Dictionary_Showcase.csv")
print("read in data dic")
# get field names from the column names - this splits the string so you're getting the bit after the first decimal point
# eg f.400 becomes 400
# this for loop will then pump out those field names to a variable called colnames_fields
colnames_fields = c()
for(i in 1:length(colnames)){
  name = stringr::str_split(colnames[i],'\\.')[[1]][2]
  colnames_fields <<- c(colnames_fields,name)
}

# make a new df with those cleaned column names
df = data.frame(colnames_fields)
colnames(df) = 'FieldID'

#change data dictionary field names to a string
data_dic$FieldID = as.character(data_dic$FieldID)
#just choose the interesting columns from the data dictionary (Field ID and Field)
new_data_dic = dplyr::select(data_dic,FieldID,Field)

#combine the phenotype data and the data dictionary - joining matched rows with same field id
combo = dplyr::left_join(df,new_data_dic,by='FieldID')

#now we want to add back in the extra information in the variable names - i.e. visit number and observation number (where one individual has a couple of different data points for a variable)
fields = combo$Field
new_colnames_fields = c()
for(i in 1:length(colnames)){
  instance = stringr::str_split(colnames[i],'\\.')[[1]][3]
  array = stringr::str_split(colnames[i],'\\.')[[1]][4]
  name = paste0(fields[i],'.',instance,'.',array)
  new_colnames_fields <<- c(new_colnames_fields,name)
}

#now we will add in these column names
# first we'll copy the data so we dont fidle with the original
new_bd = bd
colnames(new_bd) = new_colnames_fields
colnames(new_bd)[1] = 'EID'

#export cleaned dataset
write_tsv(new_bd,"colnames.tsv")
````
Then we split the phenotype file into chunks
````unix
# split pheno file into chunks
awk 'NR>1{print}' ukb_pheno.tsv > pheno_file
awk 'NR==1{print}' colnames.tsv > pheno_colnames
split -d -n10 pheno_file

for i in {0..9}
  do
    x=$(echo "x0"$i)
    cat pheno_colnames $x > chunk$x
  done

rm pheno_file
rm colnames
rm x0*
rm colnames.tsv
````
Now we extract ICD codes to exclude individuals with various autoimmune and inflammatory conditions which could introduce confounding
````unix
qsub /data/Wolfson-UKBB-Dobson/MEL_PD/ukb_convert_pheno.sh
````
And finally we combined the results to make one nice workable phenotype file
````unix
cd /data/Wolfson-UKBB-Dobson/ukb_pheno_0204
rm ukb_pheno_final_712
awk 'NR==1{print}' outputx00 >> ukb_pheno_final_712
for i in {1..10}
  do
    j=$(($i-1))
    echo "adding outputx0$j"
    awk 'NR>1{print}' outputx0$j >> ukb_pheno_final_712
  done
````
Combine with main dataset
````R
# load packages

library(dplyr)
library(readr)

# load datasets & merge
old_pheno = read_tsv("/data/Wolfson-UKBB-Dobson/ukb_pheno_911/ukb_pheno_final_PD_1_02")
new_pheno = read_tsv("/data/Wolfson-UKBB-Dobson/ukb_pheno_0204/ukb_pheno_final_712")
old_pheno = old_pheno %>% select(-contains("PD")) %>% filter(EID %in% new_pheno$EID)
selected_vars = old_pheno %>% left_join(new_pheno,by="EID")

# determine age at dx
selected_vars = selected_vars %>%
  mutate("dob"=as.numeric(as.Date(paste0(`Year of birth.0.0`,"-01-01")))) %>%
  mutate("dodx"=as.numeric(`Date of parkinson's disease report.0.0`)) %>%
  mutate("age_at_pd_dx"=(dodx-dob)/365)
print("Sum stats for age at dx")
summary(selected_vars$age_at_pd_dx)

# define incident cases
incident_cases = selected_vars %>% filter(PD_status==1 & age_at_pd_dx>=`Age at recruitment.0.0`) %>% mutate("PD_Dx_after_rec"=1)
non_incident_cases = selected_vars %>% filter(!PD_status==1) %>% filter(!EID %in% incident_cases$EID) %>% mutate("PD_Dx_after_rec"=0)
cases_after_rec = selected_vars %>% filter(PD_status==1) %>% filter(!EID %in% incident_cases$EID) %>% mutate("PD_Dx_after_rec"=NA)
selected_vars = bind_rows(incident_cases,non_incident_cases,cases_after_rec)

print("Printing number of incident cases")
table(selected_vars$PD_Dx_after_rec)
print("Printing number of prevalent cases")
table(selected_vars$PD_status)

#define ethnicity
selected_vars$`Ethnic background.0.0`=factor(recode(selected_vars$`Ethnic background.0.0`,African="Non-white",`Any other Asian background`="Non-white",`Other ethnic group`="Non-white",`White and Black Caribbean`="Non-white",`White and Asian`="Non-white",`Any other Black background`="Non-white",`White and Black African`="Non-white",`Any other mixed background`="Non-white",`Bangladeshi`="Non-white",`British Caribbean`="Non-white",`Pakistani`="Non-white",`Indian`="Non-white",`Caribbean`="Non-white",`Chinese`="Non-white",`British`="White",`Irish`="White",`Any other white background`="White",`Do not know`="NA",`Prefer not to answer`="NA"),levels=c("White","Non-white"))
write_tsv(selected_vars,"/data/Wolfson-UKBB-Dobson/MEL_PD/mel_fbc.tsv")
````

## Observational study
````R
#############################################
#               Load packages
#############################################

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(RNOmni)
library(stringr)


#############################################
#               read in data
#############################################
setwd("/data/Wolfson-UKBB-Dobson/MEL_PD")
df = read_tsv("mel_fbc.tsv")
new_df = read_tsv("../ukb_pheno_2410/ukb_coded_pheno_final.tsv",col_types=cols_only(
  `EID` = col_double(),
  `Vitamin D.0.0`=col_double(),
  `Glycated haemoglobin (HbA1c).0.0`=col_double(),
  `C-reactive protein.0.0`=col_double(),
  `Albumin.0.0`=col_double(),
  `Creatinine.0.0`=col_double()))

df = df %>% left_join(new_df,by="EID")

#############################################
#               make flow chart
#############################################

# whole dataset
included_dataset = df
table(included_dataset$PD_status)

# remove confounders
included_dataset = included_dataset %>% filter(FBC_exclusion=="Include")
table(included_dataset$PD_status)

# remove prevalent cases
included_dataset = included_dataset %>% filter(!is.na(PD_Dx_after_rec))
table(included_dataset$PD_status)

# remove prevalent cases
df = df %>% filter(!is.na(PD_Dx_after_rec))
# remove exclusions, and save excluded people
fbc_exclusions = df %>% filter(FBC_exclusion=="Exclude")
df = df %>% filter(FBC_exclusion=="Include")

# remove people who may have pd dx from other sources (e.g death cert)
df = df %>% filter(!(PD_status==0 & !is.na(age_at_pd_dx)))

# remove people missing count data
df = df %>% filter(!is.na(`Lymphocyte count.0.0`))
table(df$PD_status)


#############################################
#               matching
#############################################

pd = df %>% filter(PD_Dx_after_rec==1) %>% select(EID,`Age at recruitment.0.0`,Sex.0.0)
controls = df %>% filter(PD_Dx_after_rec==0) %>% select(EID,`Age at recruitment.0.0`,Sex.0.0)

set.seed(123456)
for (i in 1:nrow(pd)){
  print(paste0("matching case ",i))
  case = pd[i,]
  matched_controls = controls %>%
    filter(!EID %in% pd$EID) %>%
    filter(`Age at recruitment.0.0`== case$`Age at recruitment.0.0` & Sex.0.0 == case$Sex.0.0) %>% sample_n(size=4,replace=FALSE)
  pd <<- pd %>% bind_rows(matched_controls)
}

matched_df = df %>% filter(EID %in% pd$EID)

#############################################
#               descriptive characteristics of included participants
#############################################

demo_df = data.frame()

demog_cat = function(x){
  tbl = table(df[[x]],df$PD_Dx_after_rec)
  tbl = cbind(tbl,tbl[,1]/colSums(tbl)[1]*100,tbl[,2]/colSums(tbl)[2]*100)
  demo_df <<- rbind(demo_df,tbl)
}

demog_cat("Sex.0.0")
demog_cat("Ethnic background.0.0")
demog_cat("Country of birth (UK/elsewhere).0.0")

tbl = df %>% group_by(PD_Dx_after_rec) %>%
  summarise("Age at PD Dx" = mean(age_at_pd_dx,na.rm=TRUE),
            "SD Age at PD Dx" = sd(age_at_pd_dx,na.rm=TRUE),
            "Age at recruitment" = mean(`Age at recruitment.0.0`,na.rm=TRUE),
            "SD Age at recruitment" = sd(`Age at recruitment.0.0`,na.rm=TRUE),
            "Townsend score" = mean(`Townsend deprivation index at recruitment.0.0`,na.rm=TRUE),
            "SD Townsend score" = sd(`Townsend deprivation index at recruitment.0.0`,na.rm=TRUE))
tbl = data.frame(t(tbl))
tbl$trait = rownames(tbl)
demo_df$trait = rownames(demo_df)
write_csv(tbl,"cont_demo.csv")
write_csv(demo_df,"cat_demo.csv")



# repeat for matched analysis
demo_matched_df = data.frame()

demog_cat = function(x){
  tbl = table(matched_df[[x]],matched_df$PD_Dx_after_rec)
  tbl = cbind(tbl,tbl[,1]/colSums(tbl)[1]*100,tbl[,2]/colSums(tbl)[2]*100)
  demo_matched_df <<- rbind(demo_matched_df,tbl)
}

demog_cat("Sex.0.0")
demog_cat("Ethnic background.0.0")
demog_cat("Country of birth (UK/elsewhere).0.0")

tbl = matched_df %>% group_by(PD_Dx_after_rec) %>%
  summarise("Age at PD Dx" = mean(age_at_pd_dx,na.rm=TRUE),
            "SD Age at PD Dx" = sd(age_at_pd_dx,na.rm=TRUE),
            "Age at recruitment" = mean(`Age at recruitment.0.0`,na.rm=TRUE),
            "SD Age at recruitment" = sd(`Age at recruitment.0.0`,na.rm=TRUE),
            "Townsend score" = mean(`Townsend deprivation index at recruitment.0.0`,na.rm=TRUE),
            "SD Townsend score" = sd(`Townsend deprivation index at recruitment.0.0`,na.rm=TRUE))
tbl = data.frame(t(tbl))
tbl$trait = rownames(tbl)
demo_matched_df$trait = rownames(demo_matched_df)
write_csv(tbl,"matched_cont_demo.csv")
write_csv(demo_matched_df,"matched_cat_demo.csv")

# repeat for excluded people
demo_exclusions = data.frame()

demog_cat = function(x){
  tbl = table(fbc_exclusions[[x]],fbc_exclusions$PD_Dx_after_rec)
  tbl = cbind(tbl,tbl[,1]/colSums(tbl)[1]*100,tbl[,2]/colSums(tbl)[2]*100)
  demo_exclusions <<- rbind(demo_exclusions,tbl)
}

demog_cat("Sex.0.0")
demog_cat("Ethnic background.0.0")
demog_cat("Country of birth (UK/elsewhere).0.0")

tbl = fbc_exclusions %>% group_by(PD_Dx_after_rec) %>%
  summarise("Age at PD Dx" = mean(age_at_pd_dx,na.rm=TRUE),
            "SD Age at PD Dx" = sd(age_at_pd_dx,na.rm=TRUE),
            "Age at recruitment" = mean(`Age at recruitment.0.0`,na.rm=TRUE),
            "SD Age at recruitment" = sd(`Age at recruitment.0.0`,na.rm=TRUE),
            "Townsend score" = mean(`Townsend deprivation index at recruitment.0.0`,na.rm=TRUE),
            "SD Townsend score" = sd(`Townsend deprivation index at recruitment.0.0`,na.rm=TRUE))
tbl = data.frame(t(tbl))
tbl$trait = rownames(tbl)
demo_exclusions$trait = rownames(demo_exclusions)
write_csv(tbl,"excluded_people_cont_demo.csv")
write_csv(demo_exclusions,"excluded_people_cat_demo.csv")


#############################################
#               basic count distributions
#############################################
counts = df %>% select(`C-reactive protein.0.0`,
`Albumin.0.0`,
`Platelet count.0.0`,
`Lymphocyte count.0.0`,
`Monocyte count.0.0`,
`Eosinophill count.0.0`,
`Neutrophill count.0.0`,
`Basophill count.0.0`,
`White blood cell (leukocyte) count.0.0`,
`C-reactive protein.0.0`)

summary_counts = bind_rows(counts %>% summarise_all(median,na.rm=TRUE),
counts %>% summarise_all(mean,na.rm=TRUE),
counts %>% summarise_all(min,na.rm=TRUE),
counts %>% summarise_all(max,na.rm=TRUE),
counts %>% summarise_all(sd,na.rm=TRUE)) %>%
mutate(stat = c("median","mean","min","max","sd"))
write_csv(summary_counts,"summary_counts.csv")

#############################################
#               logit models
#############################################

coef_df = data.frame()
make_model = function(x){
  z_score = (df[[x]]-mean(df[[x]],na.rm=TRUE))/sd(df[[x]],na.rm=TRUE)
  model = glm(data=df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+z_score,family=binomial(link="logit"))
  coef_tbl = summary(model)$coefficients
  aov = anova(model,test="Chisq")
  coef_df <<- rbind(coef_df,c(coef_tbl[6,],aov$`Pr(>Chi)`[6]))
}

make_model("Platelet count.0.0")
make_model("Lymphocyte count.0.0")
make_model("Monocyte count.0.0")
make_model("Eosinophill count.0.0")
make_model("Neutrophill count.0.0")
make_model("Basophill count.0.0")
make_model("White blood cell (leukocyte) count.0.0")
make_model("C-reactive protein.0.0")
make_model("Albumin.0.0")

colnames(coef_df) = c("Beta","SE","Z","P","LR_P")
coef_df$Trait = c("Platelet count","Lymphocyte count","Monocyte count","Eosinophil count","Neutrophil count","Basophil count","Total White Cell Count","CRP","Albumin")

coef_df = coef_df %>% arrange(LR_P) %>% select(Trait,Beta,SE,LR_P) %>% rename("P"="LR_P")

# flip betas so that everything is framed as the effect of a decrease in the trait
coef_df$Beta = coef_df$Beta*(-1)

coef_df$or = exp(coef_df$Beta)
coef_df$lower_ci = exp(coef_df$Beta-1.96*coef_df$SE)
coef_df$upper_ci = exp(coef_df$Beta+1.96*coef_df$SE)
coef_df$Q = p.adjust(coef_df$P,method="fdr")
coef_df = coef_df %>% mutate("OR (95% CI)" = paste0(round(or,2)," (95% CI ",round(lower_ci,2)," - ",round(upper_ci,2),")"))

write_csv(coef_df %>% select(1,9,4,8),"fbc_model_coefficients.csv")

# plot

p1 = ggplot(coef_df,aes(Beta,Trait))+geom_vline(xintercept=0,alpha=0.2)+geom_point()+geom_errorbarh(mapping=aes(xmin=Beta-1.96*SE,xmax=Beta+1.96*SE,y=Trait),height=0.1)+theme_classic()+labs(x="Beta: log odds ratio for incident PD per 1-SD decrease in trait",y="Trait")
png("model_coefficients_fbc.png",width=6,height=8,res=300,units="in")
p1
dev.off()


#############################################
#             logit models matched
#############################################


coef_matched_df = data.frame()
make_model = function(x){
  z_score = (matched_df[[x]]-mean(matched_df[[x]],na.rm=TRUE))/sd(matched_df[[x]],na.rm=TRUE)
  model = glm(data=matched_df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+z_score,family=binomial(link="logit"))
  coef_tbl = summary(model)$coefficients
  aov = anova(model,test="Chisq")
  coef_matched_df <<- rbind(coef_matched_df,c(coef_tbl[6,],aov$`Pr(>Chi)`[6]))
}

make_model("Platelet count.0.0")
make_model("Lymphocyte count.0.0")
make_model("Monocyte count.0.0")
make_model("Eosinophill count.0.0")
make_model("Neutrophill count.0.0")
make_model("Basophill count.0.0")
make_model("White blood cell (leukocyte) count.0.0")
make_model("C-reactive protein.0.0")
make_model("Albumin.0.0")

colnames(coef_matched_df) = c("Beta","SE","Z","P","LR_P")
coef_matched_df$Trait = c("Platelet count","Lymphocyte count","Monocyte count","Eosinophil count","Neutrophil count","Basophil count","Total White Cell Count","CRP","Albumin")

coef_matched_df = coef_matched_df %>% arrange(LR_P) %>% select(Trait,Beta,SE,LR_P) %>% rename("P"="LR_P")

# flip betas so that everything is framed as the effect of a decrease in the trait
coef_matched_df$Beta = coef_matched_df$Beta*(-1)

coef_matched_df$or = exp(coef_matched_df$Beta)
coef_matched_df$lower_ci = exp(coef_matched_df$Beta-1.96*coef_matched_df$SE)
coef_matched_df$upper_ci = exp(coef_matched_df$Beta+1.96*coef_matched_df$SE)
coef_matched_df$Q = p.adjust(coef_matched_df$P,method="fdr")
coef_matched_df = coef_matched_df %>% mutate("OR (95% CI)" = paste0(round(or,2)," (95% CI ",round(lower_ci,2)," - ",round(upper_ci,2),")"))

write_csv(coef_matched_df %>% select(1,9,4,8),"matched_fbc_model_coefficients.csv")

# plot

p1 = ggplot(coef_matched_df,aes(Beta,Trait))+geom_vline(xintercept=0,alpha=0.2)+geom_point()+geom_errorbarh(mapping=aes(xmin=Beta-1.96*SE,xmax=Beta+1.96*SE,y=Trait),height=0.1)+theme_classic()+labs(x="Beta: log odds ratio for incident PD per 1-SD decrease in trait",y="Trait")
png("matched_model_coefficients_fbc.png",width=6,height=8,res=300,units="in")
p1
dev.off()

#############################################
#             sensitivity analyses
#############################################

# recode alcohol and smoking

df$alcohol = recode(df$`Alcohol drinker status.0.0`,"Never"="Never","Previous"="Ever","Current"="Ever","Prefer not to answer"="Prefer not to answer")
df$smoking = recode(df$`Smoking status.0.0`,"Never"="Never","Previous"="Ever","Current"="Ever","Prefer not to answer"="Prefer not to answer")
df = df %>% mutate("alcohol"=na_if(alcohol,"Prefer not to answer")) %>%
  mutate("smoking"=na_if(smoking,"Prefer not to answer"))
df = df %>% filter(!is.na(alcohol))


# 1. control for BMI, smoking, alcohol
coef_df = data.frame()
make_model = function(x){
  z_score = (df[[x]]-mean(df[[x]],na.rm=TRUE))/sd(df[[x]],na.rm=TRUE)
  model = glm(data=df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+ `Body mass index (BMI).0.0`+ alcohol + smoking + z_score,family=binomial(link="logit"))
  coef_tbl = summary(model)$coefficients
  aov = anova(model,test="Chisq")
  coef_df <<- rbind(coef_df,c(coef_tbl[9,],aov$`Pr(>Chi)`[9]))
}


colnames(coef_df) = c("Beta","SE","Z","P","LR_P")
coef_df$Trait = c("Platelet count","Lymphocyte count","Monocyte count","Eosinophil count","Neutrophil count","Basophil count","Total White Cell Count","CRP","Albumin")

coef_df = coef_df %>% arrange(LR_P) %>% select(Trait,Beta,SE,LR_P) %>% rename("P"="LR_P")

# flip betas so that everything is framed as the effect of a decrease in the trait
coef_df$Beta = coef_df$Beta*(-1)

coef_df$or = exp(coef_df$Beta)
coef_df$lower_ci = exp(coef_df$Beta-1.96*coef_df$SE)
coef_df$upper_ci = exp(coef_df$Beta+1.96*coef_df$SE)
coef_df$Q = p.adjust(coef_df$P,method="fdr")
coef_df = coef_df %>% mutate("OR (95% CI)" = paste0(round(or,2)," (95% CI ",round(lower_ci,2)," - ",round(upper_ci,2),")"))

write_csv(coef_df %>% select(1,9,4,8),"extra_covars_fbc_model_coefficients.csv")

# plot

p1 = ggplot(coef_df,aes(Beta,Trait))+geom_vline(xintercept=0,alpha=0.2)+geom_point()+geom_errorbarh(mapping=aes(xmin=Beta-1.96*SE,xmax=Beta+1.96*SE,y=Trait),height=0.1)+theme_classic()+labs(x="Beta: log odds ratio for incident PD per 1-SD decrease in trait",y="Trait")
png("extracovars_model_coefficients_fbc.png",width=6,height=8,res=300,units="in")
p1
dev.off()

# 2. Exclude extreme counts

extreme_counts = df %>% filter(!is.na(`Lymphocyte count.0.0`)) %>%
mutate(lymph_z = (`Lymphocyte count.0.0`-mean(`Lymphocyte count.0.0`))/sd(`Lymphocyte count.0.0`)) %>% filter(abs(lymph_z)<3)
hist(extreme_count$lymph_z)

# recalculate z score with extre values cut out
extreme_counts = extreme_counts %>% filter(!is.na(`Lymphocyte count.0.0`)) %>%
mutate(lymph_z = (`Lymphocyte count.0.0`-mean(`Lymphocyte count.0.0`))/sd(`Lymphocyte count.0.0`))

model = glm(data=extreme_counts,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+lymph_z,family=binomial(link="logit"))
coef_tbl = summary(model)$coefficients
aov = anova(model,test="Chisq")
stats = c(coef_tbl[6,],aov$`Pr(>Chi)`[6])
or = exp(-stats[1])
upper_ci = exp(-stats[1]+1.96*stats[2])
lower_ci = exp(-stats[1]-1.96*stats[2])
paste0(round(or,2)," (95% CI ",round(lower_ci,2)," - ",round(upper_ci,2),")",", p=",round(stats[5],3))



# 3. Exclude people close to diagnosis

df = df %>% mutate(time_to_dx = age_at_pd_dx - `Age at recruitment.0.0`)
pd = df %>% filter(PD_Dx_after_rec==1)
pd$time_to_dx_z = rankNorm(pd$time_to_dx)


# repeat lympho models, excluding people within x years of diagnosis (up to 8)

coef_df = data.frame()
make_model = function(x){
  print(paste0("making models excluding people ",x," years from dx"))
  test_df = df %>% filter(!(PD_Dx_after_rec==1 & time_to_dx <= x))
  z_score = (test_df$`Lymphocyte count.0.0`-mean(test_df$`Lymphocyte count.0.0`,na.rm=TRUE))/sd(test_df$`Lymphocyte count.0.0`,na.rm=TRUE)
    model = glm(data=test_df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+ z_score,family=binomial(link="logit"))
  coef_tbl = summary(model)$coefficients
  aov = anova(model,test="Chisq")
  coef_df <<- rbind(coef_df,c(coef_tbl[6,],aov$`Pr(>Chi)`[6], table(test_df$PD_status)[1],table(test_df$PD_status)[2]))
}
sapply(c(1:8),make_model)

coef_df = cbind(c(1:8),coef_df)
colnames(coef_df) = c("Years before diagnosis filter","beta","se","z","p","lr_p","n_control","n_case")
coef_df = coef_df %>% select(1,2,3,6,7,8) %>% mutate(OR = exp(-beta)) %>% mutate(Lower_CI = exp(-beta-1.96*se)) %>% mutate(Upper_CI = exp(-beta+1.96*se)) %>% select(1,OR,Lower_CI,Upper_CI,n_case,n_control,lr_p)
write_csv(coef_df,"years_before_dx_filter.csv")


# 4 . Lymphopaenia as binary
df = df %>% mutate("lymphopaenia"=ifelse(`Lymphocyte count.0.0`>=1,0,1))
model = glm(data=df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+lymphopaenia,family=binomial(link="logit"))
coef_tbl = summary(model)$coefficients
aov = anova(model,test="Chisq")

stats = c(coef_tbl[6,],aov$`Pr(>Chi)`[6])
or = exp(stats[1])
upper_ci = exp(stats[1]+1.96*stats[2])
lower_ci = exp(stats[1]-1.96*stats[2])
paste0(round(or,2)," (95% CI ",round(lower_ci,2)," - ",round(upper_ci,2),")",", p=",round(stats[5],3))


# 5. whole cohort

df = read_tsv("mel_fbc.tsv")
new_df = read_tsv("../ukb_pheno_2410/ukb_coded_pheno_final.tsv",col_types=cols_only(
  `EID` = col_double(),
  `Vitamin D.0.0`=col_double(),
  `Glycated haemoglobin (HbA1c).0.0`=col_double(),
  `C-reactive protein.0.0`=col_double(),
  `Albumin.0.0`=col_double(),
  `Creatinine.0.0`=col_double()))

df = df %>% left_join(new_df,by="EID")

coef_df = data.frame()
make_model = function(x){
  z_score = (df[[x]]-mean(df[[x]],na.rm=TRUE))/sd(df[[x]],na.rm=TRUE)
  model = glm(data=df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+z_score,family=binomial(link="logit"))
  coef_tbl = summary(model)$coefficients
  aov = anova(model,test="Chisq")
  coef_df <<- rbind(coef_df,c(coef_tbl[6,],aov$`Pr(>Chi)`[6]))
}

make_model("Platelet count.0.0")
make_model("Lymphocyte count.0.0")
make_model("Monocyte count.0.0")
make_model("Eosinophill count.0.0")
make_model("Neutrophill count.0.0")
make_model("Basophill count.0.0")
make_model("White blood cell (leukocyte) count.0.0")
make_model("C-reactive protein.0.0")
make_model("Albumin.0.0")

colnames(coef_df) = c("Beta","SE","Z","P","LR_P")
coef_df$Trait = c("Platelet count","Lymphocyte count","Monocyte count","Eosinophil count","Neutrophil count","Basophil count","Total White Cell Count","CRP","Albumin")

coef_df = coef_df %>% arrange(LR_P) %>% select(Trait,Beta,SE,LR_P) %>% rename("P"="LR_P")

# flip betas so that everything is framed as the effect of a decrease in the trait
coef_df$Beta = coef_df$Beta*(-1)

coef_df$or = exp(coef_df$Beta)
coef_df$lower_ci = exp(coef_df$Beta-1.96*coef_df$SE)
coef_df$upper_ci = exp(coef_df$Beta+1.96*coef_df$SE)
coef_df$Q = p.adjust(coef_df$P,method="fdr")
coef_df = coef_df %>% mutate("OR (95% CI)" = paste0(round(or,2)," (95% CI ",round(lower_ci,2)," - ",round(upper_ci,2),")"))

write_csv(coef_df %>% select(1,9,4,8),"wholecohort_fbc_model_coefficients.csv")

# plot

p1 = ggplot(coef_df,aes(Beta,Trait))+geom_vline(xintercept=0,alpha=0.2)+geom_point()+geom_errorbarh(mapping=aes(xmin=Beta-1.96*SE,xmax=Beta+1.96*SE,y=Trait),height=0.1)+theme_classic()+labs(x="Beta: log odds ratio for incident PD per 1-SD decrease in trait",y="Trait")
png("whle_cohort_model_coefficients_fbc.png",width=6,height=8,res=300,units="in")
p1
dev.off()




#############################################
#             relation to age at dx
#############################################

ggplot(pd,aes(`Lymphocyte count.0.0`,time_to_dx_z))+geom_point()

coef_df = data.frame()
make_model = function(x){
  z_score = (pd[[x]]-mean(pd[[x]],na.rm=TRUE)) /sd(pd[[x]],na.rm=TRUE)
  model = lm(data=pd,time_to_dx_z~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+z_score)
  null_model = lm(data=pd %>% filter(!is.na(pd[[x]])),time_to_dx_z ~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`)
  coef_tbl = summary(model)$coefficients
  aov = anova(model,null_model)
  coef_df <<- rbind(coef_df,c(coef_tbl[6,],aov$`Pr(>F)`[2]))
}

make_model("Platelet count.0.0")
make_model("Lymphocyte count.0.0")
make_model("Monocyte count.0.0")
make_model("Eosinophill count.0.0")
make_model("Neutrophill count.0.0")
make_model("Basophill count.0.0")
make_model("White blood cell (leukocyte) count.0.0")
make_model("C-reactive protein.0.0")
make_model("Albumin.0.0")

colnames(coef_df) = c("Beta","SE","Z","P","LR_P")
coef_df$Trait = c("Platelet count","Lymphocyte count","Monocyte count","Eosinophil count","Neutrophil count","Basophil count","Total White Cell Count","CRP","Albumin")

coef_df = coef_df %>% arrange(LR_P) %>% select(Trait,Beta,SE,LR_P) %>% rename("P"="LR_P")

# flip betas so that everything is framed as the effect of a decrease in the trait
coef_df$Beta = coef_df$Beta*(-1)

coef_df$or = exp(coef_df$Beta)
coef_df$lower_ci = exp(coef_df$Beta-1.96*coef_df$SE)
coef_df$upper_ci = exp(coef_df$Beta+1.96*coef_df$SE)
coef_df$Q = p.adjust(coef_df$P,method="fdr")

write_csv(coef_df %>% select(1,2,3,4,8),"time_to_dx_t_fbc_model_coefficients.csv")
````

## Mendelian randomisation

- Download sum stats for WBC traits and cut a few columns to make them easier to work with
- Filter for P<5e-8 and INFO > 0.3

````unix
cd /data/Wolfson-UKBB-Dobson/MEL_PD
wget ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/UKBB_blood_cell_traits/lymph.assoc
awk '{FS=",";$1=$1}1' OFS="\t" lymph.assoc | awk '{if($10<5e-8 && $18>0.3) print $0}' OFS="\t" > lymph_gwas.tsv
head lymph.assoc -n1 > colnames_for_gwas

wget http://www.dropbox.com/s/5la7y38od95swcf/rf.rdata?dl=0
````

````R

#############################################
#               Load packages
#############################################

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(RNOmni)
library(stringr)
library(TwoSampleMR)

#############################################
#             MR - lymphocytes
#############################################
setwd("/data/Wolfson-UKBB-Dobson/MEL_PD/")

# read in sumstats and change col names
colnames_for_gwas = read_csv("colnames_for_gwas")
lymph = read_tsv("lymph_gwas.tsv",col_names=FALSE)
colnames(lymph) = colnames(colnames_for_gwas)

lymph = lymph %>% rename(
  SNP = ID,
  pval = P,
  effect_allele = ALT,
  other_allele = REF,
  eaf = ALT_FREQ,
  se = SE,
  beta = EFFECT)

# basic QC
lymph = lymph %>% filter(pval < 5e-08)
lymph = lymph %>% filter(!nchar(other_allele)>1)
lymph = lymph %>% filter(!nchar(effect_allele)>1)
lymph = lymph %>% filter(!nchar(other_allele)>1)
lymph = lymph %>% filter(INFO>0.3)
lymph = lymph %>% distinct(SNP,.keep_all = TRUE)
lymph = lymph %>% filter(!(CHR==6 & BP >25000000 & BP < 35000000))

# read in and merge pd meta5 sumstats
ukb_rsids = lymph %>% select(CHR,BP,SNP) %>% mutate("MarkerName"=paste0("chr",CHR,":",BP))

pd = read_table2("/data/Wolfson-UKBB-Dobson/PD/ALASTAIR_SEPT2019/META_NO_UKBB1.tbl")
pd = pd %>% filter(MarkerName %in% ukb_rsids$MarkerName) %>%
  left_join(ukb_rsids,by="MarkerName")

pd = pd %>% select(CHR,BP,SNP,Allele1,Allele2,Freq1,Effect,StdErr,`P-value`) %>%
  rename("effect_allele"=Allele1,
         "other_allele"=Allele2,
         "eaf"=Freq1,
         "beta"=Effect,
         "se"=StdErr,
         "pval"=`P-value`) %>%
  mutate("effect_allele"=toupper(effect_allele),
         "other_allele"=toupper(other_allele))
pd = pd %>% filter(!effect_allele=="D" | effect_allele=="I")
pd = pd  %>% distinct(SNP,.keep_all = TRUE)


pd = pd %>% filter(SNP %in% lymph$SNP)

# remove SNPs not in meta5
lymph = lymph %>% filter(SNP %in% pd$SNP)

# flip betas so effects are in direction of reducing lympho count
lymph = lymph %>% mutate(beta = beta * (-1) )

#clump & format
lymph_dat = clump_data(lymph)
lymph_dat = format_data(lymph_dat,type="exposure")
pd_dat = pd %>% filter(SNP %in% lymph_dat$SNP)
pd_dat = format_data(pd_dat,type="outcome")

#harmonise
combo_dat = harmonise_data(lymph_dat,pd_dat)

# steiger and moe
combo_dat$samplesize.outcome=37688+981372
combo_dat$samplesize.exposure=408112

# Load the downloaded RData object. This loads the rf object
load("rf.rdata?dl=0")

# Obtain estimates from all methods, and generate data metrics
combo_dat$units.outcome="log odds"
combo_dat$units.exposure="z score"
combo_dat$prevalence.outcome = 37688/(37688+981372)
combo_dat$ncase.outcome=37688
combo_dat$ncontrol.outcome=981372

combo_dat$id.exposure="Lymphocyte count"
combo_dat$id.outcome="PD"
combo_dat$r.outcome = get_r_from_lor(combo_dat$beta.outcome,combo_dat$eaf.outcome,combo_dat$ncase.outcome,combo_dat$ncontrol.outcome,combo_dat$prevalence.outcome)
combo_dat$r.exposure = get_r_from_pn(combo_dat$pval.exposure,combo_dat$samplesize.exposure)
combo_dat = combo_dat %>% filter(mr_keep=="TRUE")

steiger = steiger_filtering(combo_dat) %>% filter(steiger_dir == "TRUE") %>% select(SNP)
combo_dat = combo_dat %>% filter(SNP %in% steiger$SNP)


print(paste0("variance explained by all SNPs in instrument is ",sum(combo_dat$r.exposure^2)))

# f stat
r2 = sum(combo_dat$r.exposure^2)
k = nrow(combo_dat)
n = 408112
f_stat = ((n-k-1)/k) * (r2/(1-r2))
print(paste0("F stat ",f_stat))

write_csv(combo_dat,"combo_snps.csv")
res = mr_wrapper(combo_dat)

# standard mr
standard_res = mr(combo_dat)
or = exp(standard_res$b)[3]
lower_ci=exp(standard_res$b[3]-1.96*standard_res$se[3])
upper_ci=exp(standard_res$b[3]+1.96*standard_res$se[3])
paste0("OR ",round(or,2)," 95% CI (",round(lower_ci,2)," - ",round(upper_ci,2),")")
mr_pleiotropy_test(combo_dat)
mr_heterogeneity(combo_dat)

forest = mr_forest_plot(mr_singlesnp(combo_dat))
forest = forest$`Lymphocyte count.PD`+labs(x="Log(OR) for PD per 1-SD increase in lymphocyte count")+theme_classic()+theme(legend.position="none")
png("forest.png",res=300,width=6,height=8, units="in")
forest
dev.off()

scat=mr_scatter_plot(standard_res,combo_dat)
scat = scat$`Lymphocyte count.PD`+labs(x="Per-allele reduction in lymphocyte count (SD)",y="Per-allele effect on PD risk (logOR)")+theme_classic()+theme(legend.position="right")
png("scat.png",res=300,width=6,height=8, units="in")
scat
dev.off()
combo_dat

# MR-MoE - predict the performance of each method
res_moe = mr_moe(res, rf)


summ = res_moe$`Lymphocyte count.PD`$estimates %>% arrange(desc(MOE))

write_csv(summ,"moe.csv")
summ = summ %>% filter(steiger_filtered=="TRUE") %>% filter(!selection=="tophits")
summ$method2 = factor(summ$method2,levels=(summ$method2),order=TRUE)
p1=ggplot(summ,aes(b,method2,fill=MOE))+
  geom_errorbarh(mapping=(aes(xmin=b-1.96*se,xmax=b+1.96*se,y=method2)),height=0.3,col="black")+
  geom_point(size=4,shape=22,col="black")+
  scale_fill_gradient(low="yellow",high="red")+
  geom_vline(xintercept=0,alpha=0.3)+
  theme_classic()+
  labs(x="Beta (MR estimate)",y="MR method")
png("moe.png",res=300,width=6,height=8, units="in")
p1
dev.off()
sum(combo_dat$r.exposure^2)

# mr presso

set.seed(1)
presso_res = run_mr_presso(combo_dat,NbDistribution=1000)

presso_res[[1]][1]$`Main MR results` %>%
as.tbl() %>%
mutate(OR = exp(`Causal Estimate`),
SE = Sd,
lower_ci = exp(`Causal Estimate` - 1.96*SE),
upper_ci = exp(`Causal Estimate` + 1.96*SE)) %>%
select(-2,-3,-4,-5,-8)
````
