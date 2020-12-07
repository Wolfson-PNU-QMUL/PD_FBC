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


# remove prevalent cases
table(df$PD_status)
df = df %>% filter(!is.na(PD_Dx_after_rec))

#############################################
#               matched
#############################################

pd = df %>% filter(PD_Dx_after_rec==1)
controls = df %>% filter(PD_Dx_after_rec==0)

for (i in 1:nrow(pd)){
  print(paste0("matching case ",i))
  case = pd[i,]
  matched_controls = controls %>%
    filter(!EID %in% pd$EID) %>%
    filter(`Age at recruitment.0.0`== case$`Age at recruitment.0.0` & Sex.0.0 == case$Sex.0.0) %>% sample_n(size=4,replace=FALSE)
  pd <<- pd %>% bind_rows(matched_controls)
}

matched_df = pd

#############################################
#               descriptive
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



#############################################
#               logit models
#############################################

coef_df = data.frame()
make_model = function(x){
  model = glm(data=df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+df[[x]],family=binomial(link="logit"))
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
make_model("Lymphocyte percentage.0.0")
make_model("Monocyte percentage.0.0")
make_model("Eosinophill percentage.0.0")
make_model("Neutrophill percentage.0.0")
make_model("Basophill percentage.0.0")

colnames(coef_df) = c("Beta","SE","Z","P","LR_P")
coef_df$Trait = c("Platelet count","Lymphocyte count","Monocyte count","Eosinophil count","Neutrophil count","Basophil count","Total White Cell Count","CRP","Albumin","Lymphocyte percentage","Monocyte percentage","Eosinophil percentage","Neutrophil percentage","Basophil percentage")

coef_df = coef_df %>% arrange(LR_P) %>% select(Trait,Beta,SE,LR_P) %>% rename("P"="LR_P")
write_csv(coef_df,"supplemental_fbc_model_coefficients.csv")
coef_df = coef_df[-grep("percentage",coef_df$Trait),]
coef_df$Q = p.adjust(coef_df$P,method="fdr")

coef_df$or = exp(coef_df$Beta)
coef_df$lower_ci = exp(coef_df$Beta-1.96*coef_df$SE)
coef_df$upper_ci = exp(coef_df$Beta+1.96*coef_df$SE)

write_csv(coef_df,"fbc_model_coefficients.csv")

# plot

p1 = ggplot(coef_df,aes(Beta,Trait))+geom_vline(xintercept=0,alpha=0.2)+geom_point()+geom_errorbarh(mapping=aes(xmin=Beta-1.96*SE,xmax=Beta+1.96*SE,y=Trait),height=0.1)+theme_classic()+labs(x="Beta: log odds of incident PD per unit increase in trait",y="Trait")
png("model_coefficients_raw_fbc.png",width=8,height=8,res=300,units="in")
p1
dev.off()

# repeat with normalised vars

coef_df = data.frame()
normalise = function(x){
  mean = mean(x,na.rm=TRUE)
  sd = sd(x,na.rm=TRUE)
  y=(x - mean)/sd
  return(y)
}


make_model_norm = function(x){
  df$y = normalise(df[[x]])
  model = glm(data=df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+y,family=binomial(link="logit"))
  coef_tbl = summary(model)$coefficients
  aov = anova(model,test="Chisq")
  coef_df <<- rbind(coef_df,c(coef_tbl[6,],aov$`Pr(>Chi)`[6]))
}

make_model_norm("Platelet count.0.0")
make_model_norm("Lymphocyte count.0.0")
make_model_norm("Monocyte count.0.0")
make_model_norm("Eosinophill count.0.0")
make_model_norm("Neutrophill count.0.0")
make_model_norm("Basophill count.0.0")
make_model_norm("White blood cell (leukocyte) count.0.0")
make_model_norm("C-reactive protein.0.0")
make_model_norm("Albumin.0.0")
make_model_norm("Lymphocyte percentage.0.0")
make_model_norm("Monocyte percentage.0.0")
make_model_norm("Eosinophill percentage.0.0")
make_model_norm("Neutrophill percentage.0.0")
make_model_norm("Basophill percentage.0.0")

colnames(coef_df) = c("Beta","SE","Z","P","LR_P")
coef_df$Trait = c("Platelet count","Lymphocyte count","Monocyte count","Eosinophil count","Neutrophil count","Basophil count","Total White Cell Count","CRP","Albumin","Lymphocyte percentage","Monocyte percentage","Eosinophil percentage","Neutrophil percentage","Basophil percentage")

coef_df = coef_df %>% arrange(LR_P) %>% select(Trait,Beta,SE,LR_P) %>% rename("P"="LR_P")
write_csv(coef_df,"normalised_supplemental_fbc_model_coefficients.csv")
coef_df = coef_df[-grep("percentage",coef_df$Trait),]
coef_df$Q = p.adjust(coef_df$P,method="fdr")

coef_df$or = exp(coef_df$Beta)
coef_df$lower_ci = exp(coef_df$Beta-1.96*coef_df$SE)
coef_df$upper_ci = exp(coef_df$Beta+1.96*coef_df$SE)

write_csv(coef_df,"normlised_fbc_model_coefficients.csv")

# plot

p1 = ggplot(coef_df,aes(Beta,Trait))+geom_vline(xintercept=0,alpha=0.2)+geom_point()+geom_errorbarh(mapping=aes(xmin=Beta-1.96*SE,xmax=Beta+1.96*SE,y=Trait),height=0.1)+theme_classic()+labs(x="Beta: log odds of incident PD per 1-SD increase in trait",y="Trait")
png("normalised_model_coefficients_raw_fbc.png",width=8,height=8,res=300,units="in")
p1
dev.off()


#############################################
#             survival analysis
#############################################

library(survival)
library(survminer)
df = df %>% mutate("lymphopaenia"=ifelse(`Lymphocyte count.0.0`>=1,0,1))
model = glm(data=df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+lymphopaenia,family=binomial(link="logit"))
coef_tbl = summary(model)$coefficients
aov = anova(model,test="Chisq")


df = df %>% mutate("time_to_dx" = ifelse(PD_Dx_after_rec==1,age_at_pd_dx - `Age at recruitment.0.0`,12.2))
surv = Surv(time=df$time_to_dx, event=df$PD_Dx_after_rec)
surv_model = survfit(data=df,surv ~ lymphopaenia)
ggsurvplot(surv_model,pval=TRUE)

df = df %>% mutate(lymph_bin = cut2(`Lymphocyte count.0.0`,g=10))
coxph = coxph(surv~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+`Lymphocyte count.0.0`,data=df)
ggforest(coxph)

coxph = coxph(surv~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+lymphopaenia,data=df)
ggforest(coxph)

coxph = coxph(surv~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+lymph_bin,data=df)
ggforest(coxph)

#############################################
#             logit models matched
#############################################

coef_matched_df = data.frame()
make_model = function(x){
  model = glm(data=matched_df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+matched_df[[x]],family=binomial(link="logit"))
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
make_model("Lymphocyte percentage.0.0")
make_model("Monocyte percentage.0.0")
make_model("Eosinophill percentage.0.0")
make_model("Neutrophill percentage.0.0")
make_model("Basophill percentage.0.0")

colnames(coef_matched_df) = c("Beta","SE","Z","P","LR_P")
coef_matched_df$Trait = c("Platelet count","Lymphocyte count","Monocyte count","Eosinophil count","Neutrophil count","Basophil count","Total White Cell Count","CRP","Albumin","Lymphocyte percentage","Monocyte percentage","Eosinophil percentage","Neutrophil percentage","Basophil percentage")

coef_matched_df = coef_matched_df %>% arrange(LR_P) %>% select(Trait,Beta,SE,LR_P) %>% rename("P"="LR_P")
write_csv(coef_matched_df,"matched_supplemental_fbc_model_coefficients.csv")
coef_matched_df = coef_matched_df[-grep("percentage",coef_matched_df$Trait),]
coef_matched_df$Q = p.adjust(coef_matched_df$P,method="fdr")

coef_matched_df$or = exp(coef_matched_df$Beta)
coef_matched_df$lower_ci = exp(coef_matched_df$Beta-1.96*coef_matched_df$SE)
coef_matched_df$upper_ci = exp(coef_matched_df$Beta+1.96*coef_matched_df$SE)

write_csv(coef_matched_df,"matched_fbc_model_coefficients.csv")

# plot

p1 = ggplot(coef_matched_df,aes(Beta,Trait))+geom_vline(xintercept=0,alpha=0.2)+geom_point()+geom_errorbarh(mapping=aes(xmin=Beta-1.96*SE,xmax=Beta+1.96*SE,y=Trait),height=0.1)+theme_classic()+labs(x="Beta: log odds of incident PD per unit increase in trait",y="Trait")
png("matched_model_coefficients_raw_fbc.png",width=8,height=8,res=300,units="in")
p1
dev.off()

# repeat with normalised vars

coef_matched_df = data.frame()
normalise = function(x){
  mean = mean(x,na.rm=TRUE)
  sd = sd(x,na.rm=TRUE)
  y=(x - mean)/sd
  return(y)
}


make_model_norm = function(x){
  matched_df$y = normalise(matched_df[[x]])
  model = glm(data=matched_df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+y,family=binomial(link="logit"))
  coef_tbl = summary(model)$coefficients
  aov = anova(model,test="Chisq")
  coef_matched_df <<- rbind(coef_matched_df,c(coef_tbl[6,],aov$`Pr(>Chi)`[6]))
}

make_model_norm("Platelet count.0.0")
make_model_norm("Lymphocyte count.0.0")
make_model_norm("Monocyte count.0.0")
make_model_norm("Eosinophill count.0.0")
make_model_norm("Neutrophill count.0.0")
make_model_norm("Basophill count.0.0")
make_model_norm("White blood cell (leukocyte) count.0.0")
make_model_norm("C-reactive protein.0.0")
make_model_norm("Albumin.0.0")
make_model_norm("Lymphocyte percentage.0.0")
make_model_norm("Monocyte percentage.0.0")
make_model_norm("Eosinophill percentage.0.0")
make_model_norm("Neutrophill percentage.0.0")
make_model_norm("Basophill percentage.0.0")

colnames(coef_matched_df) = c("Beta","SE","Z","P","LR_P")
coef_matched_df$Trait = c("Platelet count","Lymphocyte count","Monocyte count","Eosinophil count","Neutrophil count","Basophil count","Total White Cell Count","CRP","Albumin","Lymphocyte percentage","Monocyte percentage","Eosinophil percentage","Neutrophil percentage","Basophil percentage")

coef_matched_df = coef_matched_df %>% arrange(LR_P) %>% select(Trait,Beta,SE,LR_P) %>% rename("P"="LR_P")
write_csv(coef_matched_df,"matched_normalised_supplemental_fbc_model_coefficients.csv")
coef_matched_df = coef_matched_df[-grep("percentage",coef_matched_df$Trait),]
coef_matched_df$Q = p.adjust(coef_matched_df$P,method="fdr")

coef_matched_df$or = exp(coef_matched_df$Beta)
coef_matched_df$lower_ci = exp(coef_matched_df$Beta-1.96*coef_matched_df$SE)
coef_matched_df$upper_ci = exp(coef_matched_df$Beta+1.96*coef_matched_df$SE)

write_csv(coef_matched_df,"normlised_fbc_model_coefficients.csv")

# plot

p1 = ggplot(coef_matched_df,aes(Beta,Trait))+geom_vline(xintercept=0,alpha=0.2)+geom_point()+geom_errorbarh(mapping=aes(xmin=Beta-1.96*SE,xmax=Beta+1.96*SE,y=Trait),height=0.1)+theme_classic()+labs(x="Beta: log odds of incident PD per 1-SD increase in trait",y="Trait")
png("normalised_model_coefficients_raw_fbc.png",width=8,height=8,res=300,units="in")
p1
dev.off()




#############################################
#             raw count plots
#############################################

make_plot = function(x){
  p<<-ggplot(df,aes(factor(PD_Dx_after_rec),df[[x]]))+geom_violin(fill="blue",alpha=0.3)+geom_boxplot(fill="red",alpha=0.2,width=0.2)+scale_y_log10()+labs(x="Incident PD (0=control, 1=case)",y="Cell count")
}

make_plot("Eosinophill count.0.0")
png("eos.png",height=8,width=8,res=300,units="in")
p
dev.off()
make_plot("Lymphocyte count.0.0")
png("lymphs.png",height=8,width=8,res=300,units="in")
p
dev.off()

#############################################
#             relation to age at dx
#############################################

df = df %>% mutate(time_to_dx = age_at_pd_dx - `Age at recruitment.0.0`)
pd = df %>% filter(PD_Dx_after_rec==1)
pd$time_to_dx_z = rankNorm(pd$time_to_dx)

ggplot(pd,aes(`Lymphocyte count.0.0`,time_to_dx_z))+geom_point()

coef_df = data.frame()
make_model = function(x){
  model = lm(data=pd,time_to_dx_z~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+pd[[x]])
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
make_model("Lymphocyte percentage.0.0")
make_model("Monocyte percentage.0.0")
make_model("Eosinophill percentage.0.0")
make_model("Neutrophill percentage.0.0")
make_model("Basophill percentage.0.0")

colnames(coef_df) = c("Beta","SE","Z","P","LR_P")
coef_df$Trait = c("Platelet count","Lymphocyte count","Monocyte count","Eosinophil count","Neutrophil count","Basophil count","Total White Cell Count","CRP","Albumin","Lymphocyte percentage","Monocyte percentage","Eosinophil percentage","Neutrophil percentage","Basophil percentage")

coef_df = coef_df %>% arrange(LR_P) %>% select(Trait,Beta,SE,LR_P) %>% rename("P"="LR_P")
write_csv(coef_df,"relation_to_ageatdx_supplemental_fbc_model_coefficients.csv")
coef_df = coef_df[-grep("percentage",coef_df$Trait),]
coef_df$Q = p.adjust(coef_df$P,method="fdr")

coef_df$lower_ci = coef_df$Beta-1.96*coef_df$SE
coef_df$upper_ci = coef_df$Beta+1.96*coef_df$SE
write_csv(coef_df,"relation_to_ageatdx_fbc_model_coefficients.csv")

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
  model = glm(data=df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+ `Body mass index (BMI).0.0`+ alcohol + smoking + df[[x]],family=binomial(link="logit"))
  coef_tbl = summary(model)$coefficients
  aov = anova(model,test="Chisq")
  coef_df <<- rbind(coef_df,c(coef_tbl[9,],aov$`Pr(>Chi)`[9]))
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
coef_df$Q = p.adjust(coef_df$P,method="fdr")

coef_df$or = exp(coef_df$Beta)
coef_df$lower_ci = exp(coef_df$Beta-1.96*coef_df$SE)
coef_df$upper_ci = exp(coef_df$Beta+1.96*coef_df$SE)

write_csv(coef_df,"supplemental_with_extra_covariates_fbc_model_coefficients.csv")

# plot

p1 = ggplot(coef_df,aes(Beta,Trait))+geom_vline(xintercept=0,alpha=0.2)+geom_point()+geom_errorbarh(mapping=aes(xmin=Beta-1.96*SE,xmax=Beta+1.96*SE,y=Trait),height=0.1)+theme_classic()+labs(x="Beta: log odds of incident PD per unit increase in trait",y="Trait")
png("extra_covariates_model_coefficients_raw_fbc.png",width=8,height=8,res=300,units="in")
p1
dev.off()


# 2. Exclude people close to diagnosis

df = df %>% filter(!(PD_Dx_after_rec==1 & time_to_dx <= 5))


coef_df = data.frame()
make_model = function(x){
  model = glm(data=df,PD_Dx_after_rec~`Age at recruitment.0.0`+Sex.0.0+`Townsend deprivation index at recruitment.0.0`+`Ethnic background.0.0`+ df[[x]],family=binomial(link="logit"))
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
coef_df$Q = p.adjust(coef_df$P,method="fdr")

coef_df$or = exp(coef_df$Beta)
coef_df$lower_ci = exp(coef_df$Beta-1.96*coef_df$SE)
coef_df$upper_ci = exp(coef_df$Beta+1.96*coef_df$SE)

write_csv(coef_df,"supplemental_with_time_to_dx_filter_fbc_model_coefficients.csv")

# plot

p1 = ggplot(coef_df,aes(Beta,Trait))+geom_vline(xintercept=0,alpha=0.2)+geom_point()+geom_errorbarh(mapping=aes(xmin=Beta-1.96*SE,xmax=Beta+1.96*SE,y=Trait),height=0.1)+theme_classic()+labs(x="Beta: log odds of incident PD per unit increase in trait",y="Trait")
png("time_to_dx_filter_model_coefficients_raw_fbc.png",width=8,height=8,res=300,units="in")
p1
dev.off()


#############################################
#             MR
#############################################

# read in sumstats and change col names
sumstats = read_csv("lymph.assoc")

sumstats = sumstats %>% rename(
SNP = ID,
effect_allele = ALT,
other_allele = REF,
beta = EFFECT,
pval = P,
se = SE,
eaf = ALT_FREQ)

sumstats = sumstats %>% select(-CHR_num,-GENPOS,-MLOG10P,-ALT_MINOR,-MA_FREQ,-R2,-GWSIG)


# remove SNPs above p threshold
lymph = sumstats %>% filter(pval < 5e-08)
lymph = lymph %>% filter(!nchar(other_allele)>1)
lymph = lymph %>% filter(!nchar(effect_allele)>1)
lymph = lymph %>% filter(!nchar(other_allele)>1)
lymph = lymph %>% filter(INFO>0.3)
rm(sumstats)
#remove duplicates
lymph = lymph %>% distinct(SNP,.keep_all = TRUE)
# exclude mhc
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

###################################
# mr
###################################

pd = pd %>% filter(SNP %in% lymph$SNP)


# remove SNPs not in meta5
lymph = lymph %>% filter(SNP %in% pd$SNP)

#clump & format
lymph_dat = clump_data(lymph)
lymph_gsmr = format_data(lymph,type="exposure")
lymph_dat = format_data(lymph_dat,type="exposure")

pd_gsmr = pd %>% filter(SNP %in% lymph_gsmr$SNP)
pd_gsmr = format_data(pd_gsmr,type="outcome")

pd_dat = pd %>% filter(SNP %in% lymph_dat$SNP)
pd_dat = format_data(pd_dat,type="outcome")

#harmonise
combo_dat = harmonise_data(lymph_dat,pd_dat)
gsmr_combo = harmonise_data(lymph_gsmr,pd_gsmr)

gsmr_combo = clump_data(gsmr_combo)
# steiger and moe
combo_dat$samplesize.outcome=37688+981372
combo_dat$samplesize.exposure=408112

# Load the downloaded RData object. This loads the rf object
load("/data/home/hmy117/rf.rdata")

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
png("forest.png",res=300,width=8,height=8, units="in")
forest
dev.off()

scat=mr_scatter_plot(standard_res,combo_dat)
scat = scat$`Lymphocyte count.PD`+labs(x="Per-allele effect on lymphocyte count (SD)",y="Per-allele effect on PD risk (logOR)")+theme_classic()+theme(legend.position="right")
png("scat.png",res=300,width=8,height=8, units="in")
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
png("moe.png",res=300,width=8,height=8, units="in")
p1
dev.off()
sum(combo_dat$r.exposure^2)
