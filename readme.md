

````unix
# navigate to the working directory
cd /data/Wolfson-UKBB-Dobson/ukb_pheno_0204/

# get headers
head -n1 ukb_pheno.tsv > colnames

module load R
R
````

````R
library(readr)
library(dplyr)

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

#make a new df with those cleaned column names
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

qsub /data/Wolfson-UKBB-Dobson/MEL_PD/ukb_convert_pheno.sh
````

##### Combine
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

#### Combine with main dataset
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

# Download sum stats for WBC traits and cut a few columns to make them easier to work with

````unix
cd /data/Wolfson-UKBB-Dobson/MEL_PD
wget ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/UKBB_blood_cell_traits/lymph.assoc
awk '{FS=",";$1=$1}1' OFS="\t" lymph.assoc | awk '{print $2,$7,$8,$9,$10,$12,$13}' OFS="\t" > lymph_gwas.tsv
wget ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/UKBB_blood_cell_traits/eo.assoc
awk '{FS=",";$1=$1}1' OFS="\t" eo.assoc | awk '{print $2,$7,$8,$9,$10,$12,$13}' OFS="\t" > eo_gwas.tsv
wget ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/UKBB_blood_cell_traits/mono.assoc
awk '{FS=",";$1=$1}1' OFS="\t" mono.assoc | awk '{print $2,$7,$8,$9,$10,$12,$13}' OFS="\t" > mono_gwas.tsv
````
