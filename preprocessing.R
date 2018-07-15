Merged=read.csv('Followup1_March2017.csv',header=T, as.is = c('Housing_1','Merged6_Q2_11a','DaysHomeless_1_Streets','CD4_count_1','Merged6_Q2_10a','NumTimesInject_1','P30Drug_Freq_1'))
colnames(Merged)[1] = 'ParticipantID_1'
var_name = colnames(Merged)
attach(Merged)
missing_del = var_name[colSums(sapply(Merged,is.na))>nrow(Merged)*0.9] # remove variables with more than 90% missing
unrelated_del = c('Month_CD4_1','year_CD4_1','Merged6_CD4_date','Merged6_CD4_date2','Month_ViralLoad_1','Year_ViralLoad_1',
                  'Merged6_VL_date','Merged6_VL_date2','Interviewer_Name_1','agency_name_1','ID1','ID','Mobile_ID_1',
                  'date_1','Merged6_Surv_date2','Minutes_1','Merged6_Q13_6'
)

# Delete the columns with the same response
common_del = c()
for (var in var_name){
  column = Merged[,var]
  freq = table(column)
  if(any(freq>=nrow(Merged)*1))
  common_del = c(common_del,var)
}
Merged_new = Merged[,!var_name%in%c(missing_del,unrelated_del,common_del)]

# Replace the missing value with NAs for numerical variables
var_num = c('Housing_1','DaysHomeless_1_Streets','CD4_count_1','NumTimesInject_1','P30Drug_Freq_1')
var_cat = setdiff(colnames(Merged_new),c(var_num,'ParticipantID_1'))
for(var in var_num){
  for(j in 1:nrow(Merged_new)){
    if(Merged_new[j,var]%in%c("Don't know","dont know","don't know")){
      Merged_new[j,var] = ''
    }
    else if(Merged_new[j,var]%in%c('undetectable')){
      Merged_new[j,var] = 0
    }
  }
}

# Change the type of categorical variables to factor
for(var in var_cat){
 Merged_new[,var] = factor(Merged_new[,var])
}

write.csv(x = Merged_new,file = 'Merged1.csv',row.names = F,na = '')
save(file='Merged1.Rdata',Merged_new,missing_del,unrelated_del,common_del,var_cat,var_num)


#############################################################################################
Merged=read.csv('Followup2_March2017.csv',header=T, as.is = c('Housing_2','Merged6_Q2_22a','DaysHomeless_2_Streets','CD4_count_2','Merged6_Q2_20a','NumTimesInject_2','P30Drug_Freq_2'))
colnames(Merged)[1] = 'ParticipantID_2'
var_name = colnames(Merged)
attach(Merged)
missing_del = var_name[colSums(sapply(Merged,is.na))>nrow(Merged)*0.9] # remove variables with more than 90% missing
unrelated_del = c('Month_CD4_2','year_CD4_2','Merged6_CD4_date','Merged6_CD4_date2','Month_ViralLoad_2','Year_ViralLoad_2',
                  'Merged6_VL_date','Merged6_VL_date2','Interviewer_Name_2','agency_name_2','ID1','ID','Mobile_ID_2',
                  'date_2','Merged6_Surv_date2','Minutes_2','Merged6_Q13_6'
)

# Delete the columns with the same response
common_del = c()
for (var in var_name){
  column = Merged[,var]
  freq = table(column)
  if(any(freq>=nrow(Merged)*1))
    common_del = c(common_del,var)
}

Merged_new = Merged[,!var_name%in%c(missing_del,unrelated_del,common_del)]
# Replace the missing value with NAs for numerical variables
var_num = c('Housing_2','DaysHomeless_2_Streets','CD4_count_2','NumTimesInject_2','P30Drug_Freq_2')
var_cat = setdiff(colnames(Merged_new),c(var_num,'ParticipantID_2'))
for(var in var_num){
  for(j in 1:nrow(Merged_new)){
    if(Merged_new[j,var]%in%c("Don't know","dont know","don't know")){
      Merged_new[j,var] = ''
    }
    else if(Merged_new[j,var]%in%c('undetectable')){
      Merged_new[j,var] = 0
    }
  }
  # Merged_new[,var] = as.numeric(Merged_new[,var])
}

# Change the type of categorical variables to factor
# for(var in var_cat){
#   Merged_new[,var] = factor(Merged_new[,var])  
# }

write.csv(x = Merged_new,file = 'Merged2.csv',row.names = F,na = '')
save(file='Merged2.Rdata',Merged_new,missing_del,unrelated_del,common_del,var_cat,var_num)

#################################################################################################
Merged = read.csv('Baseline_March2017.csv',header=T)
colnames(Merged)[1] = 'Participant_ID'
var_name = colnames(Merged)
missing_del = var_name[colSums(sapply(Merged,is.na))>nrow(Merged)*0.9] # remove variables with more than 90% missing
unrelated_del = c('Specify_healthcare_0', 'Specify_Other_HIVtest', 'Specify_Other_HIVlocation', 'HIV_result_0', 
                  'Month_CD4_0', 'year_CD4_0', 'CD4_Date', 'CD4_Date2', 'Month_ViralLoad_0', 'Year_ViralLoad_0', 'VL_Date',
                  'VL_Date2','Specify_Other_ARTLocation_0', 'PYOther_Drug_specify', 'P30Other_Drug_specify_0', 'Other_PhysAbu_adolescence',
                  'other_PhysAbu_pastYear', 'other_SexAbu_Child', 'Other_SexAbu_ado', 'Other_Sex_abu_adult', 'PYSexabuse', 
                  'Prevent_healthCare_specify', 'Specify_PYhealthcare_MostOften', 'Prevent_Keep_Appt_0_specify', 'Interviewer_Name_0',
                  'agency_name_0', 'ID', 'Mobile_ID_0', 'date_0', 'Surv_Date2', 'Surv_Date2', 'Minutes_0', 'Q13_6',
                  'FU_CD4_Date','FU_CD4_Date2','FU3_Q2_11a','FU3_VL_date','FU3_VL_date2','FU3_surv_date2',
                  'ParticipantID_2','FU6_VL_date','FU6_VL_date2','FU6_CD4_Date','FU6_CD4_Date2','FU6_surv_date2','FU6_Q13_6',
                  'Q1_12_1','Q2_1_1','Q2_1_2','Q2_3_1old','Q2_3_1','Q2_4_1','Q2_5_1','FU3_Surv_date2','FU6_CD4_Date','FU6_CD4_Date2',
                  'FU6_Surv_date2'
)

# Delete the columns with the same response
common_del = c()
for (var in var_name){
  column = Merged[,var]
  freq = table(column)
  if(any(freq>=nrow(Merged)*1))
    common_del = c(common_del,var)
}

BL_new = Merged[,!var_name%in%c(missing_del,unrelated_del,common_del)]
# Replace the missing value with NAs for numerical variables
var_num = c('Participant_ID', 'age_0', 'Child_und18', 'DaysHomeless_0_Streets', 'DaysHomeless_0_Shelter', 'Q2_10a1', 
          'Num_Male_Partners_0', 'Num_female_partners_0', 'Num_trans_partners_0', 'Years_Rel_0', 'Months_Rel_0',
            'Num_SexPartners', 'yearsInjecting', 'NumTimesInject_0', 'P90_jail_0', 'NumMissed_0',
          'FU6_Q2_20a','FU6_Q2_22a'
)
var_cat = setdiff(colnames(BL_new),var_num)
for(var in var_num){
  for(j in 1:nrow(BL_new)){
    if(BL_new[j,var]%in%c("Don't know","dont know","don't know", "dont", "no", "unknown")){
      BL_new[j,var] = ''
    }
    else if(BL_new[j,var]%in%c('undetectable', 'undectable')){
      BL_new[j,var] = 0
    }
  }
  #BL_new[,var] = as.numeric(BL_new[,var])
}

# Change the type of categorical variables to factor
for(var in var_cat){
  BL_new[,var] = factor(BL_new[,var])
}
BL_new$NumTimesInject_2[BL_new$NumTimesInject_2==999]=NA

write.csv(x = BL_new,file = 'Merged.csv',row.names = F,na = '')
save(file='Merged.Rdata',BL_new,missing_del,unrelated_del,common_del,var_cat,var_num)

# Viral Load: ViralLoad_0,Viral_loadLVL_1,Viral_loadLVL_1

# Substance misuse variable
Merged$P30Substanceuse_0=as.numeric((Merged$P30Drug_Freq_0>0)|(Merged$P30_binge_0>0))
Merged$P30Substanceuse_1=as.numeric((Merged$P30Drug_Freq_1>0)|(Merged$P30_binge_1>0))
Merged$P30Substanceuse_2=as.numeric((Merged$P30Drug_Freq_2>0)|(Merged$P30_binge_2>0))

# Health issues related to HIV 
health_0 = grep('P90_healthissues_0_+',colnames(Merged),perl=T,value=T)
health_0 = setdiff(health_0,'P90_healthissues_0_none')
Merged$P90_healthissues_0 = rowSums(Merged[,health_0]>0)
health_1 = grep('P90_healthissues_1_+',colnames(Merged),perl=T,value=T)
health_1 = setdiff(health_1,'P90_healthissues_1_none')
Merged$P90_healthissues_1 = rowSums(Merged[,health_1]>0)
health_2 = grep('P90_healthissues_2_+',colnames(Merged),perl=T,value=T)
health_2 = setdiff(health_2,'P90_healthissues_2_none')
Merged$P90_healthissues_2 = rowSums(Merged[,health_2]>0)

# HIVStigma
Merged$HIVStigma_0 = relationshps0_0 + Relationships1_0 + Relationships2_0 + Relationships3_0 +
  Relationships4_0 + Relationships5_0 + Relationships6_0 + Relationships7_0 + Relationships8_0 +
  Relationships9_0 + Relationships10_0

Merged$HIVStigma_1= relationshps0_1 + Relationships1_1 + Relationships2_1 + Relationships3_1 + 
  Relationships4_1 + Relationships5_1 + Relationships6_1 + Relationships7_1 + Relationships8_1 + 
  Relationships9_1 + Relationships10_1

Merged$HIVStigma_2= relationshps0_2 + Relationships1_2 + Relationships2_2 + Relationships3_2 + 
  Relationships4_2 + Relationships5_2 + Relationships6_2 + Relationships7_2 + Relationships8_2 + 
  Relationships9_2 + Relationships20_2

# Sex Risk
for (var in c('Num_Male_Partners_0','Num_female_partners_0','Num_trans_partners_0','P90SexforGoods_0','P90GivenGoods_exchSex_0')){
  Merged_new[is.na(Merged_new[,var]),var] = 0
}
attach(Merged_new)
Merged_new$Sex_Risk_0 = as.numeric((Num_Male_Partners_0+Num_female_partners_0+Num_trans_partners_0)>1)+P90SexforGoods_0+P90GivenGoods_exchSex_0
Merged_new$Sex_Risk_1 = as.numeric((Num_Male_Partners_1+Num_female_partners_1+Num_trans_partners_1)>1)+P90SexforGoods_1+P90GivenGoods_exchSex_1
Merged_new$Sex_Risk_2 = as.numeric((Num_Male_Partners_2+Num_female_partners_2+Num_trans_partners_2)>1)+P90SexforGoods_2+P90GivenGoods_exchSex_2
analysis$Sex_Risk_2 = Merged_new$Sex_Risk_2
analysis$Sex_Risk_1 = Merged_new$Sex_Risk_1
analysis$Sex_Risk_0 = Merged_new$Sex_Risk_0


# Select analysis dataset
analysis_exp = c('age_0','gender','sex_orient','race_0','education','income','Housing_0','Health_Insurance_0','Jail',
                 'Current_ART_0','P90Overdose_0','P30Injecting_0','ShareNeedle_0',
                 'TX_Assgn_0','TX_Assign_1','TX_Assign_2',
                 'Recieve_ART_0','Recieve_ART_1','Recieve_ART_2')
analysis_res = c('ViralLoad_0','Viral_loadLVL_1','Viral_loadLVL_1',
                 'P30Substanceuse_0','P30Substanceuse_1','P30Substanceuse_2',
                 'P90_healthissues_0','P90_healthissues_1','P90_healthissues_2',
                 'HIVStigma_0','HIVStigma_1','HIVStigma_2',
                 'HIV_Risk_0','HIV_Risk_1','HIV_Risk_2',
                 'P30Injecting_0','P30Injecting_1','P30Injecting_2',
                 'NumMissed_0','NumMissed_1','NumMissed_2'
                 )
analysis_data = Merged[,c(analysis_exp,analysis_res)]
save(file='analysis.Rdata',analysis_data)
