library(MESS)
load('analysis.Rdata')
analysis$black = as.numeric(analysis$race_0==2)
analysis$Latino = as.numeric(analysis$race_0==3)
analysis$heter = as.numeric(analysis$sex_orient==1)
analysis$women = as.numeric(analysis$gender==2)
analysis$highedu = as.numeric(analysis$education>=3)
analysis$highinc = as.numeric(analysis$income>=3)
analysis$stable_housing = as.numeric(analysis$Housing_0<=4)
analysis$jail_hist = as.numeric(analysis$Jail>=1)
save(file='analysis.Rdata',analysis)

# Sample decriptive summary table
load('analysis.Rdata')
library(geepack)
table(analysis$P90_healthissues_0,useNA='always')
table(analysis$P30Injecting_0,useNA='always')
table(analysis$ShareNeedle_0,useNA='always')
# For ViralLoad (No significant result)
var_exp = c('Age','women','heter','black','Latino','highedu','highinc','stable_housing','jail_hist',
            'Health_Insurance_0','CurART','Group','Participant_ID') # 'P90Overdose_0','P30Injecting_0','ShareNeedle_0' too many missings
var_res = c('VL','FU3_VL_Recent','FU6_VL_Recent')

subset_viral = analysis[,c(var_exp,var_res)]
gee_subset_viral = reshape(subset_viral,varying=var_res,
                           v.names='Viral_Load',idvar='Participant_ID',direction = 'long',times= c(0,3,6))
gee_subset_viral = gee_subset_viral[order(gee_subset_viral$Participant_ID),]

# Viral Load 
g = ggplot(data=gee_subset_viral,aes(x=as.factor(time),y=Viral_Load))
g+facet_grid(~Group,labeller=label_both)+geom_boxplot()+labs(x='month',y='Viral Load/count')
ggsave('output_new/viral_load_box_plot.png',width=5,height=5)

# summary of non-zero observations
tb=table(subset_viral$VL>0,subset_viral$Group)
sweep(tb,2,STATS = colSums(tb),FUN = '/')

tb=table(subset_viral$FU3_VL_Recent>0,subset_viral$Group)
sweep(tb,2,STATS = colSums(tb),FUN = '/')

tb=table(subset_viral$FU6_VL_Recent>0,subset_viral$Group)
sweep(tb,2,STATS = colSums(tb),FUN = '/')
model_gee_viral = geeglm(Viral_Load~(.-Participant_ID)*time,
                         family=gaussian, id=Participant_ID, corstr = 'ar1',
                         data = gee_subset_viral)

sink('output_new/viral_load_summary.txt') # Put output in a txt file
QIC(model_gee_viral) # 8.33e+11
summary(model_gee_viral)
sink()
model_gee_viral = update(model_gee_viral,Viral_Load~(.-heter:time))
QIC(model_gee_viral) # 8.33e+11
summary(model_gee_viral)
model_gee_viral = update(model_gee_viral,Viral_Load~(.-women:time))
QIC(model_gee_viral) # 8.33e+11
summary(model_gee_viral)
model_gee_viral = update(model_gee_viral,Viral_Load~(.-Group:time  ))
QIC(model_gee_viral) # 8.36e+11
summary(model_gee_viral)
model_gee_viral = update(model_gee_viral,Viral_Load~(.-stable_housing:time    ))
QIC(model_gee_viral) # 8.33e+11
summary(model_gee_viral)
model_gee_viral = update(model_gee_viral,Viral_Load~(.-black:time  ))
QIC(model_gee_viral) # 8.38e+11
summary(model_gee_viral)
model_gee_viral = update(model_gee_viral,Viral_Load~(.-Latino:time  ))
QIC(model_gee_viral) # 8.38e+11
summary(model_gee_viral)

#############################################################################################
# HIVStigma
var_exp = c('Age','women','heter','black','Latino','highedu','highinc','stable_housing','jail_hist',
            'Health_Insurance_0','CurART','Group','Participant_ID') # 'P90Overdose_0','P30Injecting_0','ShareNeedle_0' too many missings
var_res = c('HIVStigma_0','HIVStigma_1','HIVStigma_2')

subset = analysis[,c(var_exp,var_res)]
gee_subset = reshape(subset,varying=c('HIVStigma_0','HIVStigma_1','HIVStigma_2'),
                     v.names='HIVStigma',idvar='Participant_ID',direction = 'long',times= c(0,3,6))
gee_subset = gee_subset[order(gee_subset$Participant_ID),]

g = ggplot(data=gee_subset,aes(x=as.factor(time),y=HIVStigma))
g+facet_grid(~Group,labeller=label_both)+geom_boxplot(fill="#FF9999",color='gray50')+labs(x='month',y='HIV Stigma')
ggsave('output_new/HIV_Stigma.png',width=5,height=5)
f = ggplot(data=gee_subset,aes(x=time,y=HIVStigma,group=Participant_ID))
f+geom_line()+facet_grid(~Group,labeller=label_both)+stat_summary(aes(group=1),geom='line',fun.y = mean, color='blue')+labs(x='time/month',y='HIV Stigma')
# model fitting
model_gee_HIV = geeglm(HIVStigma~(.-Participant_ID+Group)*time,
                       family=gaussian, id=Participant_ID, corstr = 'ar1',
                       data = gee_subset)
QIC(model_gee_HIV) # 6286.0
summary(model_gee_HIV)
model_gee_HIV = update(model_gee_HIV,HIVStigma~(.-jail_hist:time  ))
QIC(model_gee_HIV) # 6284.6
summary(model_gee_HIV)
model_gee_HIV = update(model_gee_HIV,HIVStigma~(.-highedu:time   ))
sink('output_new/HIV_Stigma_summary.txt')
QIC(model_gee_HIV) # 6283.8
summary(model_gee_HIV)
sink()
model_gee_HIV = update(model_gee_HIV,HIVStigma~(.-stable_housing:time     ))
QIC(model_gee_HIV) # 6290.0  
summary(model_gee_HIV)
model_gee_HIV = update(model_gee_HIV,HIVStigma~(.-highinc:time    ))
QIC(model_gee_HIV) # 6291.2
summary(model_gee_HIV)
model_gee_HIV = update(model_gee_HIV,HIVStigma~(.-heter:time     ))
QIC(model_gee_HIV) # 6309.1
summary(model_gee_HIV)
model_gee_HIV = update(model_gee_HIV,HIVStigma~(.-Age:time   ))
QIC(model_gee_HIV) # 6308.5 
summary(model_gee_HIV)
model_gee_HIV = update(model_gee_HIV,HIVStigma~(.-women:time   ))
QIC(model_gee_HIV) # 6317.0
summary(model_gee_HIV)
model_gee_HIV = update(model_gee_HIV,HIVStigma~(.-Health_Insurance_0:time   ))
QIC(model_gee_HIV) # 6334.9 
summary(model_gee_HIV)

#########################################################################
# Drug Risk 
var_exp = c('Age','women','heter','black','Latino','highedu','highinc','stable_housing','jail_hist',
            'Health_Insurance_0','CurART','Group','Participant_ID') # 'P90Overdose_0','P30Injecting_0','ShareNeedle_0' too many missings
var_res = c('P30Injecting_0','P30Injecting_1','P30Injecting_2')

subset2 = analysis[,c(var_exp,var_res)]

gee_subset2 = reshape(subset2,varying=var_res,
                      v.names='P30Injecting',idvar='Participant_ID',direction = 'long',times= c(0,3,6))
gee_subset2 = gee_subset2[order(gee_subset2$Participant_ID),]

tmp0 = subset(gee_subset2,Group==0)
tb = table(tmp0$P30Injecting,tmp0$time,useNA = 'always')
tb = sweep(tb,2,colSums(tb),'/')

f = ggplot(gee_subset2,aes(fill=as.factor(P30Injecting),x=time,y=time+1))
f+geom_bar( stat="identity", position="fill")+labs(x='month',y='Percentage')+scale_fill_discrete(name='Drug Risk')+
  scale_x_continuous(breaks=c(0,3,6))+facet_grid(~Group,labeller=label_both)
ggsave('output_new/Drug_Risk.png',width=5,height=5)
model_gee_drug_risk = geeglm(P30Injecting~(.-Participant_ID)*time-highinc:time-Latino:time,
                             family=binomial, id=Participant_ID, corstr = 'ar1',
                             data = gee_subset2) #highinc:time, Latino:time ferfectly separart the data
QIC(model_gee_drug_risk) #317.02
summary(model_gee_drug_risk)
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.-Group:time))
QIC(model_gee_drug_risk) #221
summary(model_gee_drug_risk)
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.-black:time))
QIC(model_gee_drug_risk) #218.6
summary(model_gee_drug_risk)
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.-stable_housing:time))
QIC(model_gee_drug_risk) #260.04
summary(model_gee_drug_risk)
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.-Age:time))
QIC(model_gee_drug_risk) #256.34
summary(model_gee_drug_risk)
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.-Health_Insurance_0:time))
QIC(model_gee_drug_risk) #216.12
summary(model_gee_drug_risk)
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.-CurART:time))
sink('output_new/Drug_Risk_summary.txt')
QIC(model_gee_drug_risk) #211
summary(model_gee_drug_risk)
sink()
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.-heter:time))
QIC(model_gee_drug_risk) #215.3
summary(model_gee_drug_risk)
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.-women:time))
QIC(model_gee_drug_risk) #213.4
summary(model_gee_drug_risk)
###########################################################################################
# Look for the variable that causes perfect separation
n = 0
for (term in attr(model_gee_drug_risk$terms,'term.labels')){
  name = paste0('model_',term)
  formula = as.formula(('P30Injecting~(.-Participant_ID)*time'))
  if (n>1){
    formula = as.formula(paste0('P30Injecting~',as.character(model_test$formula)[3],'-time:',term))
  }
  model_test = geeglm(formula = formula,
                      family=binomial, id=Participant_ID, corstr = 'ar1',
                      data = gee_subset2)
  print(term)
  print(model_test$coefficients[1:5])
  n = n+1
} 

############################################################################################
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.-heter:time))
QIC(model_gee_drug_risk) #1334.6?
summary(model_gee_drug_risk) # perfectly separable case
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.-CurART:time))
QIC(model_gee_drug_risk) #110.8
summary(model_gee_drug_risk)
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.-Age:time))
QIC(model_gee_drug_risk) #115.6
summary(model_gee_drug_risk)
model_gee_drug_risk=update(model_gee_drug_risk,P30Injecting~(.+Group:time))
QIC(model_gee_drug_risk) #116.5
summary(model_gee_drug_risk)
###########################################################################################

# Comparing 2 phases
subset_drug_risk_2 = analysis[,c(var_exp,var_res)]
subset_drug_risk_2$P30Injecting_2 = NULL
subset_drug_risk_2 = subset_drug_risk_2[rowSums(is.na(subset_drug_risk_2))==0,]
gee_subset_drug_risk_2 = reshape(subset_drug_risk_2,varying=c('P30Injecting_0','P30Injecting_1'),
                                 v.names='P30Injecting',idvar='Participant_ID',direction = 'long',times= c(0,3))
gee_subset_drug_risk_2 = gee_subset_drug_risk_2[order(gee_subset_drug_risk_2$Participant_ID),]
model_gee_drug_risk_2 = geeglm(P30Injecting~(.-Participant_ID-Latino-heter)*time,
                               family=binomial, id=Participant_ID, corstr = 'ar1',
                               data = gee_subset_drug_risk_2)
summary(model_gee_drug_risk_2)
model_gee_drug_risk_2 = update(model_gee_drug_risk_2,P30Injecting~(.-women:time-Health_Insurance_0-Health_Insurance_0:time))
summary(model_gee_drug_risk_2) # perfectly seperable case

# CD4
var_exp = c('Age','women','heter','black','Latino','highedu','highinc','stable_housing','jail_hist',
            'Health_Insurance_0','CurART','Group','Participant_ID') # 'P90Overdose_0','P30Injecting_0','ShareNeedle_0' too many missings
var_res = c('CD4','FU3_CD4_Recent','FU6_CD4_Recent')

subset_CD4 = analysis[,c(var_exp,var_res)]
subset_CD4$Group = as.factor(subset_CD4$Group)
levels(subset_CD4$Group) = c('Control','Intervention')
subset_CD4 = subset_CD4[-32,] #delete abnormal record with cd4=4236

gee_subset_CD4 = reshape(subset_CD4,varying=var_res,
                         v.names='CD4_Load',idvar='Participant_ID',direction = 'long',times= c(0,3,6))
gee_subset_CD4 = gee_subset_CD4[order(gee_subset_CD4$Participant_ID),]

f = ggplot(data=gee_subset_CD4,aes(x=time,y=CD4_Load,group=Participant_ID))
f+geom_line()+facet_grid(~Group,labeller=label_both)+stat_summary(aes(group=1),geom='line',fun.y = mean, color='blue')+labs(x='time/month',y='CD4 Count/cells per mm3 ')
ggsave('output_new/CD4.png',width=5,height=5)
model_gee_CD4 = geeglm(CD4_Load~(.-Participant_ID)*time,
                       family=gaussian, id=Participant_ID, corstr = "ar1",
                       data = gee_subset_CD4)
QIC(model_gee_CD4) #7582271
summary(model_gee_CD4)
sink()
model_gee_CD4=update(model_gee_CD4,CD4_Load~.-Group:time) #QIC increases
QIC(model_gee_CD4) #7.58e+06
summary(model_gee_CD4)
model_gee_CD4=update(model_gee_CD4,CD4_Load~.-highedu:time) #QIC increases
QIC(model_gee_CD4) #7.59e+06
summary(model_gee_CD4)
model_gee_CD4=update(model_gee_CD4,CD4_Load~.-highinc:time) #QIC increases
QIC(model_gee_CD4) #7586075
summary(model_gee_CD4)
model_gee_CD4=update(model_gee_CD4,CD4_Load~.-stable_housing:time) #QIC increases
QIC(model_gee_CD4) #7.62e+06
summary(model_gee_CD4)
model_gee_CD4=update(model_gee_CD4,CD4_Load~.-women:time) #QIC increases
QIC(model_gee_CD4) #7.65e+06
summary(model_gee_CD4)
sink('output_new/CD4_summary.txt')
model_gee_CD4=update(model_gee_CD4,CD4_Load~.-CurART:time) #QIC increases
QIC(model_gee_CD4) #7.66e+06
summary(model_gee_CD4)
sink()
model_gee_CD4=update(model_gee_CD4,CD4_Load~.-jail_hist:time) #QIC increases
QIC(model_gee_CD4) #7.81e+06
summary(model_gee_CD4)
model_gee_CD4=update(model_gee_CD4,CD4_Load~.-Health_Insurance_0:time) #QIC increases
QIC(model_gee_CD4) #7.98e+06
summary(model_gee_CD4)
# f = ggplot(data=gee_subset_CD4_2,aes(x=time,y=CD4_Load,group=Participant_ID))
# f+geom_line(color='grey40')+facet_grid(~Group,labeller=label_both)+stat_summary(aes(group=1),geom='line',fun.y = mean, color='blue',size=0.8)+labs(x='time/month',y='CD4 Load')
# ggsave()
# 
# plot(x=c(0,3,6),y=rep(0,3),type='n',ylim=c(0,max(subset_CD4[,c('CD4','FU3_CD4_Recent','FU6_CD4_Recent')])))
# for(n in 1:nrow(subset_CD4)){
#   lines(c(CD4[n],FU3_CD4_Recent[n],FU6_CD4_Recent[n])~c(0,3,6),data=subset_CD4,
#        col=ifelse(Group[n]==1,'red','blue'))
# }

# two phase
var_exp = c('Age','women','heter','black','Latino','highedu','highinc','stable_housing','jail_hist',
            'Health_Insurance_0','CurART','Group','Participant_ID') # 'P90Overdose_0','P30Injecting_0','ShareNeedle_0' too many missings
var_res = c('CD4','FU3_CD4_Recent','FU6_CD4_Recent')

subset_CD4_2 = analysis[,c(var_exp,var_res)]
subset_CD4_2$FU6_CD4_Recent = NULL
#subset_CD4_2 = subset_CD4_2[rowSums(is.na(subset_CD4_2))==0,] #delete obervations with any NA's w/21 left
subset_CD4_2$Group = as.factor(subset_CD4_2$Group)
levels(subset_CD4_2$Group) = c('Control','Intervention')

gee_subset_CD4_2 = reshape(subset_CD4_2,varying=c('CD4','FU3_CD4_Recent'),
                           v.names='CD4_Load',idvar='Participant_ID',direction = 'long',times= c(0,3))
gee_subset_CD4_2 = gee_subset_CD4_2[order(gee_subset_CD4_2$Participant_ID),]

model.gee_CD4_2 = geeglm(CD4_Load~(.-Participant_ID-Latino)*time,
                         family=gaussian, id=Participant_ID, corstr = 'ar1',
                         data = gee_subset_CD4_2) # no Latino since no Latino in the subset

summary(model.gee_CD4_2)

# stepwise model selection through QIC

model.gee_CD4_2 = update(model.gee_CD4_2,CD4_Load~.-Age:time-highedu:time-Health_Insurance_0:time 
                         -women:time-highinc:time-heter:time-stable_housing:time-jail_hist:time  )
summary(model.gee_CD4_2)
QIC(model.gee_CD4_2)

#######################################################################################
# Sex Risk
var_exp = c('Age','women','heter','black','Latino','highedu','highinc','stable_housing','jail_hist',
            'Health_Insurance_0','CurART','Group','Participant_ID') # 'P90Overdose_0','P30Injecting_0','ShareNeedle_0' too many missings
var_res = c('Sex_Risk_0','Sex_Risk_1','Sex_Risk_2')

subset_sex_risk = analysis[,c(var_exp,var_res)]
gee_subset_sex_risk = reshape(subset_sex_risk,varying=var_res,
                              v.names='sex_risk',idvar='Participant_ID',direction = 'long',times= c(0,3,6))
gee_subset_sex_risk = gee_subset_sex_risk[order(gee_subset_sex_risk$Participant_ID),]
# plot
f = ggplot(gee_subset_sex_risk,aes(fill=as.factor(sex_risk),x=time,y=time+1))
f+geom_bar( stat="identity", position="fill")+labs(x='month',y='Percentage')+scale_fill_discrete(name='Sex Risk')+
  scale_x_continuous(breaks=c(0,3,6))+facet_grid(~Group,labeller=label_both)
ggsave('output_new/Sex_Risk.png',width=5,height=5)

model.gee_sex_risk = geeglm(cbind(sex_risk,3-sex_risk)~(.-Participant_ID)*time,
                            family=binomial, id=Participant_ID, corstr='ar1',
                            data = gee_subset_sex_risk)
QIC(model.gee_sex_risk) #285.8 
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-highedu:time))
QIC(model.gee_sex_risk) #280.2  
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-highinc:time))
QIC(model.gee_sex_risk) #276.0   
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-Health_Insurance_0:time))
QIC(model.gee_sex_risk) #277.0   
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-stable_housing:time))
QIC(model.gee_sex_risk) #276.6  
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-jail_hist:time ))
QIC(model.gee_sex_risk) #274.3    
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-Latino:time))

QIC(model.gee_sex_risk) #273.3     
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-CurART:time ))
QIC(model.gee_sex_risk) #262.6     
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-women:time ))
QIC(model.gee_sex_risk) #256.7     
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-highedu ))
QIC(model.gee_sex_risk) #247.9    
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-stable_housing ))
QIC(model.gee_sex_risk) #240.3  
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-jail_hist ))
QIC(model.gee_sex_risk) #237.8  
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-Latino ))
QIC(model.gee_sex_risk) #233.4
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-Health_Insurance_0 ))
QIC(model.gee_sex_risk) #225.7  
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-highinc ))
QIC(model.gee_sex_risk) #219 
summary(model.gee_sex_risk)
model.gee_sex_risk = update(model.gee_sex_risk,cbind(sex_risk,3-sex_risk)~(.-women ))
sink('output_new/Sex_Risk_summary.txt')
QIC(model.gee_sex_risk) #211.5
summary(model.gee_sex_risk)
sink()

########################################################################################
# Substance Abuse
var_exp = c('Age','women','heter','race_0','highedu','highinc','stable_housing','jail_hist',
            'Health_Insurance_0','CurART','Group','Participant_ID') # 'P90Overdose_0','P30Injecting_0','ShareNeedle_0' too many missings
var_res = c('P30Substanceuse_0','P30Substanceuse_1','P30Substanceuse_2')

subset_sub = analysis[,c(var_exp,var_res)]
gee_subset_sub = reshape(subset_sub,varying=var_res,
                         v.names='P30Substanceuse',idvar='Participant_ID',direction = 'long',times= c(0,3,6))
gee_subset_sub = gee_subset_sub[order(gee_subset_sub$Participant_ID),]
f = ggplot(gee_subset_sub,aes(fill=as.factor(P30Substanceuse),x=time,y=time+1))
f+geom_bar( stat="identity", position="fill")+labs(x='month',y='Percentage')+scale_fill_discrete(name='Substance overuse')+
  scale_x_continuous(breaks=c(0,3,6))+facet_grid(~Group,labeller=label_both)
ggsave('output_new/Substance_Overuse.png',width=5,height=5)
model.gee_sub = geeglm(P30Substanceuse~.-Participant_ID+time:Group-race_0,
                       family=binomial, id=Participant_ID, corstr = 'ar1',
                       data = gee_subset_sub) # 38/43 are black
summary(model.gee_sub)
foo=gee_subset_sub[gee_subset_sub$P30Substanceuse==0,] # reason why 
model_gee_sub = geeglm(P30Substanceuse~(.-Participant_ID-race_0)+
                         (Age+women+highedu+highinc+Health_Insurance_0+CurART+Group):time,
                       family=binomial, id=Participant_ID, corstr = 'ar1',
                       data = gee_subset_sub) # completely divided

for (term in attr(model_test$terms,'term.labels')){
  name = paste0('model_',term)
  formula = as.formula(paste0('P30Substanceuse~(.-Participant_ID-race_0)+time:',term))
  model_test = geeglm(formula = formula,
                      family=binomial, id=Participant_ID, corstr = 'ar1',
                      data = gee_subset_sub)
  print(term)
  print(model_test$coefficients[1:5])
  assign(name,model_test)
} 
tmp0 = subset(gee_subset_sub, Group==0)
tb = table(tmp0$P30Substanceuse,tmp0$time)
sweep(tb,2,colSums(tb),FUN='/')

tmp1 = subset(gee_subset_sub, Group==1)
tb = table(tmp1$P30Substanceuse,tmp1$time)
sweep(tb,2,colSums(tb),FUN='/')

QIC(model_gee_sub) # 623
summary(model_gee_sub)


# Baseline summary
save(file='Kaiser_analysis.Rdata',list = c(grep('gee_subset+',ls(),value=T),grep('subset+',ls(),value=T),grep('model_gee+',ls(),value=T)))
var_exp = c('women','heter','black','Latino','highedu','highinc','stable_housing','jail_hist',
            'Health_Insurance_0','CurART') 
tb = NULL
for(k in var_exp){
  new=t(as.matrix(table(analysis[,k],analysis$Group)[2,]))
  rownames(new)=k
  tb = rbind(tb,new) 
}
tb=sweep(tb,2,table(analysis$Group),'/')
