
# SEM models with: PCA from keystroke dynamics, NIH toolbox tests
# Author: Emma Ning

# Disable scientific notation
# options(scipen = 999)
library(tidyverse)
library(data.table)
library(OpenMx)


# Mac directory
dat_dir <- '' #specify directory data is saved in
dat <- read.csv(paste0(dat_dir,'')) #post-processed data file name, should be 'PCA_NoCut_6-12-24_Logged_NewNames.csv'

# Rename visits
dat <- dat %>%
  mutate(Study_Window2=ifelse(Study_Window2=='Visit 1-2','12','23'))


# get variable names
nihNames <- grep('nih',names(dat),ignore.case=T,value=T)
qidsNames <- grep('qidsc',names(dat),ignore.case=T,value=T)
trailNames <- grep('^trail',names(dat),ignore.case=T,value=T)

# Pivot wider
temp1 <- dat %>%
  select(healthCode,age,gender,Group,
         all_of(c(nihNames,qidsNames,trailNames))) %>% distinct()

temp2 <- dat %>%
  select(healthCode,Study_Window2,PC1,PC2,PC3,PC4) %>%
  pivot_wider(names_from = Study_Window2,
              values_from = c(PC1:PC4))

wide_dat <- left_join(temp2,temp1)


### Define variables ###
kp2 <- grep('PC.*12$', names(wide_dat), value=T, ignore.case=T)
kp3 <- grep('PC.*23$', names(wide_dat), value=T, ignore.case=T)
kp2 <- kp2[1:2]
kp3 <- kp3[1:2]
kpvars <- c(kp2,kp3)
nih2 <- grep('nih_t.*v2$', names(wide_dat), value=T, ignore.case=T)
nih3 <- grep('nih_t.*v3$', names(wide_dat), value=T, ignore.case=T)
nihvars <- c(nih2,nih3)
manVars <- c(kpvars, nihvars) #manifest variables
lVars <- c('cog2','cog3') #latent variables


cor(wide_dat[,manVars],use='pairwise.complete.obs') #check correlation matrix
# looking at the variances
apply(wide_dat[,manVars],2,var,na.rm=T)
# the means
apply(wide_dat[,manVars],2,mean,na.rm=T)


#### Variable transformation
wide_dat[,nihvars] <- wide_dat[,nihvars]/10 # RESCALE NIH vars

apply(wide_dat[,manVars],2,var,na.rm=T)
apply(wide_dat[,manVars],2,mean,na.rm=T)

#### KP & NIH: Base model ####
m0 <- mxModel("NIH Base model",
              # preamble
              mxData(wide_dat[,manVars], "raw"),
              type="RAM",
              manifestVars = manVars,
              latentVars = lVars,
              
              # (residual) variances
              mxPath(manVars, arrows=2, values=1, free=T,
                     labels=paste('var_',manVars,sep='')),
              
              # latent variances and covariance
              mxPath(lVars, arrows=2, connect='unique.pairs', free=T,
                     values=c(1,.5,1),labels=c('varL2','covL','varL3')),
              
              # Factor loadings
              mxPath('cog2',nih2,arrows=1,free=c(F,T,T,T),values=c(1,1,1,1),
                     labels=c('l1','l2','l3','l4')),
              
              mxPath('cog3',nih3,arrows=1,free=c(F,T,T,T),values=c(1,1,1,1),
                     labels=c('l1','l2','l3','l4')),

              # mxPath('cog3',nih3,arrows=1,free=c(F,T,T,T),values=c(1,1,1,1),
              #        labels=c('l5','l6','l7','l8')),
              
              # means
              mxPath("one", c(manVars,lVars), values=c(rep(1,4),rep(4,8),0,0), free=c(rep(T,12),F,F),
                     labels=c(paste('m_',manVars,sep=''),NA,NA)),
              
              # Manifest factor covariances
              mxPath(from=nih2, to=nih3, arrows=2, connect='single',free=T,
                     labels=c('cov_nih1','cov_nih2','cov_nih3','cov_nih4')), 
              
              mxPath(from=kp2, to=kp3, arrows=2, connect='single',free=T,
                     labels=c('cov_pc1','cov_pc2')), 
              
              mxPath(from=kp2, to=kp2, arrows=2, connect='unique.bivariate',free=T,
                     labels='cov_PC_within'),
              
              mxPath(from=kp3, to=kp3, arrows=2, connect='unique.bivariate',free=T,
                     labels='cov_PC_within'),
             
              # regressions to latent vars
              mxPath(from=kp2, 
                     to='cog2', free=c(T,T),
                     labels=c('b1','b2')), 
              
              mxPath(from=kp3, 
                     to='cog3', free=c(T,T),
                     labels=c('b3','b4'))
              
              
)

res_m0 <- mxRun(m0)
RefModels_m0 <- mxRefModels(res_m0, run=TRUE)
summary(res_m0,refModels=RefModels_m0)
mxStandardizeRAMpaths(res_m0,SE=FALSE)
mxStandardizeRAMpaths(res_m0,SE=TRUE)

# #### Test if factor loadings can be fixed between visit 2 and 3 & if covariance within visit can be fixed
# fixLoadings1 <- mxModel(m0,
#                        mxPath('cog3',nih3,arrows=1,free=c(F,T,T,T),values=c(1,1,1,1),
#                               labels=c('l1','l2','l3','l4')),
#                        
#                        mxPath(from=kp2, to=kp2, arrows=2, connect='unique.bivariate',free=T,
#                               labels='cov_PC_within'),
#                        
#                        mxPath(from=kp3, to=kp3, arrows=2, connect='unique.bivariate',free=T,
#                               labels='cov_PC_within'))
# res_fixLoadings1 <- mxRun(fixLoadings1)
# mxCompare(res_m0,res_fixLoadings1)

##### Set regressions to fixed #####
# This step is to get p-values
m1 <- mxModel(m0,
              mxPath(from=kp2, 
                     to='cog2', free=c(T,F),
                     labels=c('b1','b2')))
res_m1 <- mxRun(m1)
mxCompare(res_m0,res_m1)
# # mxCheckIdentification(res_m1)
# clpdRefModels_m1 <- mxRefModels(res_m1, run=TRUE)
# summary(res_m1, refModels=clpdRefModels_m1)

m2 <- mxModel(m0,
              mxPath(from=kp3, 
                     to='cog3', free=c(T,F),
                     labels=c('b3','b4')))
res_m2 <- mxRun(m2)
mxCompare(res_m0,res_m2)
# summary(res_m2, refModels=clpdRefModels_m1)

m3 <- mxModel(m0,
              # regressions to latent vars
              mxPath(from=kpvars[1], 
                     to='cog2', free=T,
                     labels=c('b_kp_cog')), 
              
              mxPath(from=kpvars[2], 
                     to='cog3', free=T,
                     labels=c('b_kp_cog')))
res_m3 <- mxRun(m3)
mxCompare(res_m0,res_m3)
