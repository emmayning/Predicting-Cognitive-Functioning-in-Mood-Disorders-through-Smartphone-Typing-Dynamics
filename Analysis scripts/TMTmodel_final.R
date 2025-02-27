
# SEM models with: PCA from keystroke dynamics, TMT-B
# Author: Emma Ning

# Disable scientific notation
options(scipen = 999)
library(tidyverse)
library(data.table)
library(OpenMx)


# Mac directory
dat_dir <- '' #specify directory data is saved in
dat <- read.csv(paste0(dat_dir,'')) #post-processed data file name, should be 'PCA_NoCut_2-20-2025.csv'

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

# Log transform paper trail time
wide_dat <- wide_dat %>%
  mutate(log_trail_a_time_V2 = log(trail_a_time_V2),
         log_trail_a_time_V3 = log(trail_a_time_V3),
         log_trail_b_time_V2 = log(trail_b_time_V2),
         log_trail_b_time_V3 = log(trail_b_time_V3))

### Define variables ###
kp2 <- grep('PC.*12$', names(wide_dat), value=T, ignore.case=T)
kp3 <- grep('PC.*23$', names(wide_dat), value=T, ignore.case=T)
kp2 <- kp2[1:2]
kp3 <- kp3[1:2]
kpvars <- c(kp2,kp3)
trail2 <- grep('log_trail_b.*v2$', names(wide_dat), value=T, ignore.case=T)
trail3 <- grep('log_trail_b.*v3$', names(wide_dat), value=T, ignore.case=T)
trailvars <- c(trail2,trail3) 
manVars <- c(kpvars, trailvars) #manifest variables


cor(wide_dat[,manVars],use='pairwise.complete.obs')
# looking at the variances
apply(wide_dat[,manVars],2,var,na.rm=T)
# the means
apply(wide_dat[,manVars],2,mean,na.rm=T)

# # Remove an outlier (run both ways)
# wide_dat <- wide_dat[wide_dat$healthCode!='rehvcm7ls7mqlp_giz4uahdv',]

#### KP & NIH: Base model ####
m0 <- mxModel("TMT Base model",
              # preamble
              mxData(wide_dat[,manVars], "raw"),
              type="RAM",
              manifestVars = manVars,
              
              # variances
              mxPath(manVars, arrows=2, values=1, free=T,
                     labels=paste('var_',manVars,sep='')),
              
              # means
              mxPath("one", manVars, values=c(rep(1,4),rep(3,2)), free=T,
                     labels=c(paste('m_',manVars,sep=''))),
              
              # Manifest factor covariances
              # mxPath(from=kp2, to=kp3, arrows=2, connect='single',free=T,
              #        labels=c('cov_pc1','cov_pc2')),
              
              mxPath(from=kp2, to=kp3, arrows=2, connect='single',free=T,
                     labels=c('cov_pc_across','cov_pc_across')),
              
              mxPath(from=trail2, to=trail3, arrows=2, connect='single',free=T,
                     labels='cov_tb'), 
              
              mxPath(from=kp2, to=kp2, arrows=2, connect='unique.bivariate',free=T,
                     labels='cov_PC_within'),
              
              mxPath(from=kp3, to=kp3, arrows=2, connect='unique.bivariate',free=T,
                     labels='cov_PC_within'),
              
              # regressions to latent vars
              mxPath(from=kp2, 
                     to=trail2, connect='unique.pairs', free=T,
                     labels=c('b1','b2')), 
              
              mxPath(from=kp3, 
                     to=trail3, connect='unique.pairs', free=T,
                     labels=c('b3','b4'))
              
              
)

res_m0 <- mxRun(m0)
RefModels_m0 <- mxRefModels(res_m0, run=TRUE)
summary(res_m0,refModels=RefModels_m0)
mxStandardizeRAMpaths(res_m0,SE=FALSE)
mxStandardizeRAMpaths(res_m0,SE=TRUE)

# #### Test if covariance within visit can be fixed
fixLoadings1 <- mxModel(m0,
                        # mxPath(from=kp2, to=kp2, arrows=2, connect='unique.bivariate',free=T,
                        #        labels='cov_PC_within'),
                        # 
                        # mxPath(from=kp3, to=kp3, arrows=2, connect='unique.bivariate',free=T,
                        #        labels='cov_PC_within')
                        # ,
                        mxPath(from=kp2, to=kp3, arrows=2, connect='single',free=T,
                               labels=c('cov_pc_across','cov_pc_across'))
)
res_fixLoadings1 <- mxRun(fixLoadings1)
mxCompare(res_m0,res_fixLoadings1)

#### Test if regressions across timepoints are the same
fixRegs <- mxModel(
  
  # m0,
  # # Testing regression paths for PC1
  # mxPath(from=kpvars[1],
  #        to=trail2, free=T,
  #        labels='beta_PC1'),
  # 
  # mxPath(from=kpvars[3],
  #        to=trail3, free=T,
  #        labels='beta_PC1')
  # ,
  # 
  mConfigural,
  # #                  # Testing regression paths for PC2
  mxPath(from=kpvars[2],
         to=trail2, free=T,
         labels='beta_PC2'),
  
  mxPath(from=kpvars[4],
         to=trail3, free=T,
         labels='beta_PC2')
  
)
res_fixRegs <- mxRun(fixRegs)
mxCompare(res_m0,res_fixRegs)
mxCompare(res_mConfigural,res_fixRegs)


##### Final model before multiple group
mConfigural <- mxModel(m0,
                       
                       # Testing regression paths for PC1
                       mxPath(from=kpvars[1],
                              to=trail2, free=T,
                              labels='beta_PC1'),
                       
                       mxPath(from=kpvars[3],
                              to=trail3, free=T,
                              labels='beta_PC1')
                       
                       ,
                       
                       # Testing regression paths for PC2
                       mxPath(from=kpvars[2],
                              to=trail2, free=T,
                              labels='beta_PC2'),
                       
                       mxPath(from=kpvars[4],
                              to=trail3, free=T,
                              labels='beta_PC2')
)

res_mConfigural <- mxRun(mConfigural)
summary(res_mConfigural,refModels=RefModels_m0)
mxCompare(res_m0,res_mConfigural)
mxStandardizeRAMpaths(res_mConfigural,SE=TRUE)


##### Set regressions to fixed #####
m1 <- mxModel(mConfigural,
              # Testing regression paths for PC1
              mxPath(from=kpvars[1],
                     to=trail2, free=F,
                     labels='beta_PC1'),
              
              mxPath(from=kpvars[3],
                     to=trail3, free=F,
                     labels='beta_PC1')
              # ,
              # 
              # Testing regression paths for PC1
              # mxPath(from=kpvars[2],
              #        to=trail2, free=F,
              #        labels='beta_PC2'),
              # 
              # mxPath(from=kpvars[4],
              #        to=trail3, free=F,
              #        labels='beta_PC2')
)
res_m1 <- mxRun(m1)
mxCompare(res_mConfigural,res_m1)


#### Multiple groups ####
datMD <- wide_dat[wide_dat$Group=="Mood Disorder",]
datHC <- wide_dat[wide_dat$Group=="Healthy Control",]

# look at means and variances
apply(datMD[,manVars],2,mean,na.rm=T)
apply(datHC[,manVars],2,mean,na.rm=T)

apply(datMD[,manVars],2,var,na.rm=T)
apply(datHC[,manVars],2,var,na.rm=T)


# Models
# All variables have their own means and variances
multiMD0 <- mxModel("base MD",
                    # preamble
                    mxData(datMD[,manVars], "raw"),
                    type="RAM",
                    manifestVars = manVars,
                    
                    # variances
                    mxPath(manVars, arrows=2, values=1, free=T,
                           labels=paste('var_MD_',manVars,sep='')),
                    
                    # means
                    mxPath("one", manVars, values=c(rep(1,4),rep(3,2)), free=T,
                           labels=c(paste('m_MD_',manVars,sep=''))),
                    
                    # Manifest factor covariances
                    mxPath(from=kp2, to=kp3, arrows=2, connect='single',free=T,
                           labels=c('cov_pc_across_MD','cov_pc_across_MD')),
                    
                    mxPath(from=trail2, to=trail3, arrows=2, connect='single',free=T,
                           labels='cov_tb_MD'), 
                    
                    mxPath(from=kp2, to=kp2, arrows=2, connect='unique.bivariate',free=T,
                           labels='cov_PC_within_MD'),
                    
                    mxPath(from=kp3, to=kp3, arrows=2, connect='unique.bivariate',free=T,
                           labels='cov_PC_within_MD'),
                    
                    # regressions to latent vars
                    mxPath(from=kp2, 
                           to=trail2, connect='unique.pairs', free=T,
                           labels=c('beta_PC1','beta_PC2')), 
                    
                    mxPath(from=kp3, 
                           to=trail3, connect='unique.pairs', free=T,
                           labels=c('beta_PC1','beta_PC2'))
                    
                    
)

multiHC0 <- mxModel("base HC",
                    # preamble
                    mxData(datHC[,manVars], "raw"),
                    type="RAM",
                    manifestVars = manVars,
                    
                    # variances
                    mxPath(manVars, arrows=2, values=1, free=T,
                           labels=paste('var_HC_',manVars,sep='')),
                    
                    # means
                    mxPath("one", manVars, values=c(rep(1,4),rep(3,2)), free=T,
                           labels=c(paste('m_HC_',manVars,sep=''))),
                    
                    # Manifest factor covariances
                    mxPath(from=kp2, to=kp3, arrows=2, connect='single',free=T,
                           labels=c('cov_pc_across_HC','cov_pc_across_HC')),
                    
                    mxPath(from=trail2, to=trail3, arrows=2, connect='single',free=T,
                           labels='cov_tb_HC'), 
                    
                    mxPath(from=kp2, to=kp2, arrows=2, connect='unique.bivariate',free=T,
                           labels='cov_PC_within_HC'),
                    
                    mxPath(from=kp3, to=kp3, arrows=2, connect='unique.bivariate',free=T,
                           labels='cov_PC_within_HC'),
                    
                    # regressions to latent vars
                    mxPath(from=kp2, 
                           to=trail2, connect='unique.pairs', free=T,
                           labels=c('beta_PC1','beta_PC2')), 
                    
                    mxPath(from=kp3, 
                           to=trail3, connect='unique.pairs', free=T,
                           labels=c('beta_PC1','beta_PC2'))
                    
                    
)

# container
bucketBase <- mxModel("base multigroup",
                      multiMD0, multiHC0,
                      mxFitFunctionMultigroup(c("base MD", "base HC")))

# run it
resBase <- mxRun(bucketBase)
summary(resBase)
# Ref model
refBase <- mxRefModels(resBase, run=TRUE)


# Free regression paths

# Free both
multiMD0_0 <- mxModel(multiMD0,
                      name="multiMD0_0",
                      mxPath(from=kp2[1], 
                             to=trail2, free=T,
                             labels='beta_PC1_MD'),
                      
                      mxPath(from=kp3[1], 
                             to=trail3, free=T,
                             labels='beta_PC1_MD'),
                      
                      mxPath(from=kp2[2], 
                             to=trail2, free=T,
                             labels='beta_PC2_MD'),
                      
                      mxPath(from=kp3[2], 
                             to=trail3, free=T,
                             labels='beta_PC2_MD')
)


multiHC0_0 <- mxModel(multiHC0,
                      name='multiHC0_0',
                      mxPath(from=kp2[1], 
                             to=trail2, free=T,
                             labels='beta_PC1_HC'),
                      
                      mxPath(from=kp3[1], 
                             to=trail3, free=T,
                             labels='beta_PC1_HC'),
                      
                      mxPath(from=kp2[2], 
                             to=trail2, free=T,
                             labels='beta_PC2_HC'),
                      
                      mxPath(from=kp3[2], 
                             to=trail3, free=T,
                             labels='beta_PC2_HC')
)

freeBoth <- mxModel("free both PC1 and PC2 reg",
                    multiMD0_0, multiHC0_0,
                    mxFitFunctionMultigroup(c("multiMD0_0", "multiHC0_0")))

# run model
res_freeBoth <- mxRun(freeBoth)
summary(res_freeBoth, refModels=refBase)

mxCompare(res_freeBoth, resBase)

mxStandardizeRAMpaths(res_freeBoth,SE=TRUE)


# Finally models
# Do not need estimate different PC1 and PC2 reg paths across the two groups
summary(resBase, refModels=refBase)
mxStandardizeRAMpaths(resBase,SE=TRUE)


#### Multiple groups, letting the means and variances be fixed to make sure the effect is not due to that (in case reviewers ask) ####

checkMD0 <- mxModel("checkMD0",
                    # preamble
                    mxData(datMD[,manVars], "raw"),
                    type="RAM",
                    manifestVars = manVars,
                    
                    # variances
                    mxPath(manVars, arrows=2, values=1, free=T,
                           labels=paste('var_',manVars,sep='')),
                    
                    # means
                    mxPath("one", manVars, values=c(rep(1,4),rep(3,2)), free=T,
                           labels=c(paste('m_',manVars,sep=''))),
                    
                    # Manifest factor covariances
                    mxPath(from=kp2, to=kp3, arrows=2, connect='single',free=T,
                           labels=c('cov_pc_across','cov_pc_across')),
                    
                    mxPath(from=trail2, to=trail3, arrows=2, connect='single',free=T,
                           labels='cov_tb'), 
                    
                    mxPath(from=kp2, to=kp2, arrows=2, connect='unique.bivariate',free=T,
                           labels='cov_PC_within'),
                    
                    mxPath(from=kp3, to=kp3, arrows=2, connect='unique.bivariate',free=T,
                           labels='cov_PC_within'),
                    
                    # regressions to latent vars
                    mxPath(from=kp2, 
                           to=trail2, connect='unique.pairs', free=T,
                           labels=c('beta_PC1','beta_PC2')), 
                    
                    mxPath(from=kp3, 
                           to=trail3, connect='unique.pairs', free=T,
                           labels=c('beta_PC1','beta_PC2'))
                    
)

checkHC0 <- mxModel("checkHC0",
                    # preamble
                    mxData(datHC[,manVars], "raw"),
                    type="RAM",
                    manifestVars = manVars,
                    
                    # variances
                    mxPath(manVars, arrows=2, values=1, free=T,
                           labels=paste('var_',manVars,sep='')),
                    
                    # means
                    mxPath("one", manVars, values=c(rep(1,4),rep(3,2)), free=T,
                           labels=c(paste('m_',manVars,sep=''))),
                    
                    # Manifest factor covariances
                    mxPath(from=kp2, to=kp3, arrows=2, connect='single',free=T,
                           labels=c('cov_pc_across','cov_pc_across')),
                    
                    mxPath(from=trail2, to=trail3, arrows=2, connect='single',free=T,
                           labels='cov_tb'), 
                    
                    mxPath(from=kp2, to=kp2, arrows=2, connect='unique.bivariate',free=T,
                           labels='cov_PC_within'),
                    
                    mxPath(from=kp3, to=kp3, arrows=2, connect='unique.bivariate',free=T,
                           labels='cov_PC_within'),
                    
                    # regressions to latent vars
                    mxPath(from=kp2, 
                           to=trail2, connect='unique.pairs', free=T,
                           labels=c('beta_PC1','beta_PC2')), 
                    
                    mxPath(from=kp3, 
                           to=trail3, connect='unique.pairs', free=T,
                           labels=c('beta_PC1','beta_PC2'))
                    
                    
)

# container
bucketCheck0 <- mxModel("base multigroup check",
                        checkMD0, checkHC0,
                        mxFitFunctionMultigroup(c("checkMD0", "checkHC0")))

# run it
resCheck0 <- mxRun(bucketCheck0)
summary(resCheck0)
# Ref model
refCheck0 <- mxRefModels(resCheck0, run=TRUE)


# Free regression paths

# Free PC1
checkMD0_1 <- mxModel(checkMD0,
                      name="checkMD0_1",
                      mxPath(from=kp2[1], 
                             to=trail2, free=T,
                             labels='beta_PC1_MD'),
                      
                      mxPath(from=kp3[1], 
                             to=trail3, free=T,
                             labels='beta_PC1_MD')
)

checkHC0_1 <- mxModel(checkHC0,
                      name='checkHC0_1',
                      mxPath(from=kp2[1], 
                             to=trail2, free=T,
                             labels='beta_PC1_HC'),
                      
                      mxPath(from=kp3[1], 
                             to=trail3, free=T,
                             labels='beta_PC1_HC')
)

freePC1_check <- mxModel("free PC1 reg check",
                         checkMD0_1, checkHC0_1,
                         mxFitFunctionMultigroup(c("checkMD0_1", "checkHC0_1")))




# Free PC2
checkMD0_2 <- mxModel(checkMD0,
                      name="checkMD0_2",
                      mxPath(from=kp2[2], 
                             to=trail2, free=T,
                             labels='beta_PC2_MD'),
                      
                      mxPath(from=kp3[2], 
                             to=trail3, free=T,
                             labels='beta_PC2_MD')
)

checkHC0_2 <- mxModel(checkHC0,
                      name='checkHC0_2',
                      mxPath(from=kp2[2], 
                             to=trail2, free=T,
                             labels='beta_PC2_HC'),
                      
                      mxPath(from=kp3[2], 
                             to=trail3, free=T,
                             labels='beta_PC2_HC')
)

freePC2_check <- mxModel("free PC2 reg check",
                         checkMD0_2, checkHC0_2,
                         mxFitFunctionMultigroup(c("checkMD0_2", "checkHC0_2")))




# run models

# PC1 
res_freePC1_check <- mxRun(freePC1_check)
summary(res_freePC1_check, refModels=refCheck0)

mxCompare(res_freePC1_check, resCheck0)



# PC2
res_freePC2_check <- mxRun(freePC2_check)
summary(res_freePC2_check, refModels=refCheck0)

mxCompare(res_freePC2_check, resCheck0)

