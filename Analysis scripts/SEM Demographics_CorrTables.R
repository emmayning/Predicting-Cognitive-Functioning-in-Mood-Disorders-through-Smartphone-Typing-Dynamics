
# Getting correlation tables for neuropsych tests
# Author: Alexander Demos


library(tidyverse)
library(conflicted)
library(table1)
conflicts_prefer(dplyr::recode)

Demographics <- read_csv("Demographics.csv")

# PCA and Biaffect saved Long
Visit.Level<-read.csv("PCA_NoCut_6-12-24_Logged_NewNames.csv")
Visit.Level.PCA.Long<-Visit.Level %>% select(healthCode,Study_Window2,PC1:PC4)
Visit.Level.PCA.Long$Visit<-factor(Visit.Level.PCA.Long$Study_Window2, 
                                   levels=c("Visit 1-2","Visit 2-3"),
                                   labels=c("V1_2","V2_3"))
Visit.Level.PCA.Wide<-Visit.Level.PCA.Long %>% select(-Study_Window2) %>% 
  pivot_wider(names_from=Visit, 
              values_from=PC1:PC4)

###########################################
Demographics$visit1_diag_dr<-ifelse(is.na(Demographics$visit1_diag_dr)==T, 0,Demographics$visit1_diag_dr)

Demographics$Diagnosis<-factor(Demographics$visit1_diag_dr, 
                               levels=c(0,1,2,3,4,5,6),
                               labels=c("Healthy Control", "Manic episode","Bipolar disorder","MDD single","MDD recurrent","Persistent mood disorder","Unspecified mood disorder"))
table(Demographics$Diagnosis)
Demographics$Diagnosis<-factor(Demographics$visit1_diag_dr, 
                               levels=c(0,4,2,1,3,5,6),
                               labels=c("Healthy Control","MDD recurrent","Bipolar disorder","Other mood disorder","Other mood disorder","Other mood disorder","Other mood disorder"))

table(Demographics$Diagnosis)

Demographics.Select<-Demographics %>% select(healthCode, Diagnosis, age:education,hispanic, race, trail_a_time_V2:trail_b_time_V3,
                                             NIH_T_Pattern_Procesing_V2:NIH_T_Flanker_V2,
                                             NIH_T_Pattern_Procesing_V3:NIH_T_Flanker_V3, V32_diff,V21_diff)

SEM.Demos.Wide<-left_join(Visit.Level.PCA.Wide,Demographics.Select)

SEM.Demos.Wide$gender<-factor(SEM.Demos.Wide$gender, 
                                 levels=c(1,2), 
                                 labels=c("Male","Female"))

SEM.Demos.Wide$education<-factor(SEM.Demos.Wide$education, 
         levels=c("HS or lower","AA/Some college","BA","MA","Professional/Doctoral"), 
         labels=c("HS or lower","AA/Some college","BA","MA","Professional/Doctoral"))


SEM.Demos.Wide$hispanic<-factor(SEM.Demos.Wide$hispanic, 
                            levels=c(1,2), 
                            labels=c("Hispanic/Latino","Not Hispanic/Latino"))


SEM.Demos.Wide$race<-factor(SEM.Demos.Wide$race, 
                                 levels=c(1,2,3,4,5,6,7), 
                                 labels=c("White","Black","Asian","Native Hawaiian/Pacific Islander","American Indian/Alaskan Native","Other", "Not reported"))

###################################################
### Tables ########################################
###################################################

table1(~  age + education + gender +hispanic+race+V32_diff+V21_diff| Diagnosis, data=SEM.Demos.Wide)

table1(~  V21_diff+V32_diff, data=SEM.Demos.Wide)

table1(~  NIH_T_Pattern_Procesing_V2+NIH_T_CardSort_V2 + NIH_T_ListSort_V2 + NIH_T_Flanker_V2| Diagnosis, data=SEM.Demos.Wide)


###################################################
### Boxplots ######################################
###################################################
SEM.Demos.Wide.NIH<-SEM.Demos.Wide %>% select(healthCode, Diagnosis,
                                              NIH_T_Pattern_Procesing_V2:NIH_T_Flanker_V2,
                                              NIH_T_Pattern_Procesing_V3:NIH_T_Flanker_V3)


SEM.Demos.Wide.NIH <- SEM.Demos.Wide.NIH %>%
  rename_with(~ gsub("NIH_T_", "", .), starts_with("NIH_T_"))  %>%
  rename_with(~ gsub("Pattern_Procesing", "Pattern", .))

SEM.Demos.Long.NIH <- SEM.Demos.Wide.NIH %>%
  pivot_longer(cols = Pattern_V2:Flanker_V3, 
               names_to = c(".value", "Visit"), 
               names_pattern = "(.*)_(V[23])") %>%
  mutate(Visit = recode(Visit, V2 = 2, V3 = 3)) %>% 
  pivot_longer(cols = Pattern:Flanker,
               names_to ="Measurement", 
               values_to = "TScore")

SEM.Demos.Long.NIH$Measurement<-factor(SEM.Demos.Long.NIH$Measurement, 
                                       levels=c("Pattern",  "CardSort", "ListSort", "Flanker" ),
                                       labels=c("Pattern",  "CardSort", "ListSort", "Flanker" ))


Box.NIH<-ggplot(data = SEM.Demos.Long.NIH, aes(x = as.factor(Visit), y=TScore), group=Diagnosis) +
  facet_grid(~Measurement)+
  coord_cartesian(ylim=c(0,100))+
  geom_violin(aes(fill=Measurement), trim = F,adjust = .75, alpha=.75)+
  geom_boxplot(aes(fill=Measurement),width=.2)+
  xlab("Visit")+theme_bw()+theme(legend.position = "none")
Box.NIH

ggsave(filename="NIHbox.png", plot=Box.NIH,dpi=1200, width=6.5, height=2.5)


###################################################
### ICC  ##########################################
###################################################
SEM.Demos.Long.NIH.ICC<-SEM.Demos.Long.NIH %>% pivot_wider(names_from=Visit, values_from = TScore)
psych::ICC(SEM.Demos.Long.NIH.ICC[,c(4,5)])

SEM.Demos.Wide.TMT<-SEM.Demos.Wide %>% select(healthCode, Diagnosis,  trail_b_time_V2,trail_b_time_V3)
psych::ICC(SEM.Demos.Wide.TMT[,c(3,4)])
###################################################
### Boxplots ######################################
###################################################
SEM.Demos.Long.TMT<- SEM.Demos.Wide.TMT %>%
  pivot_longer(cols = trail_b_time_V2:trail_b_time_V3, 
               names_to = c(".value", "Visit"), 
               names_pattern = "(.*)_(V[23])") %>%
  mutate(Visit = recode(Visit, V2 = 2, V3 = 3)) %>% 
  pivot_longer(cols = trail_b_time,
               names_to ="Measurement", 
               values_to = "Time")


Box.TMT<-ggplot(data = SEM.Demos.Long.TMT, aes(x = as.factor(Visit), y=log(Time)), group=Diagnosis) +
  #facet_grid(~Measurement)+
  #coord_cartesian(ylim=c(0,300))+
  geom_violin(aes(fill=Measurement), trim = F,adjust = .75, alpha=.75)+
  geom_boxplot(aes(fill=Measurement),width=.2)+
  xlab("Visit")+ylab("Log Time (Seconds)")+theme_bw()+theme(legend.position = "none")
Box.TMT

ggsave(filename="Trialsbox.png",plot=Box.TMT, dpi=1200, width=3, height=2.5)



###############################################################
### Correlation tables ########################################
###############################################################

# Cor and means table
library(apaTables)
SEM.Demos.Wide.NIH.TMT<-left_join(SEM.Demos.Wide.NIH, SEM.Demos.Wide.TMT) %>% left_join(Visit.Level.PCA.Wide)
SEM.Demos.Wide.NIH.TMT$trail_b_time_V2<-log(SEM.Demos.Wide.NIH.TMT$trail_b_time_V2)
SEM.Demos.Wide.NIH.TMT$trail_b_time_V3<-log(SEM.Demos.Wide.NIH.TMT$trail_b_time_V3)

SEM.Demos.Wide.NIH.TMT<- SEM.Demos.Wide.NIH.TMT %>% select(Pattern_V2:Flanker_V2, trail_b_time_V2,PC1_V1_2,PC2_V1_2,PC3_V1_2,PC4_V1_2,
                                                           Pattern_V3:Flanker_V3, trail_b_time_V3,PC1_V2_3,PC2_V2_3,PC3_V2_3,PC4_V2_3)

############################
apa.cor.table(
  SEM.Demos.Wide.NIH.TMT,
  filename = "Cor_final.doc",
  table.number = NA,
  show.conf.interval = FALSE,
  show.sig.stars = F,
  landscape = FALSE
)

