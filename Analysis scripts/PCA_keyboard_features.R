

# PCA analysis of keyboard features
# Author: Alexander Demos

library(tidyverse)
library(conflicted)
library(psych)
library(nFactors)
library(missMDA)
library(FactoMineR)
library(factoextra)
library(ggrepel)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::lead)


#######################
### PCA ###############
#######################

conflicts_prefer(nFactors::parallel)

Visit.Level.L<-read.csv("Visit_Level_AAAAX_6-10-24_KitchenSink_Long_NoCut.csv")

###### BY Window

Visit.Level.L1<-Visit.Level.L %>% filter(Study_Window2=="Visit 1-2")
Visit.Level.Lf<-Visit.Level.L1 %>% select(-c(X:Study_Window2, mean_energy:focus_RMSSD))

PCA.Data.clean <- Visit.Level.Lf %>%
  select_if(~ !anyNA(.))

logtrans<-function(x){
  z = log(x+1)
}

PCA.Data.clean<-PCA.Data.clean %>% 
  mutate(across(c(n_aaaaa_total:Mean_win,MAD_win,median_n_aaaaa_day:Mean_day,mean_MAD_day,sessions_n_win,
                  med_sessions_n_day,mad_sessions_n_day,med_All_KP_day:mad_Jerk_perc_day), logtrans))

PCA.Data.clean<-PCA.Data.clean %>% rename(Acc.C.D.1=med_Jerk_perc_day,
                                          Acc.V.D.2=mad_Jerk_perc_day,
                                          Acc.V.W.3=RMS_Jerk_perc_session,
                                          Speed.C.D.1=median_Q25_day,
                                          Speed.C.D.2=median_Q95_day,
                                          Speed.C.D.3=median_Q75_day,
                                          Speed.C.D.4=median_Q50_day,
                                          Speed.C.D.5=Mean_day,
                                          Speed.C.W.6=Q25_win,
                                          Speed.C.W.7=Q75_win,
                                          Speed.C.W.8=Q95_win,
                                          Speed.C.W.9=Mean_win,
                                          Speed.C.W.10=Q50_win,
                                          Speed.V.D.11=MAD_RMS_IKD_day,
                                          Speed.V.D.12=mean_MAD_day,
                                          Speed.V.D.13=mad_Q25_day,
                                          Speed.V.D.14=mad_Q95_day,
                                          Speed.V.D.15=mad_Q50_day,
                                          Speed.V.D.16=mad_Q75_day,
                                          Speed.V.D.17=median_RMS_IKD_day,
                                          Speed.V.W.18=RMSSD_IKD_win,
                                          Speed.V.W.19=MAD_win,
                                          Time.C.W.1=hour_mean_direct_win,
                                          Time.C.W.2=hour_med_direct_win,
                                          Time.V.D.3=rmssd_hour_mean_direct_day,
                                          Time.V.D.4=rmssd_hour_med_direct_day,
                                          Time.V.D.5=rmssd_hour_r_length_day,
                                          Time.V.W.6=hour_r_length_win,
                                          Usage.C.D.1=med_sessions_n_day,
                                          Usage.C.D.2=med_All_KP_day,
                                          Usage.C.D.3=median_n_aaaaa_day,
                                          Usage.C.W.4=sessions_n_win,
                                          Usage.C.W.5=n_aaaaa_total,
                                          Usage.C.W.6=days,
                                          Usage.V.D.7=mad_sessions_n_day,
                                          Usage.V.D.8=rmssd_sessions_n_day,
                                          Usage.V.D.9=mad_n_aaaaa_day,
                                          Usage.V.D.10=mad_All_KP_day,
                                          Usage.V.W.11=RMS_Total_KP)


CM<-cor(PCA.Data.clean, method = "pearson", use="pairwise.complete")

ev <- eigen(CM) # get eigenvalues
ap <- parallel(subject=nrow(PCA.Data.clean),var=ncol(PCA.Data.clean),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

# Run PCA
PCA.Fit<-principal(CM,nfactors = 4, residuals = FALSE, rotate="none")
Loadings <- factor.scores(PCA.Data.clean,PCA.Fit, method = "Thurstone")

PC1<-Loadings$scores[,1]
PC2<-Loadings$scores[,3]
PC3<-Loadings$scores[,2]
PC4<-Loadings$scores[,4]


#Cor
Visit.Level.Keep<-Visit.Level.L1 %>% select(c(healthCode:Study_Window2, mean_energy:focus_RMSSD))
Merged.1<-cbind(Visit.Level.Keep,PC1,PC2,PC3,PC4)

############################################################
##### Figure 1 #############################################
############################################################

result <- PCA(PCA.Data.clean,scale.unit = TRUE, ncp = 4)

### Show All variables
p <- fviz_pca_var(result, ggtheme = theme_bw(), repel = TRUE)

p2 <- fviz_pca_var(result,
                   ggtheme = theme_bw(), , repel = FALSE, label = "none") +
  ggtitle("") +
  theme(
    axis.title.x = element_text(face = "bold", hjust = 1),  # x-axis label bold and right-justified
    axis.title.y = element_text(face = "bold", hjust = 1),  # y-axis label bold and right-justified
    text = element_text(size = 8)  # Adjust text size as needed
  ) # Thinner vector lines

# Extract the data from the PCA plot for ggrepel
p_data <- ggplot_build(p)$data[[1]]
p2<-p2 + 
  geom_text_repel(data = p_data, aes(x = x, y = y, label = label), size = 3, 
                  force = 10, box.padding = 0.5,max.overlaps =20)
p2

ggsave(filename = "PCA1.png", plot=p2, dpi=1200, width =3,height = 3, scale=2)




###########################################################################
##################################################### Visit 2 #############
###########################################################################

Visit.Level.L2<-Visit.Level.L %>% filter(Study_Window2=="Visit 2-3")
Visit.Level.L2f<-Visit.Level.L2 %>% select(-c(X:Study_Window2, mean_energy:focus_RMSSD))

PCA.Data.clean.2 <- Visit.Level.L2f %>%
  select_if(~ !anyNA(.))


PCA.Data.clean.2<-PCA.Data.clean.2 %>% 
  mutate(across(c(n_aaaaa_total:Mean_win,MAD_win,median_n_aaaaa_day:Mean_day,mean_MAD_day,sessions_n_win,
                  med_sessions_n_day,mad_sessions_n_day,med_All_KP_day:mad_Jerk_perc_day), logtrans))


PCA.Data.clean.2<-PCA.Data.clean.2 %>% rename(Acc.C.D.1=med_Jerk_perc_day,
                                          Acc.V.D.2=mad_Jerk_perc_day,
                                          Acc.V.W.3=RMS_Jerk_perc_session,
                                          Speed.C.D.1=median_Q25_day,
                                          Speed.C.D.2=median_Q95_day,
                                          Speed.C.D.3=median_Q75_day,
                                          Speed.C.D.4=median_Q50_day,
                                          Speed.C.D.5=Mean_day,
                                          Speed.C.W.6=Q25_win,
                                          Speed.C.W.7=Q75_win,
                                          Speed.C.W.8=Q95_win,
                                          Speed.C.W.9=Mean_win,
                                          Speed.C.W.10=Q50_win,
                                          Speed.V.D.11=MAD_RMS_IKD_day,
                                          Speed.V.D.12=mean_MAD_day,
                                          Speed.V.D.13=mad_Q25_day,
                                          Speed.V.D.14=mad_Q95_day,
                                          Speed.V.D.15=mad_Q50_day,
                                          Speed.V.D.16=mad_Q75_day,
                                          Speed.V.D.17=median_RMS_IKD_day,
                                          Speed.V.W.18=RMSSD_IKD_win,
                                          Speed.V.W.19=MAD_win,
                                          Time.C.W.1=hour_mean_direct_win,
                                          Time.C.W.2=hour_med_direct_win,
                                          Time.V.D.3=rmssd_hour_mean_direct_day,
                                          Time.V.D.4=rmssd_hour_med_direct_day,
                                          Time.V.D.5=rmssd_hour_r_length_day,
                                          Time.V.W.6=hour_r_length_win,
                                          Usage.C.D.1=med_sessions_n_day,
                                          Usage.C.D.2=med_All_KP_day,
                                          Usage.C.D.3=median_n_aaaaa_day,
                                          Usage.C.W.4=sessions_n_win,
                                          Usage.C.W.5=n_aaaaa_total,
                                          Usage.C.W.6=days,
                                          Usage.V.D.7=mad_sessions_n_day,
                                          Usage.V.D.8=rmssd_sessions_n_day,
                                          Usage.V.D.9=mad_n_aaaaa_day,
                                          Usage.V.D.10=mad_All_KP_day,
                                          Usage.V.W.11=RMS_Total_KP)


CM2<-cor(PCA.Data.clean.2, method = "pearson", use="pairwise.complete")

ev2 <- eigen(CM2) # get eigenvalues
ap2 <- parallel(subject=nrow(PCA.Data.clean.2),var=ncol(PCA.Data.clean.2),
               rep=100,cent=.05)
nS2 <- nScree(x=ev2$values, aparallel=ap2$eigen$qevpea)
plotnScree(nS2)

# PCA
PCA.Fit2<-principal(CM,nfactors = 4, residuals = FALSE, rotate="none")

Loadings2 <- factor.scores(PCA.Data.clean.2,PCA.Fit2, method = "Thurstone")

PC1a<-Loadings2$scores[,1]
PC2a<-Loadings2$scores[,3]
PC3a<-Loadings2$scores[,2]
PC4a<-Loadings2$scores[,4]

####################################################
#### Make PCA figure 2 for sup methods
result2 <- PCA(PCA.Data.clean.2,scale.unit = TRUE, ncp = 4)

p.2 <- fviz_pca_var(result2, ggtheme = theme_bw(), repel = TRUE)

p2a <- fviz_pca_var(result2,
                   ggtheme = theme_bw(), , repel = FALSE, label = "none") +
  ggtitle("") +
  theme(
    axis.title.x = element_text(face = "bold", hjust = 1),  # x-axis label bold and right-justified
    axis.title.y = element_text(face = "bold", hjust = 1),  # y-axis label bold and right-justified
    text = element_text(size = 8)  # Adjust text size as needed
  ) # Thinner vector lines


# Extract the data from the PCA plot for ggrepel
p_data.2 <- ggplot_build(p.2)$data[[1]]

# Add ggrepel to avoid label overlap
p2a<-p2a + 
  geom_text_repel(data = p_data.2, aes(x = x, y = y, label = label), size = 3, 
                  force = 10, box.padding = 0.5,max.overlaps =20)
p2a

ggsave(filename = "PCA2.png", plot=p2, dpi=1200, width =3,height = 3, scale=2)


##################################################################
#Save data
Visit.Level.Keep.2<-Visit.Level.L2 %>% select(c(healthCode:Study_Window2, mean_energy:focus_RMSSD))
Merged.2<-cbind(Visit.Level.Keep.2,PC1a,PC2a,PC3a,PC4a)
Merged<-rbind(Merged.1, Merged.2)

write.csv(Merged, "PCA_NoCut_6-12-24_Logged_NewNames.csv")
