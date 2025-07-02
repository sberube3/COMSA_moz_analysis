library(tidyverse)
library(ggplot2)
library(viridis)
library(ggforce)
library(readxl)
library(gtsummary)
library(officer)
library(dplyr)
#install.packages("ggforce")
#install.packages("readxl")

#setwd ("C:/Users/acarcelan/Dropbox/COMSA - Serosurveillance/Data/Analysis")
setwd("~/Dropbox/COMSA - Serosurveillance/Data/Analysis")

## ST: The code below is to make a combined Figure 1-4
## Figure 1: Weighted provincial seroprevalence
## Figure 2: Weighted provincial seroprevalence by age group
## Figure 3: Weighted provincial seroprevalence by sex
## Figure 4: Weighted provincial seroprevalence by cluster geography (urban/rural)

########################################
########################################
#######################################
###REVISED (15 Mar 2024)###
########################################
#######################################
#######################################

#######OVERALL FOREST#################
all_data<- read.csv("overal_seroprev_results.csv")
all_data[1:10,]

all_data_pos<- all_data[,c("Antigen", "X.1")]
all_data_pos[1:10,]

all_data_pos_filtered<- all_data_pos%>%
  filter(Antigen%in%c("b","ll", "ul"))

antigen_names<- all_data_pos%>%
  filter(X.1=="positive")%>%
  select(Antigen)

varnames<- read_excel("variablenames.xlsx")
colnames(varnames)[2] ="Antigen"
colnames(varnames)[7] ="Antigen_new"
varnamesnew = subset(varnames, select=c(Antigen, Antigen_new))

colnames(all_data_pos_filtered)<- c("quantity","value")

all_data_pos_filtered<- all_data_pos_filtered%>%
  mutate(Antigen = rep(varnamesnew$Antigen_new,each=3))

all_data_pos_wide_CIs<- spread(all_data_pos_filtered,quantity, value)
all_data_pos_wide_CIs[all_data_pos_wide_CIs == ''] <- 0

all_data_pos_wide_CIs<- all_data_pos_wide_CIs%>%
  mutate(Antigen_Type= ifelse(Antigen%in% c("Measles virus (wMev)",	"Rubella virus (wRuv)",	"Tetanus (Tet tox)",	"Diphtheria (Dip tox)"),"VPD",
                              ifelse(Antigen%in% c("Dengue (dengns1-4)",	"Dengue (dengns1-3)",	"Dengue (dengns1-2)",	"Dengue (dengns1-1)",	"Chikungunya (chike1)"),"arbovirus",
                                     ifelse(Antigen %in% c("SARS-CoV-2 (sars2rbd)",	"SARS-CoV-2 (sars2np)"),"COVID",
                                            ifelse(Antigen %in% c("Giardia lamblia (vsp5)",	"Giardia lamblia (vsp3)",	"Cryptosporidium parvum (cp23)",	"Cryptosporidium parvum (cp17)"),"enteric",
                                                   ifelse(Antigen%in% c("P. vivax (pvrbp2b)",	"P. vivax (pvmsp119)",	"P. vivax (pvdbprii)",	"P. ovale (pomsp119)",	"P. malariae (pmmsp119)",	"P. falciparum (pfmsp119)",	"P. falciparum (pfama1)",	"P. falciparum (glurpr2)",	"P. falciparum (gexp18)",	"P. falciparum (etramp5ag1)",	"P. falciparum (csp)", "P. falciparum (rh42)"),"malaria", "NTD"))))))

all_data_pos_wide_CIs<-all_data_pos_wide_CIs%>%
  filter(Antigen_Type!="arbovirus")

all_data_pos_wide_CIs$LCI<-as.numeric(all_data_pos_wide_CIs$ll)
all_data_pos_wide_CIs$UCI<-as.numeric(all_data_pos_wide_CIs$ul)
all_data_pos_wide_CIs$proportion<-as.numeric(all_data_pos_wide_CIs$b)

ggplot(data=all_data_pos_wide_CIs) +
  geom_pointrange(aes(y=Antigen, x=proportion, xmin=LCI, xmax=UCI), size=0.5,stroke = 0.5,position=position_dodge(width = 0.5))+
  scale_x_continuous(limits = c(0, 1),breaks=seq(0,1,b=.1)) + 
  ylab("Antigen") + xlab("Seroprevalence") +
  facet_grid(vars(Antigen_Type), scales="free", space="free") + 
  theme_bw() + 
  theme(strip.text.y = element_text(angle=0))+
  theme(strip.background = element_rect(fill = "white", color = "black"))+
  theme(plot.title = element_text(size=12,color="black", hjust = 0.0))+
  theme(legend.position="none")+
  coord_cartesian(xlim = c(0, 1)) -> plot_overall

#ggsave("Seroprevalence_overall_15Mar2024.png",height = 12,width = 10)
#ggsave("Seroprevalence_overall_15Mar2024.pdf",height = 9,width = 10)

##############
#Age Forest #
#############
#age_data<- read.csv("testage.csv")
age_data<- read.csv("Age_results.csv")
age_data[1:10,]

age_data_pos<- age_data[,c("Antigen", "X.3","X.4","X.5")]
age_data_pos[1:10,]

age_data_pos_filtered<- age_data_pos%>%
  filter(Antigen%in%c("b","ll", "ul"))

colnames(age_data_pos_filtered)<- c("quantity","0-4years", "5-17years", "18+years")

antigen_names<- age_data_pos%>%
  filter(X.3=="positive")%>%
  dplyr::select(Antigen)

varnames<- read_excel("variablenames.xlsx")
colnames(varnames)[2] ="Antigen"
colnames(varnames)[7] ="Antigen_new"
varnamesnew = subset(varnames, select=c(Antigen, Antigen_new))

age_data_pos_filtered<- age_data_pos_filtered%>%
  mutate(Antigen = rep(varnamesnew$Antigen_new,each=3))

age_data_pos_long_age<- gather(age_data_pos_filtered, Age_Group, value, 2:4)

age_data_pos_wide_CIs<- spread(age_data_pos_long_age,quantity, value)

age_data_pos_wide_CIs<- age_data_pos_wide_CIs%>%
  mutate(Antigen_Type= ifelse(Antigen%in% c("Measles virus (wMev)",	"Rubella virus (wRuv)",	"Tetanus (Tet tox)",	"Diphtheria (Dip tox)"),"VPD",
                              ifelse(Antigen%in% c("Dengue (dengns1-4)",	"Dengue (dengns1-3)",	"Dengue (dengns1-2)",	"Dengue (dengns1-1)",	"Chikungunya (chike1)"),"arbovirus",
                                     ifelse(Antigen %in% c("SARS-CoV-2 (sars2rbd)",	"SARS-CoV-2 (sars2np)"),"COVID",
                                            ifelse(Antigen %in% c("Giardia lamblia (vsp5)",	"Giardia lamblia (vsp3)",	"Cryptosporidium parvum (cp23)",	"Cryptosporidium parvum (cp17)"),"enteric",
                                                   ifelse(Antigen%in% c("P. vivax (pvrbp2b)",	"P. vivax (pvmsp119)",	"P. vivax (pvdbprii)",	"P. ovale (pomsp119)",	"P. malariae (pmmsp119)",	"P. falciparum (pfmsp119)",	"P. falciparum (pfama1)",	"P. falciparum (glurpr2)",	"P. falciparum (gexp18)",	"P. falciparum (etramp5ag1)",	"P. falciparum (csp)", "P. falciparum (rh42)"),"malaria", "NTD"))))))
age_data_pos_wide_CIs[age_data_pos_wide_CIs == ''] <- 0
age_data_pos_wide_CIs<-age_data_pos_wide_CIs%>%
  filter(Antigen_Type!="arbovirus")

age_data_pos_wide_CIs$Age_Group<- as.factor(age_data_pos_wide_CIs$Age_Group)

age_data_pos_wide_CIs$LCI<-as.numeric(age_data_pos_wide_CIs$ll)
age_data_pos_wide_CIs$UCI<-as.numeric(age_data_pos_wide_CIs$ul)
age_data_pos_wide_CIs$proportion<-as.numeric(age_data_pos_wide_CIs$b)
age_data_pos_wide_CIs$Age_Group<- factor(age_data_pos_wide_CIs$Age_Group, levels=c("0-4years", "5-17years", "18+years"), labels=c("0-4 years", "5-17 years", "18+ years"))


ggplot(data=age_data_pos_wide_CIs) +
  geom_pointrange(aes(y=Antigen, x=proportion, xmin=LCI, xmax=UCI, color=Age_Group), fill="white", size=.5,stroke = 0.5,position=position_dodge(width = 0.5))+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  scale_x_continuous(limits = c(0, 1),breaks=seq(0,1,b=.1)) + 
  ylab("") + xlab("Seroprevalence") +
  facet_grid(vars(Antigen_Type), scales="free", space="free") + 
  theme_bw() + 
  theme(strip.text.y = element_text(angle=0))+
  theme(strip.background = element_rect(fill = "white", color = "black"))+
  theme(plot.title = element_text(size=12,color="black", hjust = 0.0))+
  coord_cartesian(xlim = c(0, 1)) -> plot_age


age_data_pos_wide_CIs_talk_DF<- age_data_pos_wide_CIs|>
  filter(Antigen%in%c("SARS-CoV-2 (sars2rbd)",
                      "Cryptosporidium parvum (cp23)",
                      "Cryptosporidium parvum (cp17)",
                      "P. falciparum (glurpr2)",
                      "P. falciparum (pfama1)",
                      "P. falciparum (pfmsp119)",
                      "Chlamydia trachomatis (ct694)",
                      "Chlamydia trachomatis (pgp3)",
                      "Onchocerca volvulus (ov16)",
                      "Treponema palladium (rp17)",
                      "Treponema palladium (tmpa)",
                      "Measles virus (wMev)"))|>
  rename(Age= Age_Group)
ggplot(data=age_data_pos_wide_CIs_talk_DF) +
  geom_pointrange(aes(y=Antigen, x=proportion, xmin=LCI, xmax=UCI, color=Age), fill="white", size=0.5,stroke = 0.5,position=position_dodge(width = 0.5))+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  scale_x_continuous(limits = c(0, 1),breaks=seq(0,1,b=.1)) + 
  ylab("") + xlab("Seroprevalence") +
  facet_grid(vars(Antigen_Type), scales="free", space="free") + 
  theme_bw() + 
  theme(strip.text.y = element_text(angle=0))+
  theme(strip.background = element_rect(fill = "white", color = "black"))+
  theme(plot.title = element_text(size=12,color="black", hjust = 0.0))+
  coord_cartesian(xlim = c(0, 1))+
  theme(axis.text.x= element_text(size=18),
        axis.text.y= element_text(size=18),
        axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        legend.title= element_text(size=20), 
        legend.text= element_text(size=18))-> plot_age

#ggsave("Seroprevalence_byAge_15Mar2024.png",height = 12,width = 10)
#ggsave("Seroprevalence_byAge_15Mar2024.pdf",height = 9,width = 10)

###################
# Cluster Forest (urban/rural)#
##################
clust_data<- read.csv("Rural_urban_results.csv")
clust_data[1:10,]

clust_data_pos<- clust_data[,c("Antigen","X.2", "X.3")]
clust_data_pos[1:10,]

clust_data_pos_filtered<- clust_data_pos%>%
  filter(Antigen%in%c("b","ll", "ul"))

colnames(clust_data_pos_filtered)<- c("quantity","urban", "rural")

antigen_names<- clust_data_pos%>%
  filter(X.3=="positive")%>%
  dplyr::select(Antigen)

varnames<- read_excel("variablenames.xlsx")
colnames(varnames)[2] ="Antigen"
colnames(varnames)[7] ="Antigen_new"
varnamesnew = subset(varnames, select=c(Antigen, Antigen_new))

clust_data_pos_filtered<- clust_data_pos_filtered%>%
  mutate(Antigen = rep(varnamesnew$Antigen_new,each=3))

clust_data_pos_long_age<- gather(clust_data_pos_filtered,Cluster, value, 2:3)

clust_data_pos_wide_CIs<- spread(clust_data_pos_long_age,quantity, value)
clust_data_pos_wide_CIs<- clust_data_pos_wide_CIs%>%
  mutate(Antigen_Type= ifelse(Antigen%in% c("Measles virus (wMev)",	"Rubella virus (wRuv)",	"Tetanus (Tet tox)",	"Diphtheria (Dip tox)"),"VPD",
                              ifelse(Antigen%in% c("Dengue (dengns1-4)",	"Dengue (dengns1-3)",	"Dengue (dengns1-2)",	"Dengue (dengns1-1)",	"Chikungunya (chike1)"),"arbovirus",
                                     ifelse(Antigen %in% c("SARS-CoV-2 (sars2rbd)",	"SARS-CoV-2 (sars2np)"),"COVID",
                                            ifelse(Antigen %in% c("Giardia lamblia (vsp5)",	"Giardia lamblia (vsp3)",	"Cryptosporidium parvum (cp23)",	"Cryptosporidium parvum (cp17)"),"enteric",
                                                   ifelse(Antigen%in% c("P. vivax (pvrbp2b)",	"P. vivax (pvmsp119)",	"P. vivax (pvdbprii)",	"P. ovale (pomsp119)",	"P. malariae (pmmsp119)",	"P. falciparum (pfmsp119)",	"P. falciparum (pfama1)",	"P. falciparum (glurpr2)",	"P. falciparum (gexp18)",	"P. falciparum (etramp5ag1)",	"P. falciparum (csp)", "P. falciparum (rh42)"),"malaria", "NTD"))))))
clust_data_pos_wide_CIs[clust_data_pos_wide_CIs == ''] <- 0
clust_data_pos_wide_CIs<-clust_data_pos_wide_CIs%>%
  filter(Antigen_Type!="arbovirus")

clust_data_pos_wide_CIs$Sex<- as.factor(clust_data_pos_wide_CIs$Cluster)

clust_data_pos_wide_CIs$LCI<-as.numeric(clust_data_pos_wide_CIs$ll)
clust_data_pos_wide_CIs$UCI<-as.numeric(clust_data_pos_wide_CIs$ul)
clust_data_pos_wide_CIs$proportion<-as.numeric(clust_data_pos_wide_CIs$b)

ggplot(data=clust_data_pos_wide_CIs) +
  geom_pointrange(aes(y=Antigen, x=proportion, xmin=LCI, xmax=UCI, color=Sex), fill="white", size=.5,stroke = 0.5,position=position_dodge(width = 0.5))+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette="Set1")+
  scale_x_continuous(limits = c(0, 1),breaks=seq(0,1,b=.1)) + 
  ylab("Antigen") + xlab("Seroprevalence") + 
  facet_grid(vars(Antigen_Type), scales="free", space="free") + 
  theme_bw() + 
  theme(strip.text.y = element_text(angle=0))+
  theme(strip.background = element_rect(fill = "white", color = "black"))+
  theme(plot.title = element_text(size=12,color="black", hjust = 0.0))+
  labs(color="Cluster geography")+
  coord_cartesian(xlim = c(0, 1)) -> plot_ruralUrban

#ggsave("Seroprevalence_byrural_15Mar2024.png",height = 12,width = 10)
#ggsave("Seroprevalence_byrural_15Mar2024.pdf", height=9, width=10)

##############
#Sex Forest #
#############
#sex_data<- read.csv("testsex.csv")
sex_data<- read.csv("Sex_results.csv")
sex_data[1:10,]

sex_data_pos<- sex_data[,c("Antigen","X.2", "X.3")]
sex_data_pos[1:10,]

sex_data_pos_filtered<- sex_data_pos%>%
  filter(Antigen%in%c("b","ll", "ul"))

colnames(sex_data_pos_filtered)<- c("quantity","male", "female")

antigen_names<- sex_data_pos%>%
  filter(X.3=="positive")%>%
  select(Antigen)

varnames<- read_excel("variablenames.xlsx")
colnames(varnames)[2] ="Antigen"
colnames(varnames)[7] ="Antigen_new"
varnamesnew = subset(varnames, select=c(Antigen, Antigen_new))

sex_data_pos_filtered<- sex_data_pos_filtered%>%
  mutate(Antigen = rep(varnamesnew$Antigen_new,each=3))

sex_data_pos_long_age<- gather(sex_data_pos_filtered,Sex, value, 2:3)

sex_data_pos_wide_CIs<- spread(sex_data_pos_long_age,quantity, value)
sex_data_pos_wide_CIs<- sex_data_pos_wide_CIs%>%
  mutate(Antigen_Type= ifelse(Antigen%in% c("Measles virus (wMev)",	"Rubella virus (wRuv)",	"Tetanus (Tet tox)",	"Diphtheria (Dip tox)"),"VPD",
                              ifelse(Antigen%in% c("Dengue (dengns1-4)",	"Dengue (dengns1-3)",	"Dengue (dengns1-2)",	"Dengue (dengns1-1)",	"Chikungunya (chike1)"),"arbovirus",
                                     ifelse(Antigen %in% c("SARS-CoV-2 (sars2rbd)",	"SARS-CoV-2 (sars2np)"),"COVID",
                                            ifelse(Antigen %in% c("Giardia lamblia (vsp5)",	"Giardia lamblia (vsp3)",	"Cryptosporidium parvum (cp23)",	"Cryptosporidium parvum (cp17)"),"enteric",
                                                   ifelse(Antigen%in% c("P. vivax (pvrbp2b)",	"P. vivax (pvmsp119)",	"P. vivax (pvdbprii)",	"P. ovale (pomsp119)",	"P. malariae (pmmsp119)",	"P. falciparum (pfmsp119)",	"P. falciparum (pfama1)",	"P. falciparum (glurpr2)",	"P. falciparum (gexp18)",	"P. falciparum (etramp5ag1)",	"P. falciparum (csp)", "P. falciparum (rh42)"),"malaria", "NTD"))))))

sex_data_pos_wide_CIs<-sex_data_pos_wide_CIs%>%
  filter(Antigen_Type!="arbovirus")
sex_data_pos_wide_CIs[sex_data_pos_wide_CIs == ''] <- 0

sex_data_pos_wide_CIs$Sex<- as.factor(sex_data_pos_wide_CIs$Sex)

sex_data_pos_wide_CIs$LCI<-as.numeric(sex_data_pos_wide_CIs$ll)
sex_data_pos_wide_CIs$UCI<-as.numeric(sex_data_pos_wide_CIs$ul)
sex_data_pos_wide_CIs$proportion<-as.numeric(sex_data_pos_wide_CIs$b)


ggplot(data=sex_data_pos_wide_CIs) +
  geom_pointrange(aes(y=Antigen, x=proportion, xmin=LCI, xmax=UCI, color=Sex), fill="white", size=.5,stroke = 0.5,position=position_dodge(width = 0.5))+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette="Set1")+
  scale_x_continuous(limits = c(0, 1),breaks=seq(0,1,b=.1)) + 
  ylab("Antigen") + xlab("Seroprevalence") +
  facet_grid(vars(Antigen_Type), scales="free", space="free") + 
  theme_bw() + 
  theme(strip.text.y = element_text(angle=0))+
  theme(strip.background = element_rect(fill = "white", color = "black"))+
  theme(plot.title = element_text(size=12,color="black", hjust = 0.0))+
  coord_cartesian(xlim = c(0, 1)) -> plot_sex

#ggsave("Seroprevalence_bysex_15Mar2024.png",height = 12,width = 10)
#ggsave("Seroprevalence_bysex_15Mar2024.pdf",height = 9,width = 10)

## Combine panels
library(patchwork)

## Vanity things
plot_overall2 <- plot_overall + 
  theme(strip.text.y=element_blank()) + 
  ggtitle("A. Overall") + 
  scale_x_continuous(limits = c(0, 1),breaks=seq(0,1,b=.2))

plot_age2 <- plot_age + 
  theme(legend.position="bottom", axis.text.y=element_blank(), strip.text.y=element_blank(), legend.title=element_blank(), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) + 
  ggtitle("B. Age Group") + 
  scale_x_continuous(limits = c(0, 1),breaks=seq(0,1,b=.2))

plot_sex2 <- plot_sex + theme(legend.position="bottom", axis.text.y=element_blank(), strip.text.y=element_blank(), legend.title=element_blank(), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) + ylab("") + 
  ggtitle("C. Sex") + 
  scale_x_continuous(limits = c(0, 1),breaks=seq(0,1,b=.2))

plot_ruralUrban2 <- plot_ruralUrban + theme(legend.position="bottom", axis.text.y=element_blank(), legend.title=element_blank(), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) + ylab("") + 
  ggtitle("D. Cluster Geography") + 
  scale_x_continuous(limits = c(0, 1),breaks=seq(0,1,b=.2))

## Plot
#pdf("~/Dropbox/COMSA - Serosurveillance/10_Manuscripts/Main paper/Figure1_combinedPanel.pdf", width=18, height=10)
plot_overall2 + plot_age2 + plot_sex2 + plot_ruralUrban2 + plot_layout(ncol=4)
#dev.off()
