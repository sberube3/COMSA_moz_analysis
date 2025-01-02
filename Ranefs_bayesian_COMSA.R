library(tidyverse)
library(ggplot2)
library(lme4)
library(janitor)
library(readxl)
library(parameters)
library(brms)
library(cowplot)

setwd("/Users/sophieberube/Desktop/Hopkins/Hopkins/Bead_serology")

merged_data<- read.csv("merged_comsa_MFI_plus_metadata.csv")


dat_meta <- read_excel("~/Dropbox/COMSA - Serosurveillance/Data/Saki/MozCOMSA_all_26AUG2022_cutoffs.xlsx", sheet="cutoffs_ST") %>%
  clean_names() %>%
  mutate(ag_name_with_cutoff = paste0(ag_name, "_", primary_cutoff_value))

clean_antigen_names<- dat_meta$ag_name_with_cutoff[-c(24:26)]

seroprev_data<- merged_data%>%
  dplyr::select(c(clean_antigen_names, mage,cluster,mfi_gst_ls,mfi_gst_cdc,mfi_unlys))

seroprev_data$cluster<- as.factor(seroprev_data$cluster)

seroprev_data<- seroprev_data|>
  mutate(Urban= ifelse(cluster%in%c("4","64","25","6","24"),"Urban","Rural"))

adj_OR_urban<- data.frame()
for(i in 1:length(clean_antigen_names)){
  results_matrix<- matrix(NA, nrow=3000, ncol=32)
  data1<- seroprev_data
  
  colnames(data1)[i]<- "Y"
  
  brm_model<- brm(as.factor(Y) ~ as.numeric(mfi_gst_ls)+as.numeric(mfi_gst_cdc)+  as.numeric(mfi_unlys)+
                    as.factor(Urban)+
                    (1 | cluster), data = data1, family = bernoulli(link="logit"),
                  warmup=500, iter=1500, chains=3, inits="0")
  
  ranef_brm<- ranef(brm_model,summary=FALSE)
  ranef_brm_trim<- ranef_brm$cluster
  results_matrix[1:3000,1:30] <- ranef_brm_trim
  results_matrix[,31]<- colnames(seroprev_data)[i]
  results_matrix[1:3000, 32]<- "Adjusted NoArbo"
  
  sum_brm_model<- summary(brm_model)
  
  adj_OR_urban_df<- data.frame(Antigen= colnames(seroprev_data)[i],
                               OR_Urban= exp(sum_brm_model$fixed[5,1]),
                               lower95_Urban= exp(sum_brm_model$fixed[5,3]),
                               upper95_Urban= exp(sum_brm_model$fixed[5,4]))
  
  adj_OR_urban<- rbind(adj_OR_urban_df,adj_OR_urban)
  #brm_noCtrl<- brm(as.factor(Y) ~ as.numeric(mage)+
                     #(1 | cluster), data = data1, family = bernoulli(link="logit"),
                   #warmup=500, iter=1500, chains=3, inits="0")
  #ranef_brm_noCtrl<-ranef(brm_noCtrl, summary=FALSE)
  #ranef_brm_noCtrl_trim<- ranef_brm_noCtrl$cluster
  #results_matrix[3001:6000,1:30] <- ranef_brm_noCtrl_trim
  #results_matrix[3001:6000, 32]<- "Unadjusted"
  
  
  #results_df<- as.data.frame(results_matrix)
  #colnames(results_df)<- c(colnames(ranef_brm$cluster), "Antigen", "Model")
 
  print(i)
}



adj_OR_urban<- adj_OR_urban|>
  group_by(Antigen)|>
  summarize(OR= OR_Urban[1], Lower95= lower95_Urban[1], Upper95= upper95_Urban[1])
ag_full_names<- read.csv("variablenames.csv")
OR_results<- data.frame(Outcome_Antigen= character(), 
                        Regressor_Antigen= character(), 
                        Adjusted_OR= numeric(), 
                        p_val= numeric(), 
                        CI_lower=numeric(),
                        CI_upper=numeric())

seroprev_clean<- seroprev_data[, -c(1:5)]

seroprev_clean<- seroprev_clean|>
  mutate(cluster_type= ifelse(cluster%in% c(4,6,24,25,64),"Urban","Rural"))

for(i in 1:35){
  reg_df<- seroprev_clean
  
  colnames(reg_df)[i]<- "Y"
  
  reg_df$Y<- as.factor(reg_df$Y)
  
  reg_clean_df<- reg_df|>
    dplyr::select(-cluster)
  
  model1<- glm(Y~., data=reg_clean_df, family="binomial")
  sum_model1<- summary(model1)
  ci_model<- confint(model1)
  coef_sum_model1<- sum_model1$coefficients
  OR_df<- data.frame(Outcome_Antigen=rep(colnames(seroprev_clean)[i],34), 
                     Regressor_Antigen= colnames(seroprev_clean)[1:35][-i],
                     Adjusted_OR= exp(model1$coefficients[2:35]), 
                     p_val= coef_sum_model1[2:35,4],
                     CI_lower= exp(ci_model[2:35,1]),
                     CI_upper= exp(ci_model[2:35, 2]))
  
  OR_results<- rbind(OR_results, OR_df)
  
}

ag_nonVPD_names_outcome<- ag_full_names|>
  filter(name%in%c("measles", "rube", "tetanus_01", "dipth_01")==F)|>
  dplyr::select(name, varlab)|>
  mutate(Outcome_Antigen_Name= sapply(str_split(name, "_"), function(x) x[2]))|>
  rename(Outcome_Clean_Name=varlab)|>
  dplyr::select(Outcome_Antigen_Name, Outcome_Clean_Name)
  
ag_nonVPD_names_regressor<- ag_full_names|>
  filter(name%in%c("measles", "rube", "tetanus_01", "dipth_01")==F)|>
  dplyr::select(name, varlab)|>
  mutate(Regressor_Antigen_Name= sapply(str_split(name, "_"), function(x) x[2]))|>
  rename(Regressor_Clean_Name=varlab)|>
  dplyr::select(Regressor_Antigen_Name, Regressor_Clean_Name)
  
OR_results<- OR_results|>
  mutate(Antigen_Type_Outcome= ifelse(Outcome_Antigen%in% c("mfi_tettox_594", "mfi_diptox_241","mfi_wruv_255","mfi_wmev_215"),"VPD",
                              ifelse(Outcome_Antigen%in% c("mfi_dengns11_804","mfi_dengns12_347","mfi_dengns13_479","mfi_dengns14_706","mfi_chike1_585"),"arbovirus",
                                     ifelse(Outcome_Antigen %in% c("mfi_spiken_ls_644","mfi_spikerbd_ls_413"),"COVID",
                                            ifelse(Outcome_Antigen %in% c("mfi_cp17_99","mfi_cp23_800","mfi_vsp5_148","mfi_vsp3_99"),"enteric",
                                                   ifelse(Outcome_Antigen%in% c("mfi_csp_ls_1250","mfi_glurpr2_ls_212","mfi_pfama1_ls_1419","mfi_pfmsp119_ls_1557",
                                                                                "mfi_rh42_401","mfi_pmmsp119_ls_136","mfi_pomsp119_ls_144","mfi_pvdbprii_ls_133",
                                                                                "mfi_pvmsp119_ls_220","mfi_pvrbp2b_ls_151","mfi_etramp5ag1_ls_366",
                                                                                "mfi_gexp18_2_261"),"malaria", "NTD"))))))|>
  mutate(Antigen_Type_Regressor= ifelse(Regressor_Antigen%in% c("mfi_tettox_594", "mfi_diptox_241","mfi_wruv_255","mfi_wmev_215"),"VPD",
                                        ifelse(Regressor_Antigen%in% c("mfi_dengns11_804","mfi_dengns12_347","mfi_dengns13_479","mfi_dengns14_706","mfi_chike1_585"),"arbovirus",
                                               ifelse(Regressor_Antigen%in% c("mfi_spiken_ls_644","mfi_spikerbd_ls_413"),"COVID",
                                                      ifelse(Regressor_Antigen %in% c("mfi_cp17_99","mfi_cp23_800","mfi_vsp5_148","mfi_vsp3_99"),"enteric",
                                                             ifelse(Regressor_Antigen%in% c("mfi_csp_ls_1250","mfi_glurpr2_ls_212","mfi_pfama1_ls_1419",
                                                                                            "mfi_pfmsp119_ls_1557","mfi_rh42_401","mfi_pmmsp119_ls_136",
                                                                                            "mfi_pomsp119_ls_144","mfi_pvdbprii_ls_133","mfi_pvmsp119_ls_220",
                                                                                            "mfi_pvrbp2b_ls_151","mfi_etramp5ag1_ls_366",
                                                                                            "mfi_gexp18_2_261"),"malaria", "NTD"))))))|>
  mutate(Outcome_Antigen_Name= sapply(str_split(Outcome_Antigen,"_"), function(x) x[2]))|>
  left_join(ag_nonVPD_names_outcome)|>
  mutate(Regressor_Antigen_Name= sapply(str_split(Regressor_Antigen,"_"), function(x) x[2]))|>
  left_join(ag_nonVPD_names_regressor)

OR_results<- OR_results|>
  mutate(Outcome_AG_Full=ifelse(Outcome_Antigen_Name=="diptox", "Diptheria (Dip tox)",
                ifelse(Outcome_Antigen_Name=="wmev", "Measles (wMev)",
                       ifelse(Outcome_Antigen_Name=="wruv", "Rubella (wRuv)",
                              ifelse(Outcome_Antigen_Name=="tettox","Tetanus (Tet tox)",Outcome_Clean_Name)))),
         Regressor_AG_Full=ifelse(Regressor_Antigen_Name=="diptox", "Diptheria (Dip tox)",
                                 ifelse(Regressor_Antigen_Name=="wmev", "Measles (wMev)",
                                        ifelse(Regressor_Antigen_Name=="wruv", "Rubella (wRuv)",
                                               ifelse(Regressor_Antigen_Name=="tettox","Tetanus (Tet tox)",Regressor_Clean_Name)))))
  


OR_complete_results<- OR_results|>
  mutate(adj_p_val= p.adjust(p_val, "BH"))

write.csv(OR_complete_results, "Odds_Ratios_Results_Indiv_Analysis.csv")

P_val_df<- OR_results|>
  mutate(adj_p_holm=p.adjust(p_val,'bonferroni') ,
         adj_p_BH=p.adjust(p_val, "BH"),
         adj_p_fdr=p.adjust(p_val, "fdr"))|>
  filter(adj_p_BH<0.05)

library(viridis)
 OR_heatmap<- ggplot(OR_results, aes(x= Regressor_AG_Full, y= Outcome_AG_Full, fill= log(Adjusted_OR)))+
    geom_tile()+ geom_tile(data=P_val_df, aes(x=Regressor_AG_Full, y=Outcome_AG_Full), fill="transparent", color="red", size=0.5)+
   theme_bw()+scale_fill_continuous(type='viridis', breaks= c(-20,-10,0,10,20),
     labels= c("2 x 10^-9","4 x 10^-5", "1", "2 x 10^4", "4 x 10^8"))+xlab("Regressor Antigen")+
   facet_grid(Antigen_Type_Outcome~Antigen_Type_Regressor, scales='free', space='free')+
   ylab("Outcome Antigen")+theme(axis.text.x=element_text(angle=90, vjust = 0.6, face= "italic"),
                                 axis.text.y=element_text(face="italic"))+
   labs(fill="Adj OR")
 
 
 
   
  
cluster_dif_df<- data.frame()
antigens_to_check<- P_val_df[,c("Outcome_Antigen", "Regressor_Antigen", "Adjusted_OR")]

cluster_comparison_df<- data.frame()
for(i in 1:nrow(antigens_to_check)){
  analysis_df<- seroprev_data|>
    dplyr::select(antigens_to_check[i,1], antigens_to_check[i,2], cluster, mage, mfi_gst_ls, mfi_gst_cdc, mfi_unlys)
  antigens_to_check[i,]
  analysis_df[,1]<- as.factor(analysis_df[,1])
  model_check<- glmer(analysis_df[,1]~analysis_df[,2]+mage + mfi_gst_ls+mfi_gst_cdc+mfi_unlys+(1|cluster), 
                      family = 'binomial', data= analysis_df)
  
  model_check2<- glm(analysis_df[,1]~analysis_df[,2]+mage + mfi_gst_ls+mfi_gst_cdc+mfi_unlys, 
                     family = 'binomial', data= analysis_df)
  sum<- summary(model_check)
  sum2<- summary(model_check2)
  

  
  dif_df<- data.frame(Outcome_Antigen= antigens_to_check[i,1], 
                      Regressor_Antigen= antigens_to_check[i,2], 
                      Pct_change = (sum$coefficients[2,1]-sum2$coefficients[2,1])*100/sum$coefficients[2,1],
                      Multilevel_OR= sum$coefficients[2,1],
                      Flat_OR= sum2$coefficients[2,1],
                      Sig_multi= sum$coefficients[2,4], 
                      Sig_flat= sum2$coefficients[2,4])
  cluster_comparison_df<- rbind(cluster_comparison_df, dif_df)
}

cluster_list<- unique(merged_data$cluster)


bayes_results<- read.csv("ranef_bayes_draws.csv")

#unadjusted_bayes<- bayes_results%>%
  #filter(Model=="Unadjusted")


antigen_names_list<- str_split(bayes_results$Antigen, "_")
antigen_names_clean<- unlist(lapply(antigen_names_list, function(x) x[2]))


adjusted_bayes<- bayes_results%>%
mutate(Antigen_Type= ifelse(Antigen%in% c("mfi_tettox_594", "mfi_diptox_241","mfi_wruv_255","mfi_wmev_215"),"VPD",
                            ifelse(Antigen%in% c("mfi_dengns11_804","mfi_dengns13_479","mfi_dengns14_706","mfi_chike1_585"),"arbovirus",
                                   ifelse(Antigen %in% c("mfi_spiken_ls_644","mfi_spikerbd_ls_413"),"COVID",
                                          ifelse(Antigen %in% c("mfi_cp17_99","mfi_cp23_800","mfi_vsp5_148","mfi_vsp3_99"),"enteric",
                                                 ifelse(Antigen%in% c("mfi_csp_ls_1250","mfi_glurpr2_ls_212","mfi_pfama1_ls_1419","mfi_pfmsp119_ls_1557","mfi_rh42_401","mfi_pmmsp119_ls_136","mfi_pomsp119_ls_144","mfi_pvdbprii_ls_133","mfi_pvmsp119_ls_220"," mfi_pvrbp2b_ls_151"),"malaria", "NTD"))))))%>%
  mutate(Antigen_clean=antigen_names_clean)%>%
  mutate(Antigen_new_clean= ifelse(Antigen_clean=="spikerbd", "sars2rbd", ifelse(Antigen_clean=="spiken", "sars2np", 
                                                                                 ifelse(Antigen_clean=="tettox", "tetanus", 
                                                                                        ifelse(Antigen_clean=="diptox", "diptheria", 
                                                                                               ifelse(Antigen_clean=="wruv", "rubella", 
                                                                                                      ifelse(Antigen_clean=="wmev", "measles",Antigen_clean)))))))

ag_full_names<- read.csv("variablenames.csv")
district_cluster_names<- read.csv("clusters districts.csv")


ag_names_new<- ag_full_names%>%
  mutate(Antigen_split= strsplit(name, "_"))

ag_clean<- c(sapply(ag_names_new$Antigen_split[1:4], "[[", 1),
             sapply(ag_names_new$Antigen_split[-c(1:4)], "[[", 2))

ag_names_new<- ag_names_new%>%
  mutate(Ag_clean=ag_clean)%>%
  mutate(Antigen_clean= ifelse(Ag_clean== "rube", "wruv",
                               ifelse(Ag_clean=="dipth", "diptox", 
                                      ifelse(Ag_clean=="measles", "wmev",
                                             ifelse(Ag_clean=="tetanus", "tettox", Ag_clean)))))

adjusted_bayes_merge_antigens<- merge(adjusted_bayes, ag_names_new, by='Antigen_clean', all.x = T)



antigen_list<- unique(adjusted_bayes_merge_antigens$varlab)

######check this and redo the plot THEN SUM THEM

get_antigen_ranks<-function(i){
  adjusted_antigen<- adjusted_bayes_merge_antigens%>%
    filter(varlab==antigen_list[i])
    ranks_draws<- matrix(NA, nrow=3000, ncol=30)
  if(antigen_list[i]%in%c("tetanus", "diptheria","rubella","measles")){
    for(j in 1:nrow(adjusted_antigen)){
      ranks_draws[j,]<- rank(-as.numeric(adjusted_antigen[j,2:31]))
    }
  }
    else{
      for(j in 1:nrow(adjusted_antigen)){
        ranks_draws[j,]<- rank(as.numeric(adjusted_antigen[j,2:31]))
      }
    }
    ranks_draws_df<- as.data.frame(ranks_draws)
    colnames(ranks_draws_df)<- colnames(adjusted_antigen)[2:31]
    mean_ranks_df<- apply(ranks_draws_df,2,mean)
  
  mean_ranks_df<- as.data.frame(cbind(names(mean_ranks_df),mean_ranks_df,rep(antigen_list[i],30),rep("Adjusted",30)))
  colnames(mean_ranks_df)<- c("Cluster","Ranks","Antigen","Model")
  
  
  
  #unadjusted_antigen<- unadjusted_bayes%>%
    #filter(Antigen==antigen_list[i])
  
  #unadj_ranks_draws_df<- apply(unadjusted_antigen[,1:30],1, rank)
  #unadj_mean_ranks_df<- apply(unadj_ranks_draws_df,1,mean)
  
  #unadj_mean_ranks_df<- as.data.frame(cbind(names(unadj_mean_ranks_df),unadj_mean_ranks_df,rep(antigen_list[i],30),rep("Unadjusted",30)))
  #colnames(unadj_mean_ranks_df)<- c("Cluster","Ranks","Antigen","Model")
  
  return(mean_ranks_df)
  
}


ranks_results<- data.frame("Cluster"=character(),
                           "Ranks"=numeric(),
                           "Antigen"= character(), 
                           "Model"= character())  

ranks_antigen_index<-which(antigen_list%in%c("Dengue(dengns1-1)","Dengue(dengns1-2)","Dengue(dengns1-3)",
                         "Dengue(dengns1-4)","Chikungunya(chike1)",
                        "P. falciparum(gexp18)","P. falciparum(etramp5ag1)",
                        "P. vivax(pvrbp2b)")==F)

for(i in ranks_antigen_index){
  optim_results<-get_antigen_ranks(i)
  ranks_results<- rbind(ranks_results,optim_results)
}

adjusted_ranks_results<- ranks_results%>%
  mutate(Cluster_clean = substring(Cluster,2))%>%
  mutate(Cluster_Type= ifelse(Cluster_clean%in% c(4,6,24,25,64),"Urban","Rural"))%>%
  mutate(Antigen_Type= ifelse(Antigen%in% c("Tetanus (Tet tox)", "Diptheria (Dip tox)","Rubella (wRuv)","Measles (wMev)"),"VPD",
                              ifelse(Antigen%in% c("Dengue(dengns1-1)","Dengue(dengns1-3)","Dengue(dengns1-4)","Chikungunya(chike1)"),"arbovirus",
                                     ifelse(Antigen %in% c("SARS-CoV-2(sars2np)","SARS-CoV-2(sars2rbd)"),"COVID",
                                            ifelse(Antigen %in% c("Cryptosporidium parvum(cp17)","Cryptosporidium parvum(cp23)","Giardia lamblia(vsp3)","Giardia lamblia(vsp5)"),"enteric",
                                                   ifelse(Antigen%in% c("P. falciparum(csp)","P. falciparum(glurpr2)","P. falciparum(pfama1)","P. falciparum(pfmsp119)","P. falciparum(rh42)","P. malariae(pmmsp119)","P. ovale(pomsp119)","P. vivax(pvdbprii)","P. vivax(pvmsp119)","P. vivax(pvrbp2b)"),"malaria", "NTD"))))))
                                                                                                                                                                                                                                                                                                                                           

cluster_sums<- adjusted_ranks_results%>%
  group_by(Cluster_clean)%>%
  summarize(Cluster_Type= Cluster_Type[1], Rank_Sum=sum(as.numeric(Ranks)))

cluster_sums_df<- as.data.frame(cluster_sums)

cluster_sums_df$Rank_Sum<- as.numeric(cluster_sums_df$Rank_Sum)

cluster_sum_plot<- ggplot(cluster_sums_df, aes(x=Cluster_clean,y=Rank_Sum,color=Cluster_Type))+
  geom_point(size=2)+theme_bw()+scale_color_brewer(palette='Dark2')+xlab("Cluster")+ylab("Sum of Ranks")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+labs(color="Cluster Type")

cluster_sums_df<- cluster_sums_df%>%
  mutate(Model= "Adjusted", Antigen= "Overall Vulnerability", Antigen_Type= "Sum")%>%
  rename(Ranks=Rank_Sum)

adjusted_ranks_results<- adjusted_ranks_results%>%
  dplyr::select(Cluster_clean,Cluster_Type,Ranks,Model,Antigen,Antigen_Type)

adjusted_ranks_and_sums<- rbind(adjusted_ranks_results, cluster_sums_df)

colnames(district_cluster_names)<- c('district', 'Cluster_clean', 'Cluster_Name')

adjusted_ranks_and_sums_merged<- merge(adjusted_ranks_and_sums, district_cluster_names, by='Cluster_clean', all.x = T)
adjusted_ranks_rural<- adjusted_ranks_and_sums_merged%>%
  filter(Cluster_Type=="Rural")

adjusted_ranks_urban<- adjusted_ranks_and_sums_merged%>%
  filter(Cluster_Type=="Urban")
library(viridis)
library(ggh4x)

rural_arbo<-adjusted_ranks_rural%>%
  filter(Antigen_Type=="arbovirus")

rural_covid<-adjusted_ranks_rural%>%
  filter(Antigen_Type=="COVID")

rural_enteric<-adjusted_ranks_rural%>%
  filter(Antigen_Type=="enteric")

rural_malaria<-adjusted_ranks_rural%>%
  filter(Antigen_Type=="malaria")


rural_ntd<-adjusted_ranks_rural%>%
  filter(Antigen_Type=="NTD")

  rural_vpd<-adjusted_ranks_rural%>%
  filter(Antigen_Type=="VPD")
  
  rural_sum<- adjusted_ranks_rural%>%
    filter(Antigen_Type=="Sum")
  
rural_arbo$Cluster_Name<- factor(rural_arbo$Cluster_Name,
                                  levels = rural_sum$Cluster_Name[order(as.numeric(rural_sum$Ranks))])

rural_covid$Cluster_Name<- factor(rural_covid$Cluster_Name,
                                  levels = rural_sum$Cluster_Name[order(as.numeric(rural_sum$Ranks))])

rural_enteric$Cluster_Name<- factor(rural_enteric$Cluster_Name,
                                  levels = rural_sum$Cluster_Name[order(as.numeric(rural_sum$Ranks))])

rural_malaria$Cluster_Name<- factor(rural_malaria$Cluster_Name,
                                  levels = rural_sum$Cluster_Name[order(as.numeric(rural_sum$Ranks))])

rural_ntd$Cluster_Name<- factor(rural_ntd$Cluster_Name,
                                  levels = rural_sum$Cluster_Name[order(as.numeric(rural_sum$Ranks))])

rural_vpd$Cluster_Name<- factor(rural_vpd$Cluster_Name,
                                  levels = rural_sum$Cluster_Name[order(as.numeric(rural_sum$Ranks))])

rural_sum$Cluster_Name<- factor(rural_sum$Cluster_Name,
                                 levels = rural_sum$Cluster_Name[order(as.numeric(rural_sum$Ranks))])


heatmap_adjusted_rural_arbo<- ggplot(rural_arbo,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",strip.text.y = element_blank(),
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        strip.text.x=element_text(size=12), axis.text.y=element_text(size=12, face= "italics"))+ylab(NULL)+xlab(NULL)

heatmap_adjusted_rural_covid<- ggplot(rural_covid,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",strip.text.y = element_blank(),
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        strip.text.x=element_text(size=12), axis.text.y=element_text(size=12, face="italic"))+ylab(NULL)+xlab(NULL)

heatmap_adjusted_rural_enteric<- ggplot(rural_enteric,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",strip.text.y = element_blank(),
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        strip.text.x=element_blank(), axis.text.y=element_text(size=12, face="italic"))+ylab(NULL)+xlab(NULL)

heatmap_adjusted_rural_malaria<- ggplot(rural_malaria,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",strip.text.y = element_blank(),
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        strip.text.x=element_blank(), axis.text.y=element_text(size=12, face="italic"))+ylab(NULL)+xlab(NULL)

heatmap_adjusted_rural_ntd<- ggplot(rural_ntd,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",strip.text.y = element_blank(),
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        strip.text.x=element_blank(), axis.text.y=element_text(size=12, face="italic"))+ylab(NULL)+xlab(NULL)

heatmap_adjusted_rural_vpd<- ggplot(rural_vpd,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",strip.text.y = element_blank(),
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        strip.text.x=element_blank(), axis.text.y=element_text(size=12, face="italic"))+ylab(NULL)+xlab(NULL)

heatmap_adjusted_rural_sum<- ggplot(rural_sum,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",strip.text.y = element_blank(),
        axis.text.x=element_text(size=12, angle=45, vjust=0.5), strip.text.x=element_blank(), axis.text.y=element_text(size=12, face="italic"))+ylab(NULL)+xlab(NULL)

rural_grid<-plot_grid(heatmap_adjusted_rural_covid, heatmap_adjusted_rural_enteric,
          heatmap_adjusted_rural_malaria, heatmap_adjusted_rural_ntd, heatmap_adjusted_rural_vpd,heatmap_adjusted_rural_sum,
          ncol=1,align='v', rel_heights = c(1.3,1.0,1.3,1.8,1.9,1.2,1.2))

#################################################

urban_arbo<-adjusted_ranks_urban%>%
  filter(Antigen_Type=="arbovirus")

urban_covid<-adjusted_ranks_urban%>%
  filter(Antigen_Type=="COVID")

urban_enteric<-adjusted_ranks_urban%>%
  filter(Antigen_Type=="enteric")

urban_malaria<-adjusted_ranks_urban%>%
  filter(Antigen_Type=="malaria")


urban_ntd<-adjusted_ranks_urban%>%
  filter(Antigen_Type=="NTD")

urban_vpd<-adjusted_ranks_urban%>%
  filter(Antigen_Type=="VPD")

urban_sum<- adjusted_ranks_urban%>%
  filter(Antigen_Type=="Sum")

urban_arbo$Cluster_clean<- factor(urban_arbo$Cluster_clean,
                                  levels = urban_sum$Cluster_clean[order(as.numeric(urban_sum$Ranks))])

urban_covid$Cluster_Name<- factor(urban_covid$Cluster_Name,
                                   levels = urban_sum$Cluster_Name[order(as.numeric(urban_sum$Ranks))])

urban_enteric$Cluster_Name<- factor(urban_enteric$Cluster_Name,
                                     levels = urban_sum$Cluster_Name[order(as.numeric(urban_sum$Ranks))])

urban_malaria$Cluster_Name<- factor(urban_malaria$Cluster_Name,
                                     levels = urban_sum$Cluster_Name[order(as.numeric(urban_sum$Ranks))])

urban_ntd$Cluster_Name<- factor(urban_ntd$Cluster_Name,
                                 levels = urban_sum$Cluster_Name[order(as.numeric(urban_sum$Ranks))])

urban_vpd$Cluster_Name<- factor(urban_vpd$Cluster_Name,
                                 levels = urban_sum$Cluster_Name[order(as.numeric(urban_sum$Ranks))])

urban_sum$Cluster_Name<- factor(urban_sum$Cluster_Name,
                                 levels = urban_sum$Cluster_Name[order(as.numeric(urban_sum$Ranks))])

heatmap_adjusted_urban_arbo<- ggplot(urban_arbo,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),strip.text.y = element_text(size = 12), 
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=12))+ylab(NULL)+xlab(NULL)

heatmap_adjusted_urban_covid<- ggplot(urban_covid,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),strip.text.y = element_text(size = 12,face="italic"), 
        axis.text.x=element_blank(),axis.ticks.x = element_blank(), 
        strip.text.x = element_text(size=12))+ylab(NULL)+xlab(NULL)

heatmap_adjusted_urban_enteric<- ggplot(urban_enteric,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),strip.text.y = element_text(size = 12, face="italic"), 
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        strip.text.x = element_blank())+ylab(NULL)+xlab(NULL)

heatmap_adjusted_urban_malaria<- ggplot(urban_malaria,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),strip.text.y = element_text(size = 12, face="italic"), 
        axis.text.x=element_blank(),axis.ticks.x = element_blank(), 
        strip.text.x = element_blank())+ylab(NULL)+xlab(NULL)

heatmap_adjusted_urban_ntd<- ggplot(urban_ntd,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),strip.text.y = element_text(size = 12, face="italic"), 
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        strip.text.x = element_blank())+ylab(NULL)+xlab(NULL)

heatmap_adjusted_urban_vpd<- ggplot(urban_vpd,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),strip.text.y = element_text(size = 12, face="italic"),
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        strip.text.x = element_blank())+ylab(NULL)+xlab(NULL)


heatmap_adjusted_urban_sum<- ggplot(urban_sum,aes(x=Cluster_Name, y=Antigen, fill= as.numeric(Ranks)))+
  geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y", scales='free')+
  theme(legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),strip.text.y = element_text(size = 12, face="italic"), axis.text.x=element_text(size=12, angle=45, vjust = 0.5), 
        strip.text.x = element_blank())+ylab(NULL)+xlab(NULL)

urban_grid<-plot_grid(heatmap_adjusted_urban_covid, heatmap_adjusted_urban_enteric,
                      heatmap_adjusted_urban_malaria, heatmap_adjusted_urban_ntd, heatmap_adjusted_urban_vpd,
                      heatmap_adjusted_urban_sum,
                      ncol=1,align='v', rel_heights = c(1.3,1.1,1.3,1.8,1.9,1.2,0.8))

plots_adj<-plot_grid(rural_grid,urban_grid,ncol=2,align='h', rel_widths = c(1.5,0.3))

plots_adj
#unadjusted_ranks_results<- ranks_results%>%
  #filter(Model=="Unadjusted")%>%
  #filter(Antigen%in%c("mfi_dengns12_347","mfi_gexp18_2_261","mfi_etramp5ag1_ls_366","mfi_pvrbp2b_ls_151")==F)%>%
  #mutate(Cluster_Type= ifelse(Cluster%in% c(4,6,24,25,64),"Urban","Rural"))%>%
  #mutate(Antigen_Type= ifelse(Antigen%in% c("mfi_tettox_594", "mfi_diptox_241","mfi_wruv_255","mfi_wmev_215"),"VPD",
                              #ifelse(Antigen%in% c("mfi_dengns11_804","mfi_dengns13_479","mfi_dengns14_706","mfi_chike1_585"),"arbovirus",
                                     #ifelse(Antigen %in% c("mfi_spiken_ls_644","mfi_spikerbd_ls_413"),"COVID",
                                            #ifelse(Antigen %in% c("mfi_cp17_99","mfi_cp23_800","mfi_vsp5_148","mfi_vsp3_99"),"enteric",
                                                   #ifelse(Antigen%in% c("mfi_csp_ls_1250","mfi_glurpr2_ls_212","mfi_pfama1_ls_1419","mfi_pfmsp119_ls_1557","mfi_rh42_401","mfi_pmmsp119_ls_136","mfi_pomsp119_ls_144","mfi_pvdbprii_ls_133","mfi_pvmsp119_ls_220"," mfi_pvrbp2b_ls_151"),"malaria", "NTD"))))))

#unadjusted_ranks_rural<- unadjusted_ranks_results%>%
  #filter(Cluster_Type=="Rural")

#unadjusted_ranks_urban<- unadjusted_ranks_results%>%
  #filter(Cluster_Type=="Urban")

#heatmap_unadjusted_rural<- ggplot(unadjusted_ranks_rural,aes(x=Cluster, y=Antigen, fill= as.numeric(Ranks)))+
  #geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y",scales='free')+
  #theme(legend.position = "none",strip.text.y = element_blank())

#heatmap_unadjusted_urban<- ggplot(unadjusted_ranks_urban,aes(x=Cluster, y=Antigen, fill= as.numeric(Ranks)))+
  #geom_tile()+ scale_fill_viridis()+theme_classic()+facet_grid2(vars(Antigen_Type), vars(Cluster_Type), axes = "all", remove_labels = "y",scales='free')+
  #theme(
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank())+ylab(NULL)+xlab(NULL)+labs(fill = "Random Effects Ranks")



#title <- ggdraw() + 
  #draw_label(
    #"Unadjusted Models",
    #fontface = 'bold',
    #x = 0,
    #hjust = 0
  #) +
  #theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    #plot.margin = margin(0, 0, 0, 7)
  #)


#plots_unadj<-plot_grid(heatmap_unadjusted_rural,heatmap_unadjusted_urban,ncol=2,align='h', rel_widths = c(1.7,1))
#plot_grid(title,plots_adj,ncol=1,align='v',rel_heights=c(0.1,1))


