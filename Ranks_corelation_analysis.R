library(tidyverse)
library(ggplot2)
library(DescTools)

rank_sums<- read.csv("individual level ranks and sums.csv")

adj_seroprev<- read.csv("adj_seroprev_table.csv")


cluster_list<- unique(rank_sums$Cluster_clean)




antigen_list<- unique(rank_sums$Antigen)

cor_rank_df<- data.frame()

for(i in 1:length(antigen_list)){
  for(j in 1:length(antigen_list)){
    antigen_df_1<- rank_sums|>
      filter(Antigen==antigen_list[i])|>
      dplyr::select(Cluster_Type, Ranks, Antigen, Antigen_Type, Cluster_Name)
    
    antigen_df_2<- rank_sums|>
      filter(Antigen==antigen_list[j])|>
      dplyr::select(Cluster_Type, Ranks, Antigen, Antigen_Type, Cluster_Name)|>
      left_join(antigen_df_1, by =c("Cluster_Type", "Cluster_Name"))
    
    if(antigen_list[i]!=antigen_list[j]){
    
    lm_output<- lm(Ranks.x~Ranks.y+Cluster_Type, data= antigen_df_2)
    lm_summary<- summary(lm_output)
    ci_lm<- confint(lm_output)
  
    
    spear_cor<- SpearmanRho(antigen_df_2$Ranks.x,antigen_df_2$Ranks.y,
                            conf.level=0.95)
    spear_cor_p<- cor.test(antigen_df_2$Ranks.x,antigen_df_2$Ranks.y, method='spearman')
    
    spear_cor_rank_df<- data.frame(Antigen_1= antigen_list[i],
                                   Antigen_Type_1= antigen_df_1$Antigen_Type[1],
                                   Antigen_2= antigen_list[j], 
                                   Antigen_Type_2= antigen_df_2$Antigen_Type.x[1],
                                   cor_est= spear_cor_p$estimate,
                                   cor_p_val= spear_cor_p$p.value,
                                   lower_CI_cor= spear_cor[2],
                                   upper_CI_cor= spear_cor[3],
                                   lm_est= lm_output$coefficients[2], 
                                   p_val_lm= lm_summary$coefficients[2,4], 
                                   lower_CI_lm= ci_lm[2,1],
                                   upper_CI_lm= ci_lm[2,2], 
                                   rsquared_lm= lm_summary$r.squared)
    cor_rank_df<- rbind(cor_rank_df, spear_cor_rank_df)
    }
      else {
  
     spear_cor_rank_df<- data.frame(Antigen_1= antigen_list[i],
                                     Antigen_Type_1= antigen_df_1$Antigen_Type[1],
                                       Antigen_2= antigen_list[j], 
                                       Antigen_Type_2= antigen_df_2$Antigen_Type.x[1],
                                    cor_est= NA,
                                    cor_p_val= NA,
                                    lower_CI_cor= NA,
                                    upper_CI_cor= NA,
                                    lm_est= NA, 
                                    p_val_lm= NA, 
                                    lower_CI_lm= NA,
                                    upper_CI_lm= NA, 
                                    rsquared_lm= NA)
      cor_rank_df<- rbind(cor_rank_df, spear_cor_rank_df)
  }
  }  
}

cor_rank_df_full_noNA<- cor_rank_df|>
  filter(Antigen_1=="Overall Vulnerability")|>
  mutate(adj_p_val= p.adjust(cor_p_val, method="BH"))|>
  rename(Correlation = cor_est)|>
  mutate(Significance=ifelse(adj_p_val<0.05,"Sig After Adj", "Not Sig After Adj"))|>
  filter(is.na(Correlation)==F)|>
  rename(Antigen=Antigen_2)|>
  left_join(adj_seroprev, by="Antigen")



cor_rank_df_full_noNA<- cor_rank_df_full_noNA|>
  mutate(Antigen_Type_STI= ifelse(Antigen%in%c("Chlamydia trachomatis (ct694)",
                                               "Chlamydia trachomatis (pgp3)",
                                               "Treponema palladium (rp17)",
                                               "Treponema palladium (tmpa)"), 
                                  "NTD+STI", Antigen_Type_2))

cor_rank_df_full_noNA$Antigen_Type_2<- factor(cor_rank_df_full_noNA$Antigen_Type_2, 
                                    levels=c("COVID", "enteric", 
                                             "malaria", "NTD", 
                                             "VPD", "Sum"))


cor_rank_df_full_noNA$Antigen<- factor(cor_rank_df_full_noNA$Antigen, 
                                       levels= c("SARS-CoV-2(sars2np)",
                                                 "SARS-CoV-2(sars2rbd)",
                                                 "Cryptosporidium parvum(cp23)",
                                                 "Cryptosporidium parvum(cp17)",
                                                 "Giardia lamblia(vsp3)",
                                                 "Giardia lamblia(vsp5)",
                                                 "P. vivax(pvdbprii)",
                                                 "P. vivax(pvmsp119)",
                                                 "P. ovale(pomsp119)",
                                                 "P. malariae(pmmsp119)",
                                                 "P. falciparum(csp)",
                                                 "P. falciparum(pfmsp119)",
                                                 "P. falciparum(pfama1)",
                                                 "P. falciparum(glurpr2)",
                                                 "P. falciparum(rh42)",
                                                 "Brugia malayi (bm14)",
                                                 "Brugia malayi (bm33)",
                                                 "Wuchereria bancrofti (wb123)",
                                                 "Onchocerca volvulus (ov16)",
                                                 "Taenia solium(t24h)",
                                                 "Taenia solium(es33)",
                                                 "Schistosoma mansoni (sea)",
                                                 "Schistosoma mansoni (sm25)",
                                                 "Strongyloides stercoralis (nie)",
                                                 "Treponema palladium (tmpa)",
                                                 "Treponema palladium (rp17)",
                                                 "Chlamydia trachomatis (pgp3)",
                                                 "Chlamydia trachomatis (ct694)",
                                                 "Tetanus (Tet tox)",
                                                 "Rubella (wRuv)",
                                                 "Measles (wMev)",
                                                 "Diphtheria (Dip tox)"))
library(viridis)
library(ggforce)

rank_sums_df<-rank_sums|>
  filter(Antigen=="Overall Vulnerability")|>
  arrange(Ranks)

rank_sums_df$Cluster_Name<- factor(rank_sums_df$Cluster_Name, 
                                   levels= rank_sums_df$Cluster_Name[order(rank_sums_df$Ranks)])

dotplot_overall_vul<- ggplot(rank_sums_df, aes(x= Cluster_Name, y= Ranks, color=Cluster_Type))+
  geom_point()+scale_color_brewer(palette="Dark2")+theme_bw()+
  ylab("Cluster Vulnerability Score")+xlab("Cluster Name")+
  theme(axis.text.x=element_text(angle = 45, vjust=1, hjust=1), 
        legend.title=element_blank())

dotplot_overall_vul<- ggplot(rank_sums_df, aes(x= Cluster_Name, y= Ranks, color=Cluster_Type))+
  geom_point(size=3)+scale_color_brewer(palette="Dark2")+theme_bw()+
  ylab("Cluster Vulnerability Score")+xlab("Cluster Name")+
  theme(axis.text.x=element_text(angle = 45, vjust=1, hjust=1), 
        legend.title=element_blank())+
theme(axis.text.x= element_text(size=16),
      axis.text.y= element_text(size=16),
      axis.title.x=element_text(size=20), 
      axis.title.y=element_text(size=20), 
      legend.text= element_text(size=18))


cor_rank_df_full_noNA_sig<- cor_rank_df_full_noNA|>
  filter(Significance=="Sig After Adj")


forestplot_overallVul<-ggplot(cor_rank_df_full_noNA, aes(y= Antigen, x= Correlation, color=Significance,size=proportion)) +
  geom_pointrange(aes(xmin=lower_CI_cor, xmax=upper_CI_cor))+
  ylab("Antigen") + xlab("Correlation with Overall Vulnerability") +
  facet_grid(vars(Antigen_Type_STI), scales="free", space="free") + 
  theme_bw() + scale_size(rang=c(0,1))+
  theme(strip.text.y = element_text(angle=0))+
  theme(strip.background = element_rect(fill = "white", color = "black"))+
  theme(plot.title = element_text(size=12,color="black", hjust = 0.0))+
  scale_color_brewer(palette="Set1")+
  geom_vline(xintercept = 0)+
  labs(size="Seroprevalence")

forestplot_overallVul_talk<-ggplot(cor_rank_df_full_noNA_sig, aes(y= Antigen, x= Correlation, size=proportion)) +
  geom_pointrange(aes(xmin=lower_CI_cor, xmax=upper_CI_cor))+
  ylab("Antigen") + xlab("Correlation with Overall Vulnerability") +
  facet_grid(vars(Antigen_Type_STI), scales="free", space="free") + 
  theme_bw() + scale_size(rang=c(0,1))+
  theme(strip.text.y = element_text(angle=0, size=18))+
  theme(strip.background = element_rect(fill = "white", color = "black"))+
  theme(plot.title = element_text(size=12,color="black", hjust = 0.0))+
  geom_vline(xintercept = 0)+
  labs(size="Seroprevalence")+
  theme(axis.text.x= element_text(size=18),
        axis.text.y= element_text(size=18),
        axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        legend.title= element_text(size=20), 
        legend.text= element_text(size=18))

library(cowplot)

dotplot_overall_vul

forestplot_overallVul
overall_vuln<- plot_grid(dotplot_overall_vul, plot_map_Zambezia,labels=c("A","B"), ncol=1)
plot_grid( overall_vuln, forestplot_overallVul, ncol= 2, labels=c("","C"))


cor_rank_df_full_noNA_trim<- cor_rank_df_full_noNA|>
  select(Antigen_Type_2, Antigen, Association, p_val, adj_p_val, rsquared, proportion)|>
  rename(Type = Antigen_Type_2)|>
  arrange(Type)

flextable(cor_rank_df_full_noNA_trim)


rank_sums_df<- rank_sums_df|>
  mutate(Vuln_Rank= rank(as.numeric(Ranks)))
