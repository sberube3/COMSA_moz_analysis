library(tidyverse)
library(ggplot2)
library(pROC)
library(flextable)
setwd("/Users/sophieberube/Desktop/Hopkins/Hopkins/Bead_serology/COMSA_cutoff_analysis")


#############################################
#############Function for ROC analysis######
############################################

get_roc<- function(FI_matrix_with_status, ag_variable){
  

  ROC_output<- roc(FI_matrix_with_status$Status~as.numeric(ag_variable), plot=F, print.auc=T )
  
  
  ROC_output_unadjusted_performance<- data.frame(Sens=ROC_output$sensitivities,
                                               Spec=ROC_output$specificities,
                                               Threshold=ROC_output$thresholds)
  
  ROC_output_unadjusted_performance<- ROC_output_unadjusted_performance%>%
    mutate(YoudenJ= Sens+Spec-1, CDC_Lkhd= Sens/(1-Spec))
  
  
  return(list(ROC_output_unadjusted_performance, ROC_output$auc, ROC_output))
  
  
}

##############################################
######Function for Threshold Comparison#######
##############################################
get_cutoff_comparison<- function(Path, Ag, CDC_Thresh, CDC_Sen, CDC_Spe, ROC_Performance){
  
  ROC_performance_100Spec<- ROC_Performance%>%
    filter(Spec==1)%>%
    filter(Sens==max(Sens))
  
  ROC_performance_MaxYouden<- ROC_Performance%>%
    filter(YoudenJ==max(YoudenJ))
  
  ROC_summary<- data.frame(Pathogen=Path, 
                           Antigen=Ag, 
                           CDC_Cutoff= CDC_Thresh, 
                           CDC_Sens= CDC_Sen, 
                           CDC_Spec= CDC_Spe, 
                           JHU_Cutoff= round(ROC_performance_MaxYouden$Threshold,digits=3), 
                           JHU_Sens= round(ROC_performance_MaxYouden$Sens,digits=3),
                           JHU_Spec= round(ROC_performance_MaxYouden$Spec,digits=3), 
                           JHU_YoudenJ= round(ROC_performance_MaxYouden$YoudenJ,digits=3), 
                           FPR_Cutoff= round(ROC_performance_100Spec$Threshold,digits=3), 
                           FPR_Sens= round(ROC_performance_100Spec$Sens,digits=3),
                           FPR_Spec= round(ROC_performance_100Spec$Spec,digits=3), 
                           FPR_YoudenJ= round(ROC_performance_100Spec$YoudenJ,digits=3))
  ROC_summary<- ROC_summary%>%
    mutate(CDC_YoudenJ= round((CDC_Sens+CDC_Spec-1),digits=3))
  
  ROC_summary_long<- gather(ROC_summary,"Measure", "Value", 3:14)
  ROC_summary_split<- ROC_summary_long%>%
    mutate(Protocol=ifelse(substring(Measure,1,3)=="CDC", "CDC", ifelse(substring(Measure,1,3)=="JHU","MaxYJ", "MaxSpec")))%>%
    mutate(Measure=ifelse(substring(Measure,5,10000L)=="Cutoff", "Cutoff", 
                          ifelse(substring(Measure,5,10000L)=="Sens", "Sens",
                                 ifelse(substring(Measure,5,10000L)=="Spec", "Spec","YoudenJ"))))
  
  ROC_summary_full<- spread(ROC_summary_split, Protocol, Value)
  
  
  return(ROC_summary_full)
  
}

get_final_cutoff<- function(pathogen, antigen, cutoff_method, table, numpos, numneg){
  df<- data.frame(Pathogen= pathogen, 
                  Antigen= antigen, 
                  Cutoff= table[1,c(cutoff_method)], 
                  Cutoff_Method= cutoff_method, 
                  Sensitivity= table[2,c(cutoff_method)], 
                  Specificity= table[3,c(cutoff_method)],
                  Npos= numpos, 
                  Nneg= numneg)
  
}

get_AUC_plot<- function(final_cutoff, full_roc, aucVal, nPos, nNeg){
  
  auc_plot<- ggroc(full_roc)+
    theme_bw()+xlab("Specificity")+ylab("Sensitivity")+
    ggtitle(paste0(final_cutoff$Pathogen, ",", final_cutoff$Antigen))+
    annotate("text", x=0.25, y=0.8, label= paste0("AUC=",round(aucVal,digits=3)))+
    annotate("text", x=0.25, y=0.6, label = paste0("Cutoff (MFI)","=", final_cutoff$Cutoff))+
    annotate("text", x=0.25, y=0.4, label= paste0("Npos=",nPos, ", ","Nneg=", nNeg))
  return(auc_plot)
}

#############################
##########Crypto#############
############################

crypto_ID<- read.csv("crypto_ID.csv")

crypto_cp17<- crypto_ID%>%
  dplyr::select(Cp17..55., Status)%>%
  rename(Cp17=Cp17..55.)

crypto_cp23<- crypto_ID%>%
  dplyr::select(Cp23..19., Status.1)%>%
  rename(Cp23=Cp23..19., Status=Status.1)

#Cp17


cp17_roc<- get_roc(crypto_cp17, crypto_cp17$Cp17)[[1]]

cp17_table<- get_cutoff_comparison("Crypto", "Cp17", 132, 0.941, 1.00, cp17_roc)

cp17_table

cp17_cutoff<- get_final_cutoff("Crypto", "Cp17", "MaxYJ", cp17_table,
                               length(which(crypto_cp17$Status=="Pos")),
                               length(which(crypto_cp17$Status=="Neg")))


cp17_auc<- get_roc(crypto_cp17, crypto_cp17$Cp17)[[2]]
cp17_roc_obj<- get_roc(crypto_cp17, crypto_cp17$Cp17)[[3]]



cp17_plot<- get_AUC_plot(cp17_cutoff, cp17_roc_obj, cp17_auc,nPos= cp17_cutoff$Npos, nNeg=cp17_cutoff$Nneg)

#ROC unadjusted Cp23 

cp23_roc<- get_roc(crypto_cp23, crypto_cp23$Cp23)[[1]]
cp23_table<- get_cutoff_comparison("Crypto", "Cp23", 830, 1.00, 0.95, cp23_roc)
cp23_table


cp23_cutoff<- get_final_cutoff("Crypto", "Cp23", "MaxYJ", cp23_table,
                               length(which(crypto_cp23$Status=="Pos")),
                               length(which(crypto_cp23$Status=="Neg")))

cp23_auc<- get_roc(crypto_cp23, crypto_cp23$Cp23)[[2]]
cp23_roc_obj<- get_roc(crypto_cp23, crypto_cp23$Cp23)[[3]]

cp23_plot<- get_AUC_plot(cp23_cutoff, cp23_roc_obj, cp23_auc,nPos= cp23_cutoff$Npos, nNeg=cp23_cutoff$Nneg)


#############################
##########Giardia############
############################

giardia_ID<- read.csv("giardia_ID.csv")

giardia_vsp3<- giardia_ID%>%
  dplyr::select(Giardia.VSP3..42., Status)%>%
  rename(VSP3= Giardia.VSP3..42.)%>%
  filter(Status!="")

giardia_vsp5<- giardia_ID%>%
  dplyr::select(Giardia.VSP5..51.,Status.1)%>%
  rename(VSP5=Giardia.VSP5..51., Status=Status.1)%>%
  filter(Status!="")

vsp3_roc<- get_roc(giardia_vsp3, giardia_vsp3$VSP3)[[1]]

flextable(vsp3_roc[,c("Sens", "Spec", "Threshold", "YoudenJ")])

vsp3_table<- get_cutoff_comparison("Giardia", "vsp3", 99, 0.909, 0.900, vsp3_roc)
vsp3_table
vsp3_cutoff<- get_final_cutoff("Giardia", "vsp3", "MaxYJ", vsp3_table,
                               length(which(giardia_vsp3$Status=="Pos")),
                               length(which(giardia_vsp5$Status=="Neg")))

vsp3_auc<- get_roc(giardia_vsp3, giardia_vsp3$VSP3)[[2]]
vsp3_roc_obj<- get_roc(giardia_vsp3, giardia_vsp3$VSP3)[[3]]

vsp3_plot<- get_AUC_plot(vsp3_cutoff, vsp3_roc_obj, vsp3_auc, nPos= vsp3_cutoff$Npos, nNeg=vsp3_cutoff$Nneg)


vsp5_roc<- get_roc(giardia_vsp5, giardia_vsp5$VSP5)[[1]]
vsp5_table<- get_cutoff_comparison("Giardia", "vsp5", 148, 0.909, 1.00, vsp5_roc)
vsp5_table

vsp5_cutoff<- get_final_cutoff("Giardia", "vsp5", "MaxYJ", vsp5_table, 
                               length(which(giardia_vsp5$Status=="Pos")),
                               length(which(giardia_vsp5$Status=="Neg")))

vsp5_auc<- get_roc(giardia_vsp5, giardia_vsp5$VSP5)[[2]]
vsp5_roc_obj<- get_roc(giardia_vsp5, giardia_vsp5$VSP5)[[3]]

vsp5_plot<- get_AUC_plot(vsp5_cutoff, vsp5_roc_obj, vsp5_auc, nPos= vsp5_cutoff$Npos, nNeg=vsp5_cutoff$Nneg)


###########################
#######Cysticercosis#######
###########################

Cyste_FI<- read.csv('Cyste_ID.csv')

Cyste_FI_clean<- Cyste_FI%>%
  dplyr::select("T24H..13.", "Status")




#t24h match 

t24h_roc<- get_roc(Cyste_FI_clean, Cyste_FI_clean$`T24H..13.`)

flextable(t24h_roc[[1]][,c("Sens", "Spec", "Threshold", "YoudenJ")])

t24h_table<- get_cutoff_comparison("Cystecercosis", "t24h", 84, 0.975, 0.968, t24h_roc[[1]])
t24h_table

t24h_cutoff<- get_final_cutoff("Cystecercosis", "t24h", "MaxSpec", t24h_table,
                               length(which(Cyste_FI_clean$Status=="Pos")),
                               length(which(Cyste_FI_clean$Status=="Neg")))

t24h_auc<- get_roc(Cyste_FI_clean, Cyste_FI_clean$`T24H..13.`)[[2]]
t24h_roc_obj<- get_roc(Cyste_FI_clean, Cyste_FI_clean$`T24H..13.`)[[3]]

t24h_plot<- get_AUC_plot(t24h_cutoff, t24h_roc_obj, t24h_auc, nPos= t24h_cutoff$Npos, nNeg=t24h_cutoff$Nneg)




Taen_FI<- read.csv('Taeniasis_IDs.csv')

Taen_FI_clean<- Taen_FI%>%
  dplyr::select("ES33..61.", "Status")

es33_roc<- get_roc(Taen_FI_clean, Taen_FI_clean$`ES33..61.`)[[1]]
es33_table<- get_cutoff_comparison("Cystecercosis", "es33", 14, 0.938, 0.989, es33_roc)

es33_table

es33_cutoff<- get_final_cutoff("Cystecercosis", "es33", "MaxSpec", es33_table,
                               length(which(Taen_FI_clean$Status=="Pos")),
                               length(which(Taen_FI_clean$Status=="Neg")))

es33_auc<- get_roc(Taen_FI_clean, Taen_FI_clean$`ES33..61.`)[[2]]
es33_roc_obj<- get_roc(Taen_FI_clean, Taen_FI_clean$`ES33..61.`)[[3]]

es33_plot<- get_AUC_plot(es33_cutoff, es33_roc_obj, es33_auc, nPos= es33_cutoff$Npos, nNeg=es33_cutoff$Nneg)


##################################
##############Oncho##############
#################################

oncho_ID<- read.csv("Oncho_IDs.csv")

oncho_ID_status<- oncho_ID%>%
  dplyr::select("OV16..28.", "Status")

#ROC unadjusted OV16 

ov16_roc<- get_roc(oncho_ID_status, oncho_ID_status$OV16..28.)[[1]]

ov16_table<- get_cutoff_comparison("Onchocercaisis", "ov16", 422, 0.909, 0.947, ov16_roc)
ov16_table


ov16_cutoff<- get_final_cutoff("Onchocercaisis", "ov16", "MaxYJ", ov16_table,
                               length(which(oncho_ID_status$Status=="Pos")),
                               length(which(oncho_ID_status$Status=="Neg")))

ov16_auc<- get_roc(oncho_ID_status, oncho_ID_status$OV16..28.)[[2]]
ov16_roc_obj<- get_roc(oncho_ID_status, oncho_ID_status$OV16..28.)[[3]]

ov16_plot<- get_AUC_plot(ov16_cutoff, ov16_roc_obj, ov16_auc, nPos= ov16_cutoff$Npos, nNeg=ov16_cutoff$Nneg)



##########################
########Schist###########
#########################
schist_csv<- read.csv("schisto_species.csv")
schist_SEA<- read.csv("Schist_SEA_ID.csv")
schist_Sm25<- read.csv("Schist_SM25_ID.csv")

schist_csv<- schist_csv%>%
  rename(Description= Code)

schist_SEA_species<- merge(schist_SEA, schist_csv, by="Description", all.x=T)


schist_SEA_species_trim<- schist_SEA_species[,c("SEA", "Status", "Species")]

schist_sm25_species<- merge(schist_Sm25, schist_csv, by="Description", all.x=T)


schist_sm25_species_trim<- schist_sm25_species[,c("Sm25", "Status", "Species")]

schist_SEA_species_trim<- schist_SEA_species_trim%>%
  mutate(Global_staus= ifelse(Status=="Neg", "Neg", "Pos"))%>%
  mutate(Mansoni_status= ifelse(Status=="Neg", "Neg",
                                ifelse(Species=="S. mansoni", "Pos", NA)))%>%
  mutate(Inter_status= ifelse(Status=="Neg", "Neg",
                                ifelse(Species=="S. intercalatum", "Pos", NA)))%>%
  mutate(Haem_status= ifelse(Status=="Neg", "Neg",
                              ifelse(Species=="S. haematobium", "Pos", NA)))

sea_roc<- roc(schist_SEA_species_trim$Global_staus~schist_SEA_species_trim$SEA)
  
sea_unadjusted_performance<- data.frame(Sens=sea_roc$sensitivities,
                                               Spec=sea_roc$specificities,
                                               Threshold=sea_roc$thresholds)

sea_unadjusted_performance<- sea_unadjusted_performance%>%
  mutate(YoudenJ= Sens+Spec-1, CDC_Lkhd= Sens/(1-Spec))

sea_unadjusted_performance_highSens<- sea_unadjusted_performance%>%
  filter(Sens>0.7)

sea_unadjusted_performance_MaxSpec<- sea_unadjusted_performance_highSens[which(sea_unadjusted_performance_highSens$Spec== max(sea_unadjusted_performance_highSens$Spec)),]

sea_unadjusted_performance_MaxSens<- sea_unadjusted_performance_MaxSpec[which(sea_unadjusted_performance_MaxSpec$Sens==max(sea_unadjusted_performance_MaxSpec$Sens)),]

sea_global_threshold<- sea_unadjusted_performance_MaxSens$Threshold
sea_global_sens<- sea_unadjusted_performance_MaxSens$Sens
sea_global_spec<- sea_unadjusted_performance_MaxSens$Spec

sea_mans_sens<- length(which(schist_SEA_species_trim$SEA>410&schist_SEA_species_trim$Mansoni_status=="Pos"))/length(which(schist_SEA_species_trim$Mansoni_status=="Pos"))
sea_mans_spec<- length(which(schist_SEA_species_trim$SEA<=410&schist_SEA_species_trim$Mansoni_status=="Neg"))/length(which(schist_SEA_species_trim$Mansoni_status=="Neg"))

sea_inter_sens<- length(which(schist_SEA_species_trim$SEA>410&schist_SEA_species_trim$Inter_status=="Pos"))/length(which(schist_SEA_species_trim$Inter_status=="Pos"))
sea_inter_spec<- length(which(schist_SEA_species_trim$SEA<=410&schist_SEA_species_trim$Inter_status=="Neg"))/length(which(schist_SEA_species_trim$Inter_status=="Neg"))

sea_haem_sens<- length(which(schist_SEA_species_trim$SEA>410&schist_SEA_species_trim$Haem_status=="Pos"))/length(which(schist_SEA_species_trim$Haem_status=="Pos"))
sea_haem_spec<- length(which(schist_SEA_species_trim$SEA<=410&schist_SEA_species_trim$Haem_status=="Neg"))/length(which(schist_SEA_species_trim$Haem_status=="Neg"))

sea_cutoff_global<- data.frame(Pathogen= "Schistosomiasis Global", 
                               Antigen= "sea", 
                               Cutoff= sea_global_threshold, 
                               Cutoff_Method= "MaxYJ", 
                               Sensitivity= sea_global_sens, 
                               Specificity= sea_global_spec, 
                               Npos= length(which(schist_SEA_species_trim$Global_staus=="Pos")),
                               Nneg= length(which(schist_SEA_species_trim$Global_staus=="Neg")))


sea_cutoff_mansoni<- data.frame(Pathogen= "Schistosomiasis Mansoni", 
                                Antigen= "sea", 
                                Cutoff= sea_global_threshold, 
                                Cutoff_Method= "MaxYJ", 
                                Sensitivity= sea_mans_sens, 
                                Specificity= sea_mans_spec, 
                                Npos= length(which(schist_SEA_species_trim$Mansoni_status=="Pos")),
                                Nneg= length(which(schist_SEA_species_trim$Mansoni_status=="Neg")))


sea_cutoff_inter<- data.frame(Pathogen= "Schistosomiasis Intercalatum", 
                                Antigen= "sea", 
                                Cutoff= sea_global_threshold, 
                                Cutoff_Method= "MaxYJ", 
                                Sensitivity= sea_inter_sens, 
                                Specificity= sea_inter_spec,
                              Npos= length(which(schist_SEA_species_trim$Inter_status=="Pos")),
                              Nneg= length(which(schist_SEA_species_trim$Inter_status=="Neg")))


sea_cutoff_heam<- data.frame(Pathogen= "Schistosomiasis Haematobium", 
                             Antigen= "sea", 
                             Cutoff= sea_global_threshold, 
                             Cutoff_Method= "MaxSpec", 
                             Sensitivity= sea_haem_sens, 
                             Specificity= sea_haem_spec, 
                             Npos= length(which(schist_SEA_species_trim$Haem_status=="Pos")),
                             Nneg= length(which(schist_SEA_species_trim$Haem_status=="Neg")))





sea_plot<- get_AUC_plot(sea_cutoff_global, sea_roc, sea_roc$auc,sea_cutoff_inter$Npos, sea_cutoff_inter$Nneg )



sm25_roc<- roc(schist_sm25_species_trim$Status~schist_sm25_species_trim$Sm25)

sm25_unadjusted_performance<- data.frame(Sens=sm25_roc$sensitivities,
                                        Spec=sm25_roc$specificities,
                                        Threshold=sm25_roc$thresholds)

sm25_unadjusted_performance<- sm25_unadjusted_performance%>%
  mutate(YoudenJ= Sens+Spec-1, CDC_Lkhd= Sens/(1-Spec))

sm25_unadjusted_performance_highSens<- sm25_unadjusted_performance%>%
  filter(Sens>0.7)

sm25_unadjusted_performance_MaxSpec<- sm25_unadjusted_performance_highSens[which(sm25_unadjusted_performance_highSens$Spec== max(sm25_unadjusted_performance_highSens$Spec)),]

sm25_unadjusted_performance_MaxSens<- sm25_unadjusted_performance_MaxSpec[which(sm25_unadjusted_performance_MaxSpec$Sens==max(sm25_unadjusted_performance_MaxSpec$Sens)),]

sm25_threshold<- sm25_unadjusted_performance_MaxSens$Threshold
sm25_global_sens<- sm25_unadjusted_performance_MaxSens$Sens
sm25_global_spec<- sm25_unadjusted_performance_MaxSens$Spec

sm25_cutoff<- data.frame(Pathogen= "Schistosomiasis Mansoni", 
                             Antigen= "sm25", 
                             Cutoff= sm25_threshold, 
                             Cutoff_Method= "MaxSpec", 
                             Sensitivity= sm25_global_sens, 
                             Specificity= sm25_global_spec, 
                            Npos= length(which(schist_sm25_species_trim$Status=="Pos")),
                         Nneg= length(which(schist_sm25_species_trim$Status=="Neg")))

sm25_plot<- get_AUC_plot(sm25_cutoff, sm25_roc, sm25_roc$auc,sm25_cutoff$Npos, sm25_cutoff$Nneg )

#####################################
##########Trachoma##################
####################################


trachoma_FI<- read.csv('Trachoma_FI_ID.csv')

pgp_FI_with_status<- trachoma_FI%>%
  dplyr::select("Description", "pgp3..M63.", "Status")

ct694_FI_with_status<- trachoma_FI%>%
  dplyr::select("Description.1","CT694..M47.", "Status.1",)


colnames(ct694_FI_with_status)<- c("Description","CT694", "Status")



#ROC unadjusted pgp3

pgp3_roc<- get_roc(pgp_FI_with_status, pgp_FI_with_status$pgp3..M63.)[[1]]

flextable(pgp3_roc[,c("Sens", "Spec", "Threshold", "YoudenJ")])

pgp3_roc_table<- get_cutoff_comparison("Trachoma", "pgp3", 212, 0.926, 1.00, pgp3_roc)
pgp3_roc_table


pgp3_cutoff<- get_final_cutoff("Trachoma", "pgp3", "CDC", pgp3_roc_table, 
                               length(which(pgp_FI_with_status$Status=="Pos")),
                               length(which(pgp_FI_with_status$Status=="Neg")))

pgp3_auc<- get_roc(pgp_FI_with_status, pgp_FI_with_status$pgp3..M63.)[[2]]
pgp3_roc_obj<- get_roc(pgp_FI_with_status, pgp_FI_with_status$pgp3..M63.)[[3]]

pgp3_plot<- get_AUC_plot(pgp3_cutoff, pgp3_roc_obj, pgp3_auc, pgp3_cutoff$Npos, pgp3_cutoff$Nneg)


ct694_roc<- get_roc(ct694_FI_with_status, ct694_FI_with_status$CT694)[[1]]

ct694_roc_table<- get_cutoff_comparison("Trachoma", "ct694", 108, 0.958, 0.973, ct694_roc)
ct694_roc_table


ct694_cutoff<- get_final_cutoff("Trachoma", "ct694", "CDC", ct694_roc_table,
                                length(which(ct694_FI_with_status$Status=="Pos")),
                                length(which(ct694_FI_with_status$Status=="Neg")))

ct694_auc<- get_roc(ct694_FI_with_status, ct694_FI_with_status$CT694)[[2]]
ct694_roc_obj<- get_roc(ct694_FI_with_status, ct694_FI_with_status$CT694)[[3]]

ct694_plot<- get_AUC_plot(ct694_cutoff, ct694_roc_obj, ct694_auc, ct694_cutoff$Npos,  ct694_cutoff$Nneg)

#################################
########Yaws####################
################################

yaws_FI<- read.csv("Yaws_FI.csv")

yaws_FI_rp17_with_status<- yaws_FI%>%
  dplyr::select("Description", "rp17..65.", "Status")

colnames(yaws_FI_rp17_with_status)<- c("Description", "rp17", "Status1")

yaws_FI_rp17_with_status<- yaws_FI_rp17_with_status%>%
  mutate(Status=ifelse(Status1=="Neg"|Status1=="NEG", "NEG", "POS"))

yaws_FI_tmpa_with_status<- yaws_FI%>%
  dplyr::select("Description.1","TmpA..25.", "Status.1")

yaws_FI_tmpa_with_status<- yaws_FI_tmpa_with_status%>%
  mutate(Status=ifelse(Status.1=="Neg"|Status.1=="NEG", "NEG", "POS"))

rp17_roc<- get_roc(yaws_FI_rp17_with_status, yaws_FI_rp17_with_status$rp17)[[1]]

flextable(rp17_roc[,c("Sens", "Spec", "Threshold", "YoudenJ")])

rp17_roc_table<- get_cutoff_comparison("Yaws", "rp17", 86, 0.944, 0.957, rp17_roc)
rp17_roc_table


rp17_cutoff<- get_final_cutoff("Yaws", "rp17", "MaxYJ", rp17_roc_table, 
                               length(which(yaws_FI_rp17_with_status$Status=="POS")),
                               length(which(yaws_FI_rp17_with_status$Status=="NEG")))

rp17_auc<- get_roc(yaws_FI_rp17_with_status, yaws_FI_rp17_with_status$rp17)[[2]]
rp17_roc_obj<- get_roc(yaws_FI_rp17_with_status, yaws_FI_rp17_with_status$rp17)[[3]]

rp17_plot<- get_AUC_plot(rp17_cutoff, rp17_roc_obj, rp17_auc, rp17_cutoff$Npos, rp17_cutoff$Nneg)

tmpa_roc<- get_roc(yaws_FI_tmpa_with_status, yaws_FI_tmpa_with_status$TmpA..25.)[[1]]
flextable(tmpa_roc[,c("Sens", "Spec", "Threshold", "YoudenJ")])
tmpa_roc_table<- get_cutoff_comparison("Yaws", "tmpa", 235, 1.00, 0.95, tmpa_roc)
tmpa_roc_table


tmpa_cutoff<- get_final_cutoff("Yaws", "tmpa", "MaxYJ", tmpa_roc_table,
                               length(which(yaws_FI_tmpa_with_status$Status=="POS")),
                               length(which(yaws_FI_tmpa_with_status$Status=="NEG")))

tmpa_auc<- get_roc(yaws_FI_tmpa_with_status, yaws_FI_tmpa_with_status$TmpA..25.)[[2]]
tmpa_roc_obj<- get_roc(yaws_FI_tmpa_with_status, yaws_FI_tmpa_with_status$TmpA..25.)[[3]]

tmpa_plot<- get_AUC_plot(tmpa_cutoff, tmpa_roc_obj, tmpa_auc, tmpa_cutoff$Npos, tmpa_cutoff$Nneg)

################################
##########Strongy###############
################################

strongy_FI<- read.csv("Strongy_FIs.csv")

strongy_FI_with_status<- strongy_FI%>%
  dplyr::select(c("NIE..74.", "Status"))

nie_roc<- get_roc(strongy_FI_with_status, strongy_FI_with_status$NIE..74.)[[1]]

nie_roc_table<- get_cutoff_comparison("Strongyloides", "nie", 526, 0.79, 1.00, nie_roc)
nie_roc_table


nie_cutoff<- get_final_cutoff("Strongyloides", "nie", "MaxYJ", nie_roc_table, 
                              length(which(strongy_FI_with_status$Status=="Pos")),
                              length(which(strongy_FI_with_status$Status=="Neg")))


nie_auc<- get_roc(strongy_FI_with_status, strongy_FI_with_status$NIE..74.)[[2]]
nie_roc_obj<- get_roc(strongy_FI_with_status, strongy_FI_with_status$NIE..74.)[[3]]

nie_plot<- get_AUC_plot(nie_cutoff, nie_roc_obj, nie_auc, nie_cutoff$Npos, nie_cutoff$Nneg)


overall_cutoff_df<- data.frame(rbind(cp17_cutoff, cp23_cutoff,
                                     vsp3_cutoff, vsp5_cutoff,
                                     t24h_cutoff, es33_cutoff,
                                     ov16_cutoff, sea_cutoff_global,
                                     sea_cutoff_mansoni, sea_cutoff_inter,
                                     sea_cutoff_heam, sm25_cutoff,
                                     ct694_cutoff, pgp3_cutoff,
                                     rp17_cutoff, tmpa_cutoff,
                                     nie_cutoff))

write.csv(overall_cutoff_df, "NTD_Cutoffs.csv")
library(cowplot)

complete_AUC_plot_grid<- plot_grid(cp17_plot, cp23_plot, vsp3_plot,
                                   vsp5_plot, t24h_plot, es33_plot,
                                   ov16_plot, sea_plot, sm25_plot,
                                   pgp3_plot, ct694_plot, rp17_plot,
                                   tmpa_plot,nie_plot, nrow=5, labels="AUTO")

