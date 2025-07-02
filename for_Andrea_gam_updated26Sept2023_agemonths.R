library(tidyverse)
library(mgcv)
library(palmerpenguins)
library(ggplot2)
library(gridExtra)
library(grid)
library(patchwork)
detach(package:gam)

setwd ("C:/Users/acarcele/Dropbox/COMSA - Serosurveillance/Data/Analysis")
## Read in data
#dat_raw <- read.csv("clean_data_12September2023.csv")
dat_raw <- read.csv("vpd data agemos.csv")

## Pick one antigen
ANTIGEN_CHOSEN <- "measles"
# "measles", "rube", "tetanus_01", "dipth_01",
# "mfi_wb123_58", "mfi_tmpa_247", "mfi_rp17_73", "mfi_es33_22", "mfi_t24h_103", "mfi_nie_526", "mfi_sea_410", "mfi_sm25_38", "mfi_ov16_422", "mfi_ct694_108", "mfi_pgp3_212", "mfi_bm33_375", "mfi_bm14_107",
# "mfi_rh42_448", "mfi_pvrbp2b_ls_577", "mfi_pvmsp119_ls_550", "mfi_pvdbprii_ls_534", "mfi_pomsp119_ls_373", "mfi_pmmsp119_ls_370", "mfi_pfmsp119_ls_405", "mfi_pfama1_ls_143", "mfi_glurpr2_ls_160", "mfi_gexp18_2_632", "mfi_etramp5ag1_ls_168", "mfi_csp_ls_392",
# "mfi_vsp5_148", "mfi_vsp3_72", "mfi_cp23_830", "mfi_cp17_72",
# "mfi_spikerbd_ls_387", "mfi_spiken_ls_615",
# "any_deng", mfi_dengns14_706", "mfi_dengns13_479", "mfi_dengns12_347", "mfi_dengns11_804", "mfi_chike1_585"
########## DIDN'T DO: any_dengue, "mfi_pvrbp2b_ls_577" 

dat_raw %>%
  mutate(serostatus = ifelse((!!sym(ANTIGEN_CHOSEN))=="positive", 1, 0)) %>%
  select(age_months, serostatus, final_wt) -> dat_clean
#  mutate(serostatus = ifelse((!!sym(ANTIGEN_CHOSEN))=="1", 1, 0)) %>%


## Have the weights sum to the sample size, and not the population size
dat_clean %>%
  mutate(final_wt_scaled = final_wt / mean(final_wt)) -> dat_clean

## Fit the GAM model
## There are many basis functions ("bs"), here we go with "cr"=cubic regression spline
## We can use the logit or the cloglog link functions; doesn't make a difference when we looked
fit.gam.logit.cr <- mgcv::gam(serostatus~s(age_months, bs="cr"),
                              family=binomial(link="logit"), 
                              weights=final_wt_scaled,
                              data=dat_clean, 
                              method="REML")

## Ages to get 95% CIs at---change age
new_data <- data.frame(age_months = 0:59)

## Fit models to 1000 bootstrap replicates of he data, logit link
predictions.logit.cr = replicate(n=1000, {
  boot <- dat_clean[sample.int(nrow(dat_clean), replace=TRUE),]
  model <- gam(serostatus ~ s(age_months, bs="cr"), 
               family=binomial(link="logit"), 
               weights=final_wt_scaled,
               data=boot, 
               method="REML")
  predict(model, new_data, type="response")
  
}
)

## Generate point estimates and 95% CI---change age
dat_CI <- data.frame(age_months=0:59,
                  logit.lb=apply(predictions.logit.cr, 1, quantile, 0.025),
                  logit.mean = apply(predictions.logit.cr, 1, quantile, 0.5),
                  logit.ub=apply(predictions.logit.cr, 1, quantile, 0.975))

## Plot
dat_clean %>%
  group_by(age_months) %>%
  summarise(seroprev = sum(serostatus) / n(),
            seroprev.wtd = weighted.mean(serostatus, final_wt_scaled),
            N = n()) -> tmp
par(mfrow = c(1,1))            
## Graph to see output
plot(x=tmp$age_nodecim, y=tmp$seroprev, 
     cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
lines(dat_CI$age_nodecim, dat_CI$logit.mean, col="red")
lines(dat_CI$age_nodecim, dat_CI$logit.lb, col="red", lty=2)
lines(dat_CI$age_nodecim, dat_CI$logit.ub, col="red", lty=2)
title(paste0(ANTIGEN_CHOSEN))
#ggsave("dengue.png",height = 5,width = 10)

##save outputs to folder of seroprevalences
setwd ("C:/Users/acarcelan/Dropbox/COMSA - Serosurveillance/Data/Analysis/GAM")
dat_CI$ag_name <- paste0(ANTIGEN_CHOSEN)

# Create the file name
csvFileName <- paste0(ANTIGEN_CHOSEN,".csv")

# Write out
write.csv(dat_CI, file=csvFileName) 

################################################
####combined plots

setwd ("C:/Users/acarcelen/Dropbox/COMSA - Serosurveillance/Data/Analysis")
# ANTIGEN NAMES to match below graphs
# measles, rube, tet, dip
# wb123, tmpa, rp17, es33, t24h, nie, sea, sm25, ov16, ct694, pgp3, bm33, bm14
# rh42, pvrbp, pvmsp, pvdbpr, pomsp, pmmsp, pfmsp, pfama1, glurpr2, gexp, etramp, csp
# vsp5, vsp3, cp23, cp17
# rbd, sarsn
# dengue, chik

### seroprevalence
dat_clean %>%
  group_by(age_nodecim) %>%
  summarise(seroprev = sum(serostatus) / n(),
            seroprev.wtd = weighted.mean(serostatus, final_wt_scaled),
            N = n()) -> tmpcsp                                   #insert antigen name here#

###bootstrap data####
setwd ("C:/Users/acarcele/Dropbox/COMSA - Serosurveillance/Data/Analysis/GAM")
## Read in data
dat_CIcsp <- read.csv(paste(ANTIGEN_CHOSEN,".csv", sep = ""))             #insert antigen name here for dat_CI AND .csv#
            
#####FIRST set of plots
par(mar = c(1.75, 1.75, 1.75, 1.75))
par(mfrow = c(4,3))


#plot7
plot(x=tmpdengue$age_nodecim, y=tmpdengue$seroprev, 
     cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
lines(dat_CIdengue$age_nodecim, dat_CIdengue$logit.mean, col="red")
lines(dat_CIdengue$age_nodecim, dat_CIdengue$logit.lb, col="red", lty=2)
lines(dat_CIdengue$age_nodecim, dat_CIdengue$logit.ub, col="red", lty=2)
title("Any Dengue")

#plot8
plot(x=tmpchik$age_nodecim, y=tmpchik$seroprev, 
     cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
lines(dat_CIchik$age_nodecim, dat_CIchik$logit.mean, col="red")
lines(dat_CIchik$age_nodecim, dat_CIchik$logit.lb, col="red", lty=2)
lines(dat_CIchik$age_nodecim, dat_CIchik$logit.ub, col="red", lty=2)
title("Chikungunya(chike1)")

#plot20
plot(x=tmprbd$age_nodecim, y=tmprbd$seroprev, 
     cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
lines(dat_CIrbd$age_nodecim, dat_CIrbd$logit.mean, col="red")
lines(dat_CIrbd$age_nodecim, dat_CIrbd$logit.lb, col="red", lty=2)
lines(dat_CIrbd$age_nodecim, dat_CIrbd$logit.ub, col="red", lty=2)
title("SARS-CoV-2(sars2rbd)")

#plot21
plot(x=tmpsarsn$age_nodecim, y=tmpsarsn$seroprev, 
     cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
lines(dat_CIsarsn$age_nodecim, dat_CIsarsn$logit.mean, col="red")
lines(dat_CIsarsn$age_nodecim, dat_CIsarsn$logit.lb, col="red", lty=2)
lines(dat_CIsarsn$age_nodecim, dat_CIsarsn$logit.ub, col="red", lty=2)
title("SARS-CoV-2(sars2np)")

#plot19
plot(x=tmpvsp5$age_nodecim, y=tmpvsp5$seroprev, 
     cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
lines(dat_CIvsp5$age_nodecim, dat_CIvsp5$logit.mean, col="red")
lines(dat_CIvsp5$age_nodecim, dat_CIvsp5$logit.lb, col="red", lty=2)
lines(dat_CIvsp5$age_nodecim, dat_CIvsp5$logit.ub, col="red", lty=2)
title("Giardia lamblia(vsp5)")

#plot18
plot(x=tmpvsp3$age_nodecim, y=tmpvsp3$seroprev, 
     cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
lines(dat_CIvsp3$age_nodecim, dat_CIvsp3$logit.mean, col="red")
lines(dat_CIvsp3$age_nodecim, dat_CIvsp3$logit.lb, col="red", lty=2)
lines(dat_CIvsp3$age_nodecim, dat_CIvsp3$logit.ub, col="red", lty=2)
title("Giardia lamblia(vsp3)")

#plot16
plot(x=tmpcp23$age_nodecim, y=tmpcp23$seroprev, 
     cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
lines(dat_CIcp23$age_nodecim, dat_CIcp23$logit.mean, col="red")
lines(dat_CIcp23$age_nodecim, dat_CIcp23$logit.lb, col="red", lty=2)
lines(dat_CIcp23$age_nodecim, dat_CIcp23$logit.ub, col="red", lty=2)
title("Cryptosporidium parvum(cp23)")

#plot17
plot(x=tmpcp17$age_nodecim, y=tmpcp17$seroprev, 
     cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
lines(dat_CIcp17$age_nodecim, dat_CIcp17$logit.mean, col="red")
lines(dat_CIcp17$age_nodecim, dat_CIcp17$logit.lb, col="red", lty=2)
lines(dat_CIcp17$age_nodecim, dat_CIcp17$logit.ub, col="red", lty=2)
title("Cryptosporidium parvum(cp17)")

#plot3
plot(x=tmpmeasles$age_nodecim, y=tmpmeasles$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CImeasles$age_nodecim, dat_CImeasles$logit.mean, col="red")
  lines(dat_CImeasles$age_nodecim, dat_CImeasles$logit.lb, col="red", lty=2)
  lines(dat_CImeasles$age_nodecim, dat_CImeasles$logit.ub, col="red", lty=2)
  title("Measles (wMev)")
  
#plot4
  plot(x=tmprube$age_nodecim, y=tmprube$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIrube$age_nodecim, dat_CIrube$logit.mean, col="red")
  lines(dat_CIrube$age_nodecim, dat_CIrube$logit.lb, col="red", lty=2)
  lines(dat_CIrube$age_nodecim, dat_CIrube$logit.ub, col="red", lty=2)
  title("Rubella (wRuv)")  

#plot5
  plot(x=tmptet$age_nodecim, y=tmptet$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CItet$age_nodecim, dat_CItet$logit.mean, col="red")
  lines(dat_CItet$age_nodecim, dat_CItet$logit.lb, col="red", lty=2)
  lines(dat_CItet$age_nodecim, dat_CItet$logit.ub, col="red", lty=2)
  title("Tetanus (Tet tox)")  

#plot6
  plot(x=tmpdip$age_nodecim, y=tmpdip$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIdip$age_nodecim, dat_CIdip$logit.mean, col="red")
  lines(dat_CIdip$age_nodecim, dat_CIdip$logit.lb, col="red", lty=2)
  lines(dat_CIdip$age_nodecim, dat_CIdip$logit.ub, col="red", lty=2)
  title("Diptheria (Dip tox)")
  
#####SECOND set of plots  
plot.new()
  text(0.5,0.5,"NTDs",cex=2,font=2)  
  par(mar = c(1.5, 1.65, 1.25, 1))
  par(mfrow = c(5,3))
  
#plot22
  plot(x=tmpwb123$age_nodecim, y=tmpwb123$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIwb123$age_nodecim, dat_CIwb123$logit.mean, col="red")
  lines(dat_CIwb123$age_nodecim, dat_CIwb123$logit.lb, col="red", lty=2)
  lines(dat_CIwb123$age_nodecim, dat_CIwb123$logit.ub, col="red", lty=2)
  title("Wuchereria bancrofti (wb123)")

#plot23
  plot(x=tmptmpa$age_nodecim, y=tmptmpa$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CItmpa$age_nodecim, dat_CItmpa$logit.mean, col="red")
  lines(dat_CItmpa$age_nodecim, dat_CItmpa$logit.lb, col="red", lty=2)
  lines(dat_CItmpa$age_nodecim, dat_CItmpa$logit.ub, col="red", lty=2)
  title("Treponema palladium (tmpa)")
  
#plot9
  plot(x=tmpt24h$age_nodecim, y=tmpt24h$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIt24h$age_nodecim, dat_CIt24h$logit.mean, col="red")
  lines(dat_CIt24h$age_nodecim, dat_CIt24h$logit.lb, col="red", lty=2)
  lines(dat_CIt24h$age_nodecim, dat_CIt24h$logit.ub, col="red", lty=2)
  title("Taenia solium(t24h)")
  
#plot12
  plot(x=tmpsm25$age_nodecim, y=tmpsm25$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIsm25$age_nodecim, dat_CIsm25$logit.mean, col="red")
  lines(dat_CIsm25$age_nodecim, dat_CIsm25$logit.lb, col="red", lty=2)
  lines(dat_CIsm25$age_nodecim, dat_CIsm25$logit.ub, col="red", lty=2)
  title("Schistosoma mansoni (sm25)")
  
#plot11
  plot(x=tmpsea$age_nodecim, y=tmpsea$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIsea$age_nodecim, dat_CIsea$logit.mean, col="red")
  lines(dat_CIsea$age_nodecim, dat_CIsea$logit.lb, col="red", lty=2)
  lines(dat_CIsea$age_nodecim, dat_CIsea$logit.ub, col="red", lty=2)
  title("Schistosoma mansoni (sea)")
  
#plot15
  plot(x=tmprp17$age_nodecim, y=tmprp17$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIrp17$age_nodecim, dat_CIrp17$logit.mean, col="red")
  lines(dat_CIrp17$age_nodecim, dat_CIrp17$logit.lb, col="red", lty=2)
  lines(dat_CIrp17$age_nodecim, dat_CIrp17$logit.ub, col="red", lty=2)
  title("Treponema palladium (rp17)")
  
#plot13
  plot(x=tmppgp3$age_nodecim, y=tmppgp3$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIpgp3$age_nodecim, dat_CIpgp3$logit.mean, col="red")
  lines(dat_CIpgp3$age_nodecim, dat_CIpgp3$logit.lb, col="red", lty=2)
  lines(dat_CIpgp3$age_nodecim, dat_CIpgp3$logit.ub, col="red", lty=2)
  title("Chlamydia trachomatis (pgp3)")

#plot23
  plot(x=tmpov16$age_nodecim, y=tmpov16$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIov16$age_nodecim, dat_CIov16$logit.mean, col="red")
  lines(dat_CIov16$age_nodecim, dat_CIov16$logit.lb, col="red", lty=2)
  lines(dat_CIov16$age_nodecim, dat_CIov16$logit.ub, col="red", lty=2)
  title("Onchocerca volvulus (ov16)")

#plot24
  plot(x=tmpnie$age_nodecim, y=tmpnie$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CInie$age_nodecim, dat_CInie$logit.mean, col="red")
  lines(dat_CInie$age_nodecim, dat_CInie$logit.lb, col="red", lty=2)
  lines(dat_CInie$age_nodecim, dat_CInie$logit.ub, col="red", lty=2)
  title("Strongyloides stercoralis (nie)")
  
#plot25
  plot(x=tmpes33$age_nodecim, y=tmpes33$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIes33$age_nodecim, dat_CIes33$logit.mean, col="red")
  lines(dat_CIes33$age_nodecim, dat_CIes33$logit.lb, col="red", lty=2)
  lines(dat_CIes33$age_nodecim, dat_CIes33$logit.ub, col="red", lty=2)
  title("Taenia solium(es33)")
  
#plot14
  plot(x=tmpct694$age_nodecim, y=tmpct694$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIct694$age_nodecim, dat_CIct694$logit.mean, col="red")
  lines(dat_CIct694$age_nodecim, dat_CIct694$logit.lb, col="red", lty=2)
  lines(dat_CIct694$age_nodecim, dat_CIct694$logit.ub, col="red", lty=2)
  title("Chlamydia trachomatis (ct694)")
  
#plot10
  plot(x=tmpbm33$age_nodecim, y=tmpbm33$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIbm33$age_nodecim, dat_CIbm33$logit.mean, col="red")
  lines(dat_CIbm33$age_nodecim, dat_CIbm33$logit.lb, col="red", lty=2)
  lines(dat_CIbm33$age_nodecim, dat_CIbm33$logit.ub, col="red", lty=2)
  title("Brugia malayi (bm33)")

#plot10
  plot(x=tmpbm14$age_nodecim, y=tmpbm14$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIbm14$age_nodecim, dat_CIbm14$logit.mean, col="red")
  lines(dat_CIbm14$age_nodecim, dat_CIbm14$logit.lb, col="red", lty=2)
  lines(dat_CIbm14$age_nodecim, dat_CIbm14$logit.ub, col="red", lty=2)
  title("Brugia malayi (bm14)")  

  
#####THIRD set of plots  
plot.new()
  text(0.5,0.5,"malaria",cex=2,font=2)  
  par(mar = c(1.5, 1.65, 1.25, 1))
  par(mfrow = c(4,3))

#plot26
  plot(x=tmprh42$age_nodecim, y=tmprh42$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1)) 
  lines(dat_CIrh42$age_nodecim, dat_CIrh42$logit.mean, col="red")
  lines(dat_CIrh42$age_nodecim, dat_CIrh42$logit.lb, col="red", lty=2)
  lines(dat_CIrh42$age_nodecim, dat_CIrh42$logit.ub, col="red", lty=2)
  title("P. falciparum(rh42)")

#plot27
  plot(x=tmppvrbp$age_nodecim, y=tmppvrbp$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1)) 
#  lines(dat_CIpvrbp$age_nodecim, dat_CIpvrbp$logit.mean, col="red")
#  lines(dat_CIpvrbp$age_nodecim, dat_CIpvrbp$logit.lb, col="red", lty=2)
#  lines(dat_CIpvrbp$age_nodecim, dat_CIpvrbp$logit.ub, col="red", lty=2)
  title("P. vivax(pvrbp2b)*")

#plot28
  plot(x=tmppvmsp$age_nodecim, y=tmppvmsp$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1)) 
  lines(dat_CIpvmsp$age_nodecim, dat_CIpvmsp$logit.mean, col="red")
  lines(dat_CIpvmsp$age_nodecim, dat_CIpvmsp$logit.lb, col="red", lty=2)
  lines(dat_CIpvmsp$age_nodecim, dat_CIpvmsp$logit.ub, col="red", lty=2)
  title("P. vivax(pvmsp119)")

#plot29
  plot(x=tmppvdbpr$age_nodecim, y=tmppvdbpr$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1)) 
  lines(dat_CIpvdbpr$age_nodecim, dat_CIpvdbpr$logit.mean, col="red")
  lines(dat_CIpvdbpr$age_nodecim, dat_CIpvdbpr$logit.lb, col="red", lty=2)
  lines(dat_CIpvdbpr$age_nodecim, dat_CIpvdbpr$logit.ub, col="red", lty=2)
  title("P. vivax(pvdbprii)")

#plot30
  plot(x=tmppomsp$age_nodecim, y=tmppomsp$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1)) 
  lines(dat_CIpomsp$age_nodecim, dat_CIpomsp$logit.mean, col="red")
  lines(dat_CIpomsp$age_nodecim, dat_CIpomsp$logit.lb, col="red", lty=2)
  lines(dat_CIpomsp$age_nodecim, dat_CIpomsp$logit.ub, col="red", lty=2)
  title("P. ovale(pomsp119)")
  
#plot31
  plot(x=tmppmmsp$age_nodecim, y=tmppmmsp$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1)) 
  lines(dat_CIpmmsp$age_nodecim, dat_CIpmmsp$logit.mean, col="red")
  lines(dat_CIpmmsp$age_nodecim, dat_CIpmmsp$logit.lb, col="red", lty=2)
  lines(dat_CIpmmsp$age_nodecim, dat_CIpmmsp$logit.ub, col="red", lty=2)
  title("P. malariae(pmmsp119)")
  
#plot32
  plot(x=tmppfmsp$age_nodecim, y=tmppfmsp$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1)) 
  lines(dat_CIpfmsp$age_nodecim, dat_CIpfmsp$logit.mean, col="red")
  lines(dat_CIpfmsp$age_nodecim, dat_CIpfmsp$logit.lb, col="red", lty=2)
  lines(dat_CIpfmsp$age_nodecim, dat_CIpfmsp$logit.ub, col="red", lty=2)
  title("P. falciparum(pfmsp119)")
  
#plot1
  plot(x=tmppfama1$age_nodecim, y=tmppfama1$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1)) 
  lines(dat_CIpfama1$age_nodecim, dat_CIpfama1$logit.mean, col="red")
  lines(dat_CIpfama1$age_nodecim, dat_CIpfama1$logit.lb, col="red", lty=2)
  lines(dat_CIpfama1$age_nodecim, dat_CIpfama1$logit.ub, col="red", lty=2)
  title("P. falciparum(pfama1)")
  
#plot2
  plot(x=tmpglurpr2$age_nodecim, y=tmpglurpr2$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIglurpr2$age_nodecim, dat_CIglurpr2$logit.mean, col="red")
  lines(dat_CIglurpr2$age_nodecim, dat_CIglurpr2$logit.lb, col="red", lty=2)
  lines(dat_CIglurpr2$age_nodecim, dat_CIglurpr2$logit.ub, col="red", lty=2)
  title("P. falciparum(glurpr2)")
  
#plot33
  plot(x=tmpgexp$age_nodecim, y=tmpgexp$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIgexp$age_nodecim, dat_CIgexp$logit.mean, col="red")
  lines(dat_CIgexp$age_nodecim, dat_CIgexp$logit.lb, col="red", lty=2)
  lines(dat_CIgexp$age_nodecim, dat_CIgexp$logit.ub, col="red", lty=2)
  title("P. falciparum(gexp18)")
  
#plot34
  plot(x=tmpetramp$age_nodecim, y=tmpetramp$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIetramp$age_nodecim, dat_CIetramp$logit.mean, col="red")
  lines(dat_CIetramp$age_nodecim, dat_CIetramp$logit.lb, col="red", lty=2)
  lines(dat_CIetramp$age_nodecim, dat_CIetramp$logit.ub, col="red", lty=2)
  title("P. falciparum(etramp5ag1)")
  
#plot35
  plot(x=tmpcsp$age_nodecim, y=tmpcsp$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIcsp$age_nodecim, dat_CIcsp$logit.mean, col="red")
  lines(dat_CIcsp$age_nodecim, dat_CIcsp$logit.lb, col="red", lty=2)
  lines(dat_CIcsp$age_nodecim, dat_CIcsp$logit.ub, col="red", lty=2)
  title("P. falciparum(csp)")

  
#################################  
#######ASTMH PPT#################  
#####  measles, rp17,pgp3, glurpr2, csp, pfmsp
  par(pin = c(6, 4))
  par(mar = c(1.8, 1.8, 1.8, 1.8))
  par(mfrow = c(2,3))  
  
  plot(x=tmpmeasles$age_nodecim, y=tmpmeasles$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CImeasles$age_nodecim, dat_CImeasles$logit.mean, col="red")
  lines(dat_CImeasles$age_nodecim, dat_CImeasles$logit.lb, col="red", lty=2)
  lines(dat_CImeasles$age_nodecim, dat_CImeasles$logit.ub, col="red", lty=2)
  title("Measles (wMev)")
  
  plot(x=tmprp17$age_nodecim, y=tmprp17$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIrp17$age_nodecim, dat_CIrp17$logit.mean, col="red")
  lines(dat_CIrp17$age_nodecim, dat_CIrp17$logit.lb, col="red", lty=2)
  lines(dat_CIrp17$age_nodecim, dat_CIrp17$logit.ub, col="red", lty=2)
  title("Treponema palladium (rp17)")
  
  plot(x=tmppgp3$age_nodecim, y=tmppgp3$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIpgp3$age_nodecim, dat_CIpgp3$logit.mean, col="red")
  lines(dat_CIpgp3$age_nodecim, dat_CIpgp3$logit.lb, col="red", lty=2)
  lines(dat_CIpgp3$age_nodecim, dat_CIpgp3$logit.ub, col="red", lty=2)
  title("Chlamydia trachomatis (pgp3)")
  
  plot(x=tmppfmsp$age_nodecim, y=tmppfmsp$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1)) 
  lines(dat_CIpfmsp$age_nodecim, dat_CIpfmsp$logit.mean, col="red")
  lines(dat_CIpfmsp$age_nodecim, dat_CIpfmsp$logit.lb, col="red", lty=2)
  lines(dat_CIpfmsp$age_nodecim, dat_CIpfmsp$logit.ub, col="red", lty=2)
  title("P. falciparum(pfmsp119)")
  
  plot(x=tmpglurpr2$age_nodecim, y=tmpglurpr2$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIglurpr2$age_nodecim, dat_CIglurpr2$logit.mean, col="red")
  lines(dat_CIglurpr2$age_nodecim, dat_CIglurpr2$logit.lb, col="red", lty=2)
  lines(dat_CIglurpr2$age_nodecim, dat_CIglurpr2$logit.ub, col="red", lty=2)
  title("P. falciparum(glurpr2)")
  
  plot(x=tmpcsp$age_nodecim, y=tmpcsp$seroprev, 
       cex=0.5, pch=16, xlab="age", ylab="seroprevalence", xlim=c(0,50), ylim=c(-0.1,1))
  lines(dat_CIcsp$age_nodecim, dat_CIcsp$logit.mean, col="red")
  lines(dat_CIcsp$age_nodecim, dat_CIcsp$logit.lb, col="red", lty=2)
  lines(dat_CIcsp$age_nodecim, dat_CIcsp$logit.ub, col="red", lty=2)
  title("P. falciparum(csp)")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#####adjusted seroprevalence for low values#####
library(boot)
library(bbmle)

## Likelihood function
  log_lik <- function(params, y, sample_size, Se, Sp) {
    prev <- params[1]
    p_sample <- prev*Se + (1-prev)*(1-Sp)
    log_lik <- y*log(p_sample) + (sample_size-y)*log(1-p_sample)
    return(-log_lik)
  }

## Data
  y <- 20 ## Number of positives
  sample_size <- 1292
  Se <- 0.909
  Sp <- 0.947

## Estimate seroprevalence
  m0 <- mle2(log_lik, 
             start=list(params=0.5),
             lower=0,
             upper=1,
             method="L-BFGS-B",
             optimizer="nlminb",
             data=list(y=y, sample_size=sample_size, Se=Se, Sp=Sp))
  m0
  
## Point estimate of seroprevalence
  coef(m0)

## 95% confidence interval
  confint(m0)