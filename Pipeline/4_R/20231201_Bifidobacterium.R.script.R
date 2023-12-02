#######################################################
#2020/11/04
#B. animalis subsp. lactis
#OD
#Grown in anaerobic tubes, Jake's method
require("ggplot2");require("tidyverse");require("png");require("dplyr");require("grid");require("vegan");require("indicspecies");require("ggrepel");require("gridExtra")

#require("ggimage");require("grid")
#require("cowplot");
#These packages break my GPU if I try to run them twice.
#install.packages("magick")
#library(magick)
#####Don't run library(magick), it broke my Rstudio instance on Windows.

mycols=c('black','black','black')
mycols_breakout=c('black','black','black','black','black','black','black','black')
#Import ancestor data
#Import basic OD data on the 3 substrates
mydf0<-read_csv("~\\..\\OneDrive - Indiana University\\Lab.Notebook\\20190820_Bifidobacterium\\data\\20201104_OD.final_xevo_T0\\20201104_B.animalis.McK.mthd.csv")
#mydf0$carbon2<-factor(mydf0$carbon2, levels=c("F","G-F5","G-F24"))
#mydf0$carbon2<-recode_factor(mydf0$carbon2, `1` = "Fructose", `2` = "DP5", `3` = "DP24")
mydf0$carbon2<-recode(mydf0$carbon2, "F" = "DP0", "G-F5" = "DP5", "G-F24" = "DP24")
mydf0$carbon2<-factor(mydf0$carbon2, levels=c("DP0","DP5","DP24"))
mydf0$carbon<-factor(mydf0$carbon, levels=c("Monomer","Low-DP","High-DP"))

#Import OD data from nutrient breakout xpmnt
mybreakout<-read_csv("~/../OneDrive - Indiana University/Lab.Notebook/20190820_Bifidobacterium/data/nutrient.breakouts/20210512_OD.final_nutrient.breakout.ancor/20210512_B.animalis.McK.mthd.csv")
mybreakout$carbon.tech<-factor(mybreakout$carbon.tech, levels=c("G", "F", "G:F 1:1","G-F","G:F 1:5", "G-F5","G:F 1:24","G-F24" ))
#mybreakout$carbon.plotting<-recode(mybreakout$carbon$tech)
#Import OD Data from the nutrient breakout xpmnt, formatted as a 2way ANOVA: Factor 1: G-F ratio Factor 2: are the Gs and Fs monomers or in fiber molecules?
mybreakout_2way<-read_csv("~/../OneDrive - Indiana University/Lab.Notebook/20190820_Bifidobacterium/data/nutrient.breakouts/20210512_OD.final_nutrient.breakout.ancor/20210512_B.animalis.McK.mthd_2way.enabled.csv")
mybreakout_2way$GF.ratio<-factor(mybreakout_2way$GF.ratio, levels=c("1:1","1:5","1:24"))
mybreakout_2way$Saccharide.size<-factor(mybreakout_2way$Saccharide.size, levels=c("Monomer","Oligomer"))




######
#Import data from HPLC experiments
myhplc<-read_csv("~/../OneDrive - Indiana University/Lab.Notebook/20190820_Bifidobacterium/data/HPLC_Bifido/20211006_HPLC_glu.fru.lact.ace_v3.csv")
myhplc$condition_short<-factor(myhplc$condition_short,levels=c("GF","F"))
myhplc$evo_trt<-factor(myhplc$evo_trt,levels=c("F","P","X","Anc"))
myhplc$evo_trt_plotting<-recode(myhplc$evo_trt, "F" = "DP0", "P" = "DP5", "X" = "DP24", "Anc" = "Ancestor")
myhplc$evo_trt_plotting<-factor(myhplc$evo_trt_plotting,levels=c("DP0","DP5","DP24","Ancestor"))



###########Statistical test on ancestor data#################
#Does the ancestor show different growth patterns on the different carbon sources?
tubeaov0<-aov(mydf0$OD.600 ~ mydf0$carbon)
summary(tubeaov0)
TukeyHSD(tubeaov0)
#All treatments are significantly different.
#Tukey letters: A, B, C




###########Import evolved data from timepoint T151, or 1003 generations of evolution.#############
###############Import basic OD data on the 3 substrates############
mydf151<-read_csv("~\\..\\OneDrive - Indiana University\\Lab.Notebook\\20190820_Bifidobacterium\\data\\20210509_OD.final_xevo_T151_FINAL\\20210509_B.animalis.McK.mthd.csv")
#mydf151$carbon2<-factor(mydf151$carbon2, levels=c("F","G-F5","G-F24"))
mydf151$carbon<-factor(mydf151$carbon, levels=c("Monomer","Low-DP","High-DP"))
mydf151$carbon2<-recode(mydf151$carbon2, "F" = "DP0", "G-F5" = "DP5", "G-F24" = "DP24")
mydf151$carbon2<-factor(mydf151$carbon2, levels=c("DP0","DP5","DP24"))

#Do a statistical test: Grown on their "home" diets, do the evolved populations differ in their max OD on each diet?
tubeaov151<-aov(mydf151$OD.600 ~ mydf151$carbon)
summary(tubeaov151)
TukeyHSD(tubeaov151)
#F and GF7 are significantly higher than GF23. F and GF5 not significantly different from each other.
#Tukey letters: A, A, B
#Interpretation: There is no longer any difference in growth on F vs. Low-DP

#####################Do the same kind of ANOVA but for comparing absolute increase vs. ancestor####################
tubeaov151_abs.change<-aov(mydf151$abs.change.vs.anc ~ mydf151$carbon2)
summary(tubeaov151_abs.change)
TukeyHSD(tubeaov151_abs.change)
###Significant ANOVA
###Tukey letters: A, B, C

#####################Do the same kind of ANOVA but for comparing the percent increase vs. ancestor##################
tubeaov151_pct.increase<-aov(mydf151$percent.increase.vs.anc ~ mydf151$carbon2)
summary(tubeaov151_pct.increase)
TukeyHSD(tubeaov151_pct.increase)
###Significant ANOVA
###Tukey letters: A, B, B

#####
#######################Import OD data from nutrient breakout xpmnts with EVOLVED bacteria##############################
mybreakoutF<-read_csv("~/../OneDrive - Indiana University/Lab.Notebook/20190820_Bifidobacterium/data/nutrient.breakouts/20210520_OD.final_nutrient.breakout.F-evolved/20210520_nutrient.breakout.F.csv")
mybreakoutF$carbon.tech<-factor(mybreakoutF$carbon.tech, levels=c("G", "F", "G:F 1:1","G-F","G:F 1:5", "G-F5","G:F 1:24","G-F24" ))

mybreakoutP<-read_csv("~/../OneDrive - Indiana University/Lab.Notebook/20190820_Bifidobacterium/data/nutrient.breakouts/20210526_OD.final_nutrient.breakout.P-evolved/20210526_nutrient.breakout.P.csv")
mybreakoutP$carbon.tech<-factor(mybreakoutP$carbon.tech, levels=c("G", "F", "G:F 1:1","G-F","G:F 1:5", "G-F5","G:F 1:24","G-F24" ))

mybreakoutX<-read_csv("~/../OneDrive - Indiana University/Lab.Notebook/20190820_Bifidobacterium/data/nutrient.breakouts/20210531_OD.final_nutrient.breakout.X-evolved/20210531_nutrient.breakout.X.csv")
mybreakoutX$carbon.tech<-factor(mybreakoutX$carbon.tech, levels=c("G", "F", "G:F 1:1","G-F","G:F 1:5", "G-F5","G:F 1:24","G-F24" ))


#######################Import anaerobic growth curve data##############################################################
grcvdata<-read_csv("~/../OneDrive - Indiana University/Lab.Notebook/20190820_Bifidobacterium/data/anaerobic.growth.curves/grcv_data_BIfidobacterium_cases.format.csv")
grcvdata$carbon.tech.assay<-factor(grcvdata$carbon.tech.assay, levels=c("F", "G-F5","G-F24","Ancestor"))
grcvdata$carbon.tech.evo<-factor(grcvdata$carbon.tech.evo, levels=c("F", "G-F5","G-F24","Ancestor"))
grcvdata$carbon.plotting.assay<-recode(grcvdata$carbon.tech.assay, "F" = "DP0", "G-F5" = "DP5", "G-F24" = "DP24")
grcvdata$carbon.plotting.evo<-recode(grcvdata$carbon.tech.evo, "F" = "DP0", "G-F5" = "DP5", "G-F24" = "DP24")


###################Plots####################
#fig
#Make a plot of the ancestor data: OD when grown on the 3 basic substrates
odtube0 <- ggplot(mydf0, aes(x=carbon2, y=OD.600))
odtube0 +
  geom_jitter(
  aes(shape = carbon2, color = carbon2), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon2),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols) +
  scale_y_continuous(limits = c(-0.04,.9),breaks=c(0,.4,.8), expand = c(0,0)) +
  labs(x="\nDiet",y="Biomass (OD600)\n") +
  #geom_text("B",hjust=0,vjust=0) +
  annotate("text",x=.7,y=0.39,label="B",size=12,color="black")+
  annotate("text",x=1.7,y=0.8,label="A",size=12,color="black")+
  annotate("text",x=2.7,y=0.31,label="C",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

#Make the same plot but with an axis that will also fit the evolved bacteria data
odtube0_sameaxis <- ggplot(mydf0, aes(x=carbon2, y=OD.600))
odtube0_sameaxis + geom_jitter(
  aes(shape = carbon2, color = carbon2), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon2),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols) +
  scale_y_continuous(limits = c(-0.04,1.6),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  labs(x="\nCarbon source",y="OD 600\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))



###
#Combine the OD plots. This plot will display the OD of evolved bacteria in their "home" environments and compare it to the ancestral growth phenotype in each of these environments.
#fig
FOD<-mean(c(filter(mydf0, carbon=="Monomer")$OD.600))
POD<-mean(c(filter(mydf0, carbon=="Low-DP")$OD.600))
XOD<-mean(c(filter(mydf0, carbon=="High-DP")$OD.600))

anc.OD<-c(FOD,FOD,FOD,FOD,FOD,FOD,FOD,FOD,POD,POD,POD,POD,POD,POD,POD,POD,XOD,XOD,XOD,XOD,XOD,XOD,XOD,XOD)
anc.colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")

odtube151 <- ggplot(mydf151, aes(x=carbon2, y=OD.600))
odtube151 +geom_errorbar(aes(ymin=anc.OD,ymax=anc.OD),color=anc.colors,lwd = 2.25,linetype=117) +
  geom_jitter(
  aes(shape = carbon2, color = carbon2), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon2),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols) +
  scale_y_continuous(limits = c(-0.04,1.6),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  labs(x="\nDiet",y="Biomass (OD600)\n") +
  annotate("text",x=.6,y=1.4,label="A",size=12,color="black")+
  annotate("text",x=1.6,y=1.4,label="A",size=12,color="black")+
  annotate("text",x=2.6,y=0.55,label="B",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))




#########################################################Plot the evolved data as Absolute Increase Vs. Ancestor
#fig
odtube151_abs.change <- ggplot(mydf151, aes(x=carbon2, y=abs.change.vs.anc))
odtube151_abs.change +#geom_errorbar(aes(ymin=anc.OD,ymax=anc.OD),color=anc.colors,lwd = 2.25,linetype=117) +
  geom_jitter(
    aes(shape = carbon2, color = carbon2), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
    position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
    size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon2),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols) +
  scale_y_continuous(limits = c(-0.04,1.25),breaks=c(0,.4,.8,1.2), expand = c(0,0)) +
  labs(x="\nDiet",y="Evolved increase in OD600\n") +
  annotate("text",x=.6,y=1.1,label="A",size=12,color="black")+
  annotate("text",x=1.6,y=0.65,label="B",size=12,color="black")+
  annotate("text",x=2.6,y=0.35,label="C",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


#########################################################Plot the evolved data as Percent Increase Vs. Ancestor
odtube151_pct.change <- ggplot(mydf151, aes(x=carbon2, y=percent.increase.vs.anc))
odtube151_pct.change +#geom_errorbar(aes(ymin=anc.OD,ymax=anc.OD),color=anc.colors,lwd = 2.25,linetype=117) +
  geom_jitter(
    aes(shape = carbon2, color = carbon2), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
    position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
    size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon2),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols) +
  scale_y_continuous(limits = c(-4,416),breaks=c(0,100,200,300,400), expand = c(0,0)) +
  labs(x="\nDiet",y="Evolved % increase in OD 600\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))




##########################t-tests with OD600##################################
#For each carbon source, I could do two-sample t-test to ask whether growth increased after evolution.
#This is misleading, though, because there are two different levels of replication.
#We should really be doing one-sample t-tests and taking the ancestral value as fixed.

################################################################################
#Below, I do the two-sample t-tests

#Did growth on fructose increase during 866 generations of natural selection?
var.test(x=c(filter(mydf0, carbon=="Monomer")$OD.600),y=c(filter(mydf151, carbon == "Monomer")$OD.600),ratio=1,alternative = "t")
#P = .3316. Homoscedastic.
t.test(x=c(filter(mydf0, carbon=="Monomer")$OD.600),y=c(filter(mydf151, carbon == "Monomer")$OD.600),mu=0,alternative = "l",paired=F,var.equal=T)
#p = 4.543EE-16 < 0.0001

#Did growth on P95 inulin (DP = 2-8) increase during 866 generations of natural selection?
var.test(x=c(filter(mydf0, carbon=="Low-DP")$OD.600),y=c(filter(mydf151, carbon == "Low-DP")$OD.600),ratio=1,alternative = "t")
#P = .1358. Homoscedastic.
t.test(x=c(filter(mydf0, carbon=="Low-DP")$OD.600),y=c(filter(mydf151, carbon == "Low-DP")$OD.600),mu=0,alternative = "l",paired=F,var.equal=T)
#p = 5.010EE-10 < 0.0001

#Did growth on TEX inulin (DP > 23) increase during 866 generations of natural selection?
var.test(x=c(filter(mydf0, carbon=="High-DP")$OD.600),y=c(filter(mydf151, carbon == "High-DP")$OD.600),ratio=1,alternative = "t")
#P = .04689. Heteroscedastic.
t.test(x=c(filter(mydf0, carbon=="High-DP")$OD.600),y=c(filter(mydf151, carbon == "High-DP")$OD.600),mu=0,alternative = "l",paired=F,var.equal=F)
#p = 0.0004431 < 0.001

mypvec151<-c(t.test(x=c(filter(mydf0, carbon=="Monomer")$OD.600),y=c(filter(mydf151, carbon == "Monomer")$OD.600),mu=0,alternative = "l",paired=F,var.equal=T)$p.value, t.test(x=c(filter(mydf0, carbon=="Low-DP")$OD.600),y=c(filter(mydf151, carbon == "Low-DP")$OD.600),mu=0,alternative = "l",paired=F,var.equal=T)$p.value, t.test(x=c(filter(mydf0, carbon=="High-DP")$OD.600),y=c(filter(mydf151, carbon == "High-DP")$OD.600),mu=0,alternative = "l",paired=F,var.equal=F)$p.value)

mypvec151.adjust<-p.adjust(p=mypvec151, method = "BH")#cite BH&Y 2009
mypvec151.adjust

#################################################################################
#Now, I will do the one-sample t-tests
#Did growth on fructose increase during 1000 generations of natural selection?
shapiro.test(c(filter(mydf151, carbon == "Monomer")$OD.600))
#P = .4783. Data are not significantly non-Gaussian.
t.test(x=c(filter(mydf151, carbon == "Monomer")$OD.600),mu=mean(c(filter(mydf0, carbon=="Monomer")$OD.600)),alternative = 'g')
fig2p1<-t.test(x=c(filter(mydf151, carbon == "Monomer")$OD.600),mu=mean(c(filter(mydf0, carbon=="Monomer")$OD.600)),alternative = 'g')$p.value
#p < 0.0001
#Growth increased by 1.02 OD units, an increase of 310%


#Did growth on low-DP inulin increase during 1000 generations of natural selection?
shapiro.test(c(filter(mydf151, carbon == "Low-DP")$OD.600))
#P = .5522. Data are not significantly non-Gaussian.
t.test(x=c(filter(mydf151, carbon == "Low-DP")$OD.600),mu=mean(c(filter(mydf0, carbon=="Low-DP")$OD.600)),alternative = 'g')
fig2p2<-t.test(x=c(filter(mydf151, carbon == "Low-DP")$OD.600),mu=mean(c(filter(mydf0, carbon=="Low-DP")$OD.600)),alternative = 'g')$p.value
#p < 0.0001
#Growth increased by 0.57 OD units, an increase of 74%

#Did growth on high-DP inulin increase during 1000 generations of natural selection?
shapiro.test(c(filter(mydf151, carbon == "High-DP")$OD.600))
#P = .1954. Data are not significantly non-Gaussian.
t.test(x=c(filter(mydf151, carbon == "High-DP")$OD.600),mu=mean(c(filter(mydf0, carbon=="High-DP")$OD.600)),alternative = 'g')
fig2p3<-t.test(x=c(filter(mydf151, carbon == "High-DP")$OD.600),mu=mean(c(filter(mydf0, carbon=="High-DP")$OD.600)),alternative = 'g')$p.value
#p  < 0.001
#Growth increased by 0.23 OD units, an increase of 91%

#fig2pvec<-c(fig2p1,fig2p2,fig2p3)
#fig2pvec

mypvec151.onesample<-c(t.test(x=c(filter(mydf151, carbon == "Monomer")$OD.600),mu=mean(c(filter(mydf0, carbon=="Monomer")$OD.600)),alternative = 'g')$p.value,t.test(x=c(filter(mydf151, carbon == "Low-DP")$OD.600),mu=mean(c(filter(mydf0, carbon=="Low-DP")$OD.600)),alternative = 'g')$p.value,t.test(x=c(filter(mydf151, carbon == "High-DP")$OD.600),mu=mean(c(filter(mydf0, carbon=="High-DP")$OD.600)),alternative = 'g')$p.value)
mypvec151.onesample
mypvec151.onesample.adjust<-p.adjust(p=mypvec151.onesample, method = "BH")#cite BH&Y 2009
mypvec151.onesample.adjust

#################################################################################
####I compared changes ("deltas") of growth on each evolved home diet compared to the ancestor, both as absolute values and as percentages. But we can do ANOVA on the absolute deltas to statistically quantify which evolvd populations had the greatest increase in growth
home.delta.aov<-aov(mydf151$OD.600.delta ~ mydf151$carbon)
summary(home.delta.aov)
TukeyHSD(home.delta.aov)

########################################################################################
##########################Nutrient breakout experiment#################################


####Within evolution treatments (i.e., anc, F evolved, P evolved, X evolved), assess whether growth is different on the 8 different nutrient compositions.

breakoutaov<-aov(mybreakout$OD.600 ~ mybreakout$carbon.tech)
summary(breakoutaov)
TukeyHSD(breakoutaov)


#breakout2way<-aov(mybreakout_2way$OD.600 ~ mybreakout_2way$Saccharide.size + mybreakout_2way$GF.ratio + mybreakout_2way$Saccharide.size*mybreakout_2way$GF.ratio)
#summary(breakout2way)
#TukeyHSD(breakout2way)


##It looks like both the G:F ratio and the molecule size "matter"
#But we can quantify this with multiple regression

#breakoutmodel<-lm(formula = mybreakout$OD.600 ~ mybreakout$molecule.size + mybreakout$proportion.F + mybreakout$molecule.size*mybreakout$proportion.F)
breakoutmodel<-lm(formula = mybreakout$OD.600 ~ mybreakout$molecule.size + mybreakout$proportion.F)
summary(breakoutmodel)


breakoutaovF<-aov(mybreakoutF$OD.600 ~ mybreakoutF$carbon.tech)
summary(breakoutaovF)
TukeyHSD(breakoutaovF)

breakoutaovP<-aov(mybreakoutP$OD.600 ~ mybreakoutP$carbon.tech)
summary(breakoutaovP)
TukeyHSD(breakoutaovP)

breakoutaovX<-aov(mybreakoutX$OD.600 ~ mybreakoutX$carbon.tech)
summary(breakoutaovX)
TukeyHSD(breakoutaovX)

mymonomers<-filter(mybreakout,molecule.size==1)
mymonomersF<-filter(mybreakoutF,molecule.size==1)
mymonomersP<-filter(mybreakoutP,molecule.size==1)
mymonomersX<-filter(mybreakoutX,molecule.size==1)

monomersaov<-aov(mymonomers$OD.600 ~ mymonomers$carbon.tech)
summary(monomersaov)
TukeyHSD(monomersaov)

monomersaovF<-aov(mymonomersF$OD.600 ~ mymonomersF$carbon.tech)
summary(monomersaovF)
TukeyHSD(monomersaovF)

monomersaovP<-aov(mymonomersP$OD.600 ~ mymonomersP$carbon.tech)
summary(monomersaovP)
TukeyHSD(monomersaovP)

monomersaovX<-aov(mymonomersX$OD.600 ~ mymonomersX$carbon.tech)
summary(monomersaovX)
TukeyHSD(monomersaovX)

##It looks like both the G:F ratio and the molecule size "matter"
#But we can quantify this with multiple regression
breakoutmodelF<-lm(formula = mybreakoutF$OD.600 ~ mybreakoutF$molecule.size + mybreakoutF$proportion.F + mybreakoutF$molecule.size*mybreakoutF$proportion.F)
summary(breakoutmodelF)

breakoutmodelP<-lm(formula = mybreakoutP$OD.600 ~ mybreakoutP$molecule.size + mybreakoutP$proportion.F)
summary(breakoutmodelP)

breakoutmodelX<-lm(formula = mybreakoutX$OD.600 ~ mybreakoutX$molecule.size + mybreakoutX$proportion.F)
summary(breakoutmodelX)

monomermodel<-lm(formula = mymonomers$OD.600 ~ mymonomers$proportion.F)
summary(monomermodel)

monomermodelF<-lm(formula = mymonomersF$OD.600 ~ mymonomersF$proportion.F)
summary(monomermodelF)

monomermodelP<-lm(formula = mymonomersP$OD.600 ~ mymonomersP$proportion.F)
summary(monomermodelP)

monomermodelX<-lm(formula = mymonomersX$OD.600 ~ mymonomersX$proportion.F)
summary(monomermodelX)

mymonomervector<-c(mean(filter(mymonomers,carbon.tech=="F")$OD.600),mean(filter(mymonomers,carbon.tech=="G:F 1:24")$OD.600),mean(filter(mymonomers,carbon.tech=="G:F 1:5")$OD.600),mean(filter(mymonomers,carbon.tech=="G:F 1:1")$OD.600),mean(filter(mymonomers,carbon.tech=="G")$OD.600))

mymonomervectorF<-c(mean(filter(mymonomersF,carbon.tech=="F")$OD.600),mean(filter(mymonomersF,carbon.tech=="G:F 1:24")$OD.600),mean(filter(mymonomersF,carbon.tech=="G:F 1:5")$OD.600),mean(filter(mymonomersF,carbon.tech=="G:F 1:1")$OD.600),mean(filter(mymonomersF,carbon.tech=="G")$OD.600))

mymonomervectorP<-c(mean(filter(mymonomersP,carbon.tech=="F")$OD.600),mean(filter(mymonomersP,carbon.tech=="G:F 1:24")$OD.600),mean(filter(mymonomersP,carbon.tech=="G:F 1:5")$OD.600),mean(filter(mymonomersP,carbon.tech=="G:F 1:1")$OD.600),mean(filter(mymonomersP,carbon.tech=="G")$OD.600))

mymonomervectorX<-c(mean(filter(mymonomersX,carbon.tech=="F")$OD.600),mean(filter(mymonomersX,carbon.tech=="G:F 1:24")$OD.600),mean(filter(mymonomersX,carbon.tech=="G:F 1:5")$OD.600),mean(filter(mymonomersX,carbon.tech=="G:F 1:1")$OD.600),mean(filter(mymonomersX,carbon.tech=="G")$OD.600))

mymonomerFprops<-c(1,24/25,5/6,1/2,0)

monomervectormodel<-lm(formula = mymonomervector ~ mymonomerFprops)
summary(monomervectormodel)

monomervectormodelF<-lm(formula = mymonomervectorF ~ mymonomerFprops)
summary(monomervectormodelF)

monomervectormodelP<-lm(formula = mymonomervectorP ~ mymonomerFprops)
summary(monomervectormodelP)

monomervectormodelX<-lm(formula = mymonomervectorX ~ mymonomerFprops)
summary(monomervectormodelX)


#1
#Did F-evolved bacteria have higher growth on DP5 than the ancestor?
t1bmu<-mean(filter(mybreakout, carbon.tech == "G-F5")$OD.600)
t.test(c(filter(mybreakoutF, carbon.tech=="G-F5")$OD.600),mu=t1bmu,alternative = 'g')
t1bp<-t.test(c(filter(mybreakoutF, carbon.tech=="G-F5")$OD.600),mu=t1bmu,alternative = 'g')$p.value

#2
#Did F-evolved bacteria have higher growth on DP24 than the ancestor?
t2bmu<-mean(filter(mybreakout, carbon.tech == "G-F24")$OD.600)
t.test(filter(mybreakoutF, carbon.tech=="G-F24")$OD.600,mu=t2bmu,alternative = 'g')
t2bp<-t.test(filter(mybreakoutF, carbon.tech=="G-F24")$OD.600,mu=t2bmu,alternative = 'g')$p.value

#For the F-evolved, I need to correct for multiple testing.
bFpvec<-c(t1bp,t2bp)
bFpvec.adjust<-p.adjust(p=bFpvec, method = "BH")#cite BH&Y 2009
bFpvec.adjust
#FDR: 0.0362450269 0.0003116769



#3, #4, #5
#I need to test whether DP5 and DP24 evolved populations increased growth on fructose. Also, was DP24 growth significantly more than DP5 growth?
#3: DP5 vs. ancestor, growth on F
t3bmu<-mean(filter(mybreakout, carbon.tech == "F")$OD.600)
t.test(x=c(filter(mybreakoutP, carbon.tech=="F")$OD.600),mu=t3bmu,alternative = 'g')
t3bp<-t.test(x=c(filter(mybreakoutP, carbon.tech=="F")$OD.600),mu=t3bmu,alternative = 'g')$p.value

#4: DP24 vs. ancor, growth on F
t4bmu<-mean(filter(mybreakout, carbon.tech == "F")$OD.600)
t.test(filter(mybreakoutX, carbon.tech=="F")$OD.600,mu=t4bmu,alternative = 'g')
t4bp<-t.test(filter(mybreakoutX, carbon.tech=="F")$OD.600,mu=t4bmu,alternative = 'g')$p.value

#5: Is there a difference between DP5 grwth and DP24 grwth on Fructose?
var.test(x=c(filter(mybreakoutP, carbon.tech=="F")$OD.600),y=c(filter(mybreakoutX, carbon.tech=="F")$OD.600),ratoi=1,alternative = 't')
#P = .02853. Heteroscedastic.
t.test(x=c(filter(mybreakoutP, carbon.tech=="F")$OD.600),y=c(filter(mybreakoutX, carbon.tech=="F")$OD.600),mu=0,alternative='t',paired=F,var.equal = F)
t5bp<-t.test(x=c(filter(mybreakoutP, carbon.tech=="F")$OD.600),y=c(filter(mybreakoutX, carbon.tech=="F")$OD.600),mu=0,alternative='t',paired=F,var.equal = F)$p.value


#Regarding test #5: 	Finally, do a t-test: was growth higher for DP24 on F than DP5 on F? And what would that suggest??? Well, it supports the breakout idea that it's import that is key, rather than monomer metabolism. If G-F24 has higher growth on fructose, it could be because G-F24 populations are evolving to import degraded inulin dregs, while DP5 populations are evolution to import actual DP = 5 FOSs. But, it was only marginally signif at P = 0.0856

binulinpvec<-c(t3bp,t4bp,t5bp)
binulinpvec.adjust<-p.adjust(p=binulinpvec, method = "BH")#cite BH&Y 2009
binulinpvec.adjust
#FDR: 8.878737e-05 7.849274e-08 8.557475e-02


###################Plots for the breakout experiment###########################
#Ancestor plot
odbreakout <- ggplot(mybreakout, aes(x=carbon.tech, y=OD.600))
odbreakout + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,0,0,0,0,0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.6),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  labs(x="\nCarbon source",y="OD 600\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


anc.breakout.colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")

odbreakoutF <- ggplot(mybreakoutF, aes(x=carbon.tech, y=OD.600))
odbreakoutF + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  #geom_errorbar(aes(ymax=anc.breakout.vector,ymin=anc.breakout.vector),color=anc.breakout.colors,lwd = 2.2,linetype=117) + 
  #geom_errorbar(data=mybreakoutF, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=anc.breakout.colors,lwd=2.2,linetype=117) +
  geom_errorbar(data=mybreakoutF, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=anc.breakout.colors,lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0,0,0,0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.83),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  labs(x="\nCarbon source",y="OD 600\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

odbreakoutP <- ggplot(mybreakoutP, aes(x=carbon.tech, y=OD.600))
odbreakoutP + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  geom_errorbar(data=mybreakoutP, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=anc.breakout.colors,lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0,0,0,0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.83),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  labs(x="\nCarbon source",y="OD 600\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

odbreakoutX <- ggplot(mybreakoutX, aes(x=carbon.tech, y=OD.600))
odbreakoutX + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  geom_errorbar(data=mybreakoutX, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=anc.breakout.colors,lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0,0,0,0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.83),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  labs(x="\nCarbon source",y="OD 600\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


#####################################################################################################################################################################################################
########################################Breakout plots with sucrose and G:F 1:1 removed for the purpose of simplifying the figure.###################################################################
breakout_nosucrose<-filter(mybreakout,carbon.tech != "G-F" & carbon.tech != "G:F 1:1")
breakoutF_nosucrose<-filter(mybreakoutF,carbon.tech != "G-F" & carbon.tech != "G:F 1:1")
breakoutP_nosucrose<-filter(mybreakoutP,carbon.tech != "G-F" & carbon.tech != "G:F 1:1")
breakoutX_nosucrose<-filter(mybreakoutX,carbon.tech != "G-F" & carbon.tech != "G:F 1:1")


bleg<-readPNG("~/GitHub/Bifidobacterium/Figures/20211103_185311.PNG")
breakoutlegend<-rasterGrob(bleg, interpolate = TRUE)


#Ancestor plot, Fig 2.5. Needs to be converted to 2-panel figure
odbreakout_nosucrose <- ggplot(breakout_nosucrose, aes(x=carbon.tech, y=OD.600))
odbreakout_nosucrose + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 2.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,0,0,15,0,15)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.6),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G"="G","F"="F","G:F 1:5"="1:5","G-F5"="1:5","G:F 1:24"="1:24","G-F24"="1:24"))+
  labs(x="\nG:F ratio",y="Biomass (OD600)\n") +
  annotate("text",x=0.68,y=1.5,label="A",size=12,color="black")+
  annotate("text",x=1.7,y=0.55,label="D",size=12,color="black")+
  annotate("text",x=2.7,y=1.33,label="B",size=12,color="black")+
  annotate("text",x=3.7,y=1.07,label="C",size=12,color="black")+
  annotate("text",x=4.7,y=1.07,label="C",size=12,color="black")+
  annotate("text",x=5.65,y=0.4,label="E",size=12,color="black")+
  annotation_custom(grob = breakoutlegend, xmin=5.4, xmax=6.4, ymin=1.4, ymax=1.55)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))



#bp1draw <- ggdraw() +
#draw_image(bleg, x=1, y=1, scale=1)+
#  draw_plot(bp1)
#insert_element(p=breakoutlegend, left=0.5,right=1.5,top=0.9,bottom=0.5)


anc.breakout.colors_nosucrose<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")

odbreakoutF_nosucrose <- ggplot(breakoutF_nosucrose, aes(x=carbon.tech, y=OD.600))
odbreakoutF_nosucrose + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  #geom_errorbar(aes(ymax=anc.breakout.vector,ymin=anc.breakout.vector),color=anc.breakout.colors,lwd = 2.2,linetype=117) + 
  #geom_errorbar(data=breakoutF_nosucrose, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=anc.breakout.colors,lwd=2.2,linetype=117) +
  geom_errorbar(data=breakoutF_nosucrose, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=anc.breakout.colors_nosucrose,lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0,15,0,15)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.83),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G"="G","F"="F","G:F 1:5"="1:5","G-F5"="1:5","G:F 1:24"="1:24","G-F24"="1:24"))+
  labs(x="\nG:F ratio",y="Biomass (OD 600)\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

odbreakoutP_nosucrose <- ggplot(breakoutP_nosucrose, aes(x=carbon.tech, y=OD.600))
odbreakoutP_nosucrose + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  geom_errorbar(data=breakoutP_nosucrose, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=anc.breakout.colors_nosucrose,lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0,15,0,15)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.83),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G"="G","F"="F","G:F 1:5"="1:5","G-F5"="1:5","G:F 1:24"="1:24","G-F24"="1:24"))+
  labs(x="\nG:F ratio",y="Biomass (OD 600)\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

odbreakoutX_nosucrose <- ggplot(breakoutX_nosucrose, aes(x=carbon.tech, y=OD.600))
odbreakoutX_nosucrose + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  geom_errorbar(data=breakoutX_nosucrose, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=anc.breakout.colors_nosucrose,lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0,15,0,15)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.83),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G"="G","F"="F","G:F 1:5"="1:5","G-F5"="1:5","G:F 1:24"="1:24","G-F24"="1:24"))+
  labs(x="\nG:F ratio",y="Biomass (OD 600)\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


##Ancestor plot with 2 way ANOVA framework
breakout.2way.colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")

dummy.1<-mean(c(filter(mybreakout_2way, GF.ratio=="1:1", Saccharide.size=="Monosaccharide")$OD.600))
dummy.2<-mean(c(filter(mybreakout_2way, GF.ratio=="1:5", Saccharide.size=="Monosaccharide")$OD.600))
dummy.3<-mean(c(filter(mybreakout_2way, GF.ratio=="1:24", Saccharide.size=="Monosaccharide")$OD.600))
dummy.4<-mean(c(filter(mybreakout_2way, GF.ratio=="1:1", Saccharide.size=="Oligosaccharide")$OD.600))
dummy.5<-mean(c(filter(mybreakout_2way, GF.ratio=="1:5", Saccharide.size=="Oligosaccharide")$OD.600))
dummy.6<-mean(c(filter(mybreakout_2way, GF.ratio=="1:24", Saccharide.size=="Oligosaccharide")$OD.600))

anc.means.OD<-c(dummy.2,dummy.2,dummy.2,dummy.2,dummy.2,dummy.2,dummy.2,dummy.2,dummy.3,dummy.3,dummy.3,dummy.3,dummy.3,dummy.3,dummy.3,dummy.3,dummy.1,dummy.1,dummy.1,dummy.1,dummy.1,dummy.1,dummy.1,dummy.1,dummy.4,dummy.4,dummy.4,dummy.4,dummy.4,dummy.4,dummy.4,dummy.4,dummy.5,dummy.5,dummy.5,dummy.5,dummy.5,dummy.5,dummy.5,dummy.5,dummy.6,dummy.6,dummy.6,dummy.6,dummy.6,dummy.6,dummy.6,dummy.6)

###I don't understand the order in which it plots the error bars. For some reason it goes 3-->1-->2-->4-->5-->6 
###To fix the order, you have to put the values into your original spreadsheet

enaught <- ggplot(mybreakout_2way, aes(x=GF.ratio, y=OD.600))
enaught +
  #geom_errorbar(aes(ymax=anc.means.OD,ymin=anc.means.OD),color=breakout.2way.colors,lwd = 2.2,linetype=117) +
  geom_jitter(aes(shape = Saccharide.size, color = Saccharide.size),
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.35),
              size = 9.5, stroke = 1.8
  ) +
  stat_summary(
    aes(color = Saccharide.size),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 1.8, shape=95, stroke=2,
    position = position_dodge(0.35),
    show.legend=FALSE
  ) +
  scale_color_manual(values = c("black","black"), labels = c("Monomer","Oligomer")) +
  scale_shape_manual(values = c(0,15), labels = c("Monomer", "Oligomer")) +#12 is a square with a vertical cross inside it
  labs(x="\nG:F ratio",y="OD 600") +
  scale_y_continuous(limits = c(-0.04,1.64),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),legend.key.size = unit(2, "lines"),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=42),legend.title = element_text(size=34), axis.line.x.bottom = element_blank(), axis.line.y.left = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit(.3, "cm"), axis.ticks.x = element_blank())




#####Fig. 2.5
#Now make fig. 2.5, with 2 panels
fig2.5a_breakout<-filter(mybreakout,carbon.tech=="G" | carbon.tech=="F")
fig2.5a_plot <- ggplot(fig2.5a_breakout, aes(x=carbon.tech, y=OD.600))
fig2.5apanel <- fig2.5a_plot + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,0)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.7),breaks=c(.5,1,1.5), expand = c(0,0)) +
  annotate("text",x=0.5,y=1.65,label="a",size=12,color="black")+
  scale_x_discrete(labels = c("G"="Glucose","F"="Fructose"))+
  labs(x="",y="Biomass (OD600)\n") +
  annotate("text",x=0.75,y=1.5,label="A",size=12,color="black")+
  annotate("text",x=1.75,y=0.55,label="B",size=12,color="black")+
  #annotation_custom(grob = breakoutlegend, xmin=2., xmax=2.5, ymin=1.45, ymax=1.55)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))
fig2.5apanel

fig2.5bbreakout<-filter(mybreakout,carbon.tech == "G-F5" | carbon.tech == "G:F 1:5" | carbon.tech == "G-F24" | carbon.tech == "G:F 1:24")
fig2.5bcolors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")
#F plot
fig2.5bplot <- ggplot(fig2.5bbreakout, aes(x=carbon.tech, y=OD.600))
fig2.5bpanel<-fig2.5bplot + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,15,0,15)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.7),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  annotate("text",x=0.5,y=1.64,label="b",size=12,color="black")+
  scale_x_discrete(labels = c("G-F5"="1:5","G:F 1:5"="1:5","G-F24"="1:24","G:F 1:24"="1:24"))+
  labs(x="\nG:F ratio",y="Biomass (OD600)\n") +
  annotate("text",x=0.75,y=1.25,label="A",size=12,color="black")+
  annotate("text",x=1.75,y=0.95,label="B",size=12,color="black")+
  annotate("text",x=2.75,y=0.95,label="B",size=12,color="black")+
  annotate("text",x=3.75,y=0.4,label="C",size=12,color="black")+
  #geom_errorbar(data=fig2.7breakout, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=fig2.7colors,lwd=2.2,linetype=117)+
  annotation_custom(grob = breakoutlegend, xmin=3.4, xmax=4.4, ymin=1.3, ymax=1.55)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))
fig2.5bpanel

grid.arrange(fig2.5apanel,fig2.5bpanel,nrow=2)

######Fig. B.7
#This is the ancestor on sucrose vs. 1:1
B.7_breakout<-filter(mybreakout,carbon.tech == "G-F" | carbon.tech == "G:F 1:1")
#
B.7_plot <- ggplot(B.7_breakout, aes(x=carbon.tech, y=OD.600))
B.7_plot + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,15)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(.46,1.6),breaks=c(.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G-F"="Sucrose","G:F 1:1"="Monomers"))+
  labs(x="",y="Biomass (OD600)\n") +
  annotate("text",x=0.68,y=1.5,label="A",size=12,color="black")+
  annotate("text",x=1.7,y=1.35,label="B",size=12,color="black")+
  annotation_custom(grob = breakoutlegend, xmin=2., xmax=2.5, ymin=1.45, ymax=1.55)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))



#####Fig. 2.6
fig2.6breakout<-filter(mybreakoutF,carbon.tech == "G" | carbon.tech == "F")
fig2.6colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")
#F plot
fig2.6plot <- ggplot(fig2.6breakout, aes(x=carbon.tech, y=OD.600))
fig2.6plot + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,0)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.85),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G"="Glucose","F"="Fructose"))+
  labs(x="",y="Biomass (OD600)\n") +
  annotate("text",x=0.68,y=1.5,label="",size=12,color="black")+
  annotate("text",x=1.7,y=0.55,label="",size=12,color="black")+
  geom_errorbar(data=Fig2.6breakout, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=fig2.6colors,lwd=2.2,linetype=117)+
  #annotation_custom(grob = breakoutlegend, xmin=5.4, xmax=6.4, ymin=1.4, ymax=1.55)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

#####Fig. B.8
#This is F evolved bacteria on DP5, DP24, and the corresponding deconstructed monomers.
figb.8breakout<-filter(mybreakoutF,carbon.tech == "G-F5" | carbon.tech == "G:F 1:5" | carbon.tech == "G-F24" | carbon.tech == "G:F 1:24")
figb.8colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")
#F plot
figb.8plot <- ggplot(figb.8breakout, aes(x=carbon.tech, y=OD.600))
figb.8plot + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,15,0,15)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.66),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G-F5"="1:5","G:F 1:5"="1:5","G-F24"="1:24","G:F 1:24"="1:24"))+
  labs(x="\nG:F ratio",y="Biomass (OD600)\n") +
  annotate("text",x=0.68,y=1.5,label="A",size=12,color="black")+
  annotate("text",x=1.7,y=1.04,label="B",size=12,color="black")+
  annotate("text",x=2.7,y=1.5,label="A",size=12,color="black")+
  annotate("text",x=3.67,y=0.55,label="C",size=12,color="black")+
  geom_errorbar(data=figb.8breakout, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=figb.8colors,lwd=2.2,linetype=117)+
  annotation_custom(grob = breakoutlegend, xmin=3.4, xmax=4.4, ymin=1.4, ymax=1.55)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))



#####Fig. 2.7
#This is DP5 evolved popns on DP5, DP24, and corresponding monomers
fig2.7breakout<-filter(mybreakoutP,carbon.tech == "G-F5" | carbon.tech == "G:F 1:5" | carbon.tech == "G-F24" | carbon.tech == "G:F 1:24")
fig2.7colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")
#F plot
fig2.7plot <- ggplot(fig2.7breakout, aes(x=carbon.tech, y=OD.600))
fig2.7plot + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,15,0,15)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.66),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G-F5"="1:5","G:F 1:5"="1:5","G-F24"="1:24","G:F 1:24"="1:24"))+
  labs(x="\nG:F ratio",y="Biomass (OD600)\n") +
  annotate("text",x=0.68,y=1.25,label="A",size=12,color="black")+
  annotate("text",x=1.68,y=1.25,label="A",size=12,color="black")+
  annotate("text",x=2.68,y=1.25,label="A",size=12,color="black")+
  annotate("text",x=3.67,y=0.5,label="B",size=12,color="black")+
  geom_errorbar(data=fig2.7breakout, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=fig2.7colors,lwd=2.2,linetype=117)+
  annotation_custom(grob = breakoutlegend, xmin=3.4, xmax=4.4, ymin=1.35, ymax=1.6)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


######Fig. B.9
#This is the DP5 on sucrose vs. 1:1
B.9breakout<-filter(mybreakoutP,carbon.tech == "G-F" | carbon.tech == "G:F 1:1")
figb.9colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")
#
B.9_plot <- ggplot(B.9breakout, aes(x=carbon.tech, y=OD.600))
B.9_plot + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,15)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(.46,1.7),breaks=c(.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G-F"="Sucrose","G:F 1:1"="Monomers"))+
  labs(x="",y="Biomass (OD600)\n") +
  annotate("text",x=0.68,y=1.5,label="",size=12,color="black")+
  annotate("text",x=1.7,y=1.35,label="",size=12,color="black")+
  annotation_custom(grob = breakoutlegend, xmin=2.0, xmax=2.6, ymin=1.52, ymax=1.65)+
  geom_errorbar(data=B.9breakout, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=figb.9colors,lwd=2.2,linetype=117)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


######Fig. B.10
#This is the DP24-evolved on sucrose vs. 1:1
B.10breakout<-filter(mybreakoutX,carbon.tech == "G-F" | carbon.tech == "G:F 1:1")
figb.10colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")
#
B.10_plot <- ggplot(B.10breakout, aes(x=carbon.tech, y=OD.600))
B.10_plot + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,15)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(.46,1.7),breaks=c(.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G-F"="Sucrose","G:F 1:1"="Monomers"))+
  labs(x="",y="Biomass (OD600)\n") +
  annotate("text",x=0.68,y=1.5,label="",size=12,color="black")+
  annotate("text",x=1.7,y=1.35,label="",size=12,color="black")+
  annotation_custom(grob = breakoutlegend, xmin=2, xmax=2.6, ymin=1.52, ymax=1.65)+
  geom_errorbar(data=B.10breakout, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=figb.10colors,lwd=2.2,linetype=117)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))




#####Fig. 2.8
#This is DP24 evolved popns on DP5, DP24, and corresponding monomers
fig2.8breakout<-filter(mybreakoutX,carbon.tech == "G-F5" | carbon.tech == "G:F 1:5" | carbon.tech == "G-F24" | carbon.tech == "G:F 1:24")
fig2.8colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")
#F plot
fig2.8plot <- ggplot(fig2.8breakout, aes(x=carbon.tech, y=OD.600))
fig2.8plot + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,15,0,15)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.66),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G-F5"="1:5","G:F 1:5"="1:5","G-F24"="1:24","G:F 1:24"="1:24"))+
  labs(x="\nG:F ratio",y="Biomass (OD600)\n") +
  annotate("text",x=0.68,y=1.25,label="A",size=12,color="black")+
  annotate("text",x=1.68,y=1.25,label="A",size=12,color="black")+
  annotate("text",x=2.68,y=1.25,label="A",size=12,color="black")+
  annotate("text",x=3.67,y=0.5,label="B",size=12,color="black")+
  geom_errorbar(data=fig2.8breakout, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=fig2.8colors,lwd=2.2,linetype=117)+
  annotation_custom(grob = breakoutlegend, xmin=3.4, xmax=4.4, ymin=1.35, ymax=1.6)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


#####Fig. B.11
#This is DP5-evolved. Comparing growth on glu vs. fru.
figb.11breakout<-filter(mybreakoutP,carbon.tech == "G" | carbon.tech == "F")
figb.11colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")
#F plot
figb.11plot <- ggplot(figb.11breakout, aes(x=carbon.tech, y=OD.600))
figb.11plot + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,0)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.85),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G"="Glucose","F"="Fructose"))+
  labs(x="",y="Biomass (OD600)\n") +
  annotate("text",x=0.7,y=1.55,label="A",size=12,color="black")+
  annotate("text",x=1.7,y=1.15,label="B",size=12,color="black")+
  geom_errorbar(data=figb.11breakout, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=figb.11colors,lwd=2.2,linetype=117)+
  #annotation_custom(grob = breakoutlegend, xmin=5.4, xmax=6.4, ymin=1.4, ymax=1.55)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


#####Fig. B.12
#This is DP24-evolved. Comparing growth on glu vs. fru.
figb.12breakout<-filter(mybreakoutX,carbon.tech == "G" | carbon.tech == "F")
figb.12colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")
#F plot
figb.12plot <- ggplot(figb.12breakout, aes(x=carbon.tech, y=OD.600))
figb.12plot + geom_jitter(
  aes(shape = carbon.tech, color = carbon.tech), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  scale_shape_manual(values = c(0,0)) +#12 is a square with a vertical cross inside it
  #scale_color_manual(values = c("black","black","black","blue","black","blue")) +
  stat_summary(
    aes(color = carbon.tech),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,1.85),breaks=c(0,.5,1,1.5), expand = c(0,0)) +
  scale_x_discrete(labels = c("G"="Glucose","F"="Fructose"))+
  labs(x="",y="Biomass (OD600)\n") +
  annotate("text",x=0.7,y=1.55,label="A",size=12,color="black")+
  annotate("text",x=1.7,y=1.25,label="B",size=12,color="black")+
  geom_errorbar(data=figb.12breakout, mapping = aes(x=carbon.tech,ymin=anc.value,ymax=anc.value),color=figb.12colors,lwd=2.2,linetype=117)+
  #annotation_custom(grob = breakoutlegend, xmin=5.4, xmax=6.4, ymin=1.4, ymax=1.55)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


##################################Growth curves##################################################
#Do statistical tests here
grcvdataF<-filter(grcvdata, carbon.tech.assay=="F", carbon.tech.evo!="Ancestor")
grcvFaov<-aov(grcvdataF$umax ~ grcvdataF$carbon.tech.evo)
summary(grcvFaov)
TukeyHSD(grcvFaov)

grcvdataP<-filter(grcvdata, carbon.tech.assay=="G-F5", carbon.tech.evo!="Ancestor")
grcvPaov<-aov(grcvdataP$umax ~ grcvdataP$carbon.tech.evo)
summary(grcvPaov)
TukeyHSD(grcvPaov)

grcvdataX<-filter(grcvdata, carbon.tech.assay=="G-F24", carbon.tech.evo!="Ancestor")
grcvXaov<-aov(grcvdataX$umax ~ grcvdataX$carbon.tech.evo)
summary(grcvXaov)
TukeyHSD(grcvXaov)


#Now, I will do the one-sample t-tests
#Did growth rate on fructose increase for the difference populations?
#Calculate the ancestor value
ancgrrF<-mean(c(filter(grcvdata, carbon.plotting.assay=="DP0",carbon.plotting.evo=="Ancestor")$umax))

#Now do 3 one-sample t-tests
t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP0", carbon.plotting.evo=="DP0")$umax),mu=ancgrrF,alternative = 'g')
grrFp1<-t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP0", carbon.plotting.evo=="DP0")$umax),mu=ancgrrF,alternative = 'g')$p.value
#p < 0.05

t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP0", carbon.plotting.evo=="DP5")$umax),mu=ancgrrF,alternative = 'g')
grrFp2<-t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP0", carbon.plotting.evo=="DP5")$umax),mu=ancgrrF,alternative = 'g')$p.value
#p < 0.05

t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP0", carbon.plotting.evo=="DP24")$umax),mu=ancgrrF,alternative = 'g')
grrFp3<-t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP0", carbon.plotting.evo=="DP24")$umax),mu=ancgrrF,alternative = 'g')$p.value
#p < 0.05

grrFpvec<-c(grrFp1,grrFp2,grrFp3)
grrFpvec.adjust<-p.adjust(p=grrFpvec, method = "BH")#cite BH&Y 2009
grrFpvec.adjust
#FDR values: 


#Did growth rate on DP5 increase for the difference populations?
#Calculate the ancestor value
ancgrrP<-mean(c(filter(grcvdata, carbon.plotting.assay=="DP5",carbon.plotting.evo=="Ancestor")$umax))

#Now do 3 one-sample t-tests
t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP5", carbon.plotting.evo=="DP0")$umax),mu=ancgrrP,alternative = 'g')
grrPp1<-t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP5", carbon.plotting.evo=="DP0")$umax),mu=ancgrrP,alternative = 'g')$p.value
#p < 0.05

t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP5", carbon.plotting.evo=="DP5")$umax),mu=ancgrrP,alternative = 'g')
grrPp2<-t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP5", carbon.plotting.evo=="DP5")$umax),mu=ancgrrP,alternative = 'g')$p.value
#p < 0.05

t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP5", carbon.plotting.evo=="DP24")$umax),mu=ancgrrP,alternative = 'g')
grrPp3<-t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP5", carbon.plotting.evo=="DP24")$umax),mu=ancgrrP,alternative = 'g')$p.value
#p < 0.05

grrPpvec<-c(grrPp1,grrPp2,grrPp3)
grrPpvec.adjust<-p.adjust(p=grrPpvec, method = "BH")#cite BH&Y 2009
grrPpvec.adjust
#FDR values: 


grrPglcU<-read_csv("~/../OneDrive - Indiana University/Lab.Notebook/20190820_Bifidobacterium/data/anaerobic.growth.curves/umax.on.P_F-evolved_glcU.vs.lacY.csv")
#Did growth on fructose increase during 866 generations of natural selection?
var.test(x=c(filter(grrPglcU, importer=="glcU")$umax_on_P),y=c(filter(grrPglcU, importer=="lacY")$umax_on_P),ratio=1,alternative = "t")
#P = .05035. Homoscedastic.
t.test(x=c(filter(grrPglcU, importer=="glcU")$umax_on_P),y=c(filter(grrPglcU, importer=="lacY")$umax_on_P),mu=0,alternative = "t",paired=F,var.equal=T)
#p = 0.004283


#Did growth rate on DP24 increase for the difference populations?
#Calculate the ancestor value
ancgrrX<-mean(c(filter(grcvdata, carbon.plotting.assay=="DP24",carbon.plotting.evo=="Ancestor")$umax))

#Now do 3 one-sample t-tests
t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP24", carbon.plotting.evo=="DP0")$umax),mu=ancgrrX,alternative = 'g')
grrXp1<-t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP24", carbon.plotting.evo=="DP0")$umax),mu=ancgrrX,alternative = 'g')$p.value
#p < 0.05

t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP24", carbon.plotting.evo=="DP5")$umax),mu=ancgrrX,alternative = 'g')
grrXp2<-t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP24", carbon.plotting.evo=="DP5")$umax),mu=ancgrrX,alternative = 'g')$p.value
#p < 0.05

t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP24", carbon.plotting.evo=="DP24")$umax),mu=ancgrrX,alternative = 'g')
grrXp3<-t.test(x=c(filter(grcvdata, carbon.plotting.assay == "DP24", carbon.plotting.evo=="DP24")$umax),mu=ancgrrX,alternative = 'g')$p.value
#p < 0.05

grrXpvec<-c(grrXp1,grrXp2,grrXp3)
grrXpvec.adjust<-p.adjust(p=grrXpvec, method = "BH")#cite BH&Y 2009
grrXpvec.adjust
#FDR values: 0.0006389239 0.0028371206 0.0028371206

grr9onesamplespvec<-c(grrFp1,grrFp2,grrFp3,grrPp1,grrPp2,grrPp3,grrXp1,grrXp2,grrXp3)
grr9onesamplespvec.adjust<-p.adjust(p=grr9onesamplespvec, method = "BH")#cite BH&Y 2009
grr9onesamplespvec.adjust
#FDR values: 3.833544e-04 4.255681e-03 3.624473e-04 6.305623e-03 1.240404e-06 1.043232e-04 8.668466e-03 8.946585e-03 3.624473e-04
#In every case, evolved populations exhibited higher maximum growth rate than the ancestor (FDR < 0.009)

########################Growth rates plots
grcv.colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")

grcvdataF<-filter(grcvdata, carbon.tech.assay=="F", carbon.tech.evo!="Ancestor")

grcvplot_assay.F <- ggplot(grcvdataF, aes(x=carbon.plotting.evo, y=umax))
grcvplot_assay.F + geom_jitter(
  aes(shape = carbon.tech.evo, color = carbon.tech.evo),
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  #geom_errorbar(data=grcvdataF, mapping = aes(x=carbon.tech.evo,ymin=umax.anc,ymax=umax.anc),color=grcv.colors,lwd=2.2,linetype=117) +
  geom_hline(yintercept = mean(filter(grcvdata, carbon.tech.assay=="F", carbon.tech.evo=="Ancestor")$umax), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon.tech.evo),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.01,.3),breaks=c(0,.1,.2), expand = c(0,0)) +
  labs(x="\nEvolution diet",y="Max growth rate on F assay diet") +
  annotate("text",x=0.7,y=0.16,label="A",size=12,color="black")+
  annotate("text",x=1.7,y=0.1,label="B",size=12,color="black")+
  annotate("text",x=2.7,y=0.1,label="B",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

grcvdataP<-filter(grcvdata, carbon.tech.assay=="G-F5", carbon.tech.evo!="Ancestor")

grcvplot_assay.P <- ggplot(grcvdataP, aes(x=carbon.plotting.evo, y=umax))
grcvplot_assay.P + geom_jitter(
  aes(shape = carbon.tech.evo, color = carbon.tech.evo),
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  #geom_errorbar(data=grcvdataF, mapping = aes(x=carbon.tech.evo,ymin=umax.anc,ymax=umax.anc),color=grcv.colors,lwd=2.2,linetype=117) +
  geom_hline(yintercept = mean(filter(grcvdata, carbon.tech.assay=="G-F5", carbon.tech.evo=="Ancestor")$umax), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon.tech.evo),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.01,.3),breaks=c(0,.1,.2), expand = c(0,0)) +
  labs(x="\nEvolution diet",y="Max growth rate on DP5 assay diet") +
  annotate("text",x=0.7,y=0.13,label="B",size=12,color="black")+
  annotate("text",x=1.7,y=0.22,label="A",size=12,color="black")+
  annotate("text",x=2.7,y=0.13,label="B",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


grcv.colors.temp<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")

grcvdataX<-filter(grcvdata, carbon.tech.assay=="G-F24", carbon.tech.evo!="Ancestor")

grcvplot_assay.X <- ggplot(grcvdataX, aes(x=carbon.plotting.evo, y=umax))
grcvplot_assay.X + geom_jitter(
  aes(shape = carbon.tech.evo, color = carbon.tech.evo),
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  #geom_errorbar(data=grcvdataF, mapping = aes(x=carbon.tech.evo,ymin=umax.anc,ymax=umax.anc),color=grcv.colors,lwd=2.2,linetype=117) +
  geom_hline(yintercept = mean(filter(grcvdata, carbon.tech.assay=="G-F24", carbon.tech.evo=="Ancestor")$umax), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = carbon.tech.evo),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.01,.08),breaks=c(0,.05), expand = c(0,0)) +
  labs(x="\nEvolution diet",y="Max growth rate on DP24 assay diet") +
  annotate("text",x=0.7,y=0.03,label="A,B",size=12,color="black")+
  annotate("text",x=1.7,y=0.02,label="B",size=12,color="black")+
  annotate("text",x=2.7,y=0.04,label="A",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


################################################################################
#########Mutation composition analysis##########################################

################################################################################
#First, let's ask whether the number of mutns in significant genes differs among the diet treatments
mutpopdf<- read.csv("~/GitHub/Bifidobacterium/Pipeline/2_Excel_/cases_signif.muts.per.pop.csv")
mutpopaov<-aov(mutpopdf$num.muts ~ mutpopdf$Evolution_treatment)
summary(mutpopaov)
#Not signif. Same number of mutns for all the strains
#TukeyHSD(mutpopaov)


gxp.all.raw <- read.csv("~/GitHub/Bifidobacterium/Pipeline/3_Python/Outputs/Bifidobacterium_labels_gxp.csv")#<-- can do ordination analysis with ALL mutations in ALL genes
# <- read.csv("~/GitHub/Bifidobacterium/Pipeline/3_Python/Bifidobacterium_labels_gxp_fixed.only.csv")#<-- can do ordination with all FIXED mutations in ALL genes
gxp.signif.raw <- read.csv("~/GitHub/Bifidobacterium/Pipeline/3_Python/Outputs/Bifidobacterium_labels_gxp_signif.genes.only.csv")
gxp.signif.raw$sample <- as.factor(gxp.signif.raw$sample)
gxp.signif.raw$treatment <- as.factor(gxp.signif.raw$treatment)
levels(gxp.signif.raw$treatment)[levels(gxp.signif.raw$treatment)=="F"] <- "DP0"
levels(gxp.signif.raw$treatment)[levels(gxp.signif.raw$treatment)=="P"] <- "DP5"
levels(gxp.signif.raw$treatment)[levels(gxp.signif.raw$treatment)=="X"] <- "DP24"
gxp.signif <- as_tibble(gxp.signif.raw)
gxp.signif <- as.matrix(gxp.signif[,3:ncol(gxp.signif)])

gxp.all.raw$sample <- as.factor(gxp.all.raw$sample)
gxp.all.raw$treatment <- as.factor(gxp.all.raw$treatment)
levels(gxp.all.raw$treatment)[levels(gxp.all.raw$treatment)=="F"] <- "DP0"
levels(gxp.all.raw$treatment)[levels(gxp.all.raw$treatment)=="P"] <- "DP5"
levels(gxp.all.raw$treatment)[levels(gxp.all.raw$treatment)=="X"] <- "DP24"
gxp.all <- as_tibble(gxp.all.raw)
gxp.all <- as.matrix(gxp.all[,3:ncol(gxp.all)])



gxp.signif.adonis <- adonis(gxp.signif ~ gxp.signif.raw$treatment, method = "bray", permutations = 9999)
#############NEED TO UPDATE THE BELOW INFO
#                          Df  SumsOfSqs MeanSqs F.Model    R2  Pr(>F)    
#gxp.signif.raw$treatment  2    1.6067 0.80334   3.596 0.25511  1e-04 ***
#Residuals                21    4.6914 0.22340         0.74489           
#Total                    23    6.2981                 1.00000      
#

#:


# Create a distance matrix
gxp.signif.dist <- vegdist(gxp.signif, method = 'bray', upper = TRUE, diag = TRUE)
gxp.signif.dist <- as_tibble(data.matrix(gxp.signif.dist))

# Run PCoA and quantify explained variance
pcoa.eig <- cmdscale(gxp.signif.dist, eig = TRUE, k = 3)
explainvar1 <- round(pcoa.eig$eig[1] / sum(pcoa.eig$eig), 3) * 100 # 29.8 %
explainvar2 <- round(pcoa.eig$eig[2] / sum(pcoa.eig$eig), 3) * 100 # 17.1 %
explainvar3 <- round(pcoa.eig$eig[3] / sum(pcoa.eig$eig), 3) * 100 # 16.3 %
sum.eig <- sum(explainvar1, explainvar2, explainvar3) # 54 %

# Add sample and treatment IDs
gxp.signif.pcoa <- as.data.frame(pcoa.eig[1])
gxp.signif.pcoa$treatment <- gxp.signif.raw$treatment
row.names(gxp.signif.pcoa) <- gxp.signif.raw$sample
gxp.signif.pcoa$sample <- gxp.signif.raw$sample
names(gxp.signif.pcoa)[1:3] <- c('PCo1', 'PCo2', 'PCo3')



################################################################################
#################Plotting PCoA##################################################
################################################################################
png(filename="~/GitHub/Bifidobacterium/Pipeline/3_Python/pcoa.test8_signif.only.hulls.png",
    width = 1200, height = 1200, res = 96*2) 

plot.new()
par(mar = c(7, 7, 5, 7))

plot(gxp.signif.pcoa[ ,1], gxp.signif.pcoa[ ,2],
     ylim = c(-.75, .75), xlim = c(-0.75, 0.75),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     pch = 22, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1,
     axes = FALSE)

# Add axes
axis(side = 2, lwd.ticks = 2, cex.axis = 1.25, las = 1,
     labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"), at = c(-1, -0.5, 0, 0.5, 1))
#axis(side = 4, lwd.ticks = 2, cex.axis = 1.25, las = 1,
#     at=c(-1, -0.5, 0, 0.5, 1), labels = F)
axis(side = 1, lwd.ticks = 2, cex.axis = 1.25, las = 1,
     labels = c("-0.5", "0.0", "0.5"), at = c(-0.5, 0, 0.5))
#axis(side = 3, lwd.ticks = 2, cex.axis = 1.25, las = 1,
#     at = c(-0.5, 0, 0.5), labels = F)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add axis labels
mtext(expression(paste("PCo 1 (29.8 %)", sep = "")), side = 1,
      outer = TRUE, cex = 1.5, line = -3.0)
mtext(expression(paste("PCo 2 (17.1 %)", sep = "")), side = 2,
      outer = TRUE, cex = 1.5, line = -3.0)

# Plot points by strain
points(essent.F[ ,1], essent.F[ ,2], pch = 21,
       cex = 1, col = "blue", bg = "blue", lwd = 2)

points(essent.P[ ,1], essent.P[ ,2], pch = 21,
       cex = 1, col = "red", bg = "red", lwd = 2)   

points(essent.X[ ,1], essent.X[ ,2], pch = 21,
       cex = 1, col = "dark grey", bg = "dark grey", lwd = 2)  

# Add ellipses col = c("blue","red","dark grey"), 
ordihull(cbind(gxp.signif.pcoa$PCo1, gxp.signif.pcoa$PCo2), groups=gxp.signif.pcoa$treatment, lwd=2, lty=0, draw = "polygon", col = c("blue","red","dark grey"), label = TRUE)

#ordiellipse(cbind(gxp.signif.pcoa$PCo1, gxp.signif.pcoa$PCo2), groups=gxp.signif.pcoa$treatment,kind = "sd", conf = 0.95, lwd = 2, lty = 3, draw = "polygon", col = c("blue","red","dark grey"), label = TRUE)

# Subset PCoA scores by strain
essent.F <- gxp.signif.pcoa[which(gxp.signif.pcoa$treatment == "DP0"), ]
essent.P <- gxp.signif.pcoa[which(gxp.signif.pcoa$treatment == "DP5"), ]
essent.X <- gxp.signif.pcoa[which(gxp.signif.pcoa$treatment == "DP24"), ]

 

# Add P-value associated with PERMANOVA
#mtext(expression(~italic("P")~"= 0.0001"), line = -1.75, cex = 1.0, at = -0.4)

# Add treatment labels
#mtext("F", line = -19, cex = 1.2, at = -0.3)
#mtext("G-F5", line = -19, cex = 1.2, at = 0.3)
#mtext("G-F24", line = -19, cex = 1.2, at = 0.3)

# Close Plot Device
dev.off()
graphics.off()

# Show Plot
img <- readPNG("~/GitHub/Bifidobacterium/Pipeline/3_Python/pcoa.test8_signif.only.hulls.png")
grid.raster(img)



################################################################################
#########################Let's look at axis 3##################################
png(filename="~/GitHub/Bifidobacterium/Pipeline/3_Python/pcoa.test9_signif.only.PCo3.png",
    width = 1200, height = 1200, res = 96*2) 

plot.new()
par(mar = c(7, 7, 5, 7))

plot(gxp.signif.pcoa[ ,1], gxp.signif.pcoa[ ,3],
     ylim = c(-.75, .75), xlim = c(-0.75, 0.75),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     pch = 22, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1,
     axes = FALSE)

# Add axes
axis(side = 2, lwd.ticks = 2, cex.axis = 1.25, las = 1,
     labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"), at = c(-1, -0.5, 0, 0.5, 1))
#axis(side = 4, lwd.ticks = 2, cex.axis = 1.25, las = 1,
#     at=c(-1, -0.5, 0, 0.5, 1), labels = F)
axis(side = 1, lwd.ticks = 2, cex.axis = 1.25, las = 1,
     labels = c("-0.5", "0.0", "0.5"), at = c(-0.5, 0, 0.5))
#axis(side = 3, lwd.ticks = 2, cex.axis = 1.25, las = 1,
#     at = c(-0.5, 0, 0.5), labels = F)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add axis labels
mtext(expression(paste("PCo 1 (29.8 %)", sep = "")), side = 1,
      outer = TRUE, cex = 1.5, line = -3.0)
mtext(expression(paste("PCo 3 (16.3 %)", sep = "")), side = 2,
      outer = TRUE, cex = 1.5, line = -3.0)

# Plot points by strain
points(essent.F[ ,1], essent.F[ ,3], pch = 21,
       cex = 1, col = "blue", bg = "blue", lwd = 2)

points(essent.P[ ,1], essent.P[ ,3], pch = 21,
       cex = 1, col = "red", bg = "red", lwd = 2)   

points(essent.X[ ,1], essent.X[ ,3], pch = 21,
       cex = 1, col = "dark grey", bg = "dark grey", lwd = 2)  

# Add ellipses col = c("blue","red","dark grey"), 
ordihull(cbind(gxp.signif.pcoa$PCo1, gxp.signif.pcoa$PCo3), groups=gxp.signif.pcoa$treatment,conf = 0.95, lwd = 2, lty = 0, draw = "polygon", col = c("blue","red","dark grey"), label = TRUE)

#ordiellipse(cbind(gxp.signif.pcoa$PCo1, gxp.signif.pcoa$PCo2), groups=gxp.signif.pcoa$treatment,kind = "sd", conf = 0.95, lwd = 2, lty = 3, draw = "polygon", col = c("blue","red","dark grey"), label = TRUE)

# Subset PCoA scores by strain
essent.F <- gxp.signif.pcoa[which(gxp.signif.pcoa$treatment == "DP0"), ]
essent.P <- gxp.signif.pcoa[which(gxp.signif.pcoa$treatment == "DP5"), ]
essent.X <- gxp.signif.pcoa[which(gxp.signif.pcoa$treatment == "DP24"), ]



# Add P-value associated with PERMANOVA
#mtext(expression(~italic("P")~"= 0.0001"), line = -1.75, cex = 1.0, at = -0.4)

# Add treatment labels
#mtext("F", line = -19, cex = 1.2, at = -0.3)
#mtext("G-F5", line = -19, cex = 1.2, at = 0.3)
#mtext("G-F24", line = -19, cex = 1.2, at = 0.3)

# Close Plot Device
dev.off()
graphics.off()

# Show Plot
img <- readPNG("~/GitHub/Bifidobacterium/Pipeline/3_Python/pcoa.test9_signif.only.PCo3.png")
grid.raster(img)


################################################################################
################################################################################
# Analysis of variation among populations (beta diversity) within each diet treatment.
# Are there differences in composition variance between the treatments?
################################################################################
gd <- vegdist(gxp.signif, method = 'bray', upper = TRUE, diag = TRUE)
compare.variance <- betadisper(d=gd, group = c("F","F","F","F","F","F","F","F","G-F5","G-F5","G-F5","G-F5","G-F5","G-F5","G-F5","G-F5","G-F24","G-F24","G-F24","G-F24","G-F24","G-F24","G-F24","G-F24"))
c.v.result <- anova(compare.variance)
c.v.result
TukeyHSD(compare.variance)
###Significant difference in the amount of variance, i.e. beta diversity, among the treatments. F has more beta diversity than G-F5 or G-F24 do


################################################################################
#################Indicator gene analysis########################################
trts.i<- c("F","F","F","F","F","F","F","F","G-F5","G-F5","G-F5","G-F5","G-F5","G-F5","G-F5","G-F5","G-F24","G-F24","G-F24","G-F24","G-F24","G-F24","G-F24","G-F24")
gxp.i<-gxp.signif
indvals = multipatt(x=gxp.i, cluster=trts.i,func="IndVal.g", control = how(nperm=999))
summary(indvals, indvalcomp = TRUE)
#Component A: "specificity" or "positive predictive value" of the species as indicator of the site group. It is the sample estimate of the probability that the surveyed site belongs to the target site group given the fact that the species has been found.
#Component B: "Fidelity" or "sensitivity" of the species as indicator of the target site group. It is is sample estimate of the probability of finding the species in sites belonging to the site group.

gxpcomb<- combinespecies(gxp.i, max.order = 2)$XC
dim(gxpcomb)
indvalcomb = multipatt(x=gxpcomb, cluster=trts.i, duleg = TRUE,control = how(nperm=999))
summary(indvalcomb, indvalcomp = TRUE)

####This is using the whole matrix instead of specifically the significant genes:
gxpa.i<-gxp.all
indval=multipatt(x=gxpa.i,cluster=trts.i,func="IndVal.g",control=how(nperm=999))
summary(indval,indvalcomp = TRUE)

gxpacomb<- combinespecies(gxpa.i, max.order = 2)$XC
dim(gxpacomb)
#indvalacomb = multipatt(x=gxpacomb, cluster=trts.i, duleg = TRUE,control = how(nperm=999))
#summary(indvalacomb, indvalcomp = TRUE)
#I just commented out the 2 abopve lines so that the whole script runs faster when I run it to declare my variables
################################################################################
#################HPLC experiments###############################################

################1###############################################################
#First, we'll show whether or not the ancestor BB-12 preferentially imports glucose over fructose when both G and F are present at 38.86 mM
var.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="Anc")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="Anc")$flux.fru.pos),ratio=1,alternative = "t")
#P = .2942. Homoscedastic.
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="Anc")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="Anc")$flux.fru.pos),mu=0,alternative = "t",paired=T,var.equal=T)
#P = 0.007261. There is a significant difference between the extracellular glucose flux and the extracellular fructose flux in the ancestor bacterium.
#wilcox.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="Anc")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="Anc")$flux.fru.pos),alternative = "t", mu = 0, paired = T, exact = NULL, correct = T, conf.int = F, conf.level = 0.95)

#Also, I recognized that there is no such argument as 'var.equal' when using a PAIRED t-test, but I'm leaving it in for explicative value (for "expliciticity"???? xD). It is also good to do the variance test itself as a sanity check.

###############2###############################################################
#The ancestor preferentially imports glucose. Do the evolved bacteria (evolved on F, DP5, and DP24) still prefer glucose after evolution? One might expect that, for example, bacteria evolved on fructose no longer prefer glucose.
#We'll test the fructose-evolved populations first.
var.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="F")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="F")$flux.fru.pos),ratio=1,alternative = "t")
#P = .4000. Homoscedastic.
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="F")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="F")$flux.fru.pos),mu=0,alternative = "t",paired=T,var.equal=T)
#P = 0.0002978. There is still a significant difference between the extracellular glucose flux and the extracellular fructose flux in the F-evolved bacteria.
#FDR #0.0002978366

#We'll test the DP5-evolved populations next.
var.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="P")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="P")$flux.fru.pos),ratio=1,alternative = "t")
#P = .004368. Heteroscedastic.
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="P")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="P")$flux.fru.pos),mu=0,alternative = "t",paired=T,var.equal=F)
#P = 3.665EE-05. There is still a significant difference between the extracellular glucose flux and the extracellular fructose flux in the DP5-evolved bacteria.
#FDR 0.0001099467

#We'll test the DP24-evolved populations last.
var.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="X")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="X")$flux.fru.pos),ratio=1,alternative = "t")
#P = .5905. Homoscedastic.
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="X")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="X")$flux.fru.pos),mu=0,alternative = "t",paired=T,var.equal=T)
#P = 8.857EE-05. There is still a significant difference between the extracellular glucose flux and the extracellular fructose flux in the DP24-evolved bacteria.
#FDR 0.0001328539

#Adjust P-values for multiple tests (FDR method)
hplcpvec2<-c(t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="F")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="F")$flux.fru.pos),mu=0,alternative = "t",paired=T,var.equal=T)$p.value,t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="P")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="P")$flux.fru.pos),mu=0,alternative = "t",paired=T,var.equal=F)$p.value,t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="X")$flux.glu.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="X")$flux.fru.pos),mu=0,alternative = "t",paired=T,var.equal=T)$p.value)
hplcpvec2.adjust<-p.adjust(p=hplcpvec2, method = "BH")#cite BH&Y 2009
hplcpvec2.adjust
#0.0002978366 0.0001099467 0.0001328539


############3###################################################################
#What about in the presence of fructose alone?
#In that case, do the evolved bacteria have increased fructose flux relative to the ancestor?
#Also, do fructose flux rates increase when G is removed, suggestive of an effect of metabolite repression by glucose?

#First, let's test whether the ancestor has a higher F flux on F alone (fructose at 39 mM) compared to GF (glucose 39 mM + fructose 39mM)
var.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="Anc")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="Anc")$flux.fru.pos),ratio=1,alternative = "t")
#P = .05057. Homoscedastic.
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="Anc")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="Anc")$flux.fru.pos),mu=0,alternative = "g",paired=F,var.equal=T)
hplcp3.1<-t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="Anc")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="Anc")$flux.fru.pos),mu=0,alternative = "g",paired=F,var.equal=T)$p.value
#Not paired because these are biological replicates (tube level), not evolution replicates (genotype level).
#P = 0.04039.
#FDR: 0.040389743 

#What about the evolved bacteria? They have evolved to perhaps have measurable levels of fructose import---can we detect metabolite repression for them, i.e. can we see that fructose flux is higher in the presence of fructose 39 mM alone (F) than in the presence of glucose 39 mM + fructose 39 mM (GF)?

#F-evolved populations:
var.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="F")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="F")$flux.fru.pos),ratio=1,alternative = "t")
#P = .4977. Homoscedastic.
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="F")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="F")$flux.fru.pos),mu=0,alternative = "g",paired=T,var.equal=T)
hplcp3.2<-t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="F")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="F")$flux.fru.pos),mu=0,alternative = "g",paired=T,var.equal=T)$p.value
#P < .05. F import is higher in the absence of G (10.97 vs. 4.28). Suggests there is meaningful metabolite repression.
#FDR: 0.001102758 

#DP5-evolved populations:
var.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="P")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="P")$flux.fru.pos),ratio=1,alternative = "t")
#P = .0009092. Heteroscedastic.
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="P")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="P")$flux.fru.pos),mu=0,alternative = "g",paired=T,var.equal=F)
hplcp3.3<-t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="P")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="P")$flux.fru.pos),mu=0,alternative = "g",paired=T,var.equal=F)$p.value
#P = < .05. F import is higher in the absence of G (5.83 vs. 0.41). Suggests there is meaningful metabolite repression.
#FDR:  0.002433313 

#DP24-evolved populations:
var.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="X")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="X")$flux.fru.pos),ratio=1,alternative = "t")
#P = .8434. Homoscedastic.
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="X")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="X")$flux.fru.pos),mu=0,alternative = "g",paired=T,var.equal=T)
hplcp3.4<-t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="X")$flux.fru.pos),y=c(filter(myhplc, condition_short=="GF", evo_trt=="X")$flux.fru.pos),mu=0,alternative = "g",paired=T,var.equal=T)$p.value
#P < .05. F import is higher in the absence of G (9.20 vs. 3.42). Suggests there is meaningful metabolite repression.
#FDR: 0.001102758

hplcpvec3a<-c(hplcp3.1,hplcp3.2,hplcp3.3,hplcp3.4)
hplcpvec3a.adjust<-p.adjust(p=hplcpvec3a, method = "BH")#cite BH&Y 2009
hplcpvec3a
hplcpvec3a.adjust

#####3b#########################################################################
#Importantly, can we detect an increase in fructose import rate of the evolved bacteria, relative to the ancestor?

#First, I need to calculate the average fructose import rate by the ancestor and use it as my mu value for 1-sample t-tests.
anc.flux.fru_F<-mean(c(filter(myhplc, condition_short=="F",evo_trt=="Anc")$flux.fru.pos))
anc.flux.fru_GF<-mean(c(filter(myhplc, condition_short=="GF",evo_trt=="Anc")$flux.fru.pos))

#Considering only GF condition.
#F evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="F")$flux.fru.pos),mu=anc.flux.fru_GF,alternative = "g")
hplcp3.5<-t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="F")$flux.fru.pos),mu=anc.flux.fru_GF,alternative = "g")$p.value
#P = 0.02826. F import after evolution is significantly higher in the presence of glucose BUT NOT AFTER fdr CORRECTION:
#FDR: 0.056513990 

#DP5 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="P")$flux.fru.pos),mu=anc.flux.fru_GF,alternative = "g")
hplcp3.6<-t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="P")$flux.fru.pos),mu=anc.flux.fru_GF,alternative = "g")$p.value
#P = 0.9979. F import after evolution is not significantly higher in the presence of glucose.
#FDR:  0.997888364 

#DP24 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="X")$flux.fru.pos),mu=anc.flux.fru_GF,alternative = "g")
hplcp3.7<-t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="X")$flux.fru.pos),mu=anc.flux.fru_GF,alternative = "g")$p.value
#P = 0.04082. F import after evolution is significantly higher in the presence of glucose. (This will prolly change after correcting for multiply testing, though.)
#FDR: 0.061233424 


#Now we'll consider the F condition.
#F evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="F")$flux.fru.pos),mu=anc.flux.fru_F,alternative = "g")
hplcp3.8<-t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="F")$flux.fru.pos),mu=anc.flux.fru_F,alternative = "g")$p.value
#P < .05. F import after evolution is significantly higher in the absence of glucose.
#FDR: 0.005412862 

#DP5 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="P")$flux.fru.pos),mu=anc.flux.fru_F,alternative = "g")
hplcp3.9<-t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="P")$flux.fru.pos),mu=anc.flux.fru_F,alternative = "g")$p.value
#P = 0.1165 F import after evolution is NOT significantly higher in the absence of glucose than that of the ancestor in the absence of glucose.
#THIS MAKES SENSE. DP5 can be imported directly. Meanwhile, DP24 pops are surviving on the dregs of degraded molecules.
#FDR: 0.139784891 

#DP24 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="X")$flux.fru.pos),mu=anc.flux.fru_F,alternative = "g")
hplcp3.10<-t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="X")$flux.fru.pos),mu=anc.flux.fru_F,alternative = "g")$p.value
#P < .05. F import after evolution is significantly higher in the absence of glucose.
#FDR: 0.004927588



hplcpvec3b<-c(hplcp3.5,hplcp3.6,hplcp3.7,hplcp3.8,hplcp3.9,hplcp3.10)
hplcpvec3b.adjust<-p.adjust(p=hplcpvec3b, method = "BH")#cite BH&Y 2009
hplcpvec3b
hplcpvec3b.adjust
#0.056513990 0.997888364 0.061233424 0.005412862 0.139784891 0.004927588
######4#########################################################################
#Do the evolved populations differ in their fructose import rates (in the presence of F alone)?
hplc4.1df<-filter(myhplc,condition_short=="F",evo_trt!="Anc")
hplc4.1aov<-aov(hplc4.1df$flux.fru.pos ~ hplc4.1df$evo_trt)
summary(hplc4.1aov)
TukeyHSD(hplc4.1aov)
#Tukey letters: F: A  DP5: B  DP24: A,B
#Interpretation: F pops were selected to import fructose. DP5 pops were selected to import short oligomers. This likely improved the function of lacY for importing oligomers, but the kinds of mutations that fixed were significantly less beneficial for importing fructose than those that fixed in the fructose-evolved populations. DP24 populations are somewhere intermediate. They're likely surviving on scraps of degraded big oligomers, which prolly includes both small oligomers and monomers.
#Padj values from Tukey HSD:
#P-F 0.0247383
#X-F 0.5942214
#X-P 0.1725517

#Why not also investigate the ANOVA outcome when we look at the GF cdtn
hplc4.2df<-filter(myhplc,condition_short=="GF",evo_trt!="Anc")
hplc4.2aov<-aov(hplc4.2df$flux.fru.pos ~ hplc4.2df$evo_trt)
summary(hplc4.2aov)
TukeyHSD(hplc4.2aov)
#Same story.
#Padj:
#.01477
#.77435
#.06358


######5#########################################################################
#Now I want to compare the fructose import rate between fructose-evolved populations that have LacY mutations to fructose-evolved pops that have GlcU mutations
var.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="F",Importer=="lacY")$flux.fru.pos),y=c(filter(myhplc, condition_short=="F", evo_trt=="F",Importer=="glcU")$flux.fru.pos),ratio=1,alternative = "t")
#P = .1623. Homoscedastic.
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="F",Importer=="lacY")$flux.fru.pos),y=c(filter(myhplc, condition_short=="F", evo_trt=="F",Importer=="glcU")$flux.fru.pos),mu=0,alternative = "t",paired=F,var.equal=T)
mean(c(filter(myhplc, condition_short=="F",evo_trt=="F",Importer=="lacY")$flux.fru.pos))
mean(c(filter(myhplc, condition_short=="F",evo_trt=="F",Importer=="glcU")$flux.fru.pos))
#P = 0.006563. F import is higher for the populations with lacY mutations than for those with glcU mutations (14.51 fmol/cell/hr vs 7.43 fmol/cell/hr)

##########6#####################################################################
#Next, I'll ask whether evolved bacteria have a higher flux of lact and ace production than the ancestor.

#Let's consider the general GF diet first.
#Lactate?
anc.flux.lact_GF<-mean(c(filter(myhplc, condition_short=="GF",evo_trt=="Anc")$flux.lact.0))


#F evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="F")$flux.lact.0),mu=anc.flux.lact_GF,alternative = "g")
#P = 0.00002354. Lact export after evolution is significantly higher in GF cdtn.
#FDR: 7.060763e-05

#DP5 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="P")$flux.lact.0),mu=anc.flux.lact_GF,alternative = "g")
#P = 0.0009305. Lact export after evolution is significantly higher in GF cdtn.
#FDR: 1.861000e-03

#DP24 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="X")$flux.lact.0),mu=anc.flux.lact_GF,alternative = "g")
#P = 7.845EE-06. Lact export after evolution is significantly higher in GF cdtn.
#FDR: 4.707231e-05

####
#Let's now consider acetate
#How about in general, i.e. with glucose around (GF condition)?
anc.flux.ace_GF<-mean(c(filter(myhplc, condition_short=="GF",evo_trt=="Anc")$flux.ace.0))

#F evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="F")$flux.ace.0),mu=anc.flux.ace_GF,alternative = "g")
#P = 4.615EE-05 Ace export after evolution is significantly higher in GF cdtn.
#FDR: 2.367139e-03 

#DP5 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="P")$flux.ace.0),mu=anc.flux.ace_GF,alternative = "g")
#P = 6.815ee-06. Ace export after evolution is significantly higher in GF cdtn.
#FDR: 5.489388e-02 

#DP24 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="X")$flux.ace.0),mu=anc.flux.ace_GF,alternative = "g")
#P = 4.027EE-06. Ace export after evolution is significantly higher in GF cdtn.
#FDR:2.367139e-03

hplcpvec6.1<-c(t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="F")$flux.lact.0),mu=anc.flux.lact_GF,alternative = "g")$p.value,t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="P")$flux.lact.0),mu=anc.flux.lact_GF,alternative = "g")$p.value,t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="X")$flux.lact.0),mu=anc.flux.lact_GF,alternative = "g")$p.value,t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="F")$flux.ace.0),mu=anc.flux.ace_F,alternative = "g")$p.value,t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="P")$flux.ace.0),mu=anc.flux.ace_F,alternative = "g")$p.value,t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="X")$flux.ace.0),mu=anc.flux.ace_F,alternative = "g")$p.value)
hplcpvec6.1.adjust<-p.adjust(p=hplcpvec6.1, method = "BH")#cite BH&Y 2009
hplcpvec6.1
hplcpvec6.1.adjust
#7.060763e-05 1.861000e-03 4.707231e-05 2.367139e-03 5.489388e-02 2.367139e-03


#####Next, we need to consider the experiment with 38.9 F as the sole carbon source.

#Does lactate production increase after evolution when we consider the fructose condition?
anc.flux.lact_F<-mean(c(filter(myhplc, condition_short=="F",evo_trt=="Anc")$flux.lact.0))

#F evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="F")$flux.lact.0),mu=anc.flux.lact_F,alternative = "g")
#P = 0.0132. Lact export after evolution is significantly higher in F cdtn.
#FDR: 1.584397e-02

#DP5 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="P")$flux.lact.0),mu=anc.flux.lact_F,alternative = "g")
#P = 0.003382. Lact export after evolution is significantly higher in F cdtn.
#FDR: 3.382046e-02

#DP24 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="X")$flux.lact.0),mu=anc.flux.lact_F,alternative = "g")
#P = 0.001789. Lact export after evolution is significantly higher in F cdtn.
#FDR: 2.684201e-03

####
#Now, does acetate production increase after evolution when we consider the fructose condition?
anc.flux.ace_F<-mean(c(filter(myhplc, condition_short=="F",evo_trt=="Anc")$flux.ace.0))

#F evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="F")$flux.ace.0),mu=anc.flux.ace_F,alternative = "g")
#P = 0.0009068. Ace export after evolution is significantly higher in F cdtn.
#FDR: 1.639248e-02 

#DP5 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="P")$flux.ace.0),mu=anc.flux.ace_F,alternative = "g")
#P = 0.01672. Ace export after evolution is significantly higher in F cdtn.
#FDR: 3.579952e-02 

#DP24 evolved vs. ancestor
t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="X")$flux.ace.0),mu=anc.flux.ace_F,alternative = "g")
#P = 0.0007976. Ace export after evolution is significantly higher in F cdtn.
#FDR: 2.793390e-03

hplcpvec6.2<-c(t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="F")$flux.ace.0),mu=anc.flux.ace_GF,alternative = "g")$p.value,t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="P")$flux.ace.0),mu=anc.flux.ace_GF,alternative = "g")$p.value,t.test(x=c(filter(myhplc, condition_short=="GF",evo_trt=="X")$flux.ace.0),mu=anc.flux.ace_GF,alternative = "g")$p.value,t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="F")$flux.lact.0),mu=anc.flux.lact_F,alternative = "g")$p.value,t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="P")$flux.lact.0),mu=anc.flux.lact_F,alternative = "g")$p.value,t.test(x=c(filter(myhplc, condition_short=="F",evo_trt=="X")$flux.lact.0),mu=anc.flux.lact_F,alternative = "g")$p.value)
hplcpvec6.2.adjust<-p.adjust(p=hplcpvec6.2, method = "BH")#cite BH&Y 2009
hplcpvec6.2
hplcpvec6.2.adjust
#9.229471e-05 2.044568e-05 2.044568e-05 



#####7##########################################################################
#Last, we'll ask whether lact and ace fluxes are differ depending on the evolution treatment. So, 4 ANOVAs (2 cdtns, 2 metabolites).

#Do the evolved populations differ in their lactate export rates in the general GF cdtn?
hplc7.1df<-filter(myhplc,condition_short=="GF",evo_trt!="Anc")
hplc7.1aov<-aov(hplc7.1df$flux.lact.0 ~ hplc7.1df$evo_trt)
summary(hplc7.1aov)
TukeyHSD(hplc7.1aov)
#Tukey letters: F: A  DP5: B  DP24: A
#Interpretation: F and DP24 pops are better at making lactate from glucose. I think it prolly means that F and DP24 are better at importing and metabolizing monomers, while DP5 have evolved to specialize on short oligomers. However, we would need to do another entire experiment to test whether lact and ace production are higher for DP5 pops than other pops when DP5 is the C src.
#Padj values from Tukey HSD:
#P-F 0.0004375
#X-F 0.3846448
#X-P 0.0102524

#Now we'll compare the evolved populations' abilities to produce lactate when the carbon source is fructose alone. Expectation would be that DP5 populations would be at even more of disadvantage now.
hplc7.2df<-filter(myhplc,condition_short=="F",evo_trt!="Anc")
hplc7.2aov<-aov(hplc7.2df$flux.lact.0 ~ hplc7.2df$evo_trt)
summary(hplc7.2aov)
TukeyHSD(hplc7.2aov)
#However, there are actually now no signif differences. ANOVA P = 0.255


#Do the evolved populations differ in their acetate export rates in the general GF cdtn?
hplc7.3df<-filter(myhplc,condition_short=="GF",evo_trt!="Anc")
hplc7.3aov<-aov(hplc7.3df$flux.ace.0 ~ hplc7.3df$evo_trt)
summary(hplc7.3aov)
TukeyHSD(hplc7.3aov)
#Tukey letters: F: A  DP5: B  DP24: A,B
#Interpretation: F pops are better at producing acetate from GF than are the other evolved populations. Could make sense, as F populations have been selected to import and use solely monomers.
#Padj values from Tukey HSD:
#P-F 0.0057101
#X-F 0.2898775
#X-P 0.1474614

#Now we'll compare the evolved populations' abilities to produce lactate when the carbon source is fructose alone. Expectation would be that non-F populations would be at even more of disadvantage now.
hplc7.4df<-filter(myhplc,condition_short=="F",evo_trt!="Anc")
hplc7.4aov<-aov(hplc7.4df$flux.ace.0 ~ hplc7.4df$evo_trt)
summary(hplc7.4aov)
TukeyHSD(hplc7.4aov)
#However, there are actually now no signif differences. ANOVA P = 0.078, and all Tukey HSD Padj values > 0.1


#####HPLC data: Making figures##################################################
hlpc.colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")

hplcdataFall<-filter(myhplc, condition_short=="F")
hplcdataGFall<-filter(myhplc, condition_short=="GF")
hplcdataF<-filter(myhplc, condition_short=="F", evo_trt!="Anc")
hplcdataF$evo_trt_plotting<-droplevels(hplcdataF$evo_trt_plotting)
hplcdataGF<-filter(myhplc, condition_short=="GF", evo_trt!="Anc")
hplcdataGF$evo_trt_plotting<-droplevels(hplcdataGF$evo_trt_plotting)

#Make Fig. 2.3
hplclactplot_F <- ggplot(hplcdataF, aes(x=evo_trt_plotting, y=flux.lact.0))
###
hplclactpanel_F <- hplclactplot_F + geom_jitter(aes(shape = evo_trt_plotting, color = evo_trt_plotting), position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.7), size = 5.5, stroke = 2) +
  #geom_errorbar(data=grcvdataF, mapping = aes(x=carbon.tech.evo,ymin=umax.anc,ymax=umax.anc),color=grcv.colors,lwd=2.2,linetype=117) +
  geom_hline(yintercept = mean(filter(myhplc, condition_short=="F", evo_trt_plotting=="Ancestor")$flux.lact.0), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = evo_trt_plotting),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.54,16),breaks=c(0,5,10,15), expand = c(0,0)) +
  #labs(x="\nEvolution diet",y="Lactate flux (fmol/cell/hr)") +
  xlab("\nEvolution diet")+
  ylab( expression(atop(paste("Lactate flux (fmol/cell * ",hr^{-1},")"))))+
  annotate("text",x=0.5,y=15.5,label="a",size=12,color="black")+
  #annotate("text",x=1.7,y=0.1,label="B",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))
hplclactpanel_F

hplcaceplot_F <- ggplot(hplcdataF, aes(x=evo_trt_plotting, y=flux.ace.0))
###
hplcacepanel_F <- hplcaceplot_F + geom_jitter(aes(shape = evo_trt_plotting, color = evo_trt_plotting), position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.7), size = 5.5, stroke = 2) +
  #geom_errorbar(data=grcvdataF, mapping = aes(x=carbon.tech.evo,ymin=umax.anc,ymax=umax.anc),color=grcv.colors,lwd=2.2,linetype=117) +
  geom_hline(yintercept = mean(filter(myhplc, condition_short=="F", evo_trt_plotting=="Ancestor")$flux.ace.0), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = evo_trt_plotting),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,28.4),breaks=c(0,5,10,15,20,25), expand = c(0,0)) +
  labs(x="\nEvolution diet",y="Acetate flux (fmol/cell/hr)") +
  xlab("\nEvolution diet")+
  ylab( expression(atop(paste("Acetate flux (fmol/cell * ",hr^{-1},")"))))+
  annotate("text",x=0.5,y=27.3,label="b",size=12,color="black")+
  #annotate("text",x=1.7,y=0.1,label="B",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))
hplcacepanel_F

grid.arrange(hplclactpanel_F,hplcacepanel_F,nrow=1)




#Now make the supp figure: metabolite production on the GF HPLC diet
#Fig. B.5
hplclactplot_GF <- ggplot(hplcdataGF, aes(x=evo_trt_plotting, y=flux.lact.0))
###
hplclactpanel_GF <- hplclactplot_GF + geom_jitter(aes(shape = evo_trt_plotting, color = evo_trt_plotting), position = position_jitterdodge(jitter.width = 0.85, dodge.width = 0.7), size = 5.5, stroke = 2) +
  #geom_errorbar(data=grcvdataF, mapping = aes(x=carbon.tech.evo,ymin=umax.anc,ymax=umax.anc),color=grcv.colors,lwd=2.2,linetype=117) +
  geom_hline(yintercept = mean(filter(myhplc, condition_short=="GF", evo_trt_plotting=="Ancestor")$flux.lact.0), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = evo_trt_plotting),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,26.5),breaks=c(0,5,10,15,20,25), expand = c(0,0)) +
  #labs(x="\nEvolution diet",y="Lactate flux (fmol/cell/hr)") +
  xlab("\nEvolution diet")+
  ylab( expression(atop(paste("Lactate flux (fmol/cell * ",hr^{-1},")"))))+
  annotate("text",x=0.5,y=25,label="a",size=12,color="black")+
  annotate("text",x=0.7,y=16.2,label="A",size=12,color="black")+
  annotate("text",x=1.7,y=8.1,label="B",size=12,color="black")+
  annotate("text",x=2.7,y=16.2,label="A",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))
hplclactpanel_GF

hplcaceplot_GF <- ggplot(hplcdataGF, aes(x=evo_trt_plotting, y=flux.ace.0))
###
hplcacepanel_GF <- hplcaceplot_GF + geom_jitter(aes(shape = evo_trt_plotting, color = evo_trt_plotting), position = position_jitterdodge(jitter.width = 0.9, dodge.width = 0.7), size = 5.5, stroke = 2) +
  #geom_errorbar(data=grcvdataF, mapping = aes(x=carbon.tech.evo,ymin=umax.anc,ymax=umax.anc),color=grcv.colors,lwd=2.2,linetype=117) +
  geom_hline(yintercept = mean(filter(myhplc, condition_short=="GF", evo_trt_plotting=="Ancestor")$flux.ace.0), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = evo_trt_plotting),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,41.5),breaks=c(0,10,20,30,40), expand = c(0,0)) +
  #labs(x="\nEvolution diet",y="Acetate flux (fmol/cell/hr)") +
  xlab("\nEvolution diet")+
  ylab( expression(atop(paste("Acetate flux (fmol/cell * ",hr^{-1},")"))))+
  annotate("text",x=0.5,y=40,label="b",size=12,color="black")+
  annotate("text",x=0.7,y=27.6,label="A",size=12,color="black")+
  annotate("text",x=1.7,y=20.6,label="B",size=12,color="black")+
  annotate("text",x=2.65,y=24.2,label="A,B",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))
hplcacepanel_GF

grid.arrange(hplclactpanel_GF,hplcacepanel_GF,nrow=1)



####################################################
#Make fig. 2.9
#Compare anc import of glu to import of fru on the GF hplc diet
#hplcdataGFanc<-filter(hplcdataGFall, evo_trt=="Anc")

mycols2.9<-c('black','black','black','black','black','black','black','black')

fig2.9df<-read.csv("~/GitHub/Bifidobacterium/datafiles/fig.2.9.csv")
fig2.9df$metabolite<-recode(fig2.9df$metabolite, "F" = "Fructose", "G" = "Glucose")
fig2.9df$metabolite<-factor(fig2.9df$metabolite, levels=c("Glucose","Fructose"))

fig2.9 <- ggplot(fig2.9df, aes(x=metabolite, y=flux.data.col))
###
fig2.9 + geom_jitter(aes(shape = metabolite, color = metabolite), position = position_jitterdodge(jitter.width = 0.9, dodge.width = 0.7), size = 5.5, stroke = 2) +
  #geom_errorbar(data=grcvdataF, mapping = aes(x=carbon.tech.evo,ymin=umax.anc,ymax=umax.anc),color=grcv.colors,lwd=2.2,linetype=117) +
  #geom_hline(yintercept = mean(filter(myhplc, condition_short=="GF", evo_trt_plotting=="Ancestor")$flux.ace.0), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = metabolite),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols2.9) +
  scale_y_continuous(limits = c(-0.04,9.5),breaks=c(0,2.5,5,7.5), expand = c(0,0)) +
  #labs(x="",y="Flux (fmol/cell/hr)") +
  xlab(" ")+
  ylab( expression(atop(paste("Flux (fmol/cell * ",hr^{-1},")"))))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

#Make fig. B.13
#It's 3 panels stacked vertically.
#Each panel is gonna be like 2.9 is.
#But then add the panel code that you used for fig. b.5

mycolsb.13<-c('black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black')
anccolsb.13<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")

figb.13adf<-read.csv("~/GitHub/Bifidobacterium/datafiles/fig.b.13a.csv")
figb.13adf$metabolite<-recode(figb.13adf$metabolite, "F" = "Fructose", "G" = "Glucose")
figb.13adf$metabolite<-factor(figb.13adf$metabolite, levels=c("Glucose","Fructose"))

figb.13a <- ggplot(figb.13adf, aes(x=metabolite, y=flux.data.col))
###
figb.13apanel <- figb.13a + geom_jitter(aes(shape = metabolite, color = metabolite), position = position_jitterdodge(jitter.width = 0.9, dodge.width = 0.7), size = 5.5, stroke = 2) +
  geom_errorbar(data=figb.13adf, mapping = aes(x=metabolite,ymin=flux.anc.data.col,ymax=flux.anc.data.col),color=anccolsb.13,lwd=2.2,linetype=117) +
  #geom_hline(yintercept = mean(filter(myhplc, condition_short=="GF", evo_trt_plotting=="Ancestor")$flux.ace.0), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = metabolite),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycolsb.13) +
  scale_y_continuous(limits = c(-0.5,20.35),breaks=c(0,5,10,15,20), expand = c(0,0)) +
  #labs(x="",y="Flux (fmol/cell/hr)") +
  xlab(" ")+
  ylab( expression(atop(paste("Flux (fmol/cell * ",hr^{-1},")"))))+
  annotate("text",x=0.5,y=19.3,label="a",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))
figb.13apanel


figb.13bdf<-read.csv("~/GitHub/Bifidobacterium/datafiles/fig.b.13b.csv")
figb.13bdf$metabolite<-recode(figb.13bdf$metabolite, "F" = "Fructose", "G" = "Glucose")
figb.13bdf$metabolite<-factor(figb.13bdf$metabolite, levels=c("Glucose","Fructose"))

figb.13b <- ggplot(figb.13bdf, aes(x=metabolite, y=flux.data.col))
###
figb.13bpanel <- figb.13b + geom_jitter(aes(shape = metabolite, color = metabolite), position = position_jitterdodge(jitter.width = 0.9, dodge.width = 0.7), size = 5.5, stroke = 2) +
  geom_errorbar(data=figb.13bdf, mapping = aes(x=metabolite,ymin=flux.anc.data.col,ymax=flux.anc.data.col),color=anccolsb.13,lwd=2.2,linetype=117) +
  #geom_hline(yintercept = mean(filter(myhplc, condition_short=="GF", evo_trt_plotting=="Ancestor")$flux.ace.0), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = metabolite),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycolsb.13) +
  scale_y_continuous(limits = c(-0.5,20.35),breaks=c(0,5,10,15,20), expand = c(0,0)) +
  annotate("text",x=0.5,y=19.3,label="b",size=12,color="black")+
  xlab(" ")+
  ylab( expression(atop(paste("Flux (fmol/cell * ",hr^{-1},")"))))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))
figb.13bpanel

figb.13cdf<-read.csv("~/GitHub/Bifidobacterium/datafiles/fig.b.13c.csv")
figb.13cdf$metabolite<-recode(figb.13cdf$metabolite, "F" = "Fructose", "G" = "Glucose")
figb.13cdf$metabolite<-factor(figb.13cdf$metabolite, levels=c("Glucose","Fructose"))

figb.13c <- ggplot(figb.13cdf, aes(x=metabolite, y=flux.data.col))
###
figb.13cpanel <- figb.13c + geom_jitter(aes(shape = metabolite, color = metabolite), position = position_jitterdodge(jitter.width = 0.9, dodge.width = 0.7), size = 5.5, stroke = 2) +
  geom_errorbar(data=figb.13cdf, mapping = aes(x=metabolite,ymin=flux.anc.data.col,ymax=flux.anc.data.col),color=anccolsb.13,lwd=2.2,linetype=117) +
  #geom_hline(yintercept = mean(filter(myhplc, condition_short=="GF", evo_trt_plotting=="Ancestor")$flux.ace.0), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = metabolite),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycolsb.13) +
  scale_y_continuous(limits = c(-0.5,20.35),breaks=c(0,5,10,15,20), expand = c(0,0)) +
  annotate("text",x=0.5,y=19.3,label="c",size=12,color="black")+
  xlab(" ")+
  ylab( expression(atop(paste("Flux (fmol/cell * ",hr^{-1},")"))))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))
figb.13cpanel

grid.arrange(figb.13apanel,figb.13bpanel,figb.13cpanel,nrow=3)


#Make fig. 2.10
#It's a simple ANOVA figure with 3 columns. It has one wide grey dashed line for the anc value

fig2.10cols<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")

fig2.10df<-filter(myhplc, condition_short=="F" & evo_trt_plotting!="Ancestor")

fig2.10plot <- ggplot(fig2.10df, aes(x=evo_trt_plotting, y=flux.fru.pos))
fig2.10plot + geom_jitter(
  aes(shape = evo_trt_plotting, color = evo_trt_plotting),
  position = position_jitterdodge(jitter.width = 1.5, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  #geom_errorbar(data=grcvdataF, mapping = aes(x=carbon.tech.evo,ymin=umax.anc,ymax=umax.anc),color=grcv.colors,lwd=2.2,linetype=117) +
  geom_hline(yintercept = mean(filter(myhplc, condition_short=="F" & evo_trt_plotting=="Ancestor")$flux.fru.pos), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = evo_trt_plotting),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,19),breaks=c(0,5,10,15), expand = c(0,0)) +
  #labs(x="\nEvolution diet",y="Fructose flux (fmol/cell/hr)") +
  xlab("\nEvolution diet")+
  ylab( expression(atop(paste("Fructose flux (fmol/cell * ",hr^{-1},")"))))+
  annotate("text",x=0.7,y=13.5,label="A",size=12,color="black")+
  annotate("text",x=1.7,y=9.0,label="B",size=12,color="black")+
  annotate("text",x=2.7,y=11.4,label="A,B",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))



#Make fig. 2.11
#It's comparing the lacY vs. glcU F-evolved pops
fig2.11cols<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")

fig2.11df<-filter(myhplc, condition_short=="F" & evo_trt_plotting=="DP0")
fig2.11df$Importer<-recode(fig2.11df$Importer, "glcU" = "GlcU", "lacY" = "LacY")
fig2.11df$Importer<-factor(fig2.11df$Importer, levels=c("LacY", "GlcU"))


fig2.11plot <- ggplot(fig2.11df, aes(x=Importer, y=flux.fru.pos))
fig2.11plot + geom_jitter(
  aes(shape = Importer, color = Importer),
  position = position_jitterdodge(jitter.width = 1, dodge.width = 0.7),
  size = 5.5, stroke = 2) +
  #geom_errorbar(data=grcvdataF, mapping = aes(x=carbon.tech.evo,ymin=umax.anc,ymax=umax.anc),color=grcv.colors,lwd=2.2,linetype=117) +
  geom_hline(yintercept = mean(filter(myhplc, condition_short=="F" & evo_trt_plotting=="Ancestor")$flux.fru.pos), color="dark grey", lwd=2.2,linetype=117) +
  scale_shape_manual(values = c(0,0)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = Importer),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 3, stroke=2, shape=16,
    position = position_dodge(1),
  ) +
  scale_color_manual(values = mycols_breakout) +
  scale_y_continuous(limits = c(-0.04,19),breaks=c(0,5,10,15), expand = c(0,0)) +
  #labs(x="\nMutated importer",y="Fructose flux (fmol/cell/hr\n")) +
  xlab("\nMutated importer")+
  ylab( expression(atop(paste("Fructose flux (fmol/cell * ",hr^{-1},")"))))+
  #annotate("text",x=0.7,y=13.5,label="A",size=12,color="black")+
  #annotate("text",x=1.7,y=9.0,label="B",size=12,color="black")+
  #annotate("text",x=2.7,y=11.4,label="A,B",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 10, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 10, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 2.5), axis.ticks.x = element_line(color = "black", size = 2.5), axis.ticks.length.x=unit(0.5, "cm"),axis.ticks.length.y=unit(0.5, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))