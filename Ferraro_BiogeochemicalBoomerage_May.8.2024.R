#Kristy Ferraro
#Code for "The Biogeochemical Boomerang: Site Fidelity Creates Nutritional Hotspots that May Promote Recurrent Calving Site Reuse"
#Last Updated May 2024

#Libraries
library(lme4)
library(lmerTest)
library(tidyverse)
library(datarium)
library(rstatix)
library(ggpubr)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#Data Load----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
fogo_data<-read.csv(".../Ferraro_BiogeochemicalBoomerage_May.8.2024.csv")
fogo_no_controls<-subset(fogo_data, treat !="C")
fogo_no_base<-subset(fogo_data, time !="Baseline")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#Functions and Codes----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
r2.mixed<-function(mF){
  mFX<-model.matrix(mF)
  VarF <- var(as.vector(fixef(mF) %*% t(mFX)))
  VarR<-sum(as.numeric(VarCorr(mF))) 
  VarResid<-attr(VarCorr(mF), "sc")^2
  fR2<-VarF/(VarF + VarR + VarResid)
  rfR2<-(VarF + VarR)/(VarF + VarR + VarResid)
  list(fR2=fR2,rfR2=rfR2)
}

standardize_with_na <- function(x) {
  na_positions <- is.na(x)
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  standardized_x <- (x - mean_x) / (2 * sd_x)
  standardized_x[na_positions] <- NA
  return(standardized_x)
}

plot_theme<-theme(aspect.ratio = 1, panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent", color=NA), plot.background = element_rect(fill = "transparent", color=NA),   legend.background = element_rect(fill = "transparent", colour = NA),  legend.box.background = element_rect(fill = "transparent", colour = NA),   legend.key = element_rect(fill = "transparent", colour = NA), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5)) 
#legend.position = "none",

graph_points<- geom_point(aes(fill = Treatment, color = Treatment), size = 1, shape = 16, position = position_jitterdodge(seed = .1))
#point_colors<-scale_fill_brewer(palette="BuPu")
point_colors<-scale_color_manual(values = c("grey50", "steelblue4"))
fill_colors<-scale_fill_brewer(palette="BuPu")

pd = position_dodge(.76)
mean_bar<-stat_summary(fun=mean, geom="crossbar", linetype = "dashed", size=.3, width = .68, position=pd, show.legend=FALSE, color = "darkred") 

Theme2<-theme(aspect.ratio = 1, panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.background = element_rect(fill = "transparent", color=NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA),  axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#Models ----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#Forage N Stock----
#...Table 1 (forage N Stock)----
forage_N_model<-lmer(N_stock_forage_gm2 ~ time*treat + (1|block), data = fogo_data)
summary(forage_N_model)
r2.mixed(forage_N_model)

#....Fig 1  (forage N Stock)----
fogo_data_plots<- fogo_data %>%
  rename(Treatment = treat) %>%
  mutate(Treatment = recode(Treatment, "C" = "Control", "T" = "Treatment")) %>%
  mutate(time = recode(time, "Year1" = "2 Weeks", "Year2" = "1 Year"))

fogo_data_plots$time <- factor(fogo_data_plots$time, levels = c("Baseline", "2 Weeks", "1 Year"))
  
pd = position_dodge(.76)
forage_N_plot<-ggplot(fogo_data_plots, aes(x= time, y=N_stock_forage_gm2, fill = Treatment))  + scale_fill_brewer(palette="BuPu")+  geom_boxplot(outlier.shape = NA) +  plot_theme + graph_points   + labs(x ="Time Point", y =bquote('Forage g N' ~m^{"-2"})) + mean_bar + point_colors
quartz()
forage_N_plot

#= = = = = = = = == = = = == = = = == = = = == 
#Biomass----
#...Table 2 (biomass)----
biomass_model<-lmer(tot_forage_biomass_gm2 ~ time*treat + (1|block), data = fogo_no_base)
summary(biomass_model)
r2.mixed(biomass_model)

#....Fig 2 (biomass)----
fogo_no_base_plot<-fogo_no_base %>%
  rename(Treatment = treat) %>%
  mutate(Treatment = recode(Treatment, "C" = "Control", "T" = "Treatment")) %>%
  mutate(time = recode(time, "Year1" = "2 Weeks", "Year2" = "1 Year"))
fogo_no_base_plot$time <- factor(fogo_no_base_plot$time, levels = c("2 Weeks", "1 Year"))
forage_biomass_plot<-ggplot(fogo_no_base_plot, aes(x=time, y=tot_forage_biomass_gm2, fill = Treatment)) +  scale_colour_brewer(palette="BuPu") +  geom_boxplot(outlier.shape = NA) +  plot_theme + graph_points  + labs(x ="Time Point", y =bquote('Forage biomass (g'~m^{"-2"}*')' )) + mean_bar + point_colors + fill_colors +
  theme(legend.position="bottom")
quartz()
forage_biomass_plot

#= = = = = = = = == = = = == = = = == = = = == 
#%N-----
#...Table S1 (%N) and Fig 2 (%N)----
forage_precN_col<-  c("block", "treat", "time",  "lichen_precN", "crow_precN", "Kleaf_precN", "lingon_precN", "moss_precN")
forage_precN_data<- fogo_data[, forage_precN_col]
forage_precN_data <- forage_precN_data %>%
  rename(Lichen = lichen_precN, Crowberry = crow_precN, "Sheep Laurel Leaf" = Kleaf_precN, Lingonberry = lingon_precN, Moss = moss_precN, Treatment = treat)  %>%
mutate(Treatment = recode(Treatment, "C" = "Control", "T" = "Treatment")) %>%
  mutate(time = recode(time, "Year1" = "2 Weeks", "Year2" = "1 Year"))
forage_precN_data$time <- factor(forage_precN_data$time, levels = c("Baseline", "2 Weeks", "1 Year"))

forage_precN_data_long<- forage_precN_data %>% 
  pivot_longer( cols = `Lichen`:`Moss`,  names_to = "pool", values_to = "precN")

#Lichen
# Model
lichen_precN_model<-lmer(lichen_precN ~ time*treat + (1|block), data = fogo_data)
summary(lichen_precN_model)
r2.mixed(lichen_precN_model)

#Plot
lichen_precN_plot<-subset(forage_precN_data_long, pool == "Lichen")
lichen_precN_plot<- lichen_precN_plot %>% ggplot( aes(x = time, y = precN, fill = Treatment  ) ) + geom_boxplot(outlier.shape = NA) + plot_theme + graph_points  + labs(x = "Time Point", y = "% Nitrogen") +  mean_bar + point_colors + fill_colors+ theme(legend.position = "none", axis.title.x=element_blank())
quartz()
lichen_precN_plot


#Crow
#Model
crow_precN_model<-lmer(crow_precN ~ time*treat + (1|block), data = fogo_data)
summary(crow_precN_model)
r2.mixed(crow_precN_model)

#Plot
crow_precN_plot<-subset(forage_precN_data_long, pool == "Crowberry")
crow_precN_plot<- crow_precN_plot %>% ggplot( aes( x = time, y = precN, fill = Treatment  ) ) +  geom_boxplot(outlier.shape = NA) + plot_theme + graph_points  +  labs(x = "Time Point",  y = "% Nitrogen")  +  mean_bar + point_colors + fill_colors + theme(legend.position = "none", axis.title.x=element_blank()) 
quartz()
crow_precN_plot

#Kalmia
#Model
Kleaf_precN_model<-lmer(Kleaf_precN ~ time*treat + (1|block), data = fogo_data)
summary(Kleaf_precN_model)
r2.mixed(Kleaf_precN_model)

#Plot
kalmia_precN_plot<-subset(forage_precN_data_long, pool == "Sheep Laurel Leaf")
kalmia_precN_plot<- kalmia_precN_plot %>% ggplot( aes( x = time, y = precN, fill = Treatment  ) ) +  geom_boxplot(outlier.shape = NA) + plot_theme + graph_points  +  labs(x = "Time Point",  y = "% Nitrogen")  +  mean_bar + point_colors + fill_colors + theme(legend.position = "none", axis.title.x=element_blank()) 
quartz()
kalmia_precN_plot

#lingon
#Model
lingon_precN_model<-lmer(lingon_precN ~ time*treat + (1|block), data = fogo_data)
summary(lingon_precN_model)
r2.mixed(lingon_precN_model)

#Plot
ling_precN_plot<-subset(forage_precN_data_long, pool == "Lingonberry")
ling_precN_plot<- ling_precN_plot %>% ggplot( aes( x = time, y = precN, fill = Treatment  ) ) +  geom_boxplot(outlier.shape = NA) + plot_theme + graph_points  +  labs(x = "Time Point",  y = "% Nitrogen")  +  mean_bar + point_colors + fill_colors  + theme(legend.position = "none", axis.title.x=element_blank()) 
quartz()
ling_precN_plot

#moss
#Model
moss_precN_model<-lmer(moss_precN ~ time*treat + (1|block), data = fogo_data)
summary(moss_precN_model)
r2.mixed(moss_precN_model)

#Plot
moss_precN_plot<-subset(forage_precN_data_long, pool == "Moss")
moss_precN_plot<- moss_precN_plot %>% ggplot( aes( x = time, y = precN, fill = Treatment  ) ) +  geom_boxplot(outlier.shape = NA) + plot_theme + graph_points  +  labs(x = "Time Point",  y = "% Nitrogen")  +  mean_bar + point_colors + fill_colors + theme(legend.position = "none", axis.title.x=element_blank()) 
quartz()
moss_precN_plot


# = = = = = = = = = = = = = = = 
#...%N Increases ----
# = = = = = = = = = = = = = = = 
prec_N_col<-  c("block", "treat", "time", "soil_precN", "FR_precN", "CR_precN", "lichen_precN", "crow_precN", "Kleaf_precN", "Kstem_precN", "lingon_precN", "moss_precN","Oleaf_precN",  "Ostem_precN", "CR_precN", "FR_precN")
precN_data<- fogo_data[, prec_N_col]

precN_data_long<- precN_data %>% 
  pivot_longer(
    cols = `soil_precN`:`Ostem_precN`, 
    names_to = "plant",
    values_to = "precN"
)

plant_avs_from_control<-list()
for (i in seq_along(prec_N_col)) {
plant_subset<-subset(precN_data_long, plant == prec_N_col[i],)
 plant<- prec_N_col[i]
plant_control<-subset(plant_subset, treat == "C",)
plant_treat<-subset(plant_subset, treat == "T",)
plant_control_Yr1<-subset(plant_control, time == "Year1",)
plant_treat_Yr1<-subset(plant_treat, time == "Year1",)
ave_control_Yr1<-mean(plant_control_Yr1$precN, na.rm = TRUE)
ave_control_Yr1 <- round(ave_control_Yr1, 2)
std_control_Yr1<-sd(plant_control_Yr1$precN, na.rm = TRUE)
std_control_Yr1 <- round(std_control_Yr1, 2)
ave_treat_Yr1<-mean(plant_treat_Yr1$precN, na.rm = TRUE)
ave_treat_Yr1 <- round(ave_treat_Yr1, 2)
std_treat_Yr1<-sd(plant_treat_Yr1$precN, na.rm = TRUE)
std_treat_Yr1 <- round(std_treat_Yr1, 2)
prec_increase_Y1<-(ave_treat_Yr1-ave_control_Yr1)/ave_control_Yr1*100
prec_increase_Y1 <- round(prec_increase_Y1, 2)
times_increase_Y1<-(ave_treat_Yr1/ave_control_Yr1)
times_increase_Y1 <- round(times_increase_Y1, 2)
plant_control_Yr2<-subset(plant_control, time == "Year2",)
plant_treat_Yr2<-subset(plant_treat, time == "Year2",)
ave_control_Yr2<-mean(plant_control_Yr2$precN, na.rm = TRUE)
ave_control_Yr2 <- round(ave_control_Yr2, 2)
std_control_Yr2<-sd(plant_control_Yr2$precN, na.rm = TRUE)
std_control_Yr2 <- round(std_control_Yr2, 2)
ave_treat_Yr2<-mean(plant_treat_Yr2$precN, na.rm = TRUE)
ave_treat_Yr2 <- round(ave_treat_Yr2, 2)
std_treat_Yr2<-sd(plant_treat_Yr2$precN, na.rm = TRUE)
std_treat_Yr2 <- round(std_treat_Yr2, 2)
prec_increase_Y2<-(ave_treat_Yr2-ave_control_Yr2)/ave_control_Yr2*100
prec_increase_Y2 <- round(prec_increase_Y2, 2)
times_increase_Y2<-(ave_treat_Yr2/ave_control_Yr2)
times_increase_Y2 <- round(times_increase_Y2, 2)
  Output<- cbind(plant,  ave_control_Yr1, std_control_Yr1, ave_treat_Yr1, std_treat_Yr1, prec_increase_Y1, times_increase_Y1, ave_control_Yr2, std_control_Yr2,  ave_treat_Yr2, std_treat_Yr2, prec_increase_Y2, times_increase_Y2)
  plant_avs_from_control<- rbind(plant_avs_from_control, Output )
}
plant_avs_from_control<-as.data.frame(plant_avs_from_control)
plant_avs_from_control




# = = = = = = = = = = = = = = = 
#N recovery ----
# = = = = = = = = = = = = = = = 
#...Table 2 (N Recov)----

#Froage Recovery
forage_recoveryNmass_model<-lmer(forage_N_rec ~ time*treat + (1|block), data = fogo_no_base)
summary(forage_recoveryNmass_model)
r2.mixed(forage_recoveryNmass_model)

forage_recoveryNprec_model<-lmer(forage_prec_N_rec ~ time*treat + (1|block), data = fogo_no_base)
summary(forage_recoveryNprec_model)
r2.mixed(forage_recoveryNprec_model)

#Ecosystem Recovery 
eco_recoveryNmass_model<-lmer(total_eco_N_rec ~ time*treat + (1|block), data = fogo_no_base)
summary(eco_recoveryNmass_model)
r2.mixed(eco_recoveryNmass_model)

eco_recoveryNprec_model<-lmer(prec_eco_N_rec ~ time*treat + (1|block), data = fogo_no_base)
summary(eco_recoveryNprec_model)
r2.mixed(eco_recoveryNprec_model)


# # # # # # # # # # # # # # # 
#...Fig 3 Barchart of recovery----

fogo_year1_C<-subset(fogo_data, time == "Year1" & treat == "C")
fogo_year1_T<-subset(fogo_data, time == "Year1" & treat == "T")
fogo_year2_C<-subset(fogo_data, time == "Year2"  & treat == "C")
fogo_year2_T<-subset(fogo_data, time == "Year2" & treat == "T")

pool<- c("Soil", "Fine Root", "Coarse Root", "Lichen", "Crowberry", "Sheep Laurel Leaf", "Sheep Laurel Stem", "Lingonberry", "Moss", "Other Leaf", "Other Stem")
Year1<- c("Year1", "Year1", "Year1", "Year1", "Year1", "Year1", "Year1", "Year1", "Year1", "Year1", "Year1")
Year2<- c("Year2", "Year2", "Year2", "Year2", "Year2", "Year2", "Year2", "Year2", "Year2", "Year2", "Year2")
AveMassRecov<-c("AveMassRecov","AveMassRecov","AveMassRecov","AveMassRecov","AveMassRecov","AveMassRecov", "AveMassRecov","AveMassRecov","AveMassRecov", "AveMassRecov","AveMassRecov" )
AvePrecRecov<-c("AvePrecRecov","AvePrecRecov","AvePrecRecov","AvePrecRecov","AvePrecRecov","AvePrecRecov", "AvePrecRecov","AvePrecRecov","AvePrecRecov", "AvePrecRecov","AvePrecRecov" )
Treatment<- c("Treatment", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment")
Control<- c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control")


mass_ave_year1_control<-c(mean(fogo_year1_C$soil_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_C$FR_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_C$CR_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_C$lichen_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_C$crow_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_C$Kleaf_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_C$Kstem_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_C$lingon_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_C$moss_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_C$Oleaf_mass_sol_rec), mean(fogo_year1_C$Ostem_mass_sol_rec, na.rm = TRUE) )
mass_ave_year1_control<-cbind(pool, Year1, AveMassRecov, mass_ave_year1_control, Control)
colnames(mass_ave_year1_control) <- c("Pool","Time","Type", "Amount", "Treatment")

mass_ave_year1_treat<-c(mean(fogo_year1_T$soil_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_T$FR_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_T$CR_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_T$lichen_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_T$crow_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_T$Kleaf_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_T$Kstem_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_T$lingon_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_T$moss_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_T$Oleaf_mass_sol_rec, na.rm = TRUE), mean(fogo_year1_T$Ostem_mass_sol_rec, na.rm = TRUE) )
mass_ave_year1_treat<-cbind(pool, Year1, AveMassRecov, mass_ave_year1_treat, Treatment)
colnames(mass_ave_year1_treat) <- c("Pool","Time","Type", "Amount", "Treatment")


mass_ave_year2_control<-c(mean(fogo_year2_C$soil_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_C$FR_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_C$CR_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_C$lichen_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_C$crow_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_C$Kleaf_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_C$Kstem_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_C$lingon_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_C$moss_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_C$Oleaf_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_C$Ostem_mass_sol_rec, na.rm = TRUE) )
mass_ave_year2_control<-cbind(pool, Year2, AveMassRecov, mass_ave_year2_control, Control)
colnames(mass_ave_year2_control) <- c("Pool","Time","Type", "Amount", "Treatment")


mass_ave_year2_treat<-c(mean(fogo_year2_T$soil_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_T$FR_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_T$CR_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_T$lichen_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_T$crow_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_T$Kleaf_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_T$Kstem_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_T$lingon_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_T$moss_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_T$Oleaf_mass_sol_rec, na.rm = TRUE), mean(fogo_year2_T$Ostem_mass_sol_rec, na.rm = TRUE) )
mass_ave_year2_treat<-cbind(pool, Year2, AveMassRecov, mass_ave_year2_treat, Treatment)
colnames(mass_ave_year2_treat) <- c("Pool","Time","Type", "Amount", "Treatment")


prec_ave_year1_control<-c(mean(fogo_year1_C$soil_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_C$FR_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_C$CR_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_C$lichen_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_C$crow_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_C$Kleaf_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_C$Kstem_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_C$lingon_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_C$moss_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_C$Oleaf_prec_sol_rec), mean(fogo_year1_C$Ostem_prec_sol_rec, na.rm = TRUE) )
prec_ave_year1_control<-cbind(pool, Year1, AvePrecRecov, prec_ave_year1_control, Control)
colnames(prec_ave_year1_control) <- c("Pool","Time","Type", "Amount", "Treatment")

prec_ave_year1_treat<-c(mean(fogo_year1_T$soil_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_T$FR_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_T$CR_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_T$lichen_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_T$crow_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_T$Kleaf_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_T$Kstem_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_T$lingon_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_T$moss_prec_sol_rec, na.rm = TRUE), mean(fogo_year1_T$Oleaf_prec_sol_rec), mean(fogo_year1_T$Ostem_prec_sol_rec, na.rm = TRUE) )
prec_ave_year1_treat<-cbind(pool, Year1, AvePrecRecov, prec_ave_year1_treat, Treatment)
colnames(prec_ave_year1_treat) <- c("Pool","Time","Type", "Amount", "Treatment")

prec_ave_year2_control<-c(mean(fogo_year2_C$soil_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_C$FR_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_C$CR_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_C$lichen_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_C$crow_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_C$Kleaf_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_C$Kstem_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_C$lingon_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_C$moss_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_C$Oleaf_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_C$Ostem_prec_sol_rec, na.rm = TRUE) )
prec_ave_year2_control<-cbind(pool, Year2, AvePrecRecov, prec_ave_year2_control, Control)
colnames(prec_ave_year2_control) <- c("Pool","Time","Type", "Amount", "Treatment")


prec_ave_year2_treat<-c(mean(fogo_year2_T$soil_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_T$FR_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_T$CR_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_T$lichen_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_T$crow_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_T$Kleaf_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_T$Kstem_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_T$lingon_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_T$moss_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_T$Oleaf_prec_sol_rec, na.rm = TRUE), mean(fogo_year2_T$Ostem_prec_sol_rec, na.rm = TRUE) )
prec_ave_year2_treat<-cbind(pool, Year2, AvePrecRecov, prec_ave_year2_treat, Treatment)
colnames(prec_ave_year2_treat) <- c("Pool","Time","Type", "Amount", "Treatment")


NRecovData<-rbind(mass_ave_year1_control, mass_ave_year1_treat,mass_ave_year2_control, mass_ave_year2_treat,prec_ave_year1_control,  prec_ave_year1_treat, prec_ave_year2_control, prec_ave_year2_treat)
NRecovData.df<-as.data.frame(NRecovData)
NRecovData.df<-NRecovData.df %>% arrange(Type, Pool, Time, Amount)
NRecovData.df<-as.data.frame(NRecovData.df)
NRecovData.df$Amount<-as.numeric(as.character(sub("," , ".", NRecovData.df$Amount)))

NRecovData.df$Pool <- factor(NRecovData.df$Pool,  levels = c("Moss", "Lichen", "Crowberry", "Lingonberry", "Sheep Laurel Leaf", "Other Leaf", "Sheep Laurel Stem", "Other Stem", "Fine Root", "Coarse Root", "Soil" ))

custom_colors <- c("#669933", "#999999", "#003333", "#660000", "#006633",  "#cccc00", "#999900", "#333300", "#996633", "#663300", "#330000")


#......All Recovery----
Treatment_NRrecov.plot<-ggplot(NRecovData.df, aes(fill=Pool, y=Amount, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = custom_colors) +
#    geom_text(aes( label = Amount), colour = "black") +
    facet_grid(Type ~ Time, scales = "free", 
             switch = "y",
             labeller = as_labeller(c(AveMassRecov = "gN recovered", AvePrecRecov = "%N recovered", Year1 = "2 Weeks", Year2 = "1 Year") ) ) + ylab(NULL) + xlab(NULL) + 
  scale_y_continuous(position = "left")+ theme(strip.background = element_blank()) +
  theme(strip.placement = "outside")  + Theme2 
quartz()

#......Forage Recovery----
forage_NRecovData.df<- subset(NRecovData.df, Pool != "Coarse Root" & Pool != "Fine Root" & Pool != "Soil" & Pool != "Sheep Laurel Stem" & Pool != "Other Stem" & Pool != "Other Leaf" )
custom_colors_forage <- c("#669933", "#999999", "#003333", "#660000", "#006633")

Forage_NRrecov.plot<-ggplot(forage_NRecovData.df, aes(fill=Pool, y=Amount, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +  scale_fill_manual(values = custom_colors_forage) +
#    geom_text(aes( label = Amount), colour = "black") +
    facet_grid(Type ~ Time, scales = "free", 
             switch = "y",
             labeller = as_labeller(c(AveMassRecov = "gN recovered", AvePrecRecov = "%N recovered", Year1 = "2 Weeks", Year2 = "1 Year") ) ) + ylab(NULL) + xlab(NULL) + 
  scale_y_continuous(position = "left")+ theme(strip.background = element_blank()) +
  theme(strip.placement = "outside")  + Theme2
quartz()




# = = = = = = = = = = = = = = = 
#Distance ----
# = = = = = = = = = = = = = = = 
dist_data<-read.csv(".../Ferraro_DistanceSamples_2024.csv")
dist_data <- na.omit(dist_data)
dist_data$C.N<-as.numeric(dist_data$C.N)

dist_data_base<-subset(dist_data, dist2 != "0")
dist_data_base <- dist_data_base %>%
  mutate(dist2 = ifelse(dist2 == 5, 0, dist2))
dist_data_base <- dist_data_base %>%
  mutate(dist = ifelse(dist == 5, 0, dist))
dist_data_base$block<-as.factor(dist_data_base$block)

#....Lichen-----
lichen_dist_base<-subset(dist_data_base, plant == "Lichen")
#...Table S2 (forage)----

#precN
precN_lichen_dist_model_base<-lmer(precN ~ dist2  + (1|block), data = lichen_dist_base)
summary(precN_lichen_dist_model_base)
r2.mixed(precN_lichen_dist_model_base)

#d15N
d15N_lichen_dist_model_base<-lmer(log(d15) ~ dist2  + (1|block), data = lichen_dist_base)
d15N_lichen_dist_model_base<-lmer(d15 ~ dist2  + (1|block), data = lichen_dist_base)
summary(d15N_lichen_dist_model_base)
r2.mixed(d15N_lichen_dist_model_base)

#CN
CN_lichen_dist_model_base<-lmer(C.N ~ dist2  + (1|block), data = lichen_dist_base)
summary(CN_lichen_dist_model_base)
r2.mixed(CN_lichen_dist_model_base)

#....Kalmia-----
kalmia_dist_base<-subset(dist_data_base, plant == "Kalmia")
#precN v base 
precN_kalmia_dist_model_base<-lmer(precN ~ dist2  + (1|block), data = kalmia_dist_base)
summary(precN_kalmia_dist_model_base)
r2.mixed(precN_kalmia_dist_model_base)

#d15N v base 
d15N_kalmia_dist_model_base<-lmer(d15 ~ dist2  + (1|block), data = kalmia_dist_base)
summary(d15N_kalmia_dist_model_base)
r2.mixed(d15N_kalmia_dist_model_base)

#CN v base 
CN_kalmia_dist_model_base<-lmer(C.N ~ dist2  + (1|block), data = kalmia_dist_base)
summary(CN_kalmia_dist_model_base)
r2.mixed(CN_kalmia_dist_model_base)



# ....Fig Supp 5----
#Kalmia
kalmia_dist_base <- kalmia_dist_base %>% rename(Distance = dist2)  
kalmia_dist_base$Distance <-as.factor(kalmia_dist_base$Distance) 

graph_points_2<- geom_point(aes(fill = Distance, color = Distance), size = 1, shape = 16, position = position_jitterdodge(seed = .1)) 
point_colors<-scale_fill_brewer(palette="BuPu")
fill_colors<-scale_fill_brewer(palette="BuPu")

Kalmia_precN<- kalmia_dist_base %>% ggplot( aes( x = Distance, y = precN, fill = Distance)) + 
  geom_boxplot()  +  scale_x_discrete(labels=c('Control', '0.25m', '0.50m', "0.75m", "1m"))  + labs(x = "", y="Sheep Laurel %N") + graph_points_2  + theme(legend.position = "none")  +  scale_colour_brewer(palette="BuPu") +  scale_fill_brewer(palette="BuPu") + plot_theme  + mean_bar
Kalmia_d15N<- kalmia_dist_base %>% ggplot( aes( x = Distance, y = d15, fill = Distance)) + 
  geom_boxplot()  +  scale_x_discrete(labels=c('Control', '0.25m', '0.50m', "0.75m", "1m")) + labs(x = "",y=bquote( "Sheep Laurel" ~ delta ^{15} ~ "N"))+ theme(legend.position = "none") + graph_points_2  + theme(legend.position = "none")  +  scale_colour_brewer(palette="BuPu") +  scale_fill_brewer(palette="BuPu") + plot_theme  + mean_bar
Kalmia_CN<- kalmia_dist_base %>% ggplot( aes( x = Distance, y = C.N, fill = Distance)) + 
  geom_boxplot()  +  scale_x_discrete(labels=c('Control', '0.25m', '0.50m', "0.75m", "1m")) + labs(x = "", y="Sheep Laurel C:N")+ theme(legend.position = "none") + graph_points_2  + theme(legend.position = "none")  +  scale_colour_brewer(palette="BuPu") +  scale_fill_brewer(palette="BuPu")+ plot_theme  + mean_bar


lichen_dist_base <- lichen_dist_base %>% rename(Distance = dist2)  
lichen_dist_base$Distance <-as.factor(lichen_dist_base$Distance) 

Lichen_precN<- lichen_dist_base %>% ggplot( aes( x = Distance, y = precN, fill = Distance)) + 
  geom_boxplot()  +  scale_x_discrete(labels=c('Control', '0.25m', '0.50m', "0.75m", "1m")) + labs(x = "", y="Lichen %N")+ theme(legend.position = "none") + graph_points_2  + theme(legend.position = "none")  +  scale_colour_brewer(palette="BuPu") +  scale_fill_brewer(palette="BuPu")+ plot_theme  + mean_bar
Lichen_d15N<- lichen_dist_base %>% ggplot( aes( x = Distance, y = d15, fill = Distance)) + 
  geom_boxplot()  +  scale_x_discrete(labels=c('Control', '0.25m', '0.50m', "0.75m", "1m")) +  labs(x = "",y=bquote( "Lichen" ~ delta ^{15} ~ "N"))+ theme(legend.position = "none") + graph_points_2  + theme(legend.position = "none")  +  scale_colour_brewer(palette="BuPu") +  scale_fill_brewer(palette="BuPu")+ plot_theme  + mean_bar
Lichen_CN<- lichen_dist_base %>% ggplot( aes( x = Distance, y = C.N, fill = Distance)) + 
  geom_boxplot()  +  scale_x_discrete(labels=c('Control', '0.25m', '0.50m', "0.75m", "1m")) + labs(x = "", y="Lichen C:N")+ theme(legend.position = "none") + graph_points_2  + theme(legend.position = "none")  +  scale_colour_brewer(palette="BuPu") +  scale_fill_brewer(palette="BuPu")+ plot_theme + mean_bar


dist_plot <- ggarrange( Lichen_precN, Lichen_d15N, Lichen_CN, Kalmia_precN, Kalmia_d15N,  Kalmia_CN, common.legend = FALSE, heights = c(3,2)) 

quartz()
dist_plot




# = = = = = = = = = = = = = = = 
##Other Supp Figs  -----
# = = = = = = = = = = = = = = = 

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#.S1 - Biomass----
biomass_cols<-c("block", "treat", "time", "FR_gplot","CR_gm2", "lichen_gm2", "crow_gm2", "Kleaf_gm2", "Kstem_gm2",  "lingon_gm2", "moss_gm2",  "Oleaf_gm2",  "Ostem_gm2")

fogo_biomass<- fogo_no_base[, biomass_cols]

fogo_biomass<- fogo_biomass %>%
  rename(Treatment = treat, Lichen = lichen_gm2, Crowberry = crow_gm2, "Sheep Laurel Leaf" = Kleaf_gm2, "Sheep Laurel Stem" = Kstem_gm2, Lingonberry = lingon_gm2, Moss = moss_gm2, "Other Leaf" = Oleaf_gm2, "Other Stem" = Ostem_gm2, "Fine Root" = FR_gplot, "Coarse Root" = CR_gm2) %>%
  mutate(Treatment = recode(Treatment, "C" = "Control", "T" = "Treatment")) %>%
  mutate(time = recode(time, "Year1" = "2 Weeks", "Year2" = "1 Year"))

fogo_biomass$time <- factor(fogo_biomass$time, levels = c("Baseline", "2 Weeks", "1 Year"))

plant_biomass_long<- fogo_biomass %>% 
  pivot_longer( cols = "Fine Root":"Other Stem", names_to = "plant", values_to = "biomass")

plant_biomass_long$plant <- factor(plant_biomass_long$plant, levels = c("Lichen", "Crowberry", "Lingonberry", "Moss", "Sheep Laurel Leaf", "Sheep Laurel Stem", "Other Leaf", "Other Stem", "Fine Root", "Coarse Root"))

all_biomass<- plant_biomass_long %>%ggplot( aes( x = time, y = biomass, fill = Treatment  ) ) +  scale_colour_brewer(palette="BuPu") +  geom_boxplot() +  plot_theme + graph_points + scale_fill_brewer(palette="BuPu") + mean_bar+ facet_wrap(as.factor(plant_biomass_long$plant), scales = "free") + labs(x ="Time Point", y =bquote('Ecosystem Pool Biomass (g' ~m^{"-2"}*')')) 

quartz()
all_biomass


#.S2 & S3 - C&N Stock Plot ----

#nitrogen
fogo_stocks_N<- fogo_data %>%
dplyr::select(block, treat, time,soil_N_gm2, FR_N_gm2, CR_N_gm2, lichen_N_gm2, crow_N_gm2, Kleaf_N_gm2, Kstem_N_gm2,lingon_N_gm2, moss_N_gm2, Oleaf_N_gm2, Ostem_N_gm2 ) 

fogo_stocks_N$soil_N_gm2<-as.numeric(fogo_stocks_N$soil_N_gm2)
fogo_stocks_N <- fogo_stocks_N %>%
  rename(Treatment = treat, Lichen = lichen_N_gm2, Crowberry = crow_N_gm2, "Sheep Laurel Leaf" = Kleaf_N_gm2, "Sheep Laurel Stem" = Kstem_N_gm2, Lingonberry = lingon_N_gm2, Moss = moss_N_gm2, "Fine Root" = FR_N_gm2, "Coarse Root" = CR_N_gm2,  Soil = soil_N_gm2, "Other Leaf" = Oleaf_N_gm2, "Other Stem" = Ostem_N_gm2 )  %>%
  mutate(Treatment = recode(Treatment, "C" = "Control", "T" = "Treatment")) %>%
  mutate(time = recode(time, "Year1" = "2 Weeks", "Year2" = "1 Year"))

fogo_stocks_N$time <- factor(fogo_stocks_N$time, levels = c("Baseline", "2 Weeks", "1 Year"))
fogo_stocks_N<-subset(fogo_stocks_N, Treatment == "Control" )
fogo_stocks_N_long<- fogo_stocks_N %>%  pivot_longer(cols = `Soil`:`Other Stem`, names_to = "plant", values_to = "gNm2")

fogo_stocks_N_long$plant <- factor(fogo_stocks_N_long$plant, levels = c("Lichen", "Crowberry", "Lingonberry", "Moss", "Sheep Laurel Leaf", "Sheep Laurel Stem", "Other Leaf", "Other Stem", "Fine Root", "Coarse Root", "Soil"))

nitrogen_stocks<- fogo_stocks_N_long %>% ggplot( aes(fill=Treatment, y=gNm2, x=time) ) + 
  geom_boxplot() +  plot_theme + graph_points + scale_fill_brewer(palette="BuPu") + mean_bar +
  labs(x = "Sampling Period",
       y = "N g/m2") +  scale_colour_brewer(palette="BuPu") +
  facet_wrap(~ plant, scales = "free") + labs(x ="Time Point", y =bquote('Nitrogen Stock (g N' ~m^{"-2"}*')')) 

quartz()



#carbon
fogo_stocks_C<- fogo_data %>%
dplyr::select(block, treat, time,soil_C_gm2, FR_C_gm2, CR_C_gm2, lichen_C_gm2, crow_C_gm2, Kleaf_C_gm2, Kstem_C_gm2,lingon_C_gm2, moss_C_gm2, Oleaf_C_gm2, Ostem_C_gm2 ) 
fogo_stocks_C$soil_C_gm2<-as.numeric(fogo_stocks_C$soil_C_gm2)

fogo_stocks_C <- fogo_stocks_C %>%
  rename(Treatment = treat, Lichen = lichen_C_gm2, Crowberry = crow_C_gm2, "Sheep Laurel Leaf" = Kleaf_C_gm2, "Sheep Laurel Stem" = Kstem_C_gm2, Lingonberry = lingon_C_gm2, Moss = moss_C_gm2, "Fine Root" = FR_C_gm2, "Coarse Root" = CR_C_gm2, Soil = soil_C_gm2, "Other Leaf" = Oleaf_C_gm2, "Other Stem" = Ostem_C_gm2 ) %>%
  mutate(Treatment = recode(Treatment, "C" = "Control", "T" = "Treatment")) %>%
  mutate(time = recode(time, "Year1" = "2 Weeks", "Year2" = "1 Year"))

fogo_stocks_C$time <- factor(fogo_stocks_C$time, levels = c("Baseline", "2 Weeks", "1 Year"))

fogo_stocks_C<-subset(fogo_stocks_C, Treatment == "Control" )

fogo_stocks_C_long<- fogo_stocks_C %>% 
  pivot_longer(
    cols = `Soil`:`Other Stem`, 
    names_to = "plant",
    values_to = "gCm2")

fogo_stocks_C_long$plant <- factor(fogo_stocks_C_long$plant, levels = c("Lichen", "Crowberry", "Lingonberry", "Moss", "Sheep Laurel Leaf", "Sheep Laurel Stem", "Other Leaf", "Other Stem", "Fine Root", "Coarse Root", "Soil"))

carbon_stocks<- fogo_stocks_C_long %>% ggplot( aes(fill=Treatment, y=gCm2, x=time) ) +  geom_boxplot() +  plot_theme + graph_points + scale_fill_brewer(palette="BuPu") + mean_bar+ labs(x = "Sampling Period",  y = "C g/m2") +  scale_colour_brewer(palette="BuPu") + facet_wrap(~ plant, scales = "free") + labs(x ="Time Point", y =bquote('Carbon Stock (g C' ~m^{"-2"}*')')) 

quartz()



#.S4 - Non-Forage %N----
non_forage_precN_col<-  c("block", "treat", "time",  "Kstem_precN", "Oleaf_precN", "Ostem_precN", "FR_precN", "CR_precN", "soil_precN")

non_forage_precN_data<- fogo_data[, non_forage_precN_col]
non_forage_precN_data <- non_forage_precN_data %>%
  rename("Sheep Laurel Stem" = Kstem_precN, "Other Leaf" = Oleaf_precN, "Other Stem" = Ostem_precN, "Fine Root" = FR_precN, "Coarse Root" = CR_precN,  Soil = soil_precN, Treatment = treat)  %>%
mutate(Treatment = recode(Treatment, "C" = "Control", "T" = "Treatment")) %>%
  mutate(time = recode(time, "Year1" = "2 Weeks", "Year2" = "1 Year"))

non_forage_precN_data$time <- factor(non_forage_precN_data$time, levels = c("Baseline", "2 Weeks", "1 Year"))

non_forage_precN_data_long<- non_forage_precN_data %>% pivot_longer( cols = "Sheep Laurel Stem" :`Soil`,  names_to = "pool", values_to = "precN")

non_forage_precN_plot<- non_forage_precN_data_long %>% ggplot( aes( x = time, y = precN, fill = Treatment  ) ) +  geom_boxplot(outlier.shape = NA) + plot_theme + graph_points + scale_fill_brewer(palette="BuPu") + labs(x = "Time Point", y = "% Nitrogen") +  scale_colour_brewer(palette="BuPu") +  mean_bar + point_colors + fill_colors + facet_wrap(~ factor(pool, levels=c('Sheep Laurel Stem', 'Other Leaf', 'Other Stem', 'Fine Root', 'Coarse Root', 'Soil')), scales = "free") 

quartz()
non_forage_precN_plot

