######################################################
##             Code for Analysis and Figure Creation           ##
##      Anderegg et al. 2024                                   ##
##      "Deep roots and complex phenology decouple tree water 
##      stress from climate and growth in a xeric oak"
######################################################


#load packages
library(RColorBrewer)
library(khroma)
library(dplyr)
library(reshape)
library(nlme)
library(lme4)
library(lmerTest)
library(MuMIn)
library(tidyverse)
require(lmodel2)
#require(rgdal) # getting retired 2023 and can no longer download
require(car)
require(quantreg)
# for downloading TerraClimate data
# remotes::install_github("mikejohnson51/AOI") # suggested!
# remotes::install_github("mikejohnson51/climateR")
library(AOI)
library(climateR)
library(terra)
library(khroma)

# load in default palette and some opacity versions
mypal <- c(brewer.pal(n=9, "Set1"), brewer.pal(n=8, "Dark2"))
palette(mypal)
pal2 <- brewer.pal(n=8, "Set2")
pal2light <- paste0(pal2,"33")

# also load to colorblind friendly palettes
# The palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

se <- function(x){
  se <- sd(x, na.rm=T)/sqrt(length(x[!is.na(x)]))
  return(se)
}


### Save Results?

save.figures <- T # whether to save figure pdfs
results.version <- "v240604" # full new version after code update
results.dir <- paste0("Results_",results.version)
if(save.figures == T) { dir.create(results.dir)}




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#######   BEGIN: LOAD DATA ###################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


popclim <- read.csv("DerivedData/PopulationClimate_230512.csv")[,-1]
  # this includes climate data from the Basin Characterisation Model (BCM)
  # all of Prahlad Papper's phenology sites
pop.terraclim <- read.csv("DerivedData/Population_TerraClimate_230512.csv")[,-1]
  # only 15 sites + dropped site from spring (lost access)

#______________________________________________________________________
#######   * Soils data ###################################################
#______________________________________________________________________
# Extracted using SoilWeb (NRCS) SSURGO

soils <- read.csv("DerivedData/BOINK_soil_data_v1 _revisednames.csv")



#______________________________________________________________________
#######   * Blue Oak herbarium records ###################################################
#______________________________________________________________________
# NOTE: These location data were extracted from the supplemental data in:
#Baldwin BG, Thornhill AH, Freyman WA, Ackerly DD, Kling MM, Morueta-Holme N, Mishler BD. 2017. Species richness and endemism in the native flora of California. American Journal of Botany 104: 487–501.

# you can download them and subset the data to Quercus douglasii from the Dryad repository here:
# https://dx.doi.org/10.6078/D16K5W

#_________ Climate extraction (skip to derived data)____________
# extracted from Baldwin et al. 2017
# qudo.raw <- read.csv("/Users/leeanderegg/Dropbox/CO2 Fertilization/herbarium data/Baldwin2017 PhylodiversityData/California_QuercusDouglasii_specimens_20180318.csv")
# 
# # remove duplicated records
# dups <- duplicated(qudo.raw %>% select(latitude, longitude))
# sum(dups) # number of duplicate records (which would incorrectly weight their locations in the model)
# qudo<- qudo.raw[!dups,] # remove the duplicates using the ! (NOT) boolian function
# 
# # clean up the collection date
# first <- str_extract(qudo$early_julian_day, "\\b[:digit:]+") # grab the first part of the date string
# # this is either month (mm/dd/yyyy) or year (yyyy-mm-dd)
# last <- str_extract(qudo$early_julian_day, "[:digit:]+$\\b")
# year <- rep(NA, times=nrow(qudo))
# year[which(nchar(first)==4)] <- first[which(nchar(first)==4)]
# year[which(nchar(last)==4)] <- last[which(nchar(last)==4)]
# qudo$year <- as.numeric(year)
# 
# # extract climate data for oak locations
# calbound = aoi_get(state = "CA")
# # download terrclimate normals 81-2010
# tcn <- getTerraClimNormals(calbound, varname=c("aet","def","pet","ppt","soil","srad","swe","tmax","tmin","vap","ws","vpd"),scenario = "19812010")
# # monthly 30yr normals for whole state, used to visualize the climate envelope of blue oak
# 
# ## summarize monthly into annual/seasonal variables
# tcn_an <- list()
# tcn_an$aet <- sum(tcn$aet) # total actual annual evap (mm)
# tcn_an$cwd <- sum(tcn$def) # total CWD (mm)
# tcn_an$cwd_gs <- mean(tcn$def[[4:9]]) # mean CWD of the growing season (Apr - Sept)
# tcn_an$pet <- sum(tcn$pet) # total PET (mm)
# tcn_an$ppt <- sum(tcn$ppt) # total precip (mm)
# tcn_an$soil_mean <- mean(tcn$soil) # mean annual soil moisture
# tcn_an$soil_min <- min(tcn$soil) # min soil moisture
# tcn_an$srad <- mean(tcn$srad) # mean shortwave
# tcn_an$swe_max <- max(tcn$swe) # max swee
# tcn_an$tmax <- max(tcn$tmax) # max T of hottest month
# tcn_an$tmin <- min(tcn$tmin) # min T of coldest month
# tcn_an$vpd_max <- max(tcn$vpd) # max VPD
# tcn_an$vpd_spr <- mean(tcn$vpd[[4:5]]) # mean VPD of spring
# tcn_an$vpd_gs <- mean(tcn$vpd[[4:9]]) # mean VPD of growing season
# 
# 
# # extract annual climate normals for oak occurance locations
# qudo$aet <- extract(tcn_an$aet, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$cwd <- extract(tcn_an$cwd, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$cwd_gs <- extract(tcn_an$cwd_gs, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$pet <- extract(tcn_an$pet, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$ppt <- extract(tcn_an$ppt, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$soil_mean <- extract(tcn_an$soil_mean, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$soil_min <- extract(tcn_an$soil_min, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$srad <- extract(tcn_an$srad, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$swe_max <- extract(tcn_an$swe_max, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$tmax <- extract(tcn_an$tmax, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$tmin <- extract(tcn_an$tmin, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$vpd_max <- extract(tcn_an$vpd_max, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$vpd_spr <- extract(tcn_an$vpd_spr, data.frame(qudo$longitude,qudo$latitude))[,2]
# qudo$vpd_gs <- extract(tcn_an$vpd_gs, data.frame(qudo$longitude,qudo$latitude))[,2]

# in case need data offline
#write.csv(qudo, paste("DerivedData/CA_QUDO_herbariumspecimens_wTerraClimate", results.version,".csv", sep=""))
#_______________________

# NOT included in Dryad repo for this paper, because just a subset of https://dx.doi.org/10.6078/D16K5W
qudo <- read.csv("DerivedData/CA_QUDO_herbariumspecimens_wTerraClimatev240216.csv")


#______________________________________________________________________
#######   * Water Potential and trait Data ###################################################
#______________________________________________________________________

# individual average data, from ~3 measurements per individual for traits and 3-8 measurements for water potentials
wp.ind <- read.csv("DerivedData/Fall18_WP_traits_growth_isos_230526.csv", row.names = 1)
# create log10-transformed variables for things that appear really right skewed
wp.ind$logAl_As <- log(wp.ind$mAl_As, base=10)
wp.ind$logml_ms <- log(wp.ind$mml_ms, base=10)
wp.ind$logleafsize <- log(wp.ind$mleafsize, base=10)
wp.ind$logLength <- log(wp.ind$mLength, base=10)
wp.ind$Site <- factor(wp.ind$Site)

# calculate the residaul of stem delD from meteoric delD
wp.ind$delD.resid <- NA
wp.ind$delD.resid[which(wp.ind$delD<0)] <- resid(lm(delD~wydelD, wp.ind))

# make percent of max BAI into percent rather than ratio
wp.ind$perc_maxBAI <- wp.ind$perc_maxBAI*100

soil.simple <- soils %>% filter(ChoiceOrder==1) %>% select(tree_correct, site_correct, parent.material.ssurgo = parent.material, soil.name, water.storage100cm = available.water.storage..0.100cm., PAW=total.plant.available.water.cm., bedrock.depth=min.bedrock.depth..cm..if.available )
soil.simple$bedrock.depth[which(soil.simple$bedrock.depth=="na")] <- NA
soil.simple$PAW[which(soil.simple$PAW=="na")] <- NA
soil.simple$bedrock.depth <- as.numeric(soil.simple$bedrock.depth)
soil.simple$PAW <- as.numeric(soil.simple$PAW)
soil.simple$bedrock.depth[which(soil.simple$site_correct=="DAL")] <- 38 #for soil types without bedrock depths, extracted from closest nearby reasonable
soil.simple$bedrock.depth[which(soil.simple$site_correct=="HOP")] <- 152 

wp.ind <- left_join(wp.ind, soil.simple, by=c("Site"="site_correct","Tree"="tree_correct"))
wp.ind$Site <- factor(wp.ind$Site)



# calculate leaf and stem hydraulic safety margins:
# From Skelton et al. 2019
  # P50_leaf = -3.88 MPa
  # P50_stem = -4.47 MPa

wp.ind$HSM_leaf <-  wp.ind$MD + 3.88
wp.ind$HSM_stem <- (.92*wp.ind$MD) + 4.47# reduce leaf wp by 8% to estimate branch wp



#______________________________________________________________________
#######   * Hotter Simulations ###################################################
#______________________________________________________________________

forcings <- read.csv("DerivedData/HotterSimulations_240627/met_blueOaks_clean.csv")
maxvpd <- read.csv("DerivedData/HotterSimulations_240627/hotter_demo_output_maxVPD_stdAl.csv")[,-1]
maxvpd.lowleaf <- read.csv("DerivedData/HotterSimulations_240627/hotter_demo_output_maxVPD_AlLow.csv")[,-1]
meanvpd <- read.csv("DerivedData/HotterSimulations_240627/hotter_demo_output_meanVPD_stdAl.csv")[,-1]
meanvpd.lowleaf <- read.csv("DerivedData/HotterSimulations_240627/hotter_demo_output_meanVPD_AlLow.csv")[,-1]
minvpd <- read.csv("DerivedData/HotterSimulations_240627/hotter_demo_output_minVPD_stdAl.csv")[,-1]
minvpd.lowleaf <- read.csv("DerivedData/HotterSimulations_240627/hotter_demo_output_minVPD_AlLow.csv")[,-1]
colnames(maxvpd)[1:6] <- paste(colnames(maxvpd)[1:6], "max", sep="_")
colnames(meanvpd)[1:6] <- paste(colnames(meanvpd)[1:6], "mean", sep="_")
colnames(minvpd)[1:6] <- paste(colnames(minvpd)[1:6], "min", sep="_")
colnames(maxvpd.lowleaf)[1:6] <- paste(colnames(maxvpd.lowleaf)[1:6], "max","ll", sep="_")
colnames(meanvpd.lowleaf)[1:6] <- paste(colnames(meanvpd.lowleaf)[1:6], "mean","ll", sep="_")
colnames(minvpd.lowleaf)[1:6] <- paste(colnames(minvpd.lowleaf)[1:6], "min","ll", sep="_")


# combine all the HOTTER Simulation outputs into one dataframe
hotsim <- full_join(full_join(full_join(full_join(full_join(full_join(forcings, maxvpd), meanvpd), minvpd),maxvpd.lowleaf),meanvpd.lowleaf),minvpd.lowleaf)

# join them with the indivdiual level observations
wp.ind <- full_join(wp.ind, hotsim[,-4]) # gotta keep mAl_As from screwing me with rounding error in hotsim
wp.ind$PLC_mean <- 100 - 100*wp.ind$K_frac_mean # calculate PLC from the K remaining, for plotting
wp.ind$PLC_mean_ll <- 100 - 100*wp.ind$K_frac_mean_ll # calculate PLC from the K remaining, for plotting

wp.ind$Site <- factor(wp.ind$Site)


### average to Sites
wp.site <- wp.ind %>% group_by(Site, SiteName,Lat.dd, Lon.dd, ParentMat, ParentMatcomb,plant_avail, aet, cwd, pet, ppt, aet.2018wy, pet.2018wy, pet.2018gs, cwd.2018wy, ppt.2018wy, ppt.2018sp, ppt.2018sm, tmx.2018wy, tmx.2018gs, tmn.2018wy, str.2018wy, soil_depth,
                                         tc.aet,tc.cwd, tc.cwd_gs,tc.pet, tc.ppt, tc.soil_mean, tc.soil_min, tc.srad, tc.tmax, tc.tmin, tc.vpd_max, tc.vpd_spr, tc.vpd_gs, 
                                         tc.aet.2018ds, tc.aet.2018gs, tc.aet.2018wy, tc.cwd.2018ds, tc.cwd.2018gs, tc.cwd.2018wy, tc.ppt.2018ds, tc.ppt.2018gs, tc.ppt.2018wy, tc.pet.2018ds, tc.pet.2018gs, tc.pet.2018wy,
                                         tc.soil_mean.2018wy, tc.soil_min.2018wy, tc.soil_mean.2018ds, tc.soil_min.2018ds, tc.soil_mean.2018gs, tc.soil_min.2018gs, tc.tmax.2018wy, tc.tmax.2018sp, tc.tmin.2018wy, tc.tmin.2018sp, tc.vpd_max.2018wy, tc.vpd.2018sp, tc.vpd.2018gs,
                                         water.storage100cm, PAW, bedrock.depth) %>% 
  summarise(mPD = mean(PD, na.rm=T), sdPD = sd(PD, na.rm=T), sePD = se(PD), rangePD = max(PD, na.rm=T) - min(PD, na.rm=T),
            mMD=mean(MD,na.rm=T), sdMD = sd(MD, na.rm=T), seMD = se(MD), rangeMD = max(MD, na.rm=T) - min(MD, na.rm=T),
            mHSM_leaf = mean(HSM_leaf, na.rm=T), mHSM_stem = mean(HSM_stem, na.rm=T),
            mE.drop = mean(E.drop, na.rm=T), sdE.drop = sd(E.drop, na.rm=T), seE.drop = se(E.drop),
            mdelD=mean(delD,na.rm=T), sddelD = sd(delD, na.rm=T), sedelD = se(delD),
            mdel18O=mean(del18O,na.rm=T), sddel18O = sd(del18O, na.rm=T), sedel18O = se(del18O),
            mlc_excess=mean(lc_excess,na.rm=T), sdlc_excess = sd(lc_excess, na.rm=T), selc_excess = se(lc_excess),
            mwinter_excessD=mean(winter_excessD,na.rm=T), sdwinter_excessD = sd(winter_excessD, na.rm=T), sewinter_excessD = se(winter_excessD),
            mwinter_excess18O=mean(winter_excess18O,na.rm=T), sdwinter_excess18O = sd(winter_excess18O, na.rm=T), sewinterexcess18O = se(winter_excess18O),
            mmaxAl_As=mean(maxAl_As, na.rm=T),sdmaxAl_As = sd(maxAl_As, na.rm=T), semaxAl_As = se(maxAl_As), maxAl_Al=max(maxAl_As, na.rm=T),
            sdmAl_As = sd(mAl_As, na.rm=T), semAl_As = se(mAl_As), mAl_As=mean(mAl_As, na.rm=T),
            sdmml_ms = sd(mml_ms, na.rm=T), semml_ms = se(mml_ms), mml_ms=mean(mml_ms, na.rm=T),
            sdmleafsize = sd(mleafsize, na.rm=T), semleafsize = se(mleafsize), mleafsize=mean(mleafsize, na.rm=T),
            sdmLMA = sd(mLMA, na.rm=T), semLMA = se(mLMA), mLMA=mean(mLMA, na.rm=T),
            sdmLDMC = sd(mLDMC, na.rm=T), semLDMC = se(mLDMC), mLDMC=mean(mLDMC, na.rm=T),
            sdmLength = sd(mLength, na.rm=T), semLength = se(mLength), mLength=mean(mLength, na.rm=T),
            mgrowth5yr = mean(growth5yr, na.rm=T), mBAI5yr = mean(BAI, na.rm=T), mperc_maxBAI=mean(perc_maxBAI, na.rm=T), seperc_maxBAI=se(perc_maxBAI), mHeight=mean(Height, na.rm=T), mDBH=mean(DBH, na.rm=T),
            NPP_mean = mean(NPP_mean, na.rm=T),GPP_mean = mean(GPP_mean, na.rm=T), K_frac_mean = mean(K_frac_mean, na.rm=T),psi_leaf_mean = mean(psi_leaf_mean, na.rm=T), psi_leaf_min = mean(psi_leaf_min, na.rm=T), PLC_mean=mean(PLC_mean, na.rm=T),
            NPP_mean_ll = mean(NPP_mean_ll, na.rm=T),GPP_mean_ll = mean(GPP_mean_ll, na.rm=T), K_frac_mean_ll = mean(K_frac_mean_ll, na.rm=T),psi_leaf_mean_ll = mean(psi_leaf_mean_ll, na.rm=T), psi_leaf_min = mean(psi_leaf_min, na.rm=T), PLC_mean_ll=mean(PLC_mean_ll, na.rm=T))



wp.site$maxAl_Al[which(wp.site$maxAl_Al== -Inf)] <- NA
wp.site$logLength <- log(wp.site$mLength, base=10)

# write out wp.ind dataset
#write.csv(wp.ind, "DerivedData/Wild_individualaverage_alltraits_20240515.csv")
### Also load in the branch level data for many of these traits for some analyses
branch.traits <- read.csv("DerivedData/FallWild_Leaftraits_branch_240508.csv")

#______________________________________________________________________
#######   * Gas Exchange Data ###################################################
#______________________________________________________________________
Midday.gasex <- read.csv("DerivedData/SummaryGasEx_fall2018_20190812.csv")


# load in the morning time series from Santa Maria
smr.ind.clean <- read.csv("DerivedData/SMR_GasEx_MorningTimeSeries.csv")





#_____________________________________________________________
##### ** FIG 1: Sites in climate space #######
#_____________________________________________________________

# let's make a specific color for blue oak
blueoak <- brewer.pal("Set1",n=3)[2]
blueoak.transp <- paste0(blueoak,"44") # and a transparent version


quartz(width=6.5, height=3.5)
par( mfrow=c(1,2))

# precip colors
#nuuk <- color("nuuk")
davos <- color("davos")
#batlow <- color("batlowW")
samplocscol <- brewer.pal(3,"Dark2")[2]
plot(tcn_an$ppt, col=rev(davos(14)[5:14]))
lines(calbound, col="black")
points(latitude~longitude,qudo, pch=16, cex=.3, col="black")
points(Lat.dd~Lon.dd, pop.terraclim, pch=16, col=samplocscol)
mtext(side=1, adj=1.5,text= "MAP(mm)", cex=.9, line=-.5)
par(mar=c(3.1,3,3,1), mgp=c(2.2,1,0))
plot(ppt~pet
     , data=qudo # data where to find those variables
     , ylim=c(130,1600)
     , xlim=c(1000,2000)
     , ylab="MAP (mm)" # x axis label
     , xlab="PET (mm)" # y axis label
     , pch=16 # point type (16=filled circle)
     , col="#22222222")#blueoak.transp) # point color, in a hexidecimal notation, which is #RRGGBB and then 2 digits for transparency so we don't overplot
points(ppt.2018wy~pet.2018wy
       , data=pop.terraclim
       , pch=3
       , cex=1.5
       , col=samplocscol)
points(ppt~pet
       , data=pop.terraclim
       , pch=16
       , cex=1.5
       , col=samplocscol)
legend(x=-200,2000, xpd=NA,legend = c("herbarium records","Sites (historical)","Sites (2018 wy)"), pch=c(16,16,3)
       , col=c("black",samplocscol,samplocscol), ncol=3, bty="n", cex=.9)


if(save.figures==T){quartz.save(file=paste0(results.dir,"/Fig_Map_Clim_v1.pdf"),type = "pdf")}





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########### Analysis: Test environmental predictors of traits/wp/growth ##########



## function for testing characteristic-climate relationships
  # only inluding 1 predictor, and selecting baised on AICc which variable is best
best.mermod <- function(trait, dataz, clim.vars) {
  results <- data.frame("mod"=rep(NA, times=length(clim.vars)+1)
                        , "n"=rep(NA, times=length(clim.vars)+1)
                        , "AICc"=rep(NA, times=length(clim.vars)+1)
                        , "AIC"=rep(NA, times=length(clim.vars)+1)
                        ,"BIC"=rep(NA, times=length(clim.vars)+1)
                        ,"df"=rep(NA, times=length(clim.vars)+1)
                        )
  null <- lmer(get(trait)~ 1 + (1|Site), dataz, REML=F)
  results$mod[1] <- "null"
  results$n[1] <- nrow(null@frame)
  results$AICc[1] <- AICc(null)
  results$AIC[1] <- summary(null)$AICtab[1]
  results$BIC[1] <- summary(null)$AICtab[2]
  results$df[1] <- summary(null)$AICtab[5]
  for (i in 1:length(clim.vars)){
    mod <- lmer(get(trait)~ get(clim.vars[i]) + (1|Site), dataz, REML=F)
    results$mod[i+1] <- clim.vars[i]
    results$n[i+1] <- nrow(mod@frame)
    results$AICc[i+1] <- AICc(mod)
    results$AIC[i+1] <- summary(mod)$AICtab[1]
    results$BIC[i+1] <- summary(mod)$AICtab[2]
    results$df[i+1] <- summary(mod)$AICtab[5]
  }
  results <- results %>% arrange(AICc)
  print(paste("deltaAIC:", round(results$AICc[1] - results$AICc[which(results$mod=="null")], 2)))
  return(results)
}

# function to calculate eta2 (among population variance) from an ANOVA
eta2 <- function(aovres, group = "Site"){
  eta2 <- aovres[[1]]$'Sum Sq'[grep(group, rownames(aovres[[1]]))]/sum(aovres[[1]]$`Sum Sq`)
  return(eta2)
}

# function to calculate omega2 (unbiased among pop variance) from an ANOVA
omega2 <- function(aovres, group = "Site"){
  omega2 <- (aovres[[1]]$'Sum Sq'[grep(group, rownames(aovres[[1]]))] - (aovres[[1]]$Df[grep(group, rownames(aovres[[1]]))]*aovres[[1]]$`Mean Sq`[grep("Residuals", rownames(aovres[[1]]))]))/(sum(aovres[[1]]$`Sum Sq`)+aovres[[1]]$`Mean Sq`[grep("Residuals", rownames(aovres[[1]]))])
  return(omega2)
}

# extract significance of pop factor in anova
aovsig <- function(aovres, group = "Site"){
  sig <- aovres[[1]]$'Pr(>F)'[grep(group, rownames(aovres[[1]]))]
  return(round(sig,3))
}


# Creat a Table of Summary Results:

## Response variables to test ('characteristics' since some are 'traits' and some are not)
chars <- c("PD","MD","E.drop","mAl_As","maxAl_As","mml_ms","mleafsize","mLMA","mLDMC","lc_excess","delD","del18O", "logLength","Height","perc_maxBAI")

## dataframe to hold results in
ClimRes <- data.frame(Trait=chars, n = rep(NA, length(chars)), eta2 = rep(NA, length(chars)), omega2 = rep(NA, length(chars)),sig = rep(NA, length(chars)), CV=rep(NA, length(chars)), bestclim = rep(NA, length(chars)), climsig = rep(NA, length(chars)), climvar=rep(NA, length(chars)))


# Calculate how much of the total variance is among sites, 
# just to know what sort of among-site signal we're dealing with

for(j in chars){
  ClimRes$n[which(ClimRes$Trait==j)] <- nrow(wp.ind) - length(which(is.na(wp.ind[,j])))
  ClimRes$eta2[which(ClimRes$Trait==j)] <- round(eta2(summary(aov(get(j)~Site , wp.ind))),2)
  ClimRes$omega2[which(ClimRes$Trait==j)] <- round(omega2(summary(aov(get(j)~Site , wp.ind))),2)
  ClimRes$CV[which(ClimRes$Trait==j)] <- round(abs(sd(wp.ind[,j], na.rm=T)/mean(wp.ind[,j], na.rm=T)),2)
  ClimRes$sig[which(ClimRes$Trait==j)] <- aovsig(summary(aov(get(j)~Site , wp.ind)))
  
  
}

### Climate Variables to Test
tc.vars.to.test <- c("cwd","aet","pet","ppt","tmin","soil_mean","vpd_max","vpd_spr",
                  "aet.2018wy","cwd.2018wy","ppt.2018wy","pet.2018wy", "soil_mean.2018wy","tmin.2018wy", "tmax.2018wy","ppt.2018gs",
                  "pet.2018ds", "soil_min.2018gs", "vpd.2018gs","tmax.2018sp", "tmin.2018sp")
vars.to.test <- c(paste("tc", tc.vars.to.test, sep="."), "water.storage100cm","PAW","bedrock.depth")

## PD
(PDtest <- best.mermod("PD", wp.ind, vars.to.test))
#OLD: 2018 gs PPT, VPD and wy AET/VPDmax plus historic AET all pretty similar but only marginal
# aet, 2018gs VPD, aet2018WY, 
PDmod <- lmer(PD~tc.aet + (1|Site), wp.ind, REML=T)
qqp(resid(PDmod)) # resids normal
qqp(ranef(PDmod)$Site[[1]]) # random effects normal


ClimRes$bestclim[which(ClimRes$Trait=="PD")] <- "30yr AET" #"GS PPT 2018"
ClimRes$climsig[which(ClimRes$Trait=="PD")] <- round(summary(PDmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="PD")] <- round(r.squaredGLMM(PDmod)[[1]],3) # R2marg = 0.159 

## MD
(MDtest <- best.mermod("MD", wp.ind, vars.to.test))
# PPT 2018 wy functionally same is null. now Null is best
MDmod <- lmer(MD~tc.ppt.2018wy + (1|Site), wp.ind, REML=T)
qqp(resid(MDmod)) # resids normal
qqp(ranef(MDmod)$Site[[1]]) # random effects normal
summary(MDmod) # p = 0.145
ClimRes$bestclim[which(ClimRes$Trait=="MD")] <- "WY PPT 2018"
ClimRes$climsig[which(ClimRes$Trait=="MD")] <- round(summary(MDmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="MD")] <- round(r.squaredGLMM(MDmod)[[1]],3) # R2marg = 0.159 

## E.drop
(E.droptest <- best.mermod("E.drop", wp.ind, vars.to.test))
  # VPD, Tmax and pet consistently best (2018 and normal equally good)
E.dropmod <- lmer(E.drop~tc.vpd.2018gs + (1|Site), wp.ind, REML=T)
qqp(resid(E.dropmod)) # resids normal
qqp(ranef(E.dropmod)$Site[[1]]) # random effects normal

ClimRes$bestclim[which(ClimRes$Trait=="E.drop")] <- "GS VPD 2018"
ClimRes$climsig[which(ClimRes$Trait=="E.drop")] <- round(summary(E.dropmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="E.drop")] <- round(r.squaredGLMM(E.dropmod)[[1]],3) 


## Al_As
(mAl_Astest <- best.mermod("mAl_As", wp.ind, vars.to.test))
(logAl_Astest <- best.mermod("logAl_As", wp.ind, vars.to.test))
# tmin best for raw and logged Al_As, and they look very similar
mAl_Asmod <- lmer(mAl_As~water.storage100cm + (1|Site), wp.ind, REML=T)
qqp(resid(mAl_Asmod)) # resids normal
qqp(ranef(mAl_Asmod)$Site[[1]]) # random effects normal
  # going with raw because resids not that bad
# logAl_Asmod <- lmer(logAl_As~tc.tmin.2018wy + (1|Site), wp.ind, REML=T)
# qqp(resid(logAl_Asmod)) # resids normal
# qqp(ranef(logAl_Asmod)$Site[[1]]) # random effects normal
ClimRes$bestclim[which(ClimRes$Trait=="mAl_As")] <- "1m Soil Water Storage"
ClimRes$climsig[which(ClimRes$Trait=="mAl_As")] <- round(summary(mAl_Asmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="mAl_As")] <- round(r.squaredGLMM(mAl_Asmod)[[1]],3) # R2marg = 0.159 

## max Al_As
(maxAl_Astest <- best.mermod("maxAl_As", wp.ind, vars.to.test))
# tmin best for maxAl_As, 
maxAl_Asmod <- lmer(maxAl_As~water.storage100cm + (1|Site), wp.ind, REML=T)
qqp(resid(maxAl_Asmod)) # resids normal
qqp(ranef(maxAl_Asmod)$Site[[1]]) # random effects normal
ClimRes$bestclim[which(ClimRes$Trait=="maxAl_As")] <- "water.storage100cm"
ClimRes$climsig[which(ClimRes$Trait=="maxAl_As")] <- round(summary(maxAl_Asmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="maxAl_As")] <- round(r.squaredGLMM(maxAl_Asmod)[[1]],3) # R2marg = 0.159 

## mml_ms
(mml_mstest <- best.mermod("mml_ms", wp.ind, vars.to.test))
(logml_mstest <- best.mermod("logml_ms", wp.ind, vars.to.test))
# tmin best for mml_ms and logml_ms
mml_msmod <- lmer(logml_ms~tc.tmin.2018wy + (1|Site), wp.ind, REML=T)
  # need to log-transform
qqp(resid(mml_msmod)) # some bad high resids unless logged
qqp(ranef(mml_msmod)$Site[[1]]) # random effects normal
ClimRes$bestclim[which(ClimRes$Trait=="mml_ms")] <- "Tmin 2018"
ClimRes$climsig[which(ClimRes$Trait=="mml_ms")] <- round(summary(mml_msmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="mml_ms")] <- round(r.squaredGLMM(mml_msmod)[[1]],3) # R2marg = 0.159 


## leaf size
(mleafsizetest <- best.mermod("mleafsize", wp.ind, vars.to.test))
  # normal ppt/cwd and 2018 ppt/cwd + 2018 pet/aet all pretty similar
mleafsizemod <- lmer(mleafsize~tc.ppt + (1|Site), wp.ind, REML=T)
qqp(resid(mleafsizemod)) # resids normal-ish
qqp(ranef(mleafsizemod)$Site[[1]]) # random effects normal

ClimRes$bestclim[which(ClimRes$Trait=="mleafsize")] <- "30yr PPT"
ClimRes$climsig[which(ClimRes$Trait=="mleafsize")] <- round(summary(mleafsizemod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="mleafsize")] <- round(r.squaredGLMM(mleafsizemod)[[1]],3) 


## LMA
(mLMAtest <- best.mermod("mLMA", wp.ind, vars.to.test))
# PAW best by a ways
mLMAmod <- lmer(mLMA~PAW + (1|Site), wp.ind, REML=T)
qqp(resid(mLMAmod)) # resids normal-ish, except value 62. Nothing about SMR-1 in the raw data that suggests a methodological outlier
qqp(ranef(mLMAmod)$Site[[1]]) # random effects normal
# best.mermod("mLMA", wp.ind[-62,], vars.to.test)
  # if you remove the large outlier 62, then the null model is best. but nothing is significant regardless
ClimRes$bestclim[which(ClimRes$Trait=="mLMA")] <- "PAW"
ClimRes$climsig[which(ClimRes$Trait=="mLMA")] <- round(summary(mLMAmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="mLMA")] <- round(r.squaredGLMM(mLMAmod)[[1]],3) 


## LDMC
(mLDMCtest <- best.mermod("mLDMC", wp.ind, vars.to.test))
# spring Tmax 2018 best
mLDMCmod <- lmer(mLDMC~tc.tmax.2018sp + (1|Site), wp.ind, REML=T)
qqp(resid(mLDMCmod)) # resids normal-ish
qqp(ranef(mLDMCmod)$Site[[1]]) # random effects normal

ClimRes$bestclim[which(ClimRes$Trait=="mLDMC")] <- "SP Tmax 2018"
ClimRes$climsig[which(ClimRes$Trait=="mLDMC")] <- round(summary(mLDMCmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="mLDMC")] <- round(r.squaredGLMM(mLDMCmod)[[1]],3) 

## lc_excess
(lc_excesstest <- best.mermod("lc_excess", wp.ind, vars.to.test))
# SP tmax 2018 best, but REALLY not sig
lc_excessmod <- lmer(lc_excess~tc.tmax.2018sp + (1|Site), wp.ind, REML=T)
qqp(resid(lc_excessmod)) # some very low small resids, but no transform really helps
qqp(ranef(lc_excessmod)$Site[[1]]) # random effects normal

ClimRes$bestclim[which(ClimRes$Trait=="lc_excess")] <- "SP Tmax 2018"
ClimRes$climsig[which(ClimRes$Trait=="lc_excess")] <- round(summary(lc_excessmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="lc_excess")] <- round(r.squaredGLMM(lc_excessmod)[[1]],3) 

## delD
(delDtest <- best.mermod("delD", wp.ind, vars.to.test))
# 2018 VPD and VPD norm both good preditors
delDmod <- lmer(delD~tc.vpd.2018gs + (1|Site), wp.ind, REML=T)
qqp(resid(delDmod)) # some very low small resids, but no transform really helps
qqp(ranef(delDmod)$Site[[1]]) # random effects normal

ClimRes$bestclim[which(ClimRes$Trait=="delD")] <- "GS VPD 2018"
ClimRes$climsig[which(ClimRes$Trait=="delD")] <- round(summary(delDmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="delD")] <- round(r.squaredGLMM(delDmod)[[1]],3) 

## del18O
(del18Otest <- best.mermod("del18O", wp.ind, vars.to.test))
# similar to delD, VPD2018gs and vpd_max good predictors
del18Omod <- lmer(del18O~tc.vpd.2018gs + (1|Site), wp.ind, REML=T)
qqp(resid(del18Omod)) # some very low small resids, but no transform really helps
qqp(ranef(del18Omod)$Site[[1]]) # random effects normal

ClimRes$bestclim[which(ClimRes$Trait=="del18O")] <- "GS VPD 2018"
ClimRes$climsig[which(ClimRes$Trait=="del18O")] <- round(summary(del18Omod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="del18O")] <- round(r.squaredGLMM(del18Omod)[[1]],3) 



## mLength
(mLengthtest <- best.mermod("logLength", wp.ind, vars.to.test))
# tmin 2018 best
mLengthmod <- lmer(logLength~tc.tmin.2018wy + (1|Site), wp.ind, REML=T)
qqp(resid(mLengthmod)) # much betterrt logged
qqp(ranef(mLengthmod)$Site[[1]]) # random effects normal

ClimRes$bestclim[which(ClimRes$Trait=="logLength")] <- "Tmin 2018"
ClimRes$climsig[which(ClimRes$Trait=="logLength")] <- round(summary(mLengthmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="logLength")] <- round(r.squaredGLMM(mLengthmod)[[1]],3) 


## Height
(Heighttest <- best.mermod("Height", wp.ind, vars.to.test))
# normal pet
Heightmod <- lmer(Height~tc.pet + (1|Site), wp.ind, REML=T)
qqp(resid(Heightmod)) # some very low small resids, but no transform really helps
qqp(ranef(Heightmod)$Site[[1]]) # random effects normal

ClimRes$bestclim[which(ClimRes$Trait=="Height")] <- "30yr PET"
ClimRes$climsig[which(ClimRes$Trait=="Height")] <- round(summary(Heightmod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="Height")] <- round(r.squaredGLMM(Heightmod)[[1]],3) 


## perc_maxBAI
(perc_maxBAItest <- best.mermod("perc_maxBAI", wp.ind, vars.to.test))
# tmin best - old:2018 Tmin and SP tmax 2018 best similar, and remarkably uncorrelated...
perc_maxBAImod <- lmer(perc_maxBAI~tc.tmin + (1|Site), wp.ind, REML=T)
qqp(resid(perc_maxBAImod)) 
qqp(ranef(perc_maxBAImod)$Site[[1]]) # random effects normal

ClimRes$bestclim[which(ClimRes$Trait=="perc_maxBAI")] <- "30yr Tmin"
ClimRes$climsig[which(ClimRes$Trait=="perc_maxBAI")] <- round(summary(perc_maxBAImod)$coefficients[2,5],3)
ClimRes$climvar[which(ClimRes$Trait=="perc_maxBAI")] <- round(r.squaredGLMM(perc_maxBAImod)[[1]],3) 


## some columns for plotting
ClimRes$Type <- factor(c("WS","WS","WS","T","T","T","T","T","T","WS","WS","WS","Gr","Gr","Gr"))
ClimRes$climsig.binary <- as.numeric(as.character(cut(ClimRes$climsig,breaks = c(0,0.05,1), labels=c(16,1))))
ClimRes$climcat <- factor(c("Supply","Supply","Demand","Soil","Soil","Temp","Supply","Soil","Temp","Temp","Demand","Demand","Temp","Demand","Temp"))
ClimRes$climcat.simp <- factor(c("Water","Water","Water","Water","Water","Temp","Water","Water","Temp","Temp","Water","Water","Temp","Water","Temp"))

ClimRes.clean <- ClimRes[which(!ClimRes$Trait %in% c("maxAl_As","lc_excess","delD","del18O","Height")), ]
if(save.figures==T){write.csv(ClimRes.clean,paste0(results.dir,"/TableS1_EnvironmentalPredictors_v1.csv"))}
#______________________________________________________________________________








#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#### . FIG 2: Predawn and Midday aren't predicted by climate ############

# ___V3: Best climate predictor, even if non-significant, plus HSM
quartz(width=6.8*2/3, height=6)
par(mfrow=c(2,2), mar=c(3,3,1.5,1),mgp=c(2,1,0), cex=1)
# plotting tweaks
site.cex <- 1.3
tree.col <- "#55555555"


plot(PD~tc.aet, wp.ind, pch=16, col=tree.col
     , ylab=expression(paste(Psi[PD], " (MPa)"))
     , xlab="30yr mean AET (mm)", ylim = c(-5,-0.5))
points(mPD~tc.aet, wp.site, pch=16, cex=site.cex)
modfit <- summary(lmer(PD~tc.aet + (1|Site), wp.ind))
#abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=2)
mtext(paste("p=", round(modfit$coefficients[2,5],2)),side = 3,line = -1,adj = 0.05)
mtext("a)", side=3, line=0.1, adj=-0.1)
abline(h=-3.88, col="grey", lty=2)
abline(h=-4.47, col="black", lty=2)


plot(MD~tc.ppt.2018wy, wp.ind, pch=16, col=tree.col
     , ylab=expression(paste(Psi[MD]," (MPa)"))
     , xlab="2018 wy PPT (mm)",  ylim = c(-5,-0.5))
points(mMD~tc.ppt.2018wy, wp.site, pch=16, cex=site.cex)
modfit <- summary(lmer(MD~tc.ppt.2018wy + (1|Site), wp.ind))
#abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=2)
mtext(paste("p=", round(modfit$coefficients[2,5],2)),side = 3,line = -1,adj = 0.05)
mtext("b)", side=3, line=0.1, adj=-0.1)
abline(h=-3.88, col="grey", lty=2)
abline(h=-4.47, col="black", lty=2)


plot(E.drop~tc.vpd.2018gs, wp.ind, pch=16, col=tree.col
     , ylab=expression(paste(Delta*Psi," (",Psi[PD]-Psi[MD],", MPa)"))
     , xlab="2018 Apr-Sep VPD (kPa)")
points(mE.drop~tc.vpd.2018gs, wp.site, pch=16, cex=site.cex)
modfit <- summary(lmer(E.drop~tc.vpd.2018gs + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
mtext(paste("p=", round(modfit$coefficients[2,5],2)),side = 3,line = -1,adj = 0.05)
mtext("c)", side=3, line=0.1, adj=-0.1)


samplocscol <- brewer.pal(3,"Dark2")[2]
plot(HSM_stem~tc.cwd.2018wy, wp.ind, pch=16, col=paste0(samplocscol, "55"), ylim=c(-2,3.5)
     , xlab="2018 wy CWD (mm)", ylab=expression(paste("HSM (",Psi[MD] - P50,", MPa)" )))
points(HSM_leaf~tc.cwd.2018wy, wp.ind, pch=16, col=tree.col)
points(mHSM_stem~tc.cwd.2018wy, wp.site, cex=site.cex, pch=16, col=samplocscol)
points(mHSM_leaf~tc.cwd.2018wy, wp.site, cex=site.cex, pch=16)
abline(h=0)
legend("bottom", bty="n", legend=c("leaf", "stem"), col=c("black",samplocscol), pch=16, ncol=2)
mtext("d)", side=3, line=0.1, adj=-0.1)

if(save.figures==T){quartz.save(file=paste0(results.dir,"/Fig2_PD_MD_climate_extended_v5.pdf"),type = "pdf")}






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###### Analysis: WP~Iso relationships ##################

# plan is to make a table with
# rows: PD, MD, deltPsi
# columns: p value, R2 marginal, R2 conditional (maybe among site R2?) for delD and del18, with and without PWD

isoresults <- data.frame("Potential" = c("PD","MD","E.drop"), "delD.p"=rep(NA, 3), "delD.R2m"=rep(NA, 3), "delD.R2c"=rep(NA, 3),"delD.R2site"=rep(NA,3), "del18O.p"=rep(NA, 3), "del18O.R2m"=rep(NA, 3), "del18O.R2c"=rep(NA, 3),"del18O.R2site"=rep(NA,3))
isoresults.np <- data.frame("Potential" = c("PD","MD","E.drop"), "delD.p"=rep(NA, 3), "delD.R2m"=rep(NA, 3), "delD.R2c"=rep(NA, 3),"delD.R2site"=rep(NA,3), "del18O.p"=rep(NA, 3), "del18O.R2m"=rep(NA, 3), "del18O.R2c"=rep(NA, 3),"del18O.R2site"=rep(NA,3))


for(i in 1:3){
  variable <- isoresults$Potential[i]
  mod <- lmer(get(variable)~delD + (1|Site), wp.ind, REML=T)
  isoresults$delD.p[i] <- round(summary(mod)$coefficients[2,5],3)
  isoresults$delD.R2m[i] <-  round(r.squaredGLMM(mod)[1],3)
  isoresults$delD.R2c[i] <-  round(r.squaredGLMM(mod)[2],3)
  isoresults$delD.R2site[i] <-  round(r.squaredGLMM(lm(get(paste0("m",variable))~mdelD, wp.site))[1],3)
  mod <- lmer(get(variable)~delD + (1|Site), wp.ind[which(wp.ind$Site != "PWD"),], REML=T)
  isoresults.np$delD.p[i] <-  round(summary(mod)$coefficients[2,5],3)
  isoresults.np$delD.R2m[i] <-  round(r.squaredGLMM(mod)[1],3)
  isoresults.np$delD.R2c[i] <-  round(r.squaredGLMM(mod)[2],3)
  isoresults.np$delD.R2site[i] <-  round(r.squaredGLMM(lm(get(paste0("m",variable))~mdelD, wp.site[which(wp.site$Site != "PWD"),]))[1],3)
}

for(i in 1:3){
  variable <- isoresults$Potential[i]
  mod <- lmer(get(variable)~del18O + (1|Site), wp.ind, REML=T)
  isoresults$del18O.p[i] <- round(summary(mod)$coefficients[2,5],3)
  isoresults$del18O.R2m[i] <-  round(r.squaredGLMM(mod)[1],3)
  isoresults$del18O.R2c[i] <-  round(r.squaredGLMM(mod)[2],3)
  isoresults$del18O.R2site[i] <-  round(r.squaredGLMM(lm(get(paste0("m",variable))~mdel18O, wp.site))[1],3)
  mod <- lmer(get(variable)~del18O + (1|Site), wp.ind[which(wp.ind$Site != "PWD"),], REML=T)
  isoresults.np$del18O.p[i] <-  round(summary(mod)$coefficients[2,5],3)
  isoresults.np$del18O.R2m[i] <-  round(r.squaredGLMM(mod)[1],3)
  isoresults.np$del18O.R2c[i] <-  round(r.squaredGLMM(mod)[2],3)
  isoresults.np$del18O.R2site[i] <-  round(r.squaredGLMM(lm(get(paste0("m",variable))~mdel18O, wp.site[which(wp.site$Site != "PWD"),]))[1],3)
  
}
isoresults.all <- rbind(isoresults, isoresults.np)
isoresults.all$type <- rep(c("all", "no PWD"), each=3)
if(save.figures==T){write.csv(isoresults.all, paste0(results.dir,"/Iso_WaterPotential_results_v1.csv"))}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### . FIG 3: Predawn f(delD) ###########


quartz(width=3.2, height=6)
par(mfrow=c(3,1), mar=c(0,3.3,0,1), oma=c(3,0,1,0), mgp=c(2,1,0), cex=1)

# Predawn
palette(paste0(mypal,"77"))
plot(PD~delD, wp.ind, col=factor(Site), pch=16, ylab=expression(paste(Psi[PD]," (MPa)" ))
     , xlab=expression(paste(delta,"D (\u2030)")), xaxt="n")
modfit <- summary(lmer(PD~delD + (1|Site), wp.ind[which(wp.ind$Site != "PWD"),]))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(DBH)
  fitted <- predict(lm(PD~delD, tmp),newdata=data.frame("delD"=c(min(tmp$delD, na.rm=T), max(tmp$delD, na.rm=T))))
  lines(x=c(min(tmp$delD, na.rm=T), max(tmp$delD, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mPD~mdelD, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("a)", side=3, line=-1, adj=0.05)

# Midday
palette(paste0(mypal,"77"))
plot(MD~delD, wp.ind, col=factor(Site), pch=16, ylab=expression(paste(Psi[MD]," (MPa)" )), xlab=expression(paste(delta,"D (\u2030)")), xaxt="n")
modfit <- summary(lmer(MD~delD + (1|Site), wp.ind[which(wp.ind$Site != "PWD"),]))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(DBH)
  fitted <- predict(lm(MD~delD, tmp),newdata=data.frame("delD"=c(min(tmp$delD, na.rm=T), max(tmp$delD, na.rm=T))))
  lines(x=c(min(tmp$delD, na.rm=T), max(tmp$delD, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mMD~mdelD, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("b)", side=3, line=-1, adj=0.05)

# Delta Psi
palette(paste0(mypal,"77"))
plot(E.drop~delD, wp.ind, col=factor(Site), pch=16, ylab=expression(paste(Delta*Psi," (",Psi[PD]-Psi[MD]," (MPa)")), xlab=expression(paste(delta,"D (\u2030)")))
modfit <- summary(lmer(E.drop~delD + (1|Site), wp.ind[which(wp.ind$Site != "PWD"),]))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(DBH)
  fitted <- predict(lm(E.drop~delD, tmp),newdata=data.frame("delD"=c(min(tmp$delD, na.rm=T), max(tmp$delD, na.rm=T))))
  lines(x=c(min(tmp$delD, na.rm=T), max(tmp$delD, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mE.drop~mdelD, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext(text=expression(paste(delta,"D (\u2030)")),side = 1, line=2)
mtext("c)", side=3, line=-1, adj=0.05)

if(save.figures==T){quartz.save(file=paste0(results.dir,"/Fig3_WP_delD_v1.pdf"),type = "pdf")}
#+++++++++++++++++++++++++++++++++++++++++++++++++




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############## . FIG 4: Gs closure with WP ################

quartz(width=6.8, height=3.2)
par(mar=c(4,4,1,1), mgp=c(2.5,1,0), mfrow=c(1,2))
plot(gs_mmol~Psi, Midday.gasex, pch=16, col=factor(Site), ylab=expression(paste(g[s]," ", (mmol*m^-2*s^-1))), xlab=expression(paste(Psi, " (MPa)")))
points(gs*1000~Psi, smr.ind.clean, pch=4, col="darkgrey")
abline(v=-3.88, col="grey", lty=2)
abline(v=-4.47, col="black", lty=2)
legend('topleft', legend=c("Midday across CA","timecourse, 1 tree"), pch=c(16,4), col=c("black","darkgrey"), cex=.8, bg = "white")


#quartz(width=3.5, height=3.5)
#par(mar=c(4,4,1,1), mgp=c(2.5,1,0))
plot(As~Psi, Midday.gasex, pch=16, col=factor(Site), ylab=expression(paste(A[s]," ", (mu*mol*m^-2*s^-1))), xlab=expression(paste(Psi, " (MPa)")))
points(A~Psi, smr.ind.clean, pch=4, col="darkgrey")
abline(v=-3.88, col="grey", lty=2)
abline(v=-4.47, col="black", lty=2)
#legend('topleft', legend=c("Midday across CA","morning closure, one tree"), pch=c(16,4), col=c("black","darkgrey"), cex=.8)
if(save.figures==T){quartz.save(file=paste0(results.dir,"/Fig4_GasEx_WP_v1.pdf"),type = "pdf")}






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########. FIG 5: Trait-climate ####################################

quartz(width=5.8, height=4.5)
par(mfrow=c(2,3), mar=c(3,3,1,1), mgp=c(2,1,0), cex=1)

# updated with soil info
plot(mAl_As~water.storage100cm, wp.ind, pch=16, col=tree.col
     , ylab=expression(paste(A[leaf]:A[stem]))
     , xlab="1m soil storage (cm)")
points(mAl_As~water.storage100cm, wp.site, pch=16, cex=site.cex)
modfit <- summary(lmer(mAl_As~water.storage100cm + (1|Site), wp.ind))
# old with only climate
# plot(mAl_As~tc.tmin.2018wy, wp.ind, pch=16, col=tree.col
#      , ylab=expression(paste(A[leaf]:A[stem]))
#      , xlab="Min Temp 2018wy (°C)")
# points(mAl_As~tc.tmin.2018wy, wp.site, pch=16, cex=site.cex)
# modfit <- summary(lmer(mAl_As~tc.tmin.2018wy + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
mtext("a)", side=3, line=0.1, adj=-0.1)

plot(mml_ms~tc.tmin.2018wy, wp.ind, pch=16, col=tree.col
     , ylab=expression(paste(Mass[leaf]:Mass[stem]))
     , xlab="Min Temp 2018wy (°C)")
points(mml_ms~tc.tmin.2018wy, wp.site, pch=16, cex=site.cex)
modfit <- summary(lmer(mml_ms~tc.tmin.2018wy + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
mtext("b)", side=3, line=0.1, adj=-0.1)

plot(mleafsize~tc.ppt, wp.ind, pch=16, col=tree.col
     , ylab="mean leaf size (cm2)"
     , xlab="30yr MAP (mm)")
points(mleafsize~tc.ppt, wp.site, pch=16, cex=site.cex)
modfit <- summary(lmer(mleafsize~tc.ppt + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
mtext("c)", side=3, line=0.1, adj=-0.1)

#new soils
plot(mLMA~PAW, wp.ind, pch=16, col=tree.col
     , ylab="LMA (g/cm2)"
     , xlab="Plant Avial Water (cm)")
points(mLMA~PAW, wp.site, pch=16, cex=site.cex)
# old, climate alone
# plot(mLMA~tc.aet, wp.ind, pch=16, col=tree.col
#      , ylab="LMA (g/cm2)"
# , xlab="30yr AET (mm)")
# points(mLMA~tc.aet, wp.site, pch=16, cex=site.cex)
modfit <- summary(lmer(mLMA~PAW + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
mtext("d)", side=3, line=0.1, adj=-0.1)

plot(mLDMC~tc.tmax.2018sp, wp.ind, pch=16, col=tree.col
     , ylab="LDMC (g/g)"
     , xlab="2018 spring Tmax (°C)")
points(mLDMC~tc.tmax.2018sp, wp.site, pch=16, cex=site.cex)
mtext("e)", side=3, line=0.1, adj=-0.1)

if(save.figures==T){quartz.save(file=paste0(results.dir,"/Fig5_Traits_climatesoil_v2.pdf"),type = "pdf")}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++







#++++++++++++++++++++++++++++++++++++++++++++++++++++
############# . FIG 6: Growth~climate responses ############

# Analysis: Growth rate correlation with each other
r.squaredGLMM(lmer(logLength~perc_maxBAI + (1|Site), wp.ind))
summary(lmer(logLength~perc_maxBAI + (1|Site), wp.ind))

# Analysis: Trait~WP relationships
# logLength was used in all the other analyses, stick to that here
summary(lmer(logLength~PD + (1|Site), wp.ind)) # marginally significant
summary(lmer(logLength~MD + (1|Site), wp.ind))
summary(lmer(logLength~E.drop + (1|Site), wp.ind)) # significant, barely

summary(lmer(perc_maxBAI~PD + (1|Site), wp.ind)) #ns
summary(lmer(perc_maxBAI~MD + (1|Site), wp.ind)) #ns
summary(lmer(perc_maxBAI~E.drop + (1|Site), wp.ind)) # ns


# perc_maxBAI
# plotting tweaks
site.cex <- 1.3
tree.col <- "#55555555"
quartz(width=6, height=5)
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2,1,0))
palette(paste0(mypal,"77"))
plot(perc_maxBAI~PD, wp.ind, col=factor(Site), pch=16, xlab=expression(paste(Psi[PD]," (MPa)" ))
     , ylab="% of max BAI")
#modfit <- summary(lmer(PD~delD + (1|Site), wp.ind[which(wp.ind$Site != "PWD"),]))
#abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(PD)
  fitted <- predict(lm(perc_maxBAI~PD, tmp),newdata=data.frame("PD"=c(min(tmp$PD, na.rm=T), max(tmp$PD, na.rm=T))))
  lines(x=c(min(tmp$PD, na.rm=T), max(tmp$PD, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mperc_maxBAI~mPD, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("a)", side=3, line=0, adj=0.05)

palette(paste0(mypal,"77"))
plot(perc_maxBAI~tc.tmin, wp.ind, pch=16, col=factor(Site),
     xlab=expression(paste("30yr ", T[min]," (°C)")),
     ylab= "% of max BAI")
modfit <- summary(lmer(perc_maxBAI~tc.tmin + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
points(mperc_maxBAI~tc.tmin, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("b)", side=3, line=0.1, adj=-0.05)


# palette(paste0(mypal,"77"))
# plot(perc_maxBAI~E.drop, wp.ind, col=factor(Site), pch=16, xlab=expression(paste(Delta*Psi," (MPa)" ))
#      , ylab="% of max BAI")
# #modfit <- summary(lmer(perc_maxBAI~E.drop + (1|Site), wp.ind[which(wp.ind$Site != "PWD"),]))
# #abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
# palette(mypal)
# for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
#   tmp <- wp.ind %>% filter(Site==i) %>% arrange(E.drop)
#   fitted <- predict(lm(perc_maxBAI~E.drop, tmp),newdata=data.frame("E.drop"=c(min(tmp$E.drop, na.rm=T), max(tmp$E.drop, na.rm=T))))
#   lines(x=c(min(tmp$E.drop, na.rm=T), max(tmp$E.drop, na.rm=T)), y=fitted, col="#55555555", lwd=2)
# }
# points(mperc_maxBAI~mE.drop, wp.site, col=factor(Site), pch=17, cex=site.cex)
# mtext("b)", side=3, line=-1, adj=0.05)

#### logLength
palette(paste0(mypal,"77"))
plot(mLength~E.drop, wp.ind, col=factor(Site), pch=16, xlab=expression(paste(Delta*Psi," (MPa)" ))
     , ylab="Stem Length (cm)", log="y")
modfit <- summary(lmer(logLength~E.drop + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=2)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(E.drop)
  fitted <- predict(lm(mLength~E.drop, tmp),newdata=data.frame("E.drop"=c(min(tmp$E.drop, na.rm=T), max(tmp$E.drop, na.rm=T))))
  lines(x=c(min(tmp$E.drop, na.rm=T), max(tmp$E.drop, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mLength~mE.drop, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("c)", side=3, line=0, adj=0.05)


palette(paste0(mypal,"77"))
plot(mLength~tc.tmin, wp.ind, pch=16, col=Site,
     xlab=expression(paste("30yr ", T[min]," (°C)")),
     ylab= "Stem Length (cm)", log="y")
modfit <- summary(lmer(logLength~tc.tmin + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2)
palette(mypal)
points(mLength~tc.tmin, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("d)", side=3, line=0.1, adj=-0.05)

if(save.figures==T){quartz.save(file=paste0(results.dir,"/Fig6_Growth-Climate_v2.pdf"),type = "pdf")}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++
############# . FIG 7: Trait~growth responses ############

quartz(width=6.8, height=4)
par(mfcol=c(2,5), mar=c(0,0,0,0), oma=c(3.2,3.5,1,1), mgp=c(2,1,0), cex=.95)

### Al_As
palette(paste0(mypal,"77"))
plot(perc_maxBAI~mAl_As, wp.ind, col=factor(Site), pch=16, xlab=expression(paste(A[l]:A[s]," (",cm^2*mm^-2,")" )),
     xaxt="n")
mtext("% of max BAI", side=2, line=2)
modfit <- summary(lmer(perc_maxBAI~mAl_As + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(mAl_As)
  fitted <- predict(lm(perc_maxBAI~mAl_As, tmp),newdata=data.frame("mAl_As"=c(min(tmp$mAl_As, na.rm=T), max(tmp$mAl_As, na.rm=T))))
  lines(x=c(min(tmp$mAl_As, na.rm=T), max(tmp$mAl_As, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mperc_maxBAI~mAl_As, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("a)", side=3, line=-1, adj=0.05)

palette(paste0(mypal,"77"))
plot(mLength~mAl_As, wp.ind, col=factor(Site), pch=16, xlab=expression(paste(A[l]:A[s]," (",cm^2*mm^-2,")" )), ylab="Stem Length (cm)", log="y")
mtext("Branch Length (cm)", side=2, line=2)
mtext(expression(paste(A[l]:A[s]," (",cm^2*mm^-2,")" )), side=1, line=2)
modfit <- summary(lmer(logLength~mAl_As + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=2)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(mAl_As)
  fitted <- predict(lm(mLength~mAl_As, tmp),newdata=data.frame("mAl_As"=c(min(tmp$mAl_As, na.rm=T), max(tmp$mAl_As, na.rm=T))))
  lines(x=c(min(tmp$mAl_As, na.rm=T), max(tmp$mAl_As, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mLength~mAl_As, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("f)", side=3, line=-1, adj=0.05)


### ml_ms
palette(paste0(mypal,"77"))
plot(perc_maxBAI~mml_ms, wp.ind, col=factor(Site), pch=16,
     xaxt="n", yaxt="n", ylab="")
modfit <- summary(lmer(perc_maxBAI~mml_ms + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(mml_ms)
  fitted <- predict(lm(perc_maxBAI~mml_ms, tmp),newdata=data.frame("mml_ms"=c(min(tmp$mml_ms, na.rm=T), max(tmp$mml_ms, na.rm=T))))
  lines(x=c(min(tmp$mml_ms, na.rm=T), max(tmp$mml_ms, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mperc_maxBAI~mml_ms, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("b)", side=3, line=-1, adj=0.05)

palette(paste0(mypal,"77"))
plot(mLength~mml_ms, wp.ind, col=factor(Site), pch=16, xlab=expression(paste(A[l]:A[s]," (",cm^2*mm^-2,")" )), ylab="Stem Length (cm)", log="y")
#mtext("Branch Length (cm)", side=2, line=2)
mtext(expression(paste(M[l]:M[s]," (",g*g^-1,")" )), side=1, line=2)
modfit <- summary(lmer(logLength~mml_ms + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(mml_ms)
  fitted <- predict(lm(mLength~mml_ms, tmp),newdata=data.frame("mml_ms"=c(min(tmp$mml_ms, na.rm=T), max(tmp$mml_ms, na.rm=T))))
  lines(x=c(min(tmp$mml_ms, na.rm=T), max(tmp$mml_ms, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mLength~mml_ms, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("g)", side=3, line=-1, adj=0.05)


### leaf size
palette(paste0(mypal,"77"))
plot(perc_maxBAI~mleafsize, wp.ind, col=factor(Site), pch=16,
     xaxt="n", yaxt="n", ylab="")
modfit <- summary(lmer(perc_maxBAI~mleafsize + (1|Site), wp.ind))
#abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(mleafsize)
  fitted <- predict(lm(perc_maxBAI~mleafsize, tmp),newdata=data.frame("mleafsize"=c(min(tmp$mleafsize, na.rm=T), max(tmp$mleafsize, na.rm=T))))
  lines(x=c(min(tmp$mleafsize, na.rm=T), max(tmp$mleafsize, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mperc_maxBAI~mleafsize, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("d)", side=3, line=-1, adj=0.05)

palette(paste0(mypal,"77"))
plot(mLength~mleafsize, wp.ind, col=factor(Site), pch=16, xlab=expression(paste(A[l]:A[s]," (",cm^2*mm^-2,")" )), ylab="Stem Length (cm)", log="y", yaxt="n")
#mtext("Branch Length (cm)", side=2, line=2)
mtext(expression(paste("Leaf Size (",cm^2,")" )), side=1, line=2)
modfit <- summary(lmer(logLength~mleafsize + (1|Site), wp.ind))
#abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(mleafsize)
  fitted <- predict(lm(mLength~mleafsize, tmp),newdata=data.frame("mleafsize"=c(min(tmp$mleafsize, na.rm=T), max(tmp$mleafsize, na.rm=T))))
  lines(x=c(min(tmp$mleafsize, na.rm=T), max(tmp$mleafsize, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mLength~mleafsize, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("h)", side=3, line=-1, adj=0.05)


### LMA
palette(paste0(mypal,"77"))
plot(perc_maxBAI~mLMA, wp.ind, col=factor(Site), pch=16,
     xaxt="n", yaxt="n", ylab="")
modfit <- summary(lmer(perc_maxBAI~mLMA + (1|Site), wp.ind))
#abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(mLMA)
  fitted <- predict(lm(perc_maxBAI~mLMA, tmp),newdata=data.frame("mLMA"=c(min(tmp$mLMA, na.rm=T), max(tmp$mLMA, na.rm=T))))
  lines(x=c(min(tmp$mLMA, na.rm=T), max(tmp$mLMA, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mperc_maxBAI~mLMA, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("d)", side=3, line=-1, adj=0.05)

palette(paste0(mypal,"77"))
plot(mLength~mLMA, wp.ind, col=factor(Site), pch=16, xlab=expression(paste(A[l]:A[s]," (",cm^2*mm^-2,")" )), ylab="Stem Length (cm)", log="y", yaxt="n")
#mtext("Branch Length (cm)", side=2, line=2)
mtext(expression(paste("LMA (",g*cm^-2,")" )), side=1, line=2)
modfit <- summary(lmer(logLength~mLMA + (1|Site), wp.ind))
#abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(mLMA)
  fitted <- predict(lm(mLength~mLMA, tmp),newdata=data.frame("mLMA"=c(min(tmp$mLMA, na.rm=T), max(tmp$mLMA, na.rm=T))))
  lines(x=c(min(tmp$mLMA, na.rm=T), max(tmp$mLMA, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mLength~mLMA, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("i)", side=3, line=-1, adj=0.05)


### LDMC
palette(paste0(mypal,"77"))
plot(perc_maxBAI~mLDMC, wp.ind, col=factor(Site), pch=16,
     xaxt="n", yaxt="n", ylab="")
modfit <- summary(lmer(perc_maxBAI~mLDMC + (1|Site), wp.ind))
#abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2)
palette(mypal)
for (i in unique(wp.ind$Site[which(wp.ind$mLDMC>0 & wp.ind$perc_maxBAI>0)])) {
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(mLDMC)
  fitted <- predict(lm(perc_maxBAI~mLDMC, tmp),newdata=data.frame("mLDMC"=c(min(tmp$mLDMC, na.rm=T), max(tmp$mLDMC, na.rm=T))))
  lines(x=c(min(tmp$mLDMC, na.rm=T), max(tmp$mLDMC, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mperc_maxBAI~mLDMC, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("e)", side=3, line=-1, adj=0.05)

palette(paste0(mypal,"77"))
plot(mLength~mLDMC, wp.ind, col=factor(Site), pch=16, xlab=expression(paste(A[l]:A[s]," (",cm^2*mm^-2,")" )), ylab="Stem Length (cm)", log="y", yaxt="n")
#mtext("Branch Length (cm)", side=2, line=2)
mtext(expression(paste("LDMC (",g[dry]*g[wet]^-1,")" )), side=1, line=2)
modfit <- summary(lmer(logLength~mLDMC + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=2)
palette(mypal)
for (i in unique(wp.ind$Site[which(wp.ind$mLDMC>0 & wp.ind$perc_maxBAI>0)])){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(mLDMC)
  fitted <- predict(lm(mLength~mLDMC, tmp),newdata=data.frame("mLDMC"=c(min(tmp$mLDMC, na.rm=T), max(tmp$mLDMC, na.rm=T))))
  lines(x=c(min(tmp$mLDMC, na.rm=T), max(tmp$mLDMC, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mLength~mLDMC, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("j)", side=3, line=-1, adj=0.05)

if(save.figures==T){quartz.save(file=paste0(results.dir,"/Fig7_Trait-Growth_v1.pdf"),type = "pdf")}




#++++++++++++++++++++++++++++++++++++++++++++++++++++
############# . FIG 8: HOTTER Simulations ############

quartz(width=6, height=5)
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2,1,0))


# Panel a) Prognostic MD WP
palette(mypal)
plot(MD~psi_leaf_mean, wp.ind, col=factor(Site), pch=16
     ,ylab=expression(paste("Observed ",Psi[MD]," (MPa)"))
     ,xlab= expression(paste("Simulated ", Psi[MD]," (MPa)" ))
     ,xlim=c(-5.5,-2), ylim=c(-5.5,-2));abline(a=0,b=1)
#points(MD~psi_leaf_mean, wp.ind, col=factor(Site), pch=16)
#points(MD~psi_leaf_max, wp.ind, col=factor(Site), pch=2)
palette(paste0(mypal,"77"))
arrows(x0 = wp.ind$psi_leaf_min, x1=wp.ind$psi_leaf_max, y0=wp.ind$MD, col=factor(wp.ind$Site), length = 0, lwd=2)
mtext("a)", side=3, line=-1, adj=0.05)

# Panel b) PLC not related to growth
palette(paste0(mypal,"77"))
plot(perc_maxBAI~PLC_mean, wp.ind, col=factor(Site), pch=16
     ,ylab="Observed BAI Growth (% of max)"
     ,xlab= "Simulated PLC (%)" )
# for (i in unique(wp.ind$Site[which(wp.ind$PLC_mean>0 & wp.ind$perc_maxBAI>0)])){
#   tmp <- wp.ind %>% filter(Site==i) %>% arrange(GPP_mean)
#   fitted <- predict(lm(perc_maxBAI~PLC_mean, tmp),newdata=data.frame("PLC_mean"=c(min(tmp$PLC_mean, na.rm=T), max(tmp$PLC_mean, na.rm=T))))
#   lines(x=c(min(tmp$PLC_mean, na.rm=T), max(tmp$PLC_mean, na.rm=T)), y=fitted, col="#55555555", lwd=2)
# }
palette(mypal)
points(mperc_maxBAI~PLC_mean, wp.site, col=factor(Site), pch=17, cex=site.cex)
#arrows(x0 = 100 - 100*wp.ind$K_frac_min, x1=100- 100*wp.ind$K_frac_max, y0=wp.ind$perc_maxBAI, col=factor(wp.ind$Site), length = 0, lwd=2)
mtext("b)", side=3, line=-1, adj=0.05)

# Panel c) GPP negatively related to growth
palette(paste0(mypal,"77"))
plot(perc_maxBAI~GPP_mean, wp.ind, col=factor(Site), pch=16
     ,ylab="Observed BAI Growth (% of max)"
     ,xlab= "Simulated GPP (gC/day)" )
summary(lmer(perc_maxBAI~GPP_mean + (1|Site), wp.ind))
modfit <- summary(lmer(perc_maxBAI~GPP_mean + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=2)
for (i in unique(wp.ind$Site[which(wp.ind$GPP_mean>0 & wp.ind$perc_maxBAI>0)])){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(GPP_mean)
  fitted <- predict(lm(perc_maxBAI~GPP_mean, tmp),newdata=data.frame("GPP_mean"=c(min(tmp$GPP_mean, na.rm=T), max(tmp$GPP_mean, na.rm=T))))
  lines(x=c(min(tmp$GPP_mean, na.rm=T), max(tmp$GPP_mean, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
palette(mypal)
points(mperc_maxBAI~GPP_mean, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("c)", side=3, line=-1, adj=0.05)

# Panel d) NPP negative, but positively related to growth
palette(paste0(mypal,"77"))
plot(perc_maxBAI~NPP_mean, wp.ind, col=factor(Site), pch=16
     ,ylab="Observed BAI Growth (% of max)"
     ,xlab= "Simulated NPP (gC/day)" )
summary(lmer(perc_maxBAI~NPP_mean + (1|Site), wp.ind))
modfit <- summary(lmer(perc_maxBAI~NPP_mean + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
for (i in unique(wp.ind$Site[which(wp.ind$NPP_mean<0 & wp.ind$perc_maxBAI>0)])){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(NPP_mean)
  fitted <- predict(lm(perc_maxBAI~NPP_mean, tmp),newdata=data.frame("NPP_mean"=c(min(tmp$NPP_mean, na.rm=T), max(tmp$NPP_mean, na.rm=T))))
  lines(x=c(min(tmp$NPP_mean, na.rm=T), max(tmp$NPP_mean, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
palette(mypal)
points(mperc_maxBAI~NPP_mean, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("d)", side=3, line=-1, adj=0.05)

if(save.figures==T){quartz.save(file=paste0(results.dir,"/Fig8_HotterSims_v1.pdf"),type = "pdf")}


# Compare correlations of full and LL simulations
  # but need to remove 'runaway' trees in full sims for apples-to-apples comparison
cor.test(wp.ind$psi_leaf_mean[which(wp.ind$PLC_mean<99)],wp.ind$MD[which(wp.ind$PLC_mean<99)])
cor.test(wp.ind$psi_leaf_mean_ll,wp.ind$MD)








########## ^ End: Main Text #################################################





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############## * Supplemental Figures and Analyses #############################


#_________________________________________________________________________________
########### . Length Correction Analysis Supplemental ###########
### branch level correlations
# chars <- c("Length","Leaf_num","Al_As","LDMC","LMA","leafsize","ml_ms")
# Mypairs(branch.traits[,chars])

# quick 
branch.traits$ml_ms.resid <- NA
branch.traits$ml_ms.resid[which(branch.traits$ml_ms>0)] <- resid(lm(ml_ms~log(Length,10), branch.traits)) #make a length corrected version


# Remove effect of branch length on Al:As by pulling out residuals
# note: it's a log-log relationship, so need to back transform prediction and calculate raw residuals
# Version with unlogged Al_As
branch.traits$Al_As.resid <- NA
branch.traits$Al_As.resid[which(branch.traits$Al_As>0)] <- resid(lm(Al_As~log(Length,10), branch.traits)) #make a length corrected version
# correct version with transformed predictions of log(Al_As)
al.aspreds <- predict(lm(log(Al_As)~log(Length), branch.traits))
al.aspreds.raw <- exp(al.aspreds)
plot(Al_As~log(Length), branch.traits)
points(al.aspreds.raw~log(branch.traits$Length[which(branch.traits$Al_As>0)]), col="red", pch=".")
branch.traits$Al_As.resid2 <- NA
branch.traits$Al_As.resid2[which(branch.traits$Al_As>0)] <- branch.traits$Al_As[which(branch.traits$Al_As>0)] - al.aspreds.raw #make a length corrected version
#plot(Al_As.resid2~Al_As.resid, branch.traits)
# turns out, they're almost identical

# Remove effect of branch length on Ml:Ms
ml_mspreds <- predict(lm(log(ml_ms)~log(Length), branch.traits))
ml_mspreds.raw <- exp(ml_mspreds)
plot(ml_ms~log(Length), branch.traits)
points(ml_mspreds.raw~log(branch.traits$Length[which(branch.traits$ml_ms>0)]), col="red", pch=".")
branch.traits$ml_ms.resid2 <- NA
branch.traits$ml_ms.resid2[which(branch.traits$ml_ms>0)] <- branch.traits$ml_ms[which(branch.traits$ml_ms>0)] - ml_mspreds.raw #make a length corrected version

# get these branch-level corrected values into the tree-average df
resid.df <- branch.traits %>% group_by(pop, tree) %>% summarise(ml_ms.resid = mean(ml_ms.resid2, na.rm=T), Al_As.resid=mean(Al_As.resid2, na.rm=T))
test <- left_join(wp.ind, resid.df, by=c("Site"="pop","Tree"="tree"))


# Al:As still negatively related to BAI when controlling for branch length
summary(lmer(perc_maxBAI~Al_As.resid + (1|Site), test)) # p=0.009
summary(lmer(perc_maxBAI~mAl_As + logLength + (1|Site), test)) #p= 0.009 or 0.01 (logLength)


summary(lmer(perc_maxBAI~ml_ms.resid + (1|Site), test)) # p=0.025 still significant and negative
summary(lmer(perc_maxBAI~logLength + mml_ms + (1|Site), wp.ind)) # no longer significant
summary(lmer(perc_maxBAI~logLength + logml_ms + (1|Site), wp.ind)) # still significant and negative

# Takehome point: Al:As and Ml:Ms are still negatively correlated with % max BAI after controlling for length

# plot(perc_maxBAI~ml_ms.resid, test, pch=16, col=factor(Site))
# modfit <- summary(lmer(perc_maxBAI~ml_ms.resid + (1|Site), test))
# abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
# 
# 
# plot(perc_maxBAI~Al_As.resid, test, pch=16, col=factor(Site))
# modfit <- summary(lmer(perc_maxBAI~Al_As.resid + (1|Site), test))
# abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
# palette(mypal)
# for (i in unique(test$Site[which(test$Al_As.resid>0 & test$perc_maxBAI>0)])){
#   tmp <- test %>% filter(Site==i) %>% arrange(Al_As.resid)
#   fitted <- predict(lm(perc_maxBAI~Al_As.resid, tmp),newdata=data.frame("Al_As.resid"=c(min(tmp$Al_As.resid, na.rm=T), max(tmp$Al_As.resid, na.rm=T))))
#   lines(x=c(min(tmp$Al_As.resid, na.rm=T), max(tmp$Al_As.resid, na.rm=T)), y=fitted, col="#55555555", lwd=2)
# }
# #points(mLength~mLDMC, wp.site, col=factor(Site), pch=17, cex=site.cex)








#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### . FIG S1: Predawn f(DBH) ###########

quartz(width=3, height=3)
par(mar=c(3.3,3.3,1,1), mgp=c(2.2,1,0))
palette(paste0(mypal,"77"))
plot(PD~DBH, wp.ind, col=factor(Site), pch=16, ylab=expression(paste(Psi[PD]," (MPa)" )), xlab="DBH (cm)")
modfit <- summary(lmer(PD~DBH + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(DBH)
  fitted <- predict(lm(PD~DBH, tmp),newdata=data.frame("DBH"=c(min(tmp$DBH, na.rm=T), max(tmp$DBH, na.rm=T))))
  lines(x=c(min(tmp$DBH, na.rm=T), max(tmp$DBH, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mPD~mDBH, wp.site, col=factor(Site), pch=17, cex=site.cex)
#mtext("DBH (cm)", side=1, line=2)

if(save.figures==T){quartz.save(file=paste0(results.dir,"/Fig_WP_DBH_v2.pdf"),type = "pdf")}




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### . FIG S2: Predawn f(del18O) ###########


quartz(width=3.2, height=6)
par(mfrow=c(3,1), mar=c(0,3.3,0,1), oma=c(3,0,1,0), mgp=c(2,1,0), cex=1)
palette(paste0(mypal,"77"))
plot(PD~del18O, wp.ind, col=factor(Site), pch=16, ylab=expression(paste(Psi[PD]," (MPa)" ))
     , xlab=expression(paste(delta,"D (\u2030)")), xaxt="n")
modfit <- summary(lmer(PD~del18O + (1|Site), wp.ind[which(wp.ind$Site != "PWD"),]))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(PD)
  fitted <- predict(lm(PD~del18O, tmp),newdata=data.frame("del18O"=c(min(tmp$del18O, na.rm=T), max(tmp$del18O, na.rm=T))))
  lines(x=c(min(tmp$del18O, na.rm=T), max(tmp$del18O, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mPD~mdel18O, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("a)", side=3, line=-1, adj=0.05)

palette(paste0(mypal,"77"))
plot(MD~del18O, wp.ind, col=factor(Site), pch=16, ylab=expression(paste(Psi[MD]," (MPa)" )), xlab=expression(paste(delta,"D (\u2030)")), xaxt="n")
modfit <- summary(lmer(MD~del18O + (1|Site), wp.ind[which(wp.ind$Site != "PWD"),]))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(del18O)
  fitted <- predict(lm(MD~del18O, tmp),newdata=data.frame("del18O"=c(min(tmp$del18O, na.rm=T), max(tmp$del18O, na.rm=T))))
  lines(x=c(min(tmp$del18O, na.rm=T), max(tmp$del18O, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mMD~mdel18O, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("b)", side=3, line=-1, adj=0.05)

palette(paste0(mypal,"77"))
plot(E.drop~del18O, wp.ind, col=factor(Site), pch=16, ylab=expression(paste(Delta*Psi," (",Psi[PD]-Psi[MD]," (MPa)")), xlab=expression(paste(delta,"D (\u2030)")))
modfit <- summary(lmer(E.drop~del18O + (1|Site), wp.ind[which(wp.ind$Site != "PWD"),]))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
palette(mypal)
for (i in levels(wp.ind$Site)[-c(1,3,7,8)]){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(del18O)
  fitted <- predict(lm(E.drop~del18O, tmp),newdata=data.frame("del18O"=c(min(tmp$del18O, na.rm=T), max(tmp$del18O, na.rm=T))))
  lines(x=c(min(tmp$del18O, na.rm=T), max(tmp$del18O, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
points(mE.drop~mdel18O, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext(text=expression(paste(delta^18,"O (\u2030)")),side = 1, line=2)
mtext("c)", side=3, line=-1, adj=0.05)

if(save.figures==T){quartz.save(file=paste0(results.dir,"/Fig_WP_del18O_v1.pdf"),type = "pdf")}
#+++++++++++++++++++++++++++++++++++++++++++++++++





#++++++++++++++++++++++++++++++++++++++++++++++++++
###########. FIG S3: delDprecip and WP #######

quartz(width=6.8, height=2.4)
par(mfrow=c(1,3), mar=c(3,3,1.5,1),mgp=c(2,1,0), cex=1)


plot(delD~wydelD, wp.ind, pch=16, col=Site,
     xlab=expression(paste(delta*D[precip]," (\u2030)")),
     ylab= expression(paste(delta*D[xylem]," (\u2030)")))
abline(lm(delD~wydelD, wp.ind))
mtext("a)", side=3, line=0.1, adj=-0.05)

plot(PD~wydelD, wp.ind, pch=16, col=Site,
     xlab=expression(paste(delta*D[precip]," (\u2030)")),
     ylab= expression(paste(Psi[PD]," (MPa)")))
points(PD~wydelD, wp.ind[which(wp.ind$Site=="PWD"),])
abline(lm(PD~wydelD, wp.ind[which(wp.ind$Site!= "PWD"),]))
mtext("b)", side=3, line=0.1, adj=-0.05)

plot(PD~delD.resid, wp.ind, pch=16, col=Site,
     xlab=expression(paste(delta*D[xylem]," residual (\u2030)")),
     ylab= expression(paste(Psi[PD]," (MPa)")))
points(PD~delD.resid, wp.ind[which(wp.ind$Site=="PWD"),])
abline(lm(PD~delD.resid, wp.ind[which(wp.ind$Site!= "PWD"),]))
mtext("c)", side=3, line=0.1, adj=-0.05)

if(save.figures==T){quartz.save(file=paste0(results.dir,"/FigS_delDprecip_v1.pdf"),type = "pdf")}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########. FIG S4: Trait-climate Only ####################################

quartz(width=5.8, height=4.5)
par(mfrow=c(2,3), mar=c(3,3,1,1), mgp=c(2,1,0), cex=1)

# old with only climate
plot(mAl_As~tc.tmin.2018wy, wp.ind, pch=16, col=tree.col
     , ylab=expression(paste(A[leaf]:A[stem]))
     , xlab="Min Temp 2018wy (°C)")
points(mAl_As~tc.tmin.2018wy, wp.site, pch=16, cex=site.cex)
modfit <- summary(lmer(mAl_As~tc.tmin.2018wy + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=2)
mtext("a)", side=3, line=0.1, adj=-0.1)

plot(mml_ms~tc.tmin.2018wy, wp.ind, pch=16, col=tree.col
     , ylab=expression(paste(Mass[leaf]:Mass[stem]))
     , xlab="Min Temp 2018wy (°C)")
points(mml_ms~tc.tmin.2018wy, wp.site, pch=16, cex=site.cex)
modfit <- summary(lmer(mml_ms~tc.tmin.2018wy + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
mtext("b)", side=3, line=0.1, adj=-0.1)

plot(mleafsize~tc.ppt, wp.ind, pch=16, col=tree.col
     , ylab="mean leaf size (cm2)"
     , xlab="30yr MAP (mm)")
points(mleafsize~tc.ppt, wp.site, pch=16, cex=site.cex)
modfit <- summary(lmer(mleafsize~tc.ppt + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
mtext("c)", side=3, line=0.1, adj=-0.1)


# old, climate alone
plot(mLMA~tc.aet, wp.ind, pch=16, col=tree.col
     , ylab="LMA (g/cm2)"
, xlab="30yr AET (mm)")
points(mLMA~tc.aet, wp.site, pch=16, cex=site.cex)
#modfit <- summary(lmer(mLMA~tc.aet + (1|Site), wp.ind))
#abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
mtext("d)", side=3, line=0.1, adj=-0.1)

plot(mLDMC~tc.tmax.2018sp, wp.ind, pch=16, col=tree.col
     , ylab="LDMC (g/g)"
     , xlab="2018 spring Tmax (°C)")
points(mLDMC~tc.tmax.2018sp, wp.site, pch=16, cex=site.cex)
mtext("e)", side=3, line=0.1, adj=-0.1)

if(save.figures==T){quartz.save(file=paste0(results.dir,"/FigS4_Traits_climateonly_v2.pdf"),type = "pdf")}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############ . FIG S5: Trait-Water Potentials ################
summary(lmer(mml_ms~PD + (1|Site), wp.ind))
summary(lmer(mLDMC~MD + (1|Site), wp.ind))

#Mypairs
#Make fancy pair plots
Mypairs <- function(Z, variable.names) {
  MyVarx <- colnames(Z)
  pairs(Z, labels = variable.names,
        cex.labels =  2,
        lower.panel = function(x, y, digits=2, prefix="", cex.cor = 7) {
          panel.cor(x, y, digits, prefix, cex.cor)}, 
        upper.panel =  function(x, y) { 
          points(x, y, 
                 pch = 16, cex = 0.8, 
                 col = gray(0.1))
          if(cor.test(y,x)$p.value <0.05){abline(lm(y~x))}
        })
  #print(P)
}
chars <- c("PD","MD","E.drop","mAl_As","mml_ms","mleafsize","mLMA","mLDMC")
quartz(width=6.8, height=6,2)
Mypairs(wp.ind[,chars], variable.names = c(expression(paste(Psi[PD])), expression(paste(Psi[MD])),expression(paste(Delta*Psi)),expression(paste(A[l]:A[s])),expression(paste(M[l]:M[s])), "leafsize", "LMA","LDMC"))
if(save.figures==T){quartz.save(file=paste0(results.dir,"/Fig_WP_Traits_corrs_v1.pdf"),type = "pdf")}




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###### Fig S6: PLC ~ MD for Leaf Area comparison #########


quartz(width=6, height=3.1)
par(mfrow=c(1,2),mar=c(3,3.5,1.5,1), mgp=c(2,1,0))

# a) PLC pretty reasonable based on MDs - true leaf area
palette(mypal)
plot(PLC_mean~MD, wp.ind, col=factor(Site), pch=16
     ,xlab=expression(paste("Observed ",Psi[MD]," (MPa)"))
     ,ylab= "Simulated PLC (%)" 
     ,ylim= c(0,100))
palette(paste0(mypal,"77"))
arrows(y0 = 100 - 100*wp.ind$K_frac_min, y1=100- 100*wp.ind$K_frac_max, x0=wp.ind$MD, col=factor(wp.ind$Site), length = 0, lwd=2)
mtext("a) Full Leaf Area", side=3, line=0, adj=0.05)
abline(v=-3.88, col="grey", lty=2)
abline(v=-4.47, col="black", lty=2)
#points(PLC_mean_ll~MD, wp.ind, col=factor(Site))

# b) PLC as f(observed MD) Low Leaf Area
palette(mypal)
plot(PLC_mean_ll~MD, wp.ind, col=factor(Site), pch=16
     ,xlab=expression(paste("Observed ",Psi[MD]," (MPa)"))
     ,ylab= "Simulated PLC (%)" 
     ,ylim=c(0,100))
palette(paste0(mypal,"77"))
arrows(y0 = 100 - 100*wp.ind$K_frac_min_ll, y1=100- 100*wp.ind$K_frac_max_ll, x0=wp.ind$MD, col=factor(wp.ind$Site), length = 0, lwd=2)
mtext("b) Half Leaf Area", side=3, line=0, adj=0.05)
abline(v=-3.88, col="grey", lty=2)
abline(v=-4.47, col="black", lty=2)


if(save.figures==T){quartz.save(file=paste0(results.dir,"/FigS_HotterSimPrognostic_v2.pdf"),type = "pdf")}



# 
# # bnew) PLC  based on PDs
# palette(mypal)
# plot(I(100-100*K_frac_mean)~PD, wp.ind, col=factor(Site), pch=16
#      ,xlab=expression(paste("Observed ",Psi[PD]," (MPa)"))
#      ,ylab= "Simulated PLC (%)" )
# palette(paste0(mypal,"77"))
# arrows(y0 = 100 - 100*wp.ind$K_frac_min, y1=100- 100*wp.ind$K_frac_max, x0=wp.ind$PD, col=factor(wp.ind$Site), length = 0, lwd=2)
# mtext("b)", side=3, line=0, adj=0.05)
# abline(v=-3.88, col="grey", lty=2)
# abline(v=-4.47, col="black", lty=2)
# 
# points(I(100-100*K_frac_mean_ll)~PD, wp.ind, col=factor(Site), pch=1
#      ,xlab=expression(paste("Observed ",Psi[MD]," (MPa)"))
#      ,ylab= "Simulated PLC (%)" )







#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###. FIG S7: Water potential comparison with low leaf are: ######

quartz(width=6, height=5)
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2,1,0))


# Panel a) Prognostic MD WP
palette(mypal)
plot(MD~psi_leaf_mean_ll, wp.ind, col=factor(Site), pch=16
     ,ylab=expression(paste("Observed ",Psi[MD]," (MPa)"))
     ,xlab= expression(paste("Simulated ", Psi[MD]," (MPa)" ))
     ,xlim=c(-5.5,-2), ylim=c(-5.5,-2));abline(a=0,b=1)
#points(MD~psi_leaf_mean, wp.ind, col=factor(Site), pch=16)
#points(MD~psi_leaf_max, wp.ind, col=factor(Site), pch=2)
palette(paste0(mypal,"77"))
arrows(x0 = wp.ind$psi_leaf_min_ll, x1=wp.ind$psi_leaf_max_ll, y0=wp.ind$MD, col=factor(wp.ind$Site), length = 0, lwd=2)
mtext("a)", side=3, line=-1, adj=0.05)

# Panel b) PLC not related to growth
palette(paste0(mypal,"77"))
plot(perc_maxBAI~PLC_mean_ll, wp.ind, col=factor(Site), pch=16
     ,ylab="Observed BAI Growth (% of max)"
     ,xlab= "Simulated PLC (%)" )
# for (i in unique(wp.ind$Site[which(wp.ind$PLC_mean>0 & wp.ind$perc_maxBAI>0)])){
#   tmp <- wp.ind %>% filter(Site==i) %>% arrange(GPP_mean)
#   fitted <- predict(lm(perc_maxBAI~PLC_mean, tmp),newdata=data.frame("PLC_mean"=c(min(tmp$PLC_mean, na.rm=T), max(tmp$PLC_mean, na.rm=T))))
#   lines(x=c(min(tmp$PLC_mean, na.rm=T), max(tmp$PLC_mean, na.rm=T)), y=fitted, col="#55555555", lwd=2)
# }
palette(mypal)
points(mperc_maxBAI~PLC_mean_ll, wp.site, col=factor(Site), pch=17, cex=site.cex)
#arrows(x0 = 100 - 100*wp.ind$K_frac_min, x1=100- 100*wp.ind$K_frac_max, y0=wp.ind$perc_maxBAI, col=factor(wp.ind$Site), length = 0, lwd=2)
mtext("b)", side=3, line=-1, adj=0.05)

# Panel c) GPP negatively related to growth
palette(paste0(mypal,"77"))
plot(perc_maxBAI~GPP_mean_ll, wp.ind, col=factor(Site), pch=16
     ,ylab="Observed BAI Growth (% of max)"
     ,xlab= "Simulated GPP (gC/day)" )
summary(lmer(perc_maxBAI~GPP_mean_ll + (1|Site), wp.ind))
modfit <- summary(lmer(perc_maxBAI~GPP_mean_ll + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=1)
for (i in unique(wp.ind$Site[which(wp.ind$GPP_mean_ll>0 & wp.ind$perc_maxBAI>0)])){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(GPP_mean_ll)
  fitted <- predict(lm(perc_maxBAI~GPP_mean_ll, tmp),newdata=data.frame("GPP_mean_ll"=c(min(tmp$GPP_mean_ll, na.rm=T), max(tmp$GPP_mean_ll, na.rm=T))))
  lines(x=c(min(tmp$GPP_mean_ll, na.rm=T), max(tmp$GPP_mean_ll, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
palette(mypal)
points(mperc_maxBAI~GPP_mean_ll, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("c)", side=3, line=-1, adj=0.05)

# Panel d) NPP negative, but positively related to growth
palette(paste0(mypal,"77"))
plot(perc_maxBAI~NPP_mean_ll, wp.ind, col=factor(Site), pch=16
     ,ylab="Observed BAI Growth (% of max)"
     ,xlab= "Simulated NPP (gC/day)" )
summary(lmer(perc_maxBAI~NPP_mean_ll + (1|Site), wp.ind))
modfit <- summary(lmer(perc_maxBAI~NPP_mean_ll + (1|Site), wp.ind))
abline(a=modfit$coefficients[1,1], b=modfit$coefficients[2,1], lwd=2, lty=2)
for (i in unique(wp.ind$Site[which(wp.ind$NPP_mean_ll<0 & wp.ind$perc_maxBAI>0)])){
  tmp <- wp.ind %>% filter(Site==i) %>% arrange(NPP_mean_ll)
  fitted <- predict(lm(perc_maxBAI~NPP_mean_ll, tmp),newdata=data.frame("NPP_mean_ll"=c(min(tmp$NPP_mean_ll, na.rm=T), max(tmp$NPP_mean_ll, na.rm=T))))
  lines(x=c(min(tmp$NPP_mean_ll, na.rm=T), max(tmp$NPP_mean_ll, na.rm=T)), y=fitted, col="#55555555", lwd=2)
}
palette(mypal)
points(mperc_maxBAI~NPP_mean_ll, wp.site, col=factor(Site), pch=17, cex=site.cex)
mtext("d)", side=3, line=-1, adj=0.05)

if(save.figures==T){quartz.save(file=paste0(results.dir,"/FigS_HotterSims_ll_v1.pdf"),type = "pdf")}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#####. FIG S8: LA comparison on NPP and PLC ####################
quartz(width=6, height=3.1)
par(mfrow=c(1,2),mar=c(3,3.5,1.5,1), mgp=c(2,1,0))

palette(mypal)
plot(NPP_mean~PLC_mean, wp.ind, col=factor(Site), pch=3
     , ylim=c(-250, 0), xlim=c(0,100)
     , ylab="Simulated NPP (gC/day)"
     , xlab="Simulated PLC (%)")
palette(paste0(mypal,"77"))
points(NPP_mean_ll~PLC_mean_ll, wp.ind, pch=21, bg=factor(Site))
mtext("a)", side=3, line=0, adj=0.05)


palette(mypal)
plot(NPP_mean~al_true, wp.ind, col=factor(Site), pch=3, ylim=c(-250, 0), xlim=c(10,140)
     , ylab="Simulated NPP (gC/day)"
     , xlab="Tree Leaf Area (m2)")
palette(paste0(mypal,"77"))
points(NPP_mean_ll~al_low, wp.ind, pch=21, bg=factor(Site))
legend('bottomleft', legend=c("full LA","50% LA"), pch=c(3,21), pt.bg="gray")
mtext("b)", side=3, line=0, adj=0.05)


if(save.figures==T){quartz.save(file=paste0(results.dir,"/FigS_HotterSim_LeafAreaSensitivity_v1.pdf"),type = "pdf")}

