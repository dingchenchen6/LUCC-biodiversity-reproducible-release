##%######################################################%##
#                                                          #
####    Models assessing differences between realms     ####
#                                                          #
##%######################################################%##


# In this script, separate models are run for subsets of the PREDICTS database
# based on which realm the sites are found in. 

# Inclusion of Realm in the model results in high multicolinearity between variables
# so subsets of the data and separate models are run instead.

# directories
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "9_RunModels_Tropical/"
if(!dir.exists(outDir)) dir.create(outDir)

sink(paste0(outDir,"log.txt"))
sink()
t.start <- Sys.time()

print(t.start)


# load libraries
library(StatisticalModels)
library(predictsFunctions)
source("Functions.R")
library(ggplot2)
library(cowplot)
library(sjPlot)
library(performance)
#library(pedometrics)
#library(car)

# read in PREDICTS data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData1.rds"))

#table(predictsSites$UI2,predictsSites$StdTmeanAnomalyRS)
with(predictsSites,
     table(UI2 == "Urban", is.na(StdTmeanAnomalyRS)))
predictsSites %>%
  dplyr::mutate(
    is_urban = UI2 == "Urban",
    is_na    = is.na(StdTmeanAnomalyRS)
  ) %>%
  dplyr::count(is_urban, is_na)


# categorise sites into temperate/tropical
# get the tropical values
predictsSites$Tropical <- NA

predictsSites[predictsSites$Latitude > -23.44 & predictsSites$Latitude < 23.44, 'Tropical'] <- "Tropical"

# label the remaining as temperate
predictsSites[is.na(predictsSites$Tropical), 'Tropical'] <- "Temperate"

# set as a factor
predictsSites$Tropical <- as.factor(predictsSites$Tropical)
# levels: Temperate Tropical


predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ] #6069

table(predictsSites$Tropical)

# Temperate  Tropical 
# 3498      2959 

# abundance data subset
table(predictsSites[!is.na(predictsSites$LogAbund), 'Tropical'])

# Temperate  Tropical 
# 2929      2805

# split by land use classes
table(predictsSites$UI2, predictsSites$Tropical)

# Temperate Tropical
# Primary vegetation        1377     1135
# Agriculture_High           548      749
# Agriculture_Low            474      525
# Secondary vegetation       717      494
# Urban                      382       56

# look into range of SCA values across landuses in tropical/temperate areas.

plot_data <- predictsSites[, c("Tropical", "UI2", "StdTmeanAnomaly")]
plot_data <- plot_data[!is.na(plot_data$UI2), ]
#table(plot_data$UI2)
plot_data$UI2 <- factor(plot_data$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High","Urban"))

ggplot(data = plot_data) +
  geom_freqpoly(aes(x = StdTmeanAnomaly, col = UI2)) +
  facet_wrap(~ Tropical) + 
  theme_bw()

plot_data <- predictsSites[, c("Tropical", "UI2", "StdTmaxAnomaly")]
plot_data <- plot_data[!is.na(plot_data$UI2), ]
#table(plot_data$UI2)
plot_data$UI2 <- factor(plot_data$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High","Urban"))

ggplot(data = plot_data) +
  geom_freqpoly(aes(x = StdTmaxAnomaly, col = UI2)) +
  facet_wrap(~ Tropical) + 
  theme_bw()

### look at correlations between climate metrics for realm subsets ###

trop <- predictsSites[predictsSites$Tropical == "Tropical", ]

temp <- predictsSites[predictsSites$Tropical == "Temperate", ]

# trop
cor(trop$avg_temp, trop$TmeanAnomaly) # -0.17
cor(trop$avg_temp, trop$StdTmeanAnomaly) # 0.03
cor(trop$TmeanAnomaly, trop$StdTmeanAnomaly) # 0.74

#temp
cor(temp$avg_temp, temp$TmeanAnomaly) # -0.04
cor(temp$avg_temp, temp$StdTmeanAnomaly) # -0.60
cor(temp$TmeanAnomaly, temp$StdTmeanAnomaly) # 0.50



##%######################################################%##
#                                                          #
####               Models for each realm                ####
#                                                          #
##%######################################################%##


# since the model including Realm and its interactions with other variables
# suffers from multicolinearity, we are running a separate model for subsets of data 
# from each realm to see if the patterns observed hold. 


# use dataset with realm info created above
head(predictsSites)

# sort dataset, remove NAs etc
predictsSites <- predictsSites[!is.na(predictsSites$UI2), ]
predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ] # 6457


levels(predictsSites$UI2)
table(predictsSites$UI2)

predictsSites$UI2 <- factor(predictsSites$UI2 , levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High","Urban"))

levels(predictsSites$UI2)
table(predictsSites$UI2)


predictsSites <- droplevels(predictsSites)


# SR subsets
predictsSites_sr_trop <- predictsSites[predictsSites$Tropical == "Tropical", ] # 2959
predictsSites_sr_temp <- predictsSites[predictsSites$Tropical == "Temperate", ] # 3498

length(unique(predictsSites_sr_trop$SS)) # 76
length(unique(predictsSites_sr_temp$SS)) # 63

table(predictsSites_sr_trop$UI2)

# Primary vegetation     Agriculture_High   Agriculture_Low  Secondary vegetation         Urban 
# 1135                  749                  525                  494                     56 

table(predictsSites_sr_temp$UI2)

# Primary vegetation     Agriculture_High   Agriculture_Low Secondary vegetation          Urban 
# 1377                  548                  474                  717                     382

# subset for abundance data
predictsSites_ab <- predictsSites[!is.na(predictsSites$LogAbund), ] # 5734

## dataset summaries ##

predictsSites_ab_trop <- predictsSites_ab[predictsSites_ab$Tropical == "Tropical", ] # 2805
predictsSites_ab_temp <- predictsSites_ab[predictsSites_ab$Tropical == "Temperate", ] # 2929
  

length(unique(predictsSites_ab_trop$SS)) # 62
length(unique(predictsSites_ab_temp$SS)) # 58

table(predictsSites_ab_trop$UI2)

# Primary vegetation Secondary vegetation   Agriculture_Low     Agriculture_High       Urban
# 1064                 458                  524                 744                    15

table(predictsSites_ab_temp$UI2)

# Primary vegetation Secondary vegetation  Agriculture_Low     Agriculture_High      Urban
# 1215                  575                330                 535                   274


##%######################################################%##
#                                                          #
####                    run models                      ####
#                                                          #
##%######################################################%##


# 1. Temperate, Abundance, mean anomaly

MeanAbundTemp <- GLMERSelect(modelData = predictsSites_ab_temp,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS"))

summary(MeanAbundTemp$model)

# LogAbund ~ UI2 + poly(StdTmeanAnomalyRS, 1) + UI2:poly(StdTmeanAnomalyRS, 1) + (1 | SS) + (1 | SSB)

# 2. Tropical, Abundance, mean anomaly
MeanAbundTrop <- GLMERSelect(modelData = predictsSites_ab_trop,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS"))

summary(MeanAbundTrop$model)
# LogAbund ~ UI2 + poly(StdTmeanAnomalyRS, 1) + UI2:poly(StdTmeanAnomalyRS, 1) + (1 | SS) + (1 | SSB)


# 3. Temperate, Richness, mean anomaly
MeanRichTemp <- GLMERSelect(modelData = predictsSites_sr_temp,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS"))

summary(MeanRichTemp$model)
#Species_richness ~ UI2 + poly(StdTmeanAnomalyRS, 1) + UI2:poly(StdTmeanAnomalyRS, 1) + (1 | SS) + (1 | SSB) + (1 | SSBS)


# 4. Tropical, Richness, mean anomaly
MeanRichTrop <- GLMERSelect(modelData = predictsSites_sr_trop,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS"))

summary(MeanRichTrop$model)
# Species_richness ~ UI2 + UI2:poly(StdTmeanAnomalyRS, 1) + poly(StdTmeanAnomalyRS,  1) + (1 | SS) + (1 | SSB) + (1 | SSBS)

save(MeanAbundTemp, file = paste0(outDir, "/MeanAbundTemp_output.rdata"))
save(MeanAbundTrop, file = paste0(outDir, "/MeanAbundTrop_output.rdata"))
save(MeanRichTemp, file = paste0(outDir, "/MeanRichTemp_output.rdata"))
save(MeanRichTrop, file = paste0(outDir, "/MeanRichTrop_output.rdata"))

#load(file = paste0(outDir, "/MeanAbundTemp_output.rdata"))
#load(file = paste0(outDir, "/MeanAbundTrop_output.rdata"))
#load(file = paste0(outDir, "/MeanRichTemp_output.rdata"))
#load(file = paste0(outDir, "/MeanRichTrop_output.rdata"))


#### save model output tables using sjPlot package ####

tab_model(MeanAbundTemp$model, MeanAbundTrop$model, transform = NULL, file = paste0(outDir, "/AbunMeanAnomTempTrop_output_table.html"))
summary(MeanAbundTemp$model)
R2GLMER(MeanAbundTemp$model)
summary(MeanAbundTrop$model)
R2GLMER(MeanAbundTrop$model)

tab_model(MeanRichTemp$model, MeanRichTrop$model, transform = NULL, file = paste0(outDir, "/RichMeanAnomTempTrop_output_table.html"))
summary(MeanRichTemp$model)
R2GLMER(MeanRichTemp$model) # use these values
# $conditional
# [1] 0.6652306
# $marginal
# [1] 0.00808889
summary(MeanRichTrop$model)
R2GLMER(MeanRichTrop$model) # use these values
# $conditional
# [1] 0.5596599
# $marginal
# [1] 0.03645811

# save model stats

# save the stats info
abtemp_stats <- as.data.frame(MeanAbundTemp$stats)
abtrop_stats <- as.data.frame(MeanAbundTrop$stats)
srtemp_stats <- as.data.frame(MeanRichTemp$stats)
srtrop_stats <- as.data.frame(MeanRichTrop$stats)

abtemp_stats$significant <- NA
abtrop_stats$significant <- NA
srtemp_stats$significant <- NA
srtrop_stats$significant <- NA

abtemp_stats$model <- "abtemp"
abtrop_stats$model <- "abtrop"
srtemp_stats$model <- "srtemp"
srtrop_stats$model <- "srtrop"

all_stats <- rbind(abtemp_stats, abtrop_stats, srtemp_stats, srtrop_stats)

# function to check significance
checksig <- function(x){
  if(x <= 0.05){ 
    res <- "Yes" 
  } else { 
    res <- "No" }
  return(res)}


# add values to table
all_stats$significant <- sapply(X = all_stats$P, FUN = checksig)



# save the stats tables
write.csv(all_stats, file = paste0(outDir, "/TropTemp_Stats.csv"), row.names = FALSE)



##%######################################################%##
#                                                          #
####                       Plots                        ####
#                                                          #
##%######################################################%##

exclQuantiles <- c(0.025,0.975)


### 1. MeanAbundTemp ###


nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAbundTemp$data$StdTmeanAnomalyRS),
                        to = max(MeanAbundTemp$data$StdTmeanAnomalyRS),
                        length.out = 200),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban"),
             levels = levels(MeanAbundTemp$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAbundTemp$data$StdTmeanAnomalyRS[
  MeanAbundTemp$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAbundTemp$data$StdTmeanAnomalyRS[
  MeanAbundTemp$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAbundTemp$data$StdTmeanAnomalyRS[
  MeanAbundTemp$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAbundTemp$data$StdTmeanAnomalyRS[
  MeanAbundTemp$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)
QUB <- quantile(x = MeanAbundTemp$data$StdTmeanAnomalyRS[
  MeanAbundTemp$data$UI2=="Urban"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAbundTemp$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban") & (nd$StdTmeanAnomalyRS < QUB[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban") & (nd$StdTmeanAnomalyRS > QUB[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
        geom_line(aes(col = UI2), size = 0.75) +
        geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
        geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
        scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
        scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
        theme_bw() + 
        ylab("Change in total abundance (%)") +
        xlab("Standardised Temperature Anomaly") +
        xlim(c(-0.5, 2)) +
        ylim(c(-100, 150)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        #legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
        ggtitle("a Non-tropical Realm")

p1

### 2. MeanAbundTrop ###


nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAbundTrop$data$StdTmeanAnomalyRS),
                        to = max(MeanAbundTrop$data$StdTmeanAnomalyRS),
                        length.out = 250),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban"),
             levels = levels(MeanAbundTrop$data$UI2)))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAbundTrop$data$StdTmeanAnomalyRS[
  MeanAbundTrop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAbundTrop$data$StdTmeanAnomalyRS[
  MeanAbundTrop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAbundTrop$data$StdTmeanAnomalyRS[
  MeanAbundTrop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAbundTrop$data$StdTmeanAnomalyRS[
  MeanAbundTrop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)
QUB <- quantile(x = MeanAbundTrop$data$StdTmeanAnomalyRS[
  MeanAbundTrop$data$UI2=="Urban"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAbundTrop$model,data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban") & (nd2$StdTmeanAnomalyRS < QUB[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban") & (nd2$StdTmeanAnomalyRS > QUB[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
  theme_bw() + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("b Tropical Realm")
p2

### 3. MeanRichTemp ###


nd3 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanRichTemp$data$StdTmeanAnomalyRS),
                        to = max(MeanRichTemp$data$StdTmeanAnomalyRS),
                        length.out = 200),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban"),
             levels = levels(MeanRichTemp$data$UI2)))

# back transform the predictors
nd3$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd3$LogAbund <- 0
nd3$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomaly==min(abs(nd3$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanRichTemp$data$StdTmeanAnomalyRS[
  MeanRichTemp$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanRichTemp$data$StdTmeanAnomalyRS[
  MeanRichTemp$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanRichTemp$data$StdTmeanAnomalyRS[
  MeanRichTemp$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanRichTemp$data$StdTmeanAnomalyRS[
  MeanRichTemp$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)
QUB <- quantile(x = MeanRichTemp$data$StdTmeanAnomalyRS[
  MeanRichTemp$data$UI2=="Urban"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanRichTemp$model,data = nd3)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmeanAnomalyRS > QAH[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban") & (nd3$StdTmeanAnomalyRS < QUB[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban") & (nd3$StdTmeanAnomalyRS > QUB[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


p3 <- ggplot(data = nd3, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
  theme_bw() + 
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  scale_y_continuous(breaks = c(-100, -50, 0, 50, 100, 150), limits = c(-100, 170)) +
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("c Non-tropical Realm")
p3


### 4. MeanRichTrop ###


nd4 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanRichTrop$data$StdTmeanAnomalyRS),
                        to = max(MeanRichTrop$data$StdTmeanAnomalyRS),
                        length.out = 250),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban"),
             levels = levels(MeanRichTrop$data$UI2)))

# back transform the predictors
nd4$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd4$LogAbund <- 0
nd4$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomaly==min(abs(nd4$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanRichTrop$data$StdTmeanAnomalyRS[
  MeanRichTrop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanRichTrop$data$StdTmeanAnomalyRS[
  MeanRichTrop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanRichTrop$data$StdTmeanAnomalyRS[
  MeanRichTrop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanRichTrop$data$StdTmeanAnomalyRS[
  MeanRichTrop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)
QUB <- quantile(x = MeanRichTrop$data$StdTmeanAnomalyRS[
  MeanRichTrop$data$UI2=="Urban"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanRichTrop$model,data = nd4)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmeanAnomalyRS > QAH[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Urban") & (nd4$StdTmeanAnomalyRS < QUB[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Urban") & (nd4$StdTmeanAnomalyRS > QUB[2])),] <- NA
# Get the median, upper and lower quants for the plot
nd4$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd4$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# p4 <- ggplot(data = nd4, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 0.75) +
#   geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00","pink")) +
#   theme_bw() + 
#   ylab("Change in Species Richness (%)") +
#   xlab("Standardised Temperature Anomaly") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 150)) + 
#   theme(aspect.ratio = 1, 
#         title = element_text(size = 8, face = "bold"),
#         axis.text = element_text(size = 7),
#         axis.title = element_text(size = 7),
#         legend.position = "none",
#         legend.text = element_text(size = 6), 
#         legend.title = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(size = 0.2),
#         panel.border = element_rect(size = 0.2), 
#         axis.ticks = element_line(size = 0.2)) +  
#   ggtitle("d Tropical Realm")
# p4

## ----------------------------------------------------------
## UI2 levels and colours (GLOBAL)
## ----------------------------------------------------------
ui2_levels <- c(
  "Primary vegetation",
  "Secondary vegetation",
  "Agriculture_Low",
  "Agriculture_High",
  "Urban"
)

ui2_cols <- c(
  "#009E73",  # Primary vegetation
  "#0072B2",  # Secondary vegetation
  "#E69F00",  # Agriculture_Low
  "#D55E00",  # Agriculture_High
  "pink"      # Urban
)

names(ui2_cols) <- ui2_levels

nd4$UI2 <- factor(nd4$UI2, levels = ui2_levels)

p4 <- ggplot(data = nd4, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  scale_fill_manual(values = ui2_cols, drop = FALSE) +
  scale_colour_manual(values = ui2_cols, drop = FALSE) +
  theme_bw() + 
  ylab("Change in Species Richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) + 
  theme(
    aspect.ratio = 1, 
    title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.position = "none",
    legend.text = element_text(size = 6), 
    legend.title = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2),
    panel.border = element_rect(size = 0.2), 
    axis.ticks = element_line(size = 0.2)
  ) +  
  ggtitle("d Tropical Realm")

## arrange Mean anomaly plots

leg <- get_legend(p4+ theme(legend.position = "bottom"))

cowplot::plot_grid(cowplot::plot_grid(p1 + theme(legend.position = "none"), 
          p2 + theme(legend.position = "none"), 
          p3+ theme(legend.position = "none"), 
          p4 + theme(legend.position = "none"),
          nrow = 2),
          leg, nrow= 2, rel_heights = c(5,1))


ggsave(filename = paste0(outDir, "Figure3_MeanAnom_TempTrop.pdf"), plot = last_plot(), width = 180, height = 170, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "Figure3_MeanAnom_TempTrop.png"), plot = last_plot(), width = 180, height = 170, units = "mm", dpi = 300)


##%######################################################%##
#                                                          #
####      Run separate models for the max anomaly       ####
#                                                          #
##%######################################################%##


# 1. Temperate, Abundance, max anomaly
MaxAbundTemp <- GLMERSelect(modelData = predictsSites_ab_temp,responseVar = "LogAbund",
                             fitFamily = "gaussian",fixedFactors = "UI2",
                             fixedTerms = list(StdTmaxAnomalyRS=1),
                             randomStruct = "(1|SS)+(1|SSB)",
                             fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                             saveVars = c("Species_richness", "Total_abundance", "SSBS"))

summary(MaxAbundTemp$model)
# LogAbund ~ UI2 + (1 | SS) + (1 | SSB)

# 2. Tropical, Abundance, max anomaly
MaxAbundTrop <- GLMERSelect(modelData = predictsSites_ab_trop,responseVar = "LogAbund",
                             fitFamily = "gaussian",fixedFactors = "UI2",
                             fixedTerms = list(StdTmaxAnomalyRS=1),
                             randomStruct = "(1|SS)+(1|SSB)",
                             fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                             saveVars = c("Species_richness", "Total_abundance", "SSBS"))

summary(MaxAbundTrop$model)
# LogAbund ~ UI2 + poly(StdTmaxAnomalyRS, 1) + UI2:poly(StdTmaxAnomalyRS,  1) + (1 | SS) + (1 |   SSB)


# 3. Temperate, Richness, max anomaly
MaxRichTemp <- GLMERSelect(modelData = predictsSites_sr_temp,responseVar = "Species_richness",
                            fitFamily = "poisson",fixedFactors = "UI2",
                            fixedTerms = list(StdTmaxAnomalyRS=1),
                            randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                            fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                            saveVars = c("Total_abundance", "SSBS"))

summary(MaxRichTemp$model)
# Species_richness ~ UI2 + UI2:poly(StdTmaxAnomalyRS, 1) + poly(StdTmaxAnomalyRS,      1) + (1 | SS) + (1 | SSB) + (1 | SSBS)

# 4. Tropical, Richness, max anomaly
MaxRichTrop <- GLMERSelect(modelData = predictsSites_sr_trop,responseVar = "Species_richness",
                            fitFamily = "poisson",fixedFactors = "UI2",
                            fixedTerms = list(StdTmaxAnomalyRS=1),
                            randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                            fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                            saveVars = c("Total_abundance", "SSBS"))

summary(MaxRichTrop$model)
# Species_richness ~ UI2 + UI2:poly(StdTmaxAnomalyRS, 1) + poly(StdTmaxAnomalyRS,1) + (1 | SS) + (1 | SSB) + (1 | SSBS)


save(MaxAbundTemp, file = paste0(outDir, "/MaxAbundTemp_output.rdata"))
save(MaxAbundTrop, file = paste0(outDir, "/MaxAbundTrop_output.rdata"))
save(MaxRichTemp, file = paste0(outDir, "/MaxRichTemp_output.rdata"))
save(MaxRichTrop, file = paste0(outDir, "/MaxRichTrop_output.rdata"))

# load(file = paste0(outDir, "/MaxAbundTemp_output.rdata"))
# load(file = paste0(outDir, "/MaxAbundTrop_output.rdata"))
# load(file = paste0(outDir, "/MaxRichTemp_output.rdata"))
# load(file = paste0(outDir, "/MaxRichTrop_output.rdata"))


#### save model output tables using sjPlot package ####

sjPlot::tab_model(MaxAbundTemp$model, MaxAbundTrop$model, transform = NULL, file = paste0(outDir, "/AbunMaxAnomTempTrop_output_table.html"))
summary(MaxAbundTemp$model)
R2GLMER(MaxAbundTemp$model)
summary(MaxAbundTrop$model)
R2GLMER(MaxAbundTrop$model)

tab_model(MaxRichTemp$model, MaxRichTrop$model, transform = NULL, file = paste0(outDir, "/RichMaxAnomTempTrop_output_table.html"))
summary(MaxRichTemp$model)
R2GLMER(MaxRichTemp$model) # use these values
# $conditional
#[1] 0.6597703
#$marginal
#[1] 0.004763959
summary(MaxRichTrop$model)
R2GLMER(MaxRichTrop$model) # use these values
# $conditional
# [1] 0.552677
# $marginal
# [1] 0.02567178

# save model stats

# save the stats info
abtemp_stats <- as.data.frame(MaxAbundTemp$stats)
abtrop_stats <- as.data.frame(MaxAbundTrop$stats)
srtemp_stats <- as.data.frame(MaxRichTemp$stats)
srtrop_stats <- as.data.frame(MaxRichTrop$stats)

abtemp_stats$significant <- NA
abtrop_stats$significant <- NA
srtemp_stats$significant <- NA
srtrop_stats$significant <- NA

abtemp_stats$model <- "abtemp"
abtrop_stats$model <- "abtrop"
srtemp_stats$model <- "srtemp"
srtrop_stats$model <- "srtrop"

all_stats <- rbind(abtemp_stats, abtrop_stats, srtemp_stats, srtrop_stats)

# function to check significance
checksig <- function(x){
  if(x <= 0.05){ 
    res <- "Yes" 
  } else { 
    res <- "No" }
  return(res)}


# add values to table
all_stats$significant <- sapply(X = all_stats$P, FUN = checksig)



# save the stats tables
write.csv(all_stats, file = paste0(outDir, "/TropTemp_Max_Stats.csv"), row.names = FALSE)




##%######################################################%##
#                                                          #
####                  Max anom Plots                    ####
#                                                          #
##%######################################################%##

exclQuantiles <- c(0.025,0.975)


### 1. MaxAbundTemp ###


nd <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAbundTemp$data$StdTmaxAnomalyRS),
                        to = max(MaxAbundTemp$data$StdTmaxAnomalyRS),
                        length.out = 250),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban"),
             levels = levels(MaxAbundTemp$data$UI2)))

# back transform the predictors
nd$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmaxAnomaly==min(abs(nd$StdTmaxAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MaxAbundTemp$data$StdTmaxAnomalyRS[
  MaxAbundTemp$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAbundTemp$data$StdTmaxAnomalyRS[
  MaxAbundTemp$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAbundTemp$data$StdTmaxAnomalyRS[
  MaxAbundTemp$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAbundTemp$data$StdTmaxAnomalyRS[
  MaxAbundTemp$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)
QUB <- quantile(x = MaxAbundTemp$data$StdTmaxAnomalyRS[
  MaxAbundTemp$data$UI2=="Urban"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MaxAbundTemp$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmaxAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmaxAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmaxAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmaxAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmaxAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmaxAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmaxAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmaxAnomalyRS > QAH[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban") & (nd$StdTmaxAnomalyRS < QUB[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban") & (nd$StdTmaxAnomalyRS > QUB[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

#"#009E73" - green
#"#0072B2" - blue
#"#E69F00" - yellow
#"#D55E00" - red

p1 <- ggplot(data = nd, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
  theme_bw() + 
  ylab("Change in Total Abundance (%)") +
  xlab("Maximum Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 100)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank()) + 
  ggtitle("a Non-tropical Realm")
p1


### 2. MaxAbundTrop ###


nd2 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAbundTrop$data$StdTmaxAnomalyRS),
                        to = max(MaxAbundTrop$data$StdTmaxAnomalyRS),
                        length.out = 250),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban"),
             levels = levels(MaxAbundTrop$data$UI2)))

# back transform the predictors
nd2$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomaly==min(abs(nd2$StdTmaxAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MaxAbundTrop$data$StdTmaxAnomalyRS[
  MaxAbundTrop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAbundTrop$data$StdTmaxAnomalyRS[
  MaxAbundTrop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAbundTrop$data$StdTmaxAnomalyRS[
  MaxAbundTrop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAbundTrop$data$StdTmaxAnomalyRS[
  MaxAbundTrop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)
QUB <- quantile(x = MaxAbundTrop$data$StdTmaxAnomalyRS[
  MaxAbundTrop$data$UI2=="Urban"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MaxAbundTrop$model,data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS > QAH[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban") & (nd2$StdTmaxAnomalyRS < QUB[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban") & (nd2$StdTmaxAnomalyRS > QUB[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


p2 <- ggplot(data = nd2, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
  theme_bw() + 
  ylab("Change in Total Abundance (%)") +
  xlab("Maximum Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 100)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank()) + 
  ggtitle("b Tropical Realm")

p2



### 3. MaxRichTemp ###


nd3 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxRichTemp$data$StdTmaxAnomalyRS),
                        to = max(MaxRichTemp$data$StdTmaxAnomalyRS),
                        length.out = 250),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban"),
             levels = levels(MaxRichTemp$data$UI2)))

# back transform the predictors
nd3$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd3$LogAbund <- 0
nd3$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomaly==min(abs(nd3$StdTmaxAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MaxRichTemp$data$StdTmaxAnomalyRS[
  MaxRichTemp$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxRichTemp$data$StdTmaxAnomalyRS[
  MaxRichTemp$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxRichTemp$data$StdTmaxAnomalyRS[
  MaxRichTemp$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxRichTemp$data$StdTmaxAnomalyRS[
  MaxRichTemp$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)
QUB <- quantile(x = MaxRichTemp$data$StdTmaxAnomalyRS[
  MaxRichTemp$data$UI2=="Urban"],
  probs = exclQuantiles)


# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MaxRichTemp$model,data = nd3)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS > QAH[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban") & (nd3$StdTmaxAnomalyRS < QUB[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban") & (nd3$StdTmaxAnomalyRS > QUB[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


p3 <- ggplot(data = nd3, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00","pink")) +
  theme_bw() + 
  ylab("Change in Species Richness (%)") +
  xlab("Maximum Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 100)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank()) + 
  ggtitle("c Non-tropical Realm")
p3

### 4. MaxRichTrop ###


nd4 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxRichTrop$data$StdTmaxAnomalyRS),
                        to = max(MaxRichTrop$data$StdTmaxAnomalyRS),
                        length.out = 250),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban"),
             levels = levels(MaxRichTrop$data$UI2)))

# back transform the predictors
nd4$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd4$LogAbund <- 0
nd4$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomaly==min(abs(nd4$StdTmaxAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MaxRichTrop$data$StdTmaxAnomalyRS[
  MaxRichTrop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxRichTrop$data$StdTmaxAnomalyRS[
  MaxRichTrop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxRichTrop$data$StdTmaxAnomalyRS[
  MaxRichTrop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxRichTrop$data$StdTmaxAnomalyRS[
  MaxRichTrop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)
QUB <- quantile(x = MaxRichTrop$data$StdTmaxAnomalyRS[
  MaxRichTrop$data$UI2=="Urban"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MaxRichTrop$model,data = nd4)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS > QAH[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Urban") & (nd4$StdTmaxAnomalyRS < QUB[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Urban") & (nd4$StdTmaxAnomalyRS > QUB[2])),] <- NA


# Get the median, upper and lower quants for the plot
nd4$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd4$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


ui2_levels <- c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban")
ui2_cols   <- c("#009E73", "#0072B2", "#E69F00", "#D55E00", "pink")
names(ui2_cols) <- ui2_levels   # 关键：命名！

# 确保 nd4$UI2 的 levels 正确
nd4$UI2 <- factor(nd4$UI2, levels = ui2_levels)
p4 <- ggplot(data = nd4, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_fill_manual(values = ui2_cols, drop = FALSE) +
  scale_colour_manual(values = ui2_cols, drop = FALSE) +
  theme_bw() + 
  ylab("Change in Species Richness (%)") +
  xlab("Maximum Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 100)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank()) + 
  ggtitle("d Tropical Realm")


p4

## arrange Mean anomaly plots

leg <- get_legend(p4+ theme(legend.position = "bottom"))

cowplot::plot_grid(cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                    p2 + theme(legend.position = "none"), 
                    p3+ theme(legend.position = "none"), 
                    p4 + theme(legend.position = "none"),
                    nrow = 2),
          leg, nrow= 2, rel_heights = c(5,1))


ggsave(file = paste0(outDir, "/Extended_Data5_MaxAnom_TempTrop.pdf"), width = 8, height = 8.5)
ggsave(filename = paste0(outDir, "Extended_Data5_MaxAnom_TempTrop.jpeg"), plot = last_plot(), width = 183, height = 200, units = "mm", dpi = 300)





##%######################################################%##
#                                                          #
####   brms realm models + diagnostics + plots (FULL)    ####
#                                                          #
##%######################################################%##

## =========================================================
## 0) UI2 order & colors / UI2顺序与配色（你指定）
## =========================================================
ui2_levels <- c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban")
ui2_cols   <- c("#009E73", "#0072B2", "#E69F00", "#D55E00", "pink")
names(ui2_cols) <- ui2_levels  # 关键：颜色与水平一一对应（命名）

## Ensure factor levels / 确保所有数据子集的UI2水平一致
predictsSites$UI2 <- factor(predictsSites$UI2, levels = ui2_levels)
predictsSites_ab_temp$UI2 <- factor(predictsSites_ab_temp$UI2, levels = ui2_levels)
predictsSites_ab_trop$UI2 <- factor(predictsSites_ab_trop$UI2, levels = ui2_levels)
predictsSites_sr_temp$UI2 <- factor(predictsSites_sr_temp$UI2, levels = ui2_levels)
predictsSites_sr_trop$UI2 <- factor(predictsSites_sr_trop$UI2, levels = ui2_levels)

## =========================================================
## 1) Libraries / 加载包
## =========================================================
library(brms)
library(ggplot2)
library(cowplot)
library(cmdstanr)

# 如果还没装过 CmdStan（只需一次）
# cmdstanr::install_cmdstan()

cmdstanr::cmdstan_version()

## =========================================================
## 2) Output dirs / 输出目录（沿用 outDir）
## =========================================================
brmsDir     <- file.path(outDir, "brms_models")
brmsFigDir  <- file.path(outDir, "brms_figures")
brmsDiagDir <- file.path(outDir, "brms_diagnostics")
brmsPredDir <- file.path(outDir, "brms_predictions")

dir.create(brmsDir,     showWarnings = FALSE, recursive = TRUE)
dir.create(brmsFigDir,  showWarnings = FALSE, recursive = TRUE)
dir.create(brmsDiagDir, showWarnings = FALSE, recursive = TRUE)
dir.create(brmsPredDir, showWarnings = FALSE, recursive = TRUE)

## =========================================================
## 3) Helper: save plot in 3 formats / 三格式保存
## =========================================================
save_plot_3formats <- function(p, out_prefix, width_mm = 180, height_mm = 170, dpi = 320){
  ggplot2::ggsave(paste0(out_prefix, ".png"), plot = p,
                  width = width_mm, height = height_mm, units = "mm", dpi = dpi)
  ggplot2::ggsave(paste0(out_prefix, ".pdf"), plot = p,
                  width = width_mm, height = height_mm, units = "mm")
  ggplot2::ggsave(paste0(out_prefix, ".svg"), plot = p,
                  width = width_mm, height = height_mm, units = "mm")
}

## =========================================================
## 4) brms control / 采样控制参数
## =========================================================
brms_ctrl <- list(adapt_delta = 0.99, max_treedepth = 15)

## =========================================================
## 5) Priors / 先验
## =========================================================
## Gaussian for LogAbund
priors_gaussian <- c(
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

## Poisson for richness
priors_poisson <- c(
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

## NegBinomial for richness (adds shape parameter)
priors_negbin <- c(
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "shape")
)

## =========================================================
## 6) Helper: overdispersion check for counts / 过度离散判断
## =========================================================
## 用一个 Poisson GLM（固定效应同款）计算 dispersion:
## dispersion = sum(Pearson residual^2) / df.residual
calc_poisson_dispersion <- function(data, x_rs){
  f_glm <- stats::as.formula(paste0("Species_richness ~ UI2 * poly(", x_rs, ", 1)"))
  m <- stats::glm(f_glm, data = data, family = stats::poisson())
  rp <- stats::residuals(m, type = "pearson")
  disp <- sum(rp^2, na.rm = TRUE) / stats::df.residual(m)
  return(disp)
}

choose_count_family <- function(data, x_rs, threshold = 1.2){
  disp <- calc_poisson_dispersion(data, x_rs)
  if(is.na(disp)) disp <- Inf
  if(disp > threshold){
    message("Overdispersion detected (dispersion = ", round(disp, 3), ") -> use NegBinomial")
    return(list(family = negbinomial(), priors = priors_negbin, dispersion = disp, chosen = "negbinomial"))
  } else {
    message("No strong overdispersion (dispersion = ", round(disp, 3), ") -> use Poisson")
    return(list(family = poisson(), priors = priors_poisson, dispersion = disp, chosen = "poisson"))
  }
}

## =========================================================
## 7) Helper: fit brms safely / 稳健拟合
## =========================================================

fit_brms_safe <- function(formula, data, family, priors, name,
                          iter = 8000, warmup = 4000, chains = 4, seed = 123,
                          backend = "cmdstanr",
                          control = list(adapt_delta = 0.99, max_treedepth = 15),
                          cores = NULL,
                          modelDir = NULL){
  
  message("Fitting brms: ", name)
  
  if (is.null(cores)) cores <- min(chains, parallel::detectCores())
  
  mod <- brms::brm(
    formula   = formula,
    data      = data,
    family    = family,
    prior     = priors,
    chains    = chains,
    iter      = iter,
    warmup    = warmup,
    cores     = cores,
    seed      = seed,
    control   = control,
    backend   = backend,
    save_pars = brms::save_pars(all = TRUE)
  )
  
  if(!is.null(modelDir)){
    saveRDS(mod, file = file.path(modelDir, paste0(name, ".rds")))
  }
  
  return(mod)
}

## =========================================================
## 8) Helper: diagnostics / 模型检验与保存
## =========================================================
save_brms_diagnostics <- function(fit, name, outDir = "brms_diagnostics"){
  
  if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
  
  ## 1. Summary（参数估计 + Rhat + ESS）
  sum_fit <- summary(fit)
  capture.output(
    sum_fit,
    file = file.path(outDir, paste0(name, "_summary.txt"))
  )
  
  ## 2. Rhat / ESS 检查（论文常用）
  rhats <- brms::rhat(fit)
  ess   <- brms::neff_ratio(fit)
  
  write.csv(
    data.frame(parameter = names(rhats), Rhat = rhats),
    file = file.path(outDir, paste0(name, "_Rhat.csv")),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(parameter = names(ess), ESS_ratio = ess),
    file = file.path(outDir, paste0(name, "_ESS_ratio.csv")),
    row.names = FALSE
  )
  
  ## 3. Trace / Density plots
  pdf(file.path(outDir, paste0(name, "_trace.pdf")), width = 8, height = 6)
  plot(fit)
  dev.off()
  
  ## 4. Posterior predictive check
  pdf(file.path(outDir, paste0(name, "_pp_check.pdf")), width = 7, height = 5)
  brms::pp_check(fit, ndraws = 100)
  dev.off()
  
  ## 5. Divergent transitions（cmdstanr 特别重要）
  sampler_diag <- brms::nuts_params(fit)
  n_div <- sum(sampler_diag$Parameter == "divergent__" &
                 sampler_diag$Value == 1)
  
  writeLines(
    paste("Number of divergent transitions:", n_div),
    con = file.path(outDir, paste0(name, "_divergent.txt"))
  )
  
  invisible(list(
    summary = sum_fit,
    n_divergent = n_div
  ))
}


## =========================================================
## 9) Prediction + plot (% change) with exclQuantiles
##    预测作图：按UI2分位数截断（超界Pred置NA）
## =========================================================
make_brms_pred_plot_percent <- function(mod, data,
                                        x_rs, x_bt,
                                        xlab, ylab, title_text,
                                        exclQuantiles = c(0.025, 0.975),
                                        xlim_use = c(-0.5, 2),
                                        ylim_use = c(-100, 150),
                                        breaks_y = seq(-100, 150, 50),
                                        ndraws = 2000,
                                        out_prefix){
  
  ## ---- prediction grid / 构建预测网格（RS尺度） ----
  nd <- expand.grid(
    x = seq(from = min(data[[x_rs]], na.rm = TRUE),
            to   = max(data[[x_rs]], na.rm = TRUE),
            length.out = 220),
    UI2 = factor(ui2_levels, levels = ui2_levels)
  )
  names(nd)[names(nd)=="x"] <- x_rs
  
  ## ---- back-transform for x-axis / 反变换用于横轴（原尺度） ----
  nd[[x_bt]] <- BackTransformCentreredPredictor(
    transformedX = nd[[x_rs]],
    originalX    = predictsSites[[x_bt]]
  )
  
  ## ---- reference: PV closest to 0 on back-transformed scale ----
  refRow <- which((nd$UI2 == "Primary vegetation") &
                    (nd[[x_bt]] == min(abs(nd[[x_bt]]), na.rm = TRUE)))
  
  ## ---- posterior expected predictions / 后验均值预测（population-level） ----
  ep <- brms::posterior_epred(mod, newdata = nd, ndraws = ndraws, re_formula = NA)
  ep_t <- t(ep)  # nrow(nd) x ndraws
  
  ## ---- relative to PV@~0 / 相对PV参考 ----
  rel <- sweep(ep_t, 2, ep_t[refRow, ], FUN = "/")
  
  ## ---- summarize / 汇总中位数与95%CI ----
  nd$PredMedian <- (apply(rel, 1, median, na.rm = TRUE) * 100) - 100
  nd$PredLower  <- (apply(rel, 1, quantile, probs = 0.025, na.rm = TRUE) * 100) - 100
  nd$PredUpper  <- (apply(rel, 1, quantile, probs = 0.975, na.rm = TRUE) * 100) - 100
  
  ## ---- NEW: quantile exclusion by UI2 / 新增：按UI2分位数截断 ----
  ## 对每个UI2，用“该模型数据”的x_rs分布计算分位数
  q_tbl <- lapply(ui2_levels, function(u){
    xx <- data[[x_rs]][data$UI2 == u]
    xx <- xx[is.finite(xx)]
    if(length(xx) < 5) return(data.frame(UI2 = u, q_lo = -Inf, q_hi = Inf))
    qq <- stats::quantile(xx, probs = exclQuantiles, na.rm = TRUE)
    data.frame(UI2 = u, q_lo = as.numeric(qq[1]), q_hi = as.numeric(qq[2]))
  })
  q_tbl <- do.call(rbind, q_tbl)
  
  nd <- merge(nd, q_tbl, by = "UI2", all.x = TRUE, sort = FALSE)
  
  ## 若 nd 的 x_rs 超出 [q_lo, q_hi] 则将Pred置NA（线段会断开）
  out_of_range <- (nd[[x_rs]] < nd$q_lo) | (nd[[x_rs]] > nd$q_hi)
  nd$PredMedian[out_of_range] <- NA
  nd$PredLower[out_of_range]  <- NA
  nd$PredUpper[out_of_range]  <- NA
  
  ## ---- save prediction data / 保存预测数据 ----
  saveRDS(nd, file = file.path(brmsPredDir, paste0(out_prefix, "_nd.rds")))
  
  ## ---- plot / 作图（严格按ui2_levels与ui2_cols） ----
  p <- ggplot(nd, aes(x = .data[[x_bt]], y = PredMedian)) +
    geom_line(aes(col = UI2), size = 0.75, na.rm = TRUE) +
    geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2),
                alpha = 0.2, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
    scale_colour_manual(values = ui2_cols, breaks = ui2_levels, drop = FALSE) +
    scale_fill_manual(values = ui2_cols, breaks = ui2_levels, drop = FALSE) +
    theme_bw() +
    ylab(ylab) +
    xlab(xlab) +
    coord_cartesian(xlim = xlim_use, ylim = ylim_use) +
    scale_y_continuous(breaks = breaks_y) +
    theme(
      aspect.ratio = 1,
      title = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.2),
      panel.border = element_rect(size = 0.2),
      axis.ticks = element_line(size = 0.2)
    ) +
    ggtitle(title_text)
  
  save_plot_3formats(p, file.path(brmsFigDir, out_prefix),
                     width_mm = 180, height_mm = 170)
  
  return(list(plot = p, nd = nd, q_tbl = q_tbl))
}

##%######################################################%##
#                                                          #
####               (A) Mean anomaly models               ####
#                                                          #
##%######################################################%##

## ---- Explicit settings per model / 每个模型显式参数 ----
exclQ_mean <- c(0.025, 0.975)
exclQ_mean <- c(0.01, 0.99)
ndraws_mean <- 2000

## A1. Temperate Abundance (Gaussian)
A1_name   <- "brms_MeanAbundTemp"
A1_iter   <- 6000
A1_warmup <- 3000
A1_chains <- 4
A1_seed   <- 101

b_MeanAbundTemp <- fit_brms_safe(
  formula = bf(LogAbund ~ UI2 * poly(StdTmeanAnomalyRS, 1) + (1|SS) + (1|SSB)),
  data    = predictsSites_ab_temp,
  family  = gaussian(),
  priors  = priors_gaussian,
  backend = "cmdstanr",
  name    = A1_name,
  iter    = A1_iter, warmup = A1_warmup, chains = A1_chains, seed = A1_seed
)
save_brms_diagnostics(b_MeanAbundTemp, A1_name)

## A2. Tropical Abundance (Gaussian)
A2_name   <- "brms_MeanAbundTrop"
A2_iter   <- 6000
A2_warmup <- 3000
A2_chains <- 4
A2_seed   <- 102

b_MeanAbundTrop <- fit_brms_safe(
  formula = bf(LogAbund ~ UI2 * poly(StdTmeanAnomalyRS, 1) + (1|SS) + (1|SSB)),
  data    = predictsSites_ab_trop,
  family  = gaussian(),
  priors  = priors_gaussian,
  backend = "cmdstanr",
  name    = A2_name,
  iter    = A2_iter, warmup = A2_warmup, chains = A2_chains, seed = A2_seed
)
save_brms_diagnostics(b_MeanAbundTrop, A2_name)

## A3. Temperate Richness (Poisson OR NegBin) — automatic check
A3_name   <- "brms_MeanRichTemp"
A3_iter   <- 8000
A3_warmup <- 4000
A3_chains <- 4
A3_seed   <- 201
A3_family_choice <- choose_count_family(predictsSites_sr_temp, x_rs = "StdTmeanAnomalyRS", threshold = 1.2)
saveRDS(A3_family_choice, file = file.path(brmsDiagDir, paste0(A3_name, "_family_choice.rds")))

b_MeanRichTemp <- fit_brms_safe(
  formula = bf(Species_richness ~ UI2 * poly(StdTmeanAnomalyRS, 1) + (1|SS) + (1|SSB) + (1|SSBS)),
  data    = predictsSites_sr_temp,
  family  = A3_family_choice$family,
  priors  = A3_family_choice$priors,
  backend = "cmdstanr",
  name    = A3_name,
  iter    = A3_iter, warmup = A3_warmup, chains = A3_chains, seed = A3_seed
)
save_brms_diagnostics(b_MeanRichTemp, A3_name)

## A4. Tropical Richness (Poisson OR NegBin) — automatic check
A4_name   <- "brms_MeanRichTrop"
A4_iter   <- 7000
A4_warmup <- 3500
A4_chains <- 4
A4_seed   <- 202
A4_family_choice <- choose_count_family(predictsSites_sr_trop, x_rs = "StdTmeanAnomalyRS", threshold = 1.2)
saveRDS(A4_family_choice, file = file.path(brmsDiagDir, paste0(A4_name, "_family_choice.rds")))

b_MeanRichTrop <- fit_brms_safe(
  formula = bf(Species_richness ~ UI2 * poly(StdTmeanAnomalyRS, 1) + (1|SS) + (1|SSB) + (1|SSBS)),
  data    = predictsSites_sr_trop,
  family  = A4_family_choice$family,
  priors  = A4_family_choice$priors,
  backend = "cmdstanr",
  name    = A4_name,
  iter    = A4_iter, warmup = A4_warmup, chains = A4_chains, seed = A4_seed
)
save_brms_diagnostics(b_MeanRichTrop, A4_name)

##%######################################################%##
#                                                          #
####               (B) Max anomaly models                ####
#                                                          #
##%######################################################%##
exclQ_max <- c(0.025, 0.975)
exclQ_max <- c(0.01, 0.99)
ndraws_max <- 2000

## B1. Temperate Abundance (Gaussian)
B1_name   <- "brms_MaxAbundTemp"
B1_iter   <- 6000
B1_warmup <- 3000
B1_chains <- 4
B1_seed   <- 103

b_MaxAbundTemp <- fit_brms_safe(
  formula = bf(LogAbund ~ UI2 * poly(StdTmaxAnomalyRS, 1) + (1|SS) + (1|SSB)),
  data    = predictsSites_ab_temp,
  family  = gaussian(),
  priors  = priors_gaussian,
  backend = "cmdstanr",
  name    = B1_name,
  iter    = B1_iter, warmup = B1_warmup, chains = B1_chains, seed = B1_seed
)
save_brms_diagnostics(b_MaxAbundTemp, B1_name)

## B2. Tropical Abundance (Gaussian)
B2_name   <- "brms_MaxAbundTrop"
B2_iter   <- 6000
B2_warmup <- 3000
B2_chains <- 4
B2_seed   <- 104

b_MaxAbundTrop <- fit_brms_safe(
  formula = bf(LogAbund ~ UI2 * poly(StdTmaxAnomalyRS, 1) + (1|SS) + (1|SSB)),
  data    = predictsSites_ab_trop,
  family  = gaussian(),
  priors  = priors_gaussian,
  backend = "cmdstanr",
  name    = B2_name,
  iter    = B2_iter, warmup = B2_warmup, chains = B2_chains, seed = B2_seed
)
save_brms_diagnostics(b_MaxAbundTrop, B2_name)

## B3. Temperate Richness (Poisson OR NegBin) — automatic check
B3_name   <- "brms_MaxRichTemp"
B3_iter   <- 7000
B3_warmup <- 3500
B3_chains <- 4
B3_seed   <- 203
B3_family_choice <- choose_count_family(predictsSites_sr_temp, x_rs = "StdTmaxAnomalyRS", threshold = 1.2)
saveRDS(B3_family_choice, file = file.path(brmsDiagDir, paste0(B3_name, "_family_choice.rds")))

b_MaxRichTemp <- fit_brms_safe(
  formula = bf(Species_richness ~ UI2 * poly(StdTmaxAnomalyRS, 1) + (1|SS) + (1|SSB) + (1|SSBS)),
  data    = predictsSites_sr_temp,
  family  = B3_family_choice$family,
  priors  = B3_family_choice$priors,
  backend = "cmdstanr",
  name    = B3_name,
  iter    = B3_iter, warmup = B3_warmup, chains = B3_chains, seed = B3_seed
)
save_brms_diagnostics(b_MaxRichTemp, B3_name)

## B4. Tropical Richness (Poisson OR NegBin) — automatic check
B4_name   <- "brms_MaxRichTrop"
B4_iter   <- 7000
B4_warmup <- 3500
B4_chains <- 4
B4_seed   <- 204
B4_family_choice <- choose_count_family(predictsSites_sr_trop, x_rs = "StdTmaxAnomalyRS", threshold = 1.2)
saveRDS(B4_family_choice, file = file.path(brmsDiagDir, paste0(B4_name, "_family_choice.rds")))

b_MaxRichTrop <- fit_brms_safe(
  formula = bf(Species_richness ~ UI2 * poly(StdTmaxAnomalyRS, 1) + (1|SS) + (1|SSB) + (1|SSBS)),
  data    = predictsSites_sr_trop,
  family  = B4_family_choice$family,
  priors  = B4_family_choice$priors,
  backend = "cmdstanr",
  name    = B4_name,
  iter    = B4_iter, warmup = B4_warmup, chains = B4_chains, seed = B4_seed
)
save_brms_diagnostics(b_MaxRichTrop, B4_name)

##%######################################################%##
#                                                          #
####            (C) Predictions + 4-panel plots          ####
#                                                          #
##%######################################################%##

## ---- Mean anomaly plots (4 panels) ----
m1 <- make_brms_pred_plot_percent(
  mod = b_MeanAbundTemp,
  data = predictsSites_ab_temp,
  x_rs = "StdTmeanAnomalyRS",
  x_bt = "StdTmeanAnomaly",
  xlab = "Standardised Temperature Anomaly",
  ylab = "Change in total abundance (%)",
  title_text = "a Non-tropical Realm (brms)",
  exclQuantiles = exclQ_mean,
  xlim_use = c(-0.5, 2),
  ylim_use = c(-100, 150),
  breaks_y = seq(-100, 150, 50),
  ndraws = ndraws_mean,
  out_prefix = "MeanAnom_Abund_Temperate_brms"
)

m2 <- make_brms_pred_plot_percent(
  mod = b_MeanAbundTrop,
  data = predictsSites_ab_trop,
  x_rs = "StdTmeanAnomalyRS",
  x_bt = "StdTmeanAnomaly",
  xlab = "Standardised Temperature Anomaly",
  ylab = "Change in total abundance (%)",
  title_text = "b Tropical Realm (brms)",
  exclQuantiles = exclQ_mean,
  xlim_use = c(-0.5, 2),
  ylim_use = c(-100, 150),
  breaks_y = seq(-100, 150, 50),
  ndraws = ndraws_mean,
  out_prefix = "MeanAnom_Abund_Tropical_brms"
)

m3 <- make_brms_pred_plot_percent(
  mod = b_MeanRichTemp,
  data = predictsSites_sr_temp,
  x_rs = "StdTmeanAnomalyRS",
  x_bt = "StdTmeanAnomaly",
  xlab = "Standardised Temperature Anomaly",
  ylab = "Change in species richness (%)",
  title_text = "c Non-tropical Realm (brms)",
  exclQuantiles = exclQ_mean,
  xlim_use = c(-0.5, 2),
  #coord_cartesian(xlim_use = c(-0.5, 2))
  ylim_use = c(-100, 170),
  breaks_y = c(-100, -50, 0, 50, 100, 150),
  ndraws = ndraws_mean,
  out_prefix = "MeanAnom_Rich_Temperate_brms"
)

m4 <- make_brms_pred_plot_percent(
  mod = b_MeanRichTrop,
  data = predictsSites_sr_trop,
  x_rs = "StdTmeanAnomalyRS",
  x_bt = "StdTmeanAnomaly",
  xlab = "Standardised Temperature Anomaly",
  ylab = "Change in species richness (%)",
  title_text = "d Tropical Realm (brms)",
  exclQuantiles = exclQ_mean,
  xlim_use = c(-0.5, 2),
  ylim_use = c(-100, 150),
  breaks_y = seq(-100, 150, 50),
  ndraws = ndraws_mean,
  out_prefix = "MeanAnom_Rich_Tropical_brms"
)

leg_mean <- cowplot::get_legend(m4$plot + theme(legend.position = "bottom"))

p_mean_brms_4 <- cowplot::plot_grid(
  cowplot::plot_grid(m1$plot + theme(legend.position = "none"),
                     m2$plot + theme(legend.position = "none"),
                     m3$plot + theme(legend.position = "none"),
                     m4$plot + theme(legend.position = "none"),
                     nrow = 2),
  leg_mean, nrow = 2, rel_heights = c(5, 1)
)

save_plot_3formats(p_mean_brms_4,
                   file.path(brmsFigDir, "Figure3_MeanAnom_TempTrop_brms_4panel"),
                   width_mm = 180, height_mm = 170)

## ---- Max anomaly plots (4 panels) ----
x1 <- make_brms_pred_plot_percent(
  mod = b_MaxAbundTemp,
  data = predictsSites_ab_temp,
  x_rs = "StdTmaxAnomalyRS",
  x_bt = "StdTmaxAnomaly",
  xlab = "Maximum Temperature Anomaly",
  ylab = "Change in Total Abundance (%)",
  title_text = "a Non-tropical Realm (brms)",
  exclQuantiles = exclQ_max,
  xlim_use = c(-0.5, 2),
  ylim_use = c(-100, 100),
  breaks_y = seq(-100, 100, 50),
  ndraws = ndraws_max,
  out_prefix = "MaxAnom_Abund_Temperate_brms"
)

x2 <- make_brms_pred_plot_percent(
  mod = b_MaxAbundTrop,
  data = predictsSites_ab_trop,
  x_rs = "StdTmaxAnomalyRS",
  x_bt = "StdTmaxAnomaly",
  xlab = "Maximum Temperature Anomaly",
  ylab = "Change in Total Abundance (%)",
  title_text = "b Tropical Realm (brms)",
  exclQuantiles = exclQ_max,
  xlim_use = c(-0.5, 2),
  ylim_use = c(-100, 100),
  breaks_y = seq(-100, 100, 50),
  ndraws = ndraws_max,
  out_prefix = "MaxAnom_Abund_Tropical_brms"
)

x3 <- make_brms_pred_plot_percent(
  mod = b_MaxRichTemp,
  data = predictsSites_sr_temp,
  x_rs = "StdTmaxAnomalyRS",
  x_bt = "StdTmaxAnomaly",
  xlab = "Maximum Temperature Anomaly",
  ylab = "Change in Species Richness (%)",
  title_text = "c Non-tropical Realm (brms)",
  exclQuantiles = exclQ_max,
  xlim_use = c(-0.5, 2),
  ylim_use = c(-100, 100),
  breaks_y = seq(-100, 100, 50),
  ndraws = ndraws_max,
  out_prefix = "MaxAnom_Rich_Temperate_brms"
)

x4 <- make_brms_pred_plot_percent(
  mod = b_MaxRichTrop,
  data = predictsSites_sr_trop,
  x_rs = "StdTmaxAnomalyRS",
  x_bt = "StdTmaxAnomaly",
  xlab = "Maximum Temperature Anomaly",
  ylab = "Change in Species Richness (%)",
  title_text = "d Tropical Realm (brms)",
  exclQuantiles = exclQ_max,
  xlim_use = c(-0.5, 2),
  ylim_use = c(-100, 100),
  breaks_y = seq(-100, 100, 50),
  ndraws = ndraws_max,
  out_prefix = "MaxAnom_Rich_Tropical_brms"
)

leg_max <- cowplot::get_legend(x4$plot + theme(legend.position = "bottom"))

p_max_brms_4 <- cowplot::plot_grid(
  cowplot::plot_grid(x1$plot + theme(legend.position = "none"),
                     x2$plot + theme(legend.position = "none"),
                     x3$plot + theme(legend.position = "none"),
                     x4$plot + theme(legend.position = "none"),
                     nrow = 2),
  leg_max, nrow = 2, rel_heights = c(5, 1)
)

save_plot_3formats(p_max_brms_4,
                   file.path(brmsFigDir, "Extended_Data5_MaxAnom_TempTrop_brms_4panel"),
                   width_mm = 183, height_mm = 200)

## =========================================================
## 10) Save all objects index / 保存所有对象索引
## =========================================================
brms_objects <- list(
  b_MeanAbundTemp = b_MeanAbundTemp,
  b_MeanAbundTrop = b_MeanAbundTrop,
  b_MeanRichTemp  = b_MeanRichTemp,
  b_MeanRichTrop  = b_MeanRichTrop,
  b_MaxAbundTemp  = b_MaxAbundTemp,
  b_MaxAbundTrop  = b_MaxAbundTrop,
  b_MaxRichTemp   = b_MaxRichTemp,
  b_MaxRichTrop   = b_MaxRichTrop
)
saveRDS(brms_objects, file.path(brmsDir, "ALL_brms_models_index.rds"))

message("DONE: brms models + diagnostics + plots (with exclQuantiles & auto family) saved in: ", outDir)








### models including realm - not used due to high VIF factors between variables


# # 1. Abundance, mean anomaly
# MeanAnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
#                                      fitFamily = "gaussian",
#                                      fixedFactors = c("UI2", "Tropical"),
#                                      fixedTerms = list(StdTmeanAnomalyRS=1),
#                                      randomStruct = "(1|SS)+(1|SSB)",
#                                      fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)",
#                                                            "UI2:Tropical",
#                                                            "Tropical:poly(StdTmeanAnomalyRS,1)",
#                                                            "Tropical:poly(StdTmeanAnomalyRS,1):UI2"),
#                                      saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
# 
# summary(MeanAnomalyModelAbund$model)
# 
# # model is rank deficient
# 
# # selected model:
# # LogAbund ~ UI2 + poly(StdTmeanAnomalyRS, 1) + Tropical +
# # UI2:Tropical +
# # Tropical:poly(StdTmeanAnomalyRS, 1):UI2 + 
# # (1 | SS) + (1 | SSB)
# 
# 
# check_collinearity(MeanAnomalyModelAbund$model)
# 
# # the 3 way interaction between land use, climate and realm has a high VIF
# 
# 
# 
# # save the model output
# save(MeanAnomalyModelAbund, file = paste0(outDir, "/MeanAnomalyModelAbund_3way.rdata"))
# 
# # 2. Richness, mean anomaly
# MeanAnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
#                                     fitFamily = "poisson",
#                                     fixedFactors = c("UI2", "Tropical"),
#                                     fixedTerms = list(StdTmeanAnomalyRS=1),
#                                     randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
#                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)",
#                                                           "UI2:Tropical",
#                                                           "Tropical:poly(StdTmeanAnomalyRS,1)",
#                                                           "Tropical:poly(StdTmeanAnomalyRS,1):UI2"),
#                                     saveVars = c("Total_abundance", "SSBS", "NH_3000"))
# 
# summary(MeanAnomalyModelRich$model)
# 
# # selected model:
# # Species_richness ~ UI2 + Tropical + poly(StdTmeanAnomalyRS, 1) +  
# # UI2:poly(StdTmeanAnomalyRS, 1) + UI2:Tropical +  
# # Tropical:poly(StdTmeanAnomalyRS, 1):UI2 + 
# # (1 | SS) +  (1 | SSB) + (1 | SSBS)
# 
# # Model failed to converge: degenerate  Hessian with 2 negative eigenvalues
# 
# check_collinearity(MeanAnomalyModelAbund$model)
# 
# 
# # the 3 way interaction between land use, climate and realm has a high VIF
# 
# 
# # save model output
# save(MeanAnomalyModelRich, file = paste0(outDir, "/MeanAnomalyModelRich_3way.rdata"))
# 
# 
# # 3. Abundance, max anomaly
# MaxAnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
#                                     fitFamily = "gaussian",
#                                     fixedFactors = c("UI2", "Tropical"),
#                                     fixedTerms = list(StdTmaxAnomalyRS=1),
#                                     randomStruct = "(1|SS)+(1|SSB)",
#                                     fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)",
#                                                           "UI2:Tropical",
#                                                           "Tropical:poly(StdTmaxAnomalyRS,1)",
#                                                           "Tropical:poly(StdTmaxAnomalyRS,1):UI2"),
#                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
# 
# summary(MaxAnomalyModelAbund$model)
# 
# # selected model:
# # LogAbund ~ UI2 + UI2:Tropical + Tropical + (1 | SS) + (1 | SSB)
# 
# # save model output
# save(MaxAnomalyModelAbund, file = paste0(outDir, "/MaxAnomalyModelAbund_3way.rdata"))
# 
# 
# # 4. Richness, max anomaly
# MaxAnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
#                                    fitFamily = "poisson",
#                                    fixedFactors = c("UI2", "Tropical"),
#                                    fixedTerms = list(StdTmaxAnomalyRS=1),
#                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
#                                    fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)",
#                                                          "UI2:Tropical",
#                                                          "Tropical:poly(StdTmaxAnomalyRS,1)",
#                                                          "Tropical:poly(StdTmaxAnomalyRS,1):UI2"),  
#                                    saveVars = c("Total_abundance", "SSBS", "NH_3000"))
# 
# summary(MaxAnomalyModelRich$model)
# 
# # selected model:
# # Species_richness ~ UI2 + Tropical + poly(StdTmaxAnomalyRS, 1) + 
# # UI2:poly(StdTmaxAnomalyRS, 1) + UI2:Tropical + Tropical:poly(StdTmaxAnomalyRS, 1):UI2 +  
# # (1 | SS) + (1 | SSB) + (1 | SSBS)
# 
# 
# # Model failed to converge with max|grad| = 0.00613945 (tol = 0.001, component 1)
# 
# 
# 
# # save model output
# save(MaxAnomalyModelRich, file = paste0(outDir, "/MaxAnomalyModelRich_3way.rdata"))
# 
# 
# #load(paste0(outDir, "MeanAnomalyModelAbund_3way.rdata"))
# #load(paste0(outDir, "/MeanAnomalyModelRich_3way.rdata"))
# #load(paste0(outDir, "/MaxAnomalyModelAbund_3way.rdata"))
# #load(paste0(outDir, "/MaxAnomalyModelRich_3way.rdata"))
# 
# 
# summary(MeanAnomalyModelAbund$model)
# summary(MeanAnomalyModelRich$model)
# summary(MaxAnomalyModelAbund$model)
# summary(MaxAnomalyModelRich$model)
# 


##%######################################################%##
#                                                          #
####                  Looking at VIFs                   ####
#                                                          #
##%######################################################%##

# 
# # remove all NAs from the dataset
# 
# predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]
# predictsSites <- predictsSites[!is.na(predictsSites$UI2), ] # 6069
# 
# 
# effs <- predictsSites[ , c("UI2", "StdTmeanAnomalyRS", "StdTmaxAnomalyRS", "NH_5000.rs", "Tropical")]
# 
# 
# x <- cor(effs[, 2:4]) %>% det() # just the continuous variables
# 
# 
# effdata <- predictsSites[ , c("LogAbund", "UI2", "StdTmeanAnomalyRS", "StdTmaxAnomalyRS", "NH_5000.rs", "Tropical")]
# effdata <- effdata[!is.na(effdata$LogAbund), ]
# 
# fitmod <- lm(LogAbund ~ UI2 + StdTmeanAnomalyRS + Tropical, data = effdata)
# stepVIF(fitmod, threshold = 10, verbose = TRUE)
# 
# 
# car::vif(MeanAnomalyModelAbund$model)
# car::vif(fitmod)
# 
# 
# res <- stepwise.vif(dataset = effdata[, c(2:3, 6)], 
#                     metrics = c("UI2", "StdTmeanAnomalyRS", "Tropical"), 
#                     vif.threshold = 5,
#                     verbose = T)
# 
# check_collinearity(MeanAnomalyModelAbund$model)
# check_collinearity(MeanAnomalyModelRich$model)
# check_collinearity(MaxAnomalyModelAbund$model)
# check_collinearity(MaxAnomalyModelRich$model)
# 


# ############## plotting ##############
# 
# 
# ## Realm/anomaly/LU interactions
# 
# 
# 
# ### 1. Abundance, mean anomaly
# 
# nd <- expand.grid(
#   StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
#                         to = max(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
#                         length.out = 100),
#   UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
#              levels = levels(MeanAnomalyModelAbund$data$UI2)),
#   Tropical = factor(c("Tropical", "Temperate"),
#                     levels = levels(MeanAnomalyModelAbund$data$Tropical)))
# 
# # back transform the predictors
# nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
#   transformedX = nd$StdTmeanAnomalyRS,
#   originalX = predictsSites$StdTmeanAnomaly)
# 
# # set richness and abundance to 0 - to be predicted
# nd$LogAbund <- 0
# nd$Species_richness <- 0
# 
# 
# # reference for % difference = primary vegetation and positive anomaly closest to 0
# refRow1 <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))) & (nd$Tropical == "Tropical"))
# refRow2 <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))) & (nd$Tropical == "Temperate"))
# refRow2 <- refRow2 -400
# 
# # adjust plot 1: mean anomaly and abundance
# 
# exclQuantiles <- c(0.025,0.975)
# 
# 
# QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelAbund$data$UI2=="Primary vegetation"],
#   probs = exclQuantiles)
# QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelAbund$data$UI2=="Secondary vegetation"],
#   probs = exclQuantiles)
# QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelAbund$data$UI2=="Agriculture_Low"],
#   probs = exclQuantiles)
# QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelAbund$data$UI2=="Agriculture_High"],
#   probs = exclQuantiles)
# 
# 
# # predict the results
# a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model,data = nd)
# 
# # back transform the abundance values
# a.preds.tmean <- exp(a.preds.tmean)-0.01
# 
# 
# a.preds.tmean_t <- a.preds.tmean[1:400, ] # tropical set
# a.preds.tmean <- a.preds.tmean[401:800, ] # temperate set
# 
# # convert to relative to reference
# a.preds.tmean_t <- sweep(x = a.preds.tmean_t,MARGIN = 2,STATS = a.preds.tmean_t[refRow1,],FUN = '/') # tropical
# a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow2,],FUN = '/') # temperate
# 
# # remove anything above and below the quantiles
# a.preds.tmean_t[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2]) & nd$Tropical == "Tropical"),] <- NA
# 
# # remove anything above and below the quantiles
# a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2]) & nd$Tropical == "Temperate")-400,] <- NA
# 
# # Get the median, upper and lower quants for the plot
# nd$PredMedian[1:400] <- ((apply(X = a.preds.tmean_t,MARGIN = 1,
#                          FUN = median,na.rm=TRUE))*100)-100
# nd$PredUpper[1:400] <- ((apply(X = a.preds.tmean_t,MARGIN = 1,
#                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd$PredLower[1:400] <- ((apply(X = a.preds.tmean_t,MARGIN = 1,
#                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# nd$PredMedian[401:800] <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                          FUN = median,na.rm=TRUE))*100)-100
# nd$PredUpper[401:800] <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd$PredLower[401:800] <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# 
# nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
# 
# plot_data <- nd[nd$Tropical == "Temperate",]
# 
# p1 <- ggplot(data = plot_data, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plot_data$PredLower, ymax = plot_data$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Change in total abundance (%)") +
#   xlab("Standardised Temperature Anomaly") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 750)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12)) +
#   ggtitle("Temperate")
# 
# 
# plot_data2 <- nd[nd$Tropical == "Tropical",]
# 
# p2 <- ggplot(data = plot_data2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plot_data2$PredLower, ymax = plot_data2$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Change in total abundance (%)") +
#   xlab("Standardised Climate Anomaly") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 100)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12), legend.position = "none") +
#   ggtitle("Tropical")
# 
# 
# legend <- get_legend(p1)
# p3 <- cowplot::plot_grid(p1+theme(legend.position = "none"), p2, legend, ncol = 3)
# 
# 
# ### 2. SR mean anomaly
# 
# nd2 <- expand.grid(
#   StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
#                         to = max(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
#                         length.out = 100),
#   UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
#              levels = levels(MeanAnomalyModelRich$data$UI2)),
#   Tropical = factor(c("Tropical", "Temperate"),
#                     levels = levels(MeanAnomalyModelRich$data$Tropical)))
# 
# # back transform the predictors
# nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
#   transformedX = nd2$StdTmeanAnomalyRS,
#   originalX = predictsSites$StdTmeanAnomaly)
# 
# # set richness and abundance to 0 - to be predicted
# nd2$LogAbund <- 0
# nd2$Species_richness <- 0
# 
# 
# # reference for % difference = primary vegetation and positive anomaly closest to 0
# 
# # reference for % difference = primary vegetation and positive anomaly closest to 0
# refRow1 <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & (nd2$Tropical == "Tropical"))
# refRow2 <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & (nd2$Tropical == "Temperate"))
# refRow2 <- refRow2 -400
# 
# 
# # quantiles of data to show on the plot
# exclQuantiles <- c(0.025,0.975)
# 
# 
# QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelRich$data$UI2=="Primary vegetation"],
#   probs = exclQuantiles)
# QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelRich$data$UI2=="Secondary vegetation"],
#   probs = exclQuantiles)
# QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelRich$data$UI2=="Agriculture_Low"],
#   probs = exclQuantiles)
# QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelRich$data$UI2=="Agriculture_High"],
#   probs = exclQuantiles)
# 
# # predict the results
# sr.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd2)
# 
# # back transform the abundance values
# sr.preds.tmean <- exp(sr.preds.tmean)
# 
# 
# sr.preds.tmean_t <- sr.preds.tmean[1:400, ] # tropical set
# sr.preds.tmean <- sr.preds.tmean[401:800, ] # temperate set
# 
# # convert to relative to reference
# sr.preds.tmean_t <- sweep(x = sr.preds.tmean_t,MARGIN = 2,STATS = sr.preds.tmean_t[refRow1,],FUN = '/')
# sr.preds.tmean <- sweep(x = sr.preds.tmean,MARGIN = 2,STATS = sr.preds.tmean[refRow2,],FUN = '/')
# 
# # remove anything above and below the quantiles
# sr.preds.tmean_t[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2]) & nd2$Tropical == "Tropical"),] <- NA
# 
# # remove anything above and below the quantiles
# sr.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2]) & nd2$Tropical == "Temperate")-400,] <- NA
# 
# # Get the median, upper and lower quants for the plot
# nd2$PredMedian[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                 FUN = median,na.rm=TRUE))*100)-100
# nd2$PredUpper[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2$PredLower[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# nd2$PredMedian[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                   FUN = median,nsr.rm=TRUE))*100)-100
# nd2$PredUpper[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                  FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2$PredLower[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                  FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# 
# nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
# 
# 
# plotdata <- nd2[nd2$Tropical == "Temperate", ]
#   
# p4 <- ggplot(data = plotdata, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plotdata$PredLower, ymax = plotdata$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Change in species richness (%)") +
#   xlab("Standardised Climate Anomaly") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 950)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12))+
#   ggtitle("Temperate")
# 
# 
# 
# plotdata2 <- nd2[nd2$Tropical == "Tropical", ]
# 
# p5 <- ggplot(data = plotdata2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plotdata2$PredLower, ymax = plotdata2$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Change in species richness (%)") +
#   xlab("Standardised Climate Anomaly") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 100)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12), legend.position = "none")+
#   ggtitle("Tropical")
# 
# legend2 <- get_legend(p4)
# p6 <- cowplot::plot_grid(p4+theme(legend.position = "none"), p5, legend, ncol = 3)
# 



### 3. Abundance, max anomaly

# 3-way interaction not selected for this set

#nd3 <- expand.grid(
#  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
#                        to = max(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
#                        length.out = 100),
#  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
#             levels = levels(MaxAnomalyModelAbund$data$UI2)),
#  Tropical = factor(c("Tropical", "Temperate"),
#                    levels = levels(MaxAnomalyModelAbund$data$Tropical)))

# back transform the predictors
#nd3$StdTmaxAnomaly <- BackTransformCentreredPredictor(
#  transformedX = nd3$StdTmaxAnomalyRS,
#  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
#nd3$LogAbund <- 0
#nd3$Species_richness <- 0


# reference for % difference = primary vegetation and positive anomaly closest to 0
#refRow <- which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomaly==min(abs(nd3$StdTmaxAnomaly))) & (nd3$Tropical == "Temperate"))

# quantiles of data to show on the plot
#exclQuantiles <- c(0.025,0.975)


#QPV <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
#  MaxAnomalyModelAbund$data$UI2=="Primary vegetation"],
#  probs = exclQuantiles)
#QSV <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
#  MaxAnomalyModelAbund$data$UI2=="Secondary vegetation"],
#  probs = exclQuantiles)
#QAL <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
#  MaxAnomalyModelAbund$data$UI2=="Agriculture_Low"],
#  probs = exclQuantiles)
#QAH <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
#  MaxAnomalyModelAbund$data$UI2=="Agriculture_High"],
#  probs = exclQuantiles)

# predict the results
#a.preds.tmean <- PredictGLMERRandIter(model = MaxAnomalyModelAbund$model,data = nd3)

# back transform the abundance values
#a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
#a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
#a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS < QPV[1])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS > QPV[2])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS < QSV[1])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS > QSV[2])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS < QAL[1])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS > QAL[2])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS < QAH[1])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
#nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                         FUN = median,na.rm=TRUE))*100)-100
#nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
#nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


#nd3$UI2 <- factor(nd3$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

#p3 <- ggplot(data = nd3, aes(x = StdTmaxAnomaly, y = PredMedian)) +
#  geom_line(aes(col = UI2), size = 1) +
#  geom_ribbon(aes(ymin = nd3$PredLower, ymax = nd3$PredUpper, fill = UI2), alpha = 0.2) +
#  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#  facet_wrap(~Tropical, ncol = 2) +
#  theme_bw() +
#  labs(fill = "% NH", col = "% NH") +
#  ylab("Total Abundance (%)") +
#  xlab("Standardised Climate Anomaly Maximum") +
#  xlim(c(-0.5, 2)) +
#  ylim(c(-100, 150)) +
#  theme(aspect.ratio = 1, text = element_text(size = 12))


# 
# 
# 
# ### 4. Species richness, max anomaly
# 
# 
# nd4 <- expand.grid(
#   StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelRich$data$StdTmaxAnomalyRS),
#                        to = max(MaxAnomalyModelRich$data$StdTmaxAnomalyRS),
#                        length.out = 100),
#   UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
#              levels = levels(MaxAnomalyModelRich$data$UI2)),
#   Tropical = factor(c("Tropical", "Temperate"),
#                     levels = levels(MaxAnomalyModelRich$data$Tropical)))
# 
# # back transform the predictors
# nd4$StdTmaxAnomaly <- BackTransformCentreredPredictor(
#   transformedX = nd4$StdTmaxAnomalyRS,
#   originalX = predictsSites$StdTmaxAnomaly)
# 
# # set richness and abundance to 0 - to be predicted
# nd4$LogAbund <- 0
# nd4$Species_richness <- 0
# 
# 
# # reference for % difference = primary vegetation and positive anomaly closest to 0
# refRow1 <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomaly==min(abs(nd4$StdTmaxAnomaly))) & (nd4$Tropical == "Tropical"))
# refRow2 <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomaly==min(abs(nd4$StdTmaxAnomaly))) & (nd4$Tropical == "Temperate"))
# refRow2 <- refRow2 -400
# 
# # quantiles of data to show on the plot
# exclQuantiles <- c(0.025,0.975)
# 
# 
# QPV <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
#   MaxAnomalyModelRich$data$UI2=="Primary vegetation"],
#   probs = exclQuantiles)
# QSV <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
#   MaxAnomalyModelRich$data$UI2=="Secondary vegetation"],
#   probs = exclQuantiles)
# QAL <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
#   MaxAnomalyModelRich$data$UI2=="Agriculture_Low"],
#   probs = exclQuantiles)
# QAH <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
#   MaxAnomalyModelRich$data$UI2=="Agriculture_High"],
#   probs = exclQuantiles)
# 
# # predict the results
# sr.preds.tmean <- PredictGLMERRandIter(model = MaxAnomalyModelRich$model,data = nd4)
# 
# # back transform the abundance values
# sr.preds.tmean <- exp(sr.preds.tmean)
# 
# 
# sr.preds.tmean_t <- sr.preds.tmean[1:400, ] # tropical set
# sr.preds.tmean <- sr.preds.tmean[401:800, ] # temperate set
# 
# # convert to relative to reference
# sr.preds.tmean_t <- sweep(x = sr.preds.tmean_t,MARGIN = 2,STATS = sr.preds.tmean_t[refRow1,],FUN = '/')
# sr.preds.tmean <- sweep(x = sr.preds.tmean,MARGIN = 2,STATS = sr.preds.tmean[refRow2,],FUN = '/')
# 
# # remove anything above and below the quantiles
# sr.preds.tmean_t[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS < QPV[1]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS > QPV[2]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS < QSV[1]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS > QSV[2]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS < QAL[1]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS > QAL[2]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS < QAH[1]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS > QAH[2]) & nd4$Tropical == "Tropical"),] <- NA
# 
# # remove anything above and below the quantiles
# sr.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS < QPV[1]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS > QPV[2]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS < QSV[1]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS > QSV[2]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS < QAL[1]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS > QAL[2]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS < QAH[1]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS > QAH[2]) & nd4$Tropical == "Temperate")-400,] <- NA
# 
# # Get the median, upper and lower quants for the plot
# nd4$PredMedian[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                  FUN = median,na.rm=TRUE))*100)-100
# nd4$PredUpper[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                 FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd4$PredLower[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                 FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# nd4$PredMedian[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                    FUN = median,na.rm=TRUE))*100)-100
# nd4$PredUpper[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd4$PredLower[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# 
# 
# nd4$UI2 <- factor(nd4$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
# 
# plot.data <- nd4[nd4$Tropical == "Temperate",]
# 
# p7 <- ggplot(data = plot.data, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plot.data$PredLower, ymax = plot.data$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Species Richness (%)") +
#   xlab("Standardised Climate \nAnomaly Maximum") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 100)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12))+
#   ggtitle("Temperate")
# 
# plot.data2 <- nd4[nd4$Tropical =="Tropical",]
# 
# p8 <- ggplot(data = plot.data2, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plot.data2$PredLower, ymax = plot.data2$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Species Richness (%)") +
#   xlab("Standardised Climate \nAnomaly Maximum") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 100)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12), legend.position = "none") + 
#   ggtitle("Tropical")
# 
# #"#009E73" - green
# #"#0072B2" - blue
# #"#E69F00" - yellow
# #"#D55E00" - red
# 
# legend3 <- get_legend(p7)
# p9 <- plot_grid(p7+theme(legend.position = "none"), p8, legend, ncol = 3)
# 
# 
# 
# # organise plots into one document
# 
# plot_grid(p3, p6,  p9, ncol = 1, 
#           labels = c("A", "B", "C"), label_size = 12, rel_widths = c(1,1,0.5))
# 
# # save plot
# ggsave(filename = paste0(outDir, "Plots_climate_LU_Tropical_ALL_allinteractions.pdf"), width = 9, height = 9, units = "in")
# 
# 
# 



t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()





