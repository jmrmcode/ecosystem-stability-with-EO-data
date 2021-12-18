library("INLA")

# import the dataset. Read METADATA.txt for more information on the dataset
data <- read.csv("~/ETresilience_HUC8subbasin.csv")

# make copies of ecoregionLevelIII
data$ecoregionLevelIIIcopy1 <- data$ecoregionLevelIII
data$ecoregionLevelIIIcopy2 <- data$ecoregionLevelIII
data$ecoregionLevelIIIcopy3 <- data$ecoregionLevelIII
data$ecoregionLevelIIIcopy4 <- data$ecoregionLevelIII
data$ecoregionLevelIIIcopy5 <- data$ecoregionLevelIII
data$ecoregionLevelIIIcopy6 <- data$ecoregionLevelIII

# import the graph
H <- inla.read.graph(filename = "~/ecoregionIII_graphSnap8")

# model formula
# resilience to water losses model code
# model for resilience to productivity losses (ET_resilienceExtBELOW) includes BelowIntensity instead of numbAboveAnomal.
Resilience_model.Formula <- ET_resilienceExtABOVE ~ 1 +
  f(ecoregionLevelIIIcopy1, numbAboveAnomal, model = 'iid')+
  #  f(ecoregionLevelIIIcopy1, BelowIntensity, model = 'iid') +
  f(ecoregionLevelIIIcopy2, mixedCover, model = 'iid') +
  f(ecoregionLevelIIIcopy3, deciduousCover, model = 'iid') +
  f(ecoregionLevelIIIcopy4, Precipitation, model = 'iid') +
  f(ecoregionLevelIIIcopy5, SnowWaterEquivalent, model = 'iid') +
  f(ecoregionLevelIIIcopy6, EvaporativeDeficit, model = 'iid') +
  f(ecoregionIII_ID, model = 'bym', graph = H)

# run the model
mod1 <- inla(Resilience_model.Formula, data = data, family = "gamma",
             control.family = list(control.link = list(model = "log")),
             control.compute = list(config = T, dic = T, cpo = T), verbose = F,
             control.predictor = list(compute = T, link = c(rep(1, nrow(data)))))


### estimate the conditional R2 based on Nakagawa et al 2017 Journal of The Royal Society Interface
# sigma2_epsilon in equation 6 (Requena-Mullor et al) based on the trigamma function 
# see Table 1 in Nakagawa et al 2017 Journal of The Royal Society Interface 
shape <- mod1$summary.hyperpar$mean[1]
sigma2_epsilon <- trigamma(shape)
sigma2_epsilon

# sigma2_u and sigma2_v in equation 6 (Requena-Mullor et al)
var.spatial <- inla.rmarginal(100000,inla.tmarginal(function(x) 1/x,
           mod1$marginals.hyper$`Precision for ecoregionIII_ID (spatial component)`))
var.iid <- inla.rmarginal(100000,inla.tmarginal(function(x) 1/x,
           mod1$marginals.hyper$`Precision for ecoregionIII_ID (iid component)`))
sigma2_u <- mean(var.spatial)
sigma2_v <- mean(var.iid)


# # sigma2_p; based on equation 11 in Johnson 2014 Methods in Ecology and Evolution
# # run lines 63 - 68 for every covariate to obtain the sigma_p values shown below
# data <- data[c(1:1703), ]
# Z <- model.matrix(~ -1 + data$numbAboveAnomal)
# #SIGMA <- diag(mod1$summary.random$ecoregionLevelIIIcopy2$mean)
# SIGMA <- inla.emarginal(function(x) 1/x, 
#          mod1$marginals.hyperpar$`Precision for ecoregionLevelIIIcopy1`)
# sigma2_slopeED <- sum(diag(Z %*% SIGMA %*% t(Z))) / 1703
# sigma2_slopeED


# sum of sigma2_p for resilience to water losses
slopes <- sum(0.06011861 + 2.531085e-06 + 0.04554138 + 0.1082932 + 0.0001518075 + 0.06043103)

# # sum of sigma2_p for resilience to productivity losses
# slopes <- sum(0.1241315 + 0.01284914 + 0.0001328778 + 0.07377675 + 0.0001475784 + 0.07946951)

# conditional R2
R2_gammaC <- (sigma2_u + slopes + sigma2_v) / (sigma2_u +  sigma2_epsilon + slopes + sigma2_v)
R2_gammaC
