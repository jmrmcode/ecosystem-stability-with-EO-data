# ecosystem-stability-with-EO-data
A compilation of R code to compute and model ecosystem stability metrics (resilience) from remotely-sensed time series of forest evapotranspiration (ET).
### Compute_ET_anomalies.R
Code for computing monthly ET anomalies (&#916;<sub>ET</sub>(*s, t*)) using [Zhang et al. (2019)](https://www.sciencedirect.com/science/article/abs/pii/S003442571830590X) [dataset](https://developers.google.com/earth-engine/datasets/catalog/CAS_IGSNRR_PML_V2) and following [Goetz et al. (2006)](https://www.sciencedirect.com/science/article/abs/pii/S0034425706000289) and [White et al. (2019)](https://zslpublications.onlinelibrary.wiley.com/doi/10.1002/rse2.148) approach where

&#916;<sub>ET</sub>(*s, t*) = (ET(*s, t*) - mean<sub>u&isin;m</sub>[ET(*s, u*)]) / sd<sub>u&isin;m</sub>[ET(*s, u*)]

ET(*s, t*) is the mean modeled ET in forests from subbasin *s* at date *t*, *m* is a month of the year and mean<sub>u&isin;m</sub>[ET(*s, u*)]) and sd<sub>u&isin;m</sub>[ET(*s, u*)] are the mean and standard deviation of ET for subbasin *s* over all dates, *u*, across the entire period (2003–2017) falling within month *m*, respectively.
### Compute_Resilience.R
Code for computing resilience to water (positive anomalies) and productivity (negative anomalies) losses using a temporal moving window algorithm based on [White et al. (2019)](https://zslpublications.onlinelibrary.wiley.com/doi/10.1002/rse2.148).
### Model_Resilience.R
Code to run varying-intercept, varying-slope models using INLA and the conditional R<sup>2</sup> based on [Nakagawa et al 2017](https://royalsocietypublishing.org/doi/10.1098/rsif.2017.0213) and [Johnson 2014](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12225).
