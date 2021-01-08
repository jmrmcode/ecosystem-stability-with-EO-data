# ecosystem-stability-with-EO-data
A compilation of R code to import EO data from GEE and compute ecosystem stability metrics
### Evapotranspiration.R
Grab [Evapotranspiration](https://developers.google.com/earth-engine/datasets/catalog/CAS_IGSNRR_PML_V2) from GEE and compute ET anomalies (&#916;<sub>ET</sub>(*w, t*)) similar to [Goetz et al (2006)](https://www.sciencedirect.com/science/article/abs/pii/S0034425706000289) and [White et al (2019)](https://zslpublications.onlinelibrary.wiley.com/doi/10.1002/rse2.148) where

&#916;<sub>ET</sub>(*w, t*) = (ET(*w, t*) - mean<sub>u&isin;m</sub>[ET(*w, u*)]) / sd<sub>u&isin;m</sub>[ET(*w, u*)]

ET(*w, t*) is the ET of watershed *w* at date *t*, *m* is a month of the year and mean<sub>u&isin;m</sub>[ET(*w, u*)]) and sd<sub>u&isin;m</sub>[ET(*w, u*)] are the mean and standard deviation of ET for watershed *w* over all dates, *u*, across the entire period (2003–2017) falling within month *m*, respectively.