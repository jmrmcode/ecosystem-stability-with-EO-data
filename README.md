# ecosystem-stability-with-EO-data
A compilation of R code to import EO data from GEE and compute ecosystem stability metrics
### Evapotranspiration.R
Grab [Evapotranspiration](https://developers.google.com/earth-engine/datasets/catalog/CAS_IGSNRR_PML_V2) from GEE and compute ET anomalies (&#916;<sub>ET</sub>(*w, t*)) similar to [Goetz et al. (2006)](https://www.sciencedirect.com/science/article/abs/pii/S0034425706000289) and [White et al. (2019)](https://zslpublications.onlinelibrary.wiley.com/doi/10.1002/rse2.148) where

&#916;<sub>ET</sub>(*w, t*) = (ET(*w, t*) - mean<sub>u&isin;m</sub>[ET(*w, u*)]) / sd<sub>u&isin;m</sub>[ET(*w, u*)]

ET(*w, t*) is the ET of watershed *w* at date *t*, *m* is a month of the year and mean<sub>u&isin;m</sub>[ET(*w, u*)]) and sd<sub>u&isin;m</sub>[ET(*w, u*)] are the mean and standard deviation of ET for watershed *w* over all dates, *u*, across the entire period (2003–2017) falling within month *m*, respectively.
ET values are the summation of its three components: Vegetation transpiration (ET<sub>c</sub>), Soil evaporation (ET<sub>s</sub>), and Interception from vegetation canopy (ET<sub>i</sub>).
### Compute_Resilience.R
Compute resistance and resilience measures based on ET anomalies obtained by running Evapotranspiration.R. Resistance is computed on the scaled ET anomalies (&#10698;<sub>ET</sub>):

Resistance<sub>ET</sub>(*w*) = (mean[&#10698;<sub>ET>1</sub>])<sup>-1</sup> where &#10698;<sub>ET</sub> = 2 + [((-4)(&#916;<sub>ET</sub> - &#916;<sub>ETmin</sub>)) / (&#916;<sub>ETmax</sub> - &#916;<sub>ETmin</sub>)]. Note that &#10698;<sub>ET</sub> &#x220A; (-2, +2).

Resilience is computed as recovery time similar to [White et al. (2019)](https://zslpublications.onlinelibrary.wiley.com/doi/10.1002/rse2.148).
