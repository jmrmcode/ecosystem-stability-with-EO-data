library("plyr")
library("dplyr")
library("rgee")
library("tidyverse")
library('abind')
library("sf")
library("geojsonio")
library('tidyr')
library("lubridate")
library('plotly')
ee_Initialize()


###*** IMPORT EVAPOTRANSPIRATION TIME SERIES

# import watersheds
waterS <- ee$FeatureCollection("USGS/WBD/2017/HUC06")

# extract watershed codes (level 3: Basin)
waterSheds <- waterS$aggregate_array("huc6")$getInfo()

ET.Data <- list()

# loop to grab Evapotranspiration (ET) time series and discard watersheds outside the CONUS
row <- 0
for (i in waterSheds) {
  
  row <- row + 1
  print(row)
  
  # select one watershed by iteration
  mask <- ee$FeatureCollection("USGS/WBD/2017/HUC06")$filter(ee$Filter$eq("huc6", i))
  # convert to geometry and extract the boundaries
  region <- mask$geometry()$bounds()
  
  # load the NLCD data set and build a forest mask
  nlcd <- ee$ImageCollection("USGS/NLCD")$
    filterDate("2016-01-01", "2016-12-31")$
    filterBounds(mask)$
    first()$
    select("landcover")$
    eq(list(41,42,43))$ 
    reduce(ee$Reducer$sum())
  
  # load ET imageCollection in forest pixels
  ET <- ee$ImageCollection("CAS/IGSNRR/PML/V2")$
    filterDate("2003-01-01", "2017-12-31")$
    filterBounds(mask)$
    map(function(x) x$reproject("EPSG:4326", NULL, 500)$select("Ec","Es","Ei")$reduce(ee$Reducer$sum()))$
  # mask non-forest
    map(function(x) x$updateMask(nlcd))

  # extract ET values and dates by watershed
  ee_nc_ET <- try(ee_extract(x = ET$select("sum"), y = region, sf = FALSE, scale = 500))
  if(class(ee_nc_ET) %in% 'try-error') {next} else {
    
    colnames(ee_nc_ET) <- substring(colnames(ee_nc_ET), regexpr("X2", colnames(ee_nc_ET)) + 1)
    ee_nc_ET$name <- "ET"
    
    # convert into time series
    time <- as.Date(colnames(ee_nc_ET[1, -ncol(ee_nc_ET)]), format = "%Y.%m.%d")
    if (ncol(ee_nc_ET)==1) {next} else {
      data <- t(ee_nc_ET[1, -ncol(ee_nc_ET)])
      ET.Data[[row]] <- data.frame(ET = data[, 1], Date = time, row.names = NULL)
      names(ET.Data)[[row]] <- i
      
      # clean memory
      rm(list=setdiff(ls(), c("ET.Data", "i", "waterSheds", "row")))
      gc(full = TRUE)
    }
  }
}

ET.Data <- compact(ET.Data)


### convert the list of watersheds into a data frame
annualET <- matrix(NA, nrow = length(lubridate::year(ET.Data[[1]][, 2])), ncol = length(ET.Data))
c <- 0
for (i in 1:length(ET.Data)) {
  c <- c + 1

  tt <- ts(data = ET.Data[[i]]$ET, start = c(1, nrow(data)), frequency = 46)
  ## raw time series
  annualET[, c] <- tt
}
annualET <- as.data.frame(annualET)
annualET <- gather(annualET[,c(1:length(ET.Data))], key = "basin", "ET")
date <- ET.Data[[1]][, 2]
annualET$basin <- rep(names(ET.Data), each = length(date))
annualET$date <- rep(date, length(names(ET.Data)))
annualET$region <- rep(substr(names(ET.Data), 1, 2), each = length(date))
annualET <- na.omit(annualET)


### compute ET monthly anomalies (ET_i - mean(month))/ sd(month)
ET_anomalies <- list()

obs1 <- NULL
obs2 <- NULL
date1 <- NULL
date2 <- NULL
r <- 0
for (i in unique(annualET$basin)) {
  r <- r + 1
  anom1 <- annualET[which(annualET$basin==i), ]
  months <- unique(month(anom1$date))
  for (m in months) {
    dm <- anom1[which(month(anom1$date)==m), c(2,3)]
    obs1 <- as.array((dm$ET - mean(dm$ET))/ sd(dm$ET) )
    obs2 <- abind(obs2, obs1)
    date1 <- as.array(dm$date)
    date2 <- as.Date(abind(date2, date1), origin = "1970-01-01")
    df <- dplyr::arrange(data.frame(ET_anomaly = obs2, Date = date2), Date)
  }
  obs1 <- NULL; obs2 <- NULL; date1 <- NULL; date2 <- NULL
  ET_anomalies [[r]] <- df
}
names(ET_anomalies) <- unique(annualET$basin)

# plot ET anomaly trajectory in a arbitrary watershed
basin <- "110800"
ts <- ET_anomalies[[which(names(ET_anomalies)==basin)]]
plot_ly(ts, x = ~Date, y = ~ET_anomaly, alpha = 1) %>% add_lines()