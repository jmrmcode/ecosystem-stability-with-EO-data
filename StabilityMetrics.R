library("abind")
library("lubridate")


# load the ET time series (one per subbasin)
load(~/ET_timeSeriesBySubbasin.Rdata)

#### COMPUTE RESILIENCE TO WATER LOSSES FROM THE TIME SERIES OF ET ANOMALIES
##  ANOMALIES ABOVE THE BASELINE

# create an empty data frame
resilET <- data.frame(resil = as.numeric(), basin = as.numeric(), numbAnomal = as.numeric(), intensityAbove = as.numeric())
rec1 <- NULL
mmonts <- list()
r <- 0
for (w in names(ET_anomalies)) {
  r <- r + 1
  resilB <- ET_anomalies[[which(names(ET_anomalies)==w)]]
  # not consider 2017 in the threshold selection
  threshold <- max(unique(floor(resilB$ET)))
  posit1 <- which(resilB$ET > threshold)
  # discard dates in 2017 which are above the threshold
  alldates <- resilB$date[posit1][which(year(resilB$date[posit1]) < 2017)]
  posit <- which(resilB$date %in% alldates)
  
  # when there is not any anomaly before 2017
  if (length(posit) == 0) {
    threshold2 <- threshold - 1 
    posit1 <- which(resilB$ET > threshold2)
    alldates <- resilB$date[posit1][which(year(resilB$date[posit1]) < 2017)]
    posit <- which(resilB$date %in% alldates)
    
    # if there are more than one anomaly > threshold
    if (length(posit) > 1) {
      
      for (p in posit) {
        resilB2 <- resilB[c(p:length(resilB$ET)), ]
        mo <- rollmean(resilB2$ET, k = 6)
        recTime <- (which(mo < 0)[1]) * 8 # in days as there are 8 days between two observations
        rec1  <- abind(rec1, as.array(recTime))
        resilET[r, 1] <- 1 / mean(rec1) # inverse of number of days
        resilET[r, 2] <- w # subbasin
        resilET[r, 3] <- length(posit) # number of anomalies
        resilET[r, 4] <- mean(resilB$ET[posit]) # intensity of anomalies
        mmonts[[r]] <- months(alldates) # months occurring anomalies
      }
      # when there is just one anomaly above the threshold
    } else {
      resilB2 <- resilB[c(posit:length(resilB$ET)), ]
      mo <- rollmean(resilB2$ET, k = 6)
      recTime <- (which(mo < 0)[1]) * 8 # in days as there are 8 days between two observations
      resilET[r, 1] <- 1 / recTime
      resilET[r, 2] <- w
      resilET[r, 3] <- length(posit)
      resilET[r, 4] <- mean(resilB$ET[posit])
      mmonts[[r]] <- months(alldates)
    }
    rm(list=setdiff(ls(), c("resilET", "ET_anomalies", "r", "mmonts")))
    rec1 <- NULL
    ########################################
    # when there are anomalies before 2017
  } else {
  
  # if there are more than one anomaly > threshold
  if (length(posit) > 1) {
    
    for (p in posit) {
      resilB2 <- resilB[c(p:length(resilB$ET)), ]
      mo <- rollmean(resilB2$ET, k = 6)
      recTime <- (which(mo < 0)[1]) * 8 # in days as there are 8 days between two observations
      rec1  <- abind(rec1, as.array(recTime))
      resilET[r, 1] <- 1 / mean(rec1)
      resilET[r, 2] <- w
      resilET[r, 3] <- length(posit)
      resilET[r, 4] <- mean(resilB$ET[posit])
      mmonts[[r]] <- months(alldates)
    }
     # when there is just one anomaly above the threshold
  } else {
    resilB2 <- resilB[c(posit:length(resilB$ET)), ]
    mo <- rollmean(resilB2$ET, k = 6)
    recTime <- (which(mo < 0)[1]) * 8 # in days as there are 8 days between two observations
    resilET[r, 1] <- 1 / recTime
    resilET[r, 2] <- w
    resilET[r, 3] <- length(posit)
    resilET[r, 4] <- mean(resilB$ET[posit])
    mmonts[[r]] <- months(alldates)
  }
    
  rm(list=setdiff(ls(), c("resilET", "ET_anomalies", "r", "mmonts")))
  rec1 <- NULL
  }
}


### ET ANOMALIES BELOW THE BASELINE
# create an empty data frame
resilET <- data.frame(resil = as.numeric(), basin = as.numeric(), numbAnomal = as.numeric(), intensityAbove = as.numeric())
rec1 <- NULL
mmonts <- list()
r <- 0
for (w in names(ET_anomalies)) {
  r <- r + 1
  resilB <- ET_anomalies[[which(names(ET_anomalies)==w)]]
  # not consider 2017 in the threshold selection
  threshold <- min(unique(ceiling(resilB$ET)))
  posit1 <- which(resilB$ET < threshold)
  # discard dates in 2017 which are above the threshold
  alldates <- resilB$date[posit1][which(year(resilB$date[posit1]) < 2017)]
  posit <- which(resilB$date %in% alldates)
  
  # when there is not any anomaly before 2017
  if (length(posit) == 0) {
    threshold2 <- threshold + 1 
    posit1 <- which(resilB$ET < threshold2)
    alldates <- resilB$date[posit1][which(year(resilB$date[posit1]) < 2017)]
    posit <- which(resilB$date %in% alldates)
    
    # if there are more than one anomaly > threshold
    if (length(posit) > 1) {
      
      for (p in posit) {
        resilB2 <- resilB[c(p:length(resilB$ET)), ]
        mo <- rollmean(resilB2$ET, k = 6)
        recTime <- (which(mo > 0)[1]) * 8 # in days as there are 8 days between two observations
        rec1  <- abind(rec1, as.array(recTime))
        resilET[r, 1] <- 1 / abs(mean(rec1)) # inverse of number of days
        resilET[r, 2] <- w # subbasin
        resilET[r, 3] <- length(posit) # number of anomalies
        resilET[r, 4] <- abs(mean(resilB$ET[posit])) # intensity of anomalies
        mmonts[[r]] <- months(alldates) # months occurring anomalies
      }
      # when there is just one anomaly above the threshold
    } else {
      resilB2 <- resilB[c(posit:length(resilB$ET)), ]
      mo <- rollmean(resilB2$ET, k = 6)
      recTime <- (which(mo > 0)[1]) * 8 # in days as there are 8 days between two observations
      resilET[r, 1] <- 1 / abs(recTime)
      resilET[r, 2] <- w
      resilET[r, 3] <- length(posit)
      resilET[r, 4] <- abs(mean(resilB$ET[posit]))
      mmonts[[r]] <- months(alldates)
    }
    rm(list=setdiff(ls(), c("resilET", "ET_anomalies", "r", "mmonts")))
    rec1 <- NULL
    ########################################
    # when there are anomalies before 2017
  } else {
    
    # if there are more than one anomaly > threshold
    if (length(posit) > 1) {
      
      for (p in posit) {
        resilB2 <- resilB[c(p:length(resilB$ET)), ]
        mo <- rollmean(resilB2$ET, k = 6)
        recTime <- (which(mo > 0)[1]) * 8 # in days as there are 8 days between two observations
        rec1  <- abind(rec1, as.array(recTime))
        resilET[r, 1] <- 1 / abs(mean(rec1))
        resilET[r, 2] <- w
        resilET[r, 3] <- length(posit)
        resilET[r, 4] <- abs(mean(resilB$ET[posit]))
        mmonts[[r]] <- months(alldates)
      }
      # when there is just one anomaly above the threshold
    } else {
      resilB2 <- resilB[c(posit:length(resilB$ET)), ]
      mo <- rollmean(resilB2$ET, k = 6)
      recTime <- (which(mo > 0)[1]) * 8 # in days as there are 8 days between two observations
      resilET[r, 1] <- 1 / abs(recTime)
      resilET[r, 2] <- w
      resilET[r, 3] <- length(posit)
      resilET[r, 4] <- abs(mean(resilB$ET[posit]))
      mmonts[[r]] <- months(alldates)
    }
    
    rm(list=setdiff(ls(), c("resilET", "ET_anomalies", "r", "mmonts")))
    rec1 <- NULL
  }
}
