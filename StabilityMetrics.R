## "ET_anomalies" object created in "Evapotranspiration.R" must be load 
## before running the code below

#### COMPUTE RESISTANCE OF ET ANOMALY TIME SERIES

library("BBmisc")
# create an empty data frame to populate below
resistET <- data.frame(resist = as.numeric(), basin = as.numeric())
r <- 0
for (i in names(ET_anomalies)) {
  r <- r + 1
  # select the i watershed
  resB <- ET_anomalies[[which(names(ET_anomalies)==i)]]
  # get the number of anomalies above the maximum SD
  threshold <- max(unique(floor(resB$ET_anomaly[which(year(resB$Date) < 2017)])))
  if (threshold > 1) {
    resBnorm <- normalize(resB$ET_anomaly, range = c(-2, +2), method = "range")
    resBnorm <- data.frame(ET = resBnorm, date = resB$Date)
    resistET[r, 1] <- 1 / mean(resBnorm[which(resBnorm$ET[which(year(resBnorm$date) < 2017)] > 1), 1])
    resistET[r, 2] <- as.integer(i)
  } else {}
  # set the number of SDs
  resistET[r, 1] <- 1 / mean(resB[which(resB$ET_anomaly[which(year(resB$Date) < 2017)] > 1), 1])
  resistET[r, 2] <- as.integer(i)
  rm(list=setdiff(ls(), c("resistET", "ET_anomalies", "r")))
}


## COMPUTE RESILIENCE OF ET ANOMALY TIME SERIES BASED ON WHITE ET AL 2019

library("lubridate")
# create an empty data frame to populate below
resilET <- data.frame(resil = as.numeric(), basin = as.numeric())
rec1 <- NULL
r <- 0
for (w in names(ET_anomalies)) {
  r <- r + 1
  resilB <- ET_anomalies[[which(names(ET_anomalies) == w)]]
  # not consider 2017 in the threshold selection
  threshold <- max(unique(floor(resilB$ET_anomaly[which(year(resilB$Date) < 2017)])))
  posit <- which(resilB$ET_anomaly > threshold)
  # discard dates in 2017 which are abovre the threshold
  alldates <- resilB$Date[posit][which(year(resilB$Date[posit]) < 2017)]
  # discard anomalies occurring in specific months
  win_fal <- c()
  months <- month(alldates)
  
  # while there are not anomalies above the threshold in the months selected above
  # keep decreasing the threshold by 1 unit each iteration 
  while (length(months[! months %in% win_fal]) == 0) {
    threshold <- threshold - 1
    posit <- which(resilB$ET_anomaly > threshold)
    alldates <- resilB$Date[posit][which(year(resilB$Date[posit]) < 2017)]
    months <- month(alldates)
  }
  # discard specific months
  months <- months[! months %in% win_fal]
  alldates <- alldates[! month(alldates) %in% win_fal]
  posit1 <- match(resilB$Date, alldates)
  posit <- which(is.na(posit1) == F)
  
  # when there are more than one anomaly > threshold
  if (length(posit) > 1) {
    
    for (p in posit) {
      resilB2 <- resilB[c(p:length(resilB$ET_anomaly)), ]
      mo <- rollmean(resilB2$ET_anomaly, k = 6)
      recTime <- (which(mo < 0)[1]) * 8 # expressed in days as there are 8 days between two observations
      rec1  <- abind(rec1, as.array(recTime))
      resilET[r, 1] <- 1 / mean(rec1)
      resilET[r, 2] <- w
    }
    
    # when there is just one anomaly above the threshold
  } else {
    resilB2 <- resilB[c(posit:length(resilB$ET_anomaly)), ]
    mo <- rollmean(resilB2$ET_anomaly, k = 6)
    recTime <- (which(mo < 0)[1]) * 8 # expressed in days as there are 8 days between two observations
    resilET[r, 1] <- 1 / recTime
    resilET[r, 2] <- w
  }
  
  rm(list=setdiff(ls(), c("resilET", "ET_anomalies", "r")))
  rec1 <- NULL
}
