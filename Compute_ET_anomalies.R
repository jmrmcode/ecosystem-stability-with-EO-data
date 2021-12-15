library("ggplot2")
library("reshape2")
library("tidyr")
library("plotly")
library("aTSA")
library("xts")
library("lubridate")
library("abind")
library("ForeCA")

# import ET time series (download first ET_timeSeriesBySubbasin.Rdata from the main branch)
load("~/ET_timeSeriesBySubbasin.Rdata")

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

# plot ET anomaly trajectory from an arbitrary watershed
basin <- "8030205"
ts <- ET_anomalies[[which(names(ET_anomalies)==basin)]]
plot_ly(ts, x = ~Date, y = ~ET_anomaly, alpha = 1) %>% add_lines()
