library(reshape)
library(sf)
library(raster)
library(ncdf4)
library(stringr)
library(dataRetrieval)
library(ggplot2)


LittleRiver <- "01671100"
JamesRiver <- "02037500"
MississippiRiver <- "07374000"



gauge <- MississippiRiver

precipPull <- function(x, basin) {
  precip_stack <- raster::brick(paste("/Users/jameskitchens/Documents/GitHub/streamflow-predictions/Data/NOAA-PSL/CPC-Unified-Gauge-Based-Analysis-Of-Daily-Precipitation-Over-CONUS-RT/precip.V1.0.",x,".nc",sep=""))
  precip_stack <- raster::rotate(precip_stack)
  meanDate <- raster::extract(precip_stack, as(sf::st_geometry(basin), Class="Spatial"), fun=mean, na.rm=TRUE)
  precip_measurements <- reshape::melt(meanDate)
  precip_measurements <- precip_measurements[,c(2,3)]
  colnames(precip_measurements) <- c("Date", "Precipitation")
  precip_measurements$Date <- as.Date(stringr::str_replace_all(stringr::str_replace_all(precip_measurements$Date, "X", ""), "\\.", "-"))
  return(precip_measurements)
}

streamflow <- dataRetrieval::readNWISdv(gauge, "00060", "2007-01-01", "2020-12-07")
streamflow <- streamflow[,c("Date", "X_00060_00003")]
colnames(streamflow) <- c("Date", "CFS")
basin <- dataRetrieval::findNLDI(nwis=gauge, find="basin")$basin
years <- seq(2007, 2020)
precip <- lapply(years, FUN=precipPull, basin=basin)
precip <- do.call(rbind.data.frame, precip)

combined <- merge(streamflow, precip, all=TRUE)

ggplot(data=combined) +
  geom_line(aes(Date, Precipitation*100000), color="green") +
  geom_line(aes(Date, CFS)) +
  scale_x_date(date_breaks="3 month") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

