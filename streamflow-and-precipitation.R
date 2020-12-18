library(dataRetrieval)
library(tidyr)
library(ggplot2)
library(raster)
library(ncdf4)


correlation.test <- function(x, timeseries, method, value) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    if (value=="correlation") {
      return(cor.test(x, timeseries, method=method)$estimate) 
    } else if (value == "p-value") {
      return(cor.test(x, timeseries, method=method)$p.value)
    }
  }
}

pixelwiseCorrelation <- function(brick, timeseries, method) {
  correlation <- raster::calc(brick, fun=function(x){correlation.test(x, timeseries=timeseries, method=method, value="correlation")})
  p.value <- raster::calc(brick, fun=function(x){correlation.test(x, timeseries=timeseries, method=method, value="p-value")})
  stats <- raster::brick(correlation, p.value)
  names(stats) <- c("correlation", "p.value")
  return(stats)
}


LittleRiver <- "01671100"
JamesRiver <- "02037500"
MississippiRiver <- "07374000"
gauge <- JamesRiver

streamflow <- dataRetrieval::readNWISdv(gauge, "00060", "2007-01-01", "2020-12-09")
streamflow <- streamflow[,c("Date", "X_00060_00003")]
colnames(streamflow) <- c("Date", "CFS")
streamflow <- streamflow %>% tidyr::complete(Date = seq(min(Date), max(Date), by="days"))

ggplot2::ggplot(data=streamflow) +
  geom_line(aes(Date, CFS), size=1.5) +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25,face="bold"))

basin <- dataRetrieval::findNLDI(nwis=gauge, find="basin")$basin

plot(basin, lwd=2)

precip_files <- Sys.glob("Data/NOAA-PSL/CPC-Unified-Gauge-Based-Analysis-Of-Daily-Precipitation-Over-CONUS-RT/*")
precip_brick <- raster::brick(lapply(precip_files, FUN=function(x) raster::brick(x)))
precip_brick <- raster::rotate(precip_brick)
precip_basin <- raster::crop(precip_brick, basin)
precip_basin <- raster::mask(precip_basin, basin)

raster::plot(precip_basin$X2020.12.07)

streamflow$Precip <- cellStats(precip_basin, stat="mean")
ggplot2::ggplot(streamflow) +
  geom_line(aes(Date, Precip), size=1.5) +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25,face="bold"))

cor.test(streamflow$CFS, streamflow$Precip, method="pearson")
pearson <- pixelwiseCorrelation(brick=precip_basin, timeseries=streamflow$CFS, method="pearson")
raster::plot(pearson)





library(raster)
library(ncdf4)
library(sf)
library(rasterVis)
library(RColorBrewer)


thumb <- raster::brick("Data/NOAA-PSL/CPC-Unified-Gauge-Based-Analysis-Of-Daily-Precipitation-Over-CONUS-RT/precip.V1.0.2020.nc")
thumb <- raster::rotate(thumb)
dev.off()
raster::plot(thumb$X2020.12.05)

shp_path <- "/Users/jameskitchens/Downloads/USA_Rivers_and_Streams-shp/"
shp_name <- "9ae73184-d43c-4ab8-940a-c8687f61952f2020328-1-r9gw71.0odx9.shp"
shp_file <- paste(shp_path, shp_name, sep="")

# read the shapefile
world_shp <- read_sf(shp_file)
east <- world_shp[which(world_shp$Miles>25),]
world_outline <- as(st_geometry(east), Class="Spatial")

mapTheme <- rasterTheme(region=brewer.pal(8,"Blues"), panel.background = list(col='grey30'))
plt <- levelplot(thumb$X2020.12.05, margin=F, par.settings=mapTheme)
plt + latticeExtra::layer(sp.lines(world_outline, col="gray30", lwd=1))






