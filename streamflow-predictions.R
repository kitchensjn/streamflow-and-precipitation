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






















streamflow_brick <- raster::brick(lapply(names(precip_brick@data), FUN=rasterWithValue, timeseries=streamflow, starting_brick=precip_brick))
combined_stack <- stack(precip_brick, streamflow_brick)
correlationStats_gridCorts <- gridcorts(rasterstack = combined_stack, method = "spearman", type = "both")
basin <- dataRetrieval::findNLDI(nwis=gauge, find="basin")$basin
correlationStats_gridCorts <- raster::crop(correlationStats_gridCorts, basin)
correlationStats_gridCorts <- raster::mask(correlationStats_gridCorts, basin)
plot(correlationStats_gridCorts)




precip_mean <- calc(precip_list, fun=mean, na.rm=TRUE)
precip_sd <- calc(precip_list, fun=sd, na.rm=TRUE)
plot(precip_mean)
plot(precip_sd)


years <- seq(2007, 2020)
precip <- lapply(years, FUN=precipPull, basin=basin)
precip <- do.call(rbind.data.frame, precip)

combined <- merge(streamflow, precip, all=TRUE)

ggplot(data=combined) +
  geom_line(aes(Date, Precipitation*100000), color="green") +
  geom_line(aes(Date, CFS)) +
  scale_x_date(date_breaks="3 month") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))









#library(reshape)
#library(sf)
#library(stringr)
#library(rlist)




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



gridcorts <- function(rasterstack, method, type=c("corel","pval","both")){
  # Values for (layers, ncell, ncol, nrow, method, crs, extent) come straight from the input raster stack
  # e.g. nlayers(rasterstack), ncell(rasterstack)... etc.
  print(paste("Start Gridcorts:",Sys.time()))
  print("Loading parameters")
  layers=nlayers(rasterstack);ncell=ncell(rasterstack);
  ncol=ncol(rasterstack);nrow=nrow(rasterstack);crs=crs(rasterstack);
  extent=extent(rasterstack);pb = txtProgressBar(min = 0, max = ncell, initial = 0)
  print("Done loading parameters")
  mtrx <- as.matrix(rasterstack,ncol=layers)
  empt <- matrix(nrow=ncell, ncol=2)
  print("Initiating loop operation")
  if (type == "corel"){
    for (i in 1:ncell){
      setTxtProgressBar(pb,i)
      if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
        empt[i,1] <- NA 
      } else 
        if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
          empt[i,1] <- NA 
        } else 
          empt[i,1] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$estimate)
    }
    print("Creating empty raster")
    corel <- raster(nrows=nrow,ncols=ncol,crs=crs)
    extent(corel) <- extent
    print("Populating correlation raster")
    values(corel) <- empt[,1]
    print(paste("Ending Gridcorts on",Sys.time()))
    corel
  } 
  else
    if (type == "pval"){
      for (i in 1:ncell){
        setTxtProgressBar(pb,i)
        if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
          empt[i,2] <- NA 
        } else 
          if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
            empt[i,2] <- NA 
          } else 
            empt[i,2] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$p.value)
      }
      pval <- raster(nrows=nrow,ncols=ncol,crs=crs)
      extent(pval) <- extent
      print("Populating significance raster")
      values(pval) <- empt[,2]
      print(paste("Ending Gridcorts on",Sys.time()))
      pval
    }
  else
    if (type == "both"){
      for (i in 1:ncell){
        setTxtProgressBar(pb,i)
        if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
          empt[i,] <- NA 
        } else 
          if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
            empt[i,] <- NA 
          } else {
            empt[i,1] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$estimate) 
            empt[i,2] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$p.value)
          }
      }
      c <- raster(nrows=nrow,ncols=ncol,crs=crs)
      p <- raster(nrows=nrow,ncols=ncol,crs=crs)
      print("Populating raster brick")
      values(c) <- empt[,1]
      values(p) <- empt[,2]
      brk <- brick(c,p)
      extent(brk) <- extent
      names(brk) <- c("Correlation","Pvalue")
      print(paste("Ending Gridcorts on",Sys.time()))
      brk
    }
}

rasterWithValue <- function(x, timeseries, starting_brick) {
  date_layer <- starting_brick[[x]]
  date <- stringr::str_replace_all(stringr::str_replace_all(x, "X", ""), "\\.", "-")
  level <- timeseries[which(timeseries$Date==date),]$CFS
  date_layer[!is.na(date_layer)] <- ifelse(length(level)==1, level, NA)
  date_layer
}
