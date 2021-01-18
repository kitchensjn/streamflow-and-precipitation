library(raster)
library(ncdf4)
library(sf)
library(rasterVis)
library(RColorBrewer)


thumb <- raster::brick("Data/NOAA-PSL/CPC-Unified-Gauge-Based-Analysis-Of-Daily-Precipitation-Over-CONUS-RT/precip.V1.0.2020.nc")
thumb <- raster::rotate(thumb)
dev.off()
raster::plot(thumb$X2020.12.05)

shp_path <- "Data/Shapefiles/USA_Rivers_and_Streams-shp/"
shp_name <- "9ae73184-d43c-4ab8-940a-c8687f61952f2020328-1-r9gw71.0odx9.shp"
shp_file <- paste(shp_path, shp_name, sep="")

# read the shapefile
world_shp <- read_sf(shp_file)
east <- world_shp[which(world_shp$Miles>25),]
world_outline <- as(st_geometry(east), Class="Spatial")

mapTheme <- rasterTheme(region=brewer.pal(8,"Blues"), panel.background = list(col='grey30'))
plt <- levelplot(thumb$X2020.12.05, margin=F, par.settings=mapTheme)
plt + latticeExtra::layer(sp.lines(world_outline, col="gray30", lwd=1))