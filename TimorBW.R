#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------
#' 
#'  Distribution and habitat use of blue whales (Balaenoptera musculus sp.) 
#'   in the Timor Trough, south of Timor-Leste during 2007-08
#'
#'
#' Manuscript: Burton et al.
#' Last Update: November 8, 2018
#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------

# Examples of github repositories and how they are structured
# https://github.com/dbarneche/fishEggSize
# https://github.com/cboettig/noise-phenomena

#' ====================================
# LIBRARIES, SETTINGS ====
#' ====================================

#'--------------------------------------------------------------------------
# Loop through the  list and load (+ install if needed) all packages
#'--------------------------------------------------------------------------

# Note: sf requires GDAL, which can be tricky to install on some OS
# For Mac, see: https://stackoverflow.com/questions/44973639/trouble-installing-sf-due-to-gdal

# https://www.fromthebottomoftheheap.net/2018/10/23/introducing-gratia/

# require(devtools)
# devtools::install_github("seananderson/ggsidekick")
# devtools::install_github("dkahle/ggmap", ref = "tidyup")

pack <- c("tidyverse", # Tidyverse
          "raster", # Raster and GIS operations
          "rgdal", # Working with geospatial data
          "rgeos", # Topology operations on geometries
          "pals", # Colour palettes (parula)
          "sf", # Simple features
          "WaveletComp", # Wavelet analysis
          "SDMTools", # Tools for processing data in SDMs
          "lubridate", # Date handling
          "rnaturalearth", # Shapefiles and Earth data
          "ggplot2", # Graphics
          "purrrlyr", # Intersection of purrr/dplyr
          "scales", # For comma format on plots
          "ggsidekick", # Nice ggplot theme by Sean Anderson 
          "rerddap", # R client for working with ERDDAP servers
          "janitor") # Data cleaning 

for (i in 1:length(pack)){
  p <- pack[i]
  if(as.character(p) %in% rownames(installed.packages())==TRUE){
    library(p, character.only = TRUE)
  }else{install.packages(p)
    library(p, character.only = TRUE)}
}

# Get the latest install of ggmap
# if(!requireNamespace("devtools")) install.packages("devtools")
# devtools::install_github("dkahle/ggmap", ref = "tidyup", force=TRUE)

# Also, Google has changed its policies. For use after July 2018, need to
# log on to Google Cloud Platform, create project, generate API key and enable billing
# https://stackoverflow.com/questions/52565472/get-map-not-passing-the-api-key-http-status-was-403-forbidden

library(ggmap)

set.seed(87) # Reproducible results

#'---------------------------------------------
# Set tibble options
#'---------------------------------------------

options(tibble.width = Inf) # All tibble columns shown
options(pillar.neg = FALSE) # No colouring negative numbers
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)

#'---------------------------------------------
# Directories
#'---------------------------------------------

homeDir = getwd()
dataDir = file.path(homeDir, "data")
figDir = file.path(homeDir, "figures") 

#'---------------------------------------------
# Time zone
#'---------------------------------------------

Sys.setenv(TZ = "Australia/West")

#' ====================================
# FUNCTIONS ====
#' ====================================

# Function to plot an input raster over study area using ggplot

make_map <- function(input.raster, 
                        col.ramp = pals::parula(100),
                        show.canyons = FALSE){
  
  # Name of raster as legend title
  
  raster.name <- deparse(substitute(input.raster))
  raster.name <- paste(toupper(substr(raster.name, 1, 1)), 
                      substr(raster.name, 2, nchar(raster.name)), sep="")
  raster.name <- gsub(pattern = "_", replacement = " ", raster.name)
  
  # Ensure that raster is in lat/lon
  
  if(!proj4string(input.raster)==as.character(CRSll)){
    input.raster <- raster::projectRaster(from = input.raster, crs = CRSll)
  }
  
  # Stops execution if basemap does not exist
  
  if(exists("gmap.timor")==FALSE) stop("Missing basemap")
  
  if(show.canyons) canyons.plot <- fortify(canyons)
  
  # Convert to data.frame for plotting
  
  input.raster <- raster::as.data.frame(input.raster, xy = TRUE)
  names(input.raster) <- c("lon", "lat", "z")
  
  ggmap(gmap.timor)+ # basemap
    
    coord_equal() + # Needs to be ### to produce map in right dimension pair
    
    geom_tile(data = input.raster, aes(x = lon, y = lat, fill = z))+ # raster
    
    {if(show.canyons) geom_path(data = canyons.plot, aes(long, lat, group=group))}+
    
    scale_fill_gradientn(colors = col.ramp,
                         na.value = 'transparent',
                         limits = range(input.raster$z),
                         breaks = pretty(input.raster$z),
                         labels = pretty(input.raster$z))+
    
    labs(fill = raster.name)+

    xlab("")+
    ylab("")+
    
    scale_x_continuous(limits = range(xval), 
                       breaks= xval,
                       labels = lab.x, expand=c(0,0))+
    
    scale_y_continuous(limits = range(-yval), 
                       breaks= rev(-yval),
                       labels = rev(lab.y),expand=c(0,0))+
    
    theme_sleek() + # ggsidekick magic happens here
    
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size = 13),
          axis.text.x = element_text(size = 13),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),
          legend.key = element_rect(fill = "transparent"),
          legend.position = c(0.1, 0.25),
          legend.key.width = unit(0.5,"cm"),
          legend.background = element_rect(fill = "transparent", size = 2),
          legend.text = element_text(size = 12, colour = "black"),
          legend.title = element_text(size = 13, face = "bold", colour = "black"),
          legend.key.size = unit(0.75,"cm"))
  
  
}

# Functions to change all character columns to factors and vice versa

chr2fac <- function(y){y %>% dplyr::mutate_if(is.character, as.factor)}
fac2chr <- function(y){y %>% dplyr::mutate_if(is.factor, as.character)}

# Function to create a SpatialLinesDataFrame from a data.frame

createLines <- function(df){

  df.L <- df[, c("longitude", "latitude")] %>%  
    as.matrix(.) %>% 
    raster::spLines(., crs=CRSll)
  
  Sldf <- SpatialLinesDataFrame(df.L, data = df)
  return(Sldf)}

# Function to generate climatologies from errdapp data

get_climatology <- function(dataset,
                            variable.name,
                            time.window = 8, # in days
                            no.years = 15,
                            region.shp = study.area,
                            rmatch = TRUE,
                            raster.dest = depth,
                            date.start = '2017-04-07', 
                            date.end = '2017-05-07'){
  
  message('Retrieving dataset details ...')
  
  # Get information on desired ERDDAP dataset.
  
  datInfo <- rerddap::info(dataset)
  
  message('Generating dates ...')
  
  # Generate date windows
  
  start.dates <- seq(lubridate::ymd(date.start), 
                     lubridate::ymd(date.end), 
                     by = paste0(time.window, ' days'))
  
  if(lubridate::ymd(date.end)-lubridate::ymd(date.start)<time.window){
    
    end.dates <- lubridate::ymd(date.start) + lubridate::days(time.window-1)
    
  }else{
  
    end.dates <- seq(lubridate::ymd(date.start) + lubridate::days(time.window-1), 
                   lubridate::ymd(date.end), 
                   by = paste0(time.window, ' days'))
  }
  
  # Adjust end dates if necessary
  
  if(length(end.dates)<length(start.dates)) end.dates <- c(end.dates, start.dates[length(start.dates)]+ days(time.window-1))
  
  start.dates <- purrr::map(.x = start.dates,
                            .f = ~rev(seq(ymd(.x), length = no.years, by = "-1 year"))) 
  
  end.dates <- purrr::map(.x = end.dates,
                          .f = ~rev(seq(ymd(.x), length = no.years, by = "-1 year"))) 
  # %>%
  #   purrr::map(.x = ., .f = ~as.list(.x))
  
  dates.list <- purrr::map2(.x = start.dates,
                            .y = end.dates,
                            .f = ~tibble(.x, .y))
  
  # Download the data
  
  message('Downloading data ...')
  
  dat <- purrr::map(.x = dates.list,
                    .f = ~split(.x, seq(nrow(.x))) %>% 
                      purrr::map(.x = ., .f = ~ as.data.frame(.x)) %>%
                      purrr::map(.x = ., .f = ~c(.$.x, .$.y)) %>% 
                      purrr::map(.x = ., .f = ~rerddap::griddap(x = datInfo,
                                                                longitude = c(extent(region.shp)[1], extent(region.shp)[2]),
                                                                latitude = c(extent(region.shp)[3], extent(region.shp)[4]),
                                                                time = .x,
                                                                fields = variable.name))) %>% 
    purrr::flatten(.) %>% 
    purrr::map(.x = ., .f = ~.x$data) %>% 
    bind_rows(.) %>% 
    as_tibble(.)
  
  message('Creating climatology ...')
  
  # Calculate climatology
  
  climg <- dat %>% 
    dplyr::group_by(lon, lat) %>% 
    dplyr::summarize_at(.vars = variable.name, .funs = mean)
  
  # Raster extent
  
  rex <- climg[[1]][, c("lon", "lat")]
  coordinates(rex) <- ~ lon + lat
  proj4string(rex) <- CRSll
  rex <- raster::raster(rex)
  
  message('Producing rasters ...')
  
  # Generate rasters (and store in stack)
  
  if(rmatch){
    
    r <- purrr::map(.x = climg,
                    .f = ~ raster::rasterize(x = .x[, c("lon", "lat")], 
                                             y = rex, 
                                             field = .x[,c(variable.name)], 
                                             fun = mean) %>% 
                      raster::projectRaster(from = ., to = raster.dest) %>% 
                      raster::crop(., sp::spTransform(region.shp, CRSutm)) %>% 
                      raster::mask(., sp::spTransform(region.shp, CRSutm))) %>% 
      raster::stack(unlist(.))
    
  }else{
    
    r <- purrr::map(.x = climg,
                    .f = ~ raster::rasterize(x = .x[, c("lon", "lat")], 
                                             y = rex, 
                                             field = .x[,c(variable.name)], 
                                             fun = mean) %>% 
                      raster::projectRaster(from = ., crs = CRSutm) %>% 
                      raster::crop(., sp::spTransform(region.shp, CRSutm)) %>% 
                      raster::mask(., sp::spTransform(region.shp, CRSutm))) %>% 
      raster::stack(unlist(.))
    
  }
  
  return(r)
  
}

#' ====================================
# DATA IMPORT ====
#' ====================================

#'---------------------------------------------
# Sightings
#'---------------------------------------------

bw <- readr::read_csv(file.path(dataDir, "timorbw_data.csv"))

#'---------------------------------------------
# GPS tracks
#'---------------------------------------------

gps <- list('2007'=NULL, '2008'=NULL)
gps$`2007` <- readr::read_tsv(file.path(dataDir, "gps_albacora07.txt"))
gps$`2008` <- readr::read_tsv(file.path(dataDir, "gps_bicuda08.txt"))

#'---------------------------------------------
# Clean column names
#'---------------------------------------------

bw <- bw %>% janitor::clean_names()
gps <- purrr::map(.x = gps, .f = function(x) janitor::clean_names(x))

#'---------------------------------------------
# Format columns
#'---------------------------------------------

bw$date <- as.Date(bw$date, "%yy-%mm-%dd") # date

bw <- chr2fac(bw) # Converts all character variables to factors
gps <- purrr::map(.x = gps, .f = function(x) chr2fac(x))

#'---------------------------------------------
# Add year column
#'---------------------------------------------

bw <- bw %>% 
  dplyr::mutate(year = lubridate::year(date))

#'---------------------------------------------
# Extract relevant data
#'---------------------------------------------

bw.pts <- bw %>% 
  dplyr::filter(species == "Blue whale") %>% 
  dplyr::filter(obs_type == "Dedicated") %>% 
  dplyr::filter(resighted == "Not a duplicate")

#'---------------------------------------------
# Distributon of records in time
#'---------------------------------------------

timeline <- bw %>% 
  dplyr::select(date) %>% 
  dplyr::mutate(top = 1)

timeline.plot <- ggplot(data = timeline, aes(date)) +
geom_ribbon(aes(ymin = 0, ymax = top, 
                colour = as.factor(date)), alpha=0.1) +
  theme(legend.position = "none")


#' =============================
# GIS SHAPEFILES ====
#' =============================

#'---------------------------------------------
# Define projections
#'---------------------------------------------

CRSll <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
CRSutm <- sp::CRS("+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

#'---------------------------------------------
# Download world polygon
#'---------------------------------------------

land <- rnaturalearth::ne_countries(type = 'countries', 
                                    scale = 'large', 
                                    continent = 'asia')

#'---------------------------------------------
# More detailed landmass for Timor leste
#'---------------------------------------------

# Extracted from GSHHS - http://www.soest.hawaii.edu/wessel/gshhg/

coast <- raster::shapefile(file.path("gis", "timor_leste.shp"))
coast_utm <- sp::spTransform(coast, CRSutm)

#'---------------------------------------------
# Define study area extent
#'---------------------------------------------

bw.ext <- as(extent(c(121, 132, -16, -7.5)), "SpatialPolygons")
sp::proj4string(bw.ext) <- CRSll

#'---------------------------------------------
# Survey area
#'---------------------------------------------

study.area <- raster::shapefile(file.path("gis", "study_area.shp"))
sp::proj4string(study.area) <- CRSll

study.area_utm <- sp::spTransform(study.area, CRSutm)

#'---------------------------------------------
# Crop land
#'---------------------------------------------

timor <- raster::crop(land, bw.ext)
timor_utm <- sp::spTransform(timor, CRSutm)

timor_sf <- sf::st_as_sf(timor)
timor_sf_utm <- sf::st_as_sf(timor_utm)  

#'---------------------------------------------
# Submarine canyons
#'---------------------------------------------
# From Harris and Whiteway (2011). Global distribution of large submarine canyons: Geomorphic differences between active and passive continental margins. Marine Geology, 285(1-4): 69-86.
# Canyon types: 1. shelf-incising and river associated; 
# 2. shelf-incising; and 3. blind.

canyons <- raster::shapefile("gis/global_canyons.shp")
canyons <- raster::crop(x = canyons, 
                        y = rgeos::gBuffer(spgeom = study.area, width = 0.2))
canyons_utm <- sp::spTransform(canyons, CRSutm)


# canyons_incising <- canyons_utm[canyons_utm@data$class==2,]

#'---------------------------------------------
# Create spatial lines for each survey day
#'---------------------------------------------

# Need a character string for each date

gps <- purrr::map(.x = gps, .f = function(x) mutate(.data = x, datechr = as.character(date)))

# Generate lines

gps.trks <- purrr::map(gps, ~.x %>% 
             split(.$datechr) %>% 
             purrr::map(., ~ createLines(.)))

# Combine all tracks from 2007/2008

gps.07 <- do.call(rbind, gps.trks$`2007`)
gps.08 <- do.call(rbind, gps.trks$`2008`)

# Necessary for saving shapefiles to disk

gps.07@data$time <- as.character(gps.07@data$time)
gps.08@data$time <- as.character(gps.08@data$time)

# Export shapefiles

# rgdal::writeOGR(obj = gps.07, dsn = file.path("gis"),
#                 layer = "gps_tracks_2007", driver="ESRI Shapefile")
# 
# rgdal::writeOGR(obj = gps.08, dsn = file.path("gis"),
#                 layer = "gps_tracks_2008", driver="ESRI Shapefile")

#' =============================
# MAPPING ====
#' =============================

#'---------------------------------------------
# Register Google API Key
#'---------------------------------------------

# ggmap::register_google(key = "AIzaSyAxvafx12hSaL5pMlJlOeXLBW0G5Bypvxw")

#'---------------------------------------------
# Download basemap
#'---------------------------------------------

gmap.timor <- ggmap::get_googlemap(center = c(lon=126,lat=-9.4),
                          zoom=8,
                          size = c(640, 640),
                          scale=2,
                          maptype="roadmap",
                          style = 'feature:road|element:all|visibility:simplified&style=feature:administrative.locality|element:labels|visibility:simplified')

# Quick visualisation
# ggmap(gmap.timor)

#'---------------------------------------------
# Fortify lines for plotting
#'---------------------------------------------

# Would normally convert to simple features, but transparency not implemented for lines yet

gps.07f <- fortify(gps.07)
gps.08f <- fortify(gps.08)

#'---------------------------------------------
# Axis labels
#'---------------------------------------------

xval <- seq(124.5,127.5,0.5)
yval <- seq(8,10.5,0.5)
lab.x<-c(paste(xval, "°E",sep=""))
lab.y<-c(paste(yval, "°S",sep=""))

#'---------------------------------------------
# ggplot map
#'---------------------------------------------

gg.timor <- ggmap(gmap.timor)+ # basemap

  # GPS tracks
  geom_path(data = gps.07f, aes(long, lat, group = group), alpha = 0.25)+
  geom_path(data = gps.08f, aes(long, lat, group = group), alpha = 0.25)+
  
  coord_equal() + # Needs to be ### to produce map in right dimension pair

  geom_point(data=bw.pts, 
             aes(longitude, latitude,
                 fill=as.factor(year)), pch=21, colour = "black", size = 2.5, alpha = 1)+
  
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"), name = "Year")+
  
  xlab("")+
  ylab("")+
  
  scale_x_continuous(limits = range(xval), 
                     breaks= xval,
                     labels = lab.x, expand=c(0,0))+
  
  scale_y_continuous(limits = range(-yval), 
                     breaks= rev(-yval),
                     labels = rev(lab.y),expand=c(0,0))+

  theme_sleek() + # ggsidekick magic happens here
  
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        legend.key = element_rect(fill = "transparent"),
        legend.position = c(0.1, 0.1),
        legend.background = element_rect(fill = "transparent", size = 2),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 13, face = "bold", colour = "black"),
        legend.key.size = unit(0.75,"cm")) +
  
  ggplot2::guides(fill = guide_legend(override.aes = list(size=6)))
  

#'---------------------------------------------
# North arrow and scale bar
#'---------------------------------------------

# Data argument HAS to match first data set of ggplot/ggmap
# Also requires coord_equal otherwise returns an error about need for Cartesian coordinates

gg.timor +  ggsn::north(data = gps.07f,
                        location = "bottomleft",
                        symbol = 12, 
                        scale = 0.15,
                        anchor = c(x = 126.95, y = -10.2)) +
  
  ggsn::scalebar(x.min = 125,
                 x.max = 127.5,
                 y.min = -10.5,
                 y.max = -9,
                 dist = 25, 
                 height = 0.02,
                 dd2km = F,
                 model="WGS84",
                 anchor = c(x = 127.25, y = -10.3),
                 st.dist = 0.05) +
  
  annotate("text", 
           label = "N", 
           x = 127.025, 
           y = -10, 
           size = 5)
  
# ggsave(filename = file.path(getwd(), "figures", "Figure1.pdf"), 
       # width = 25, 
       # height = 20, 
       # units = "cm")

#' =============================
# RASTERS ====
#' =============================

#'---------------------------------------------
# Bathymetry
#'---------------------------------------------

# From GEBCO 30s gridded dataset (see terms of use for acknowledgements)

bathy <- raster::raster(file.path("env", "gebco30sec.asc"))
proj4string(bathy) <- CRSll

bathy <- raster::projectRaster(bathy, crs = CRSutm) # Project to UTM

depth <- raster::mask(x = bathy, mask = study.area_utm)
depth <- raster::crop(x = depth, y = extent(study.area_utm))

# Quick plot

plot(depth, col = pals::parula(100))
plot(coast_utm, add = T, col = "lightgrey")
plot(study.area_utm, add = T)
plot(canyons_utm, add = T)


#'---------------------------------------------
# Seabed slope
#'---------------------------------------------

# Returns values representing the 'rise over run'
# Convert to degrees by taking ATAN ( output ) * 57.29578
  
seabed_slope <- SDMTools::slope(mat = bathy, latlon = F) # Must be in proj units, not latlon
seabed_slope <- atan(seabed_slope)*180/pi

seabed_slope <- raster::mask(x = seabed_slope, mask = study.area_utm)
seabed_slope <- raster::crop(x = seabed_slope, y = extent(study.area_utm))

# Quick plot

plot(seabed_slope, col = pals::parula(100))
plot(coast_utm, add=T, col="lightgrey")
plot(study.area_utm, add=T)

#'---------------------------------------------
# Distance to coast
#'---------------------------------------------

# Use rgeos to compute point-polygon distances
# This outputs a matrix with a column for each feature in spgeom1.
  
rdf <- raster::as.data.frame(depth, xy = T)
rdf <- rdf[complete.cases(rdf),] # Removes NAs
rdf <- SpatialPoints(rdf[,1:2], CRSutm) # Converts to SpatialPts

# Cartesian distance

dc <-  rgeos::gDistance(spgeom1 = coast_utm, 
                      spgeom2 = rdf, 
                      byid = TRUE)

# To get the nearest distance to any feature, apply min over rows
  
dist_coast <- apply(dc, 1, min)
dist_coast <- cbind(sp::coordinates(rdf), dist_coast)
dist_coast <- raster::rasterFromXYZ(xyz = dist_coast, crs = CRSutm)

plot(dist_coast, col = pals::parula(100))
plot(coast_utm, add = T, col = "lightgrey")
plot(study.area_utm, add = T)

#'---------------------------------------------
# Distance to nearest canyon
#'---------------------------------------------

# All canyons included

dcanyons <-  rgeos::gDistance(spgeom1 = canyons_utm, 
                              spgeom2 = rdf, 
                              byid = TRUE)

# To get the nearest distance to any feature, apply min over rows

dist_canyons <- apply(dcanyons, 1, min)
dist_canyons <- cbind(sp::coordinates(rdf), dist_canyons)
dist_canyons <- raster::rasterFromXYZ(xyz = dist_canyons, crs = CRSutm)

plot(dist_canyons, col = pals::parula(100))
plot(coast_utm, add = T, col = "lightgrey")
plot(study.area_utm, add=T)
plot(canyons_utm, add = T)

# Only shelf-incising canyon

dcanyons.incising <-  rgeos::gDistance(spgeom1 = canyons_utm[canyons_utm@data$class == 2,], spgeom2 = rdf, byid = TRUE)

# To get the nearest distance to any feature, apply min over rows

dist_incising <- apply(dcanyons.incising, 1, min)
dist_incising <- cbind(sp::coordinates(rdf), dist_incising)
dist_incising <- raster::rasterFromXYZ(xyz = dist_incising, crs = CRSutm)

plot(dist_incising, col = pals::parula(100))
plot(coast_utm, add = T, col = "lightgrey")
plot(study.area_utm, add=T)
plot(canyons_utm, add = T)
  
#'---------------------------------------------
# Sea Surface Temperature (SST)
#'---------------------------------------------

# https://cran.r-project.org/web/packages/rerddap/vignettes/Using_rerddap.html

# First determine the temporal variability of the system.
# MUR (Multi-scale Ultra-high Resolution) is an analyzed SST product at 0.01-degree resolution going back to 2002.

sstInfo <- rerddap::info('jplMURSST41')

# Define a location at which to assess temperature time series

buoy <- sp::SpatialPoints(coords = cbind(126.25, -9.55))
plot(study.area); points(buoy) # Quick plot

# Retrieve daily sst between 2002 and 2018

sst.buoy <- rerddap::griddap(x = sstInfo, 
                             longitude = rep(coordinates(buoy)[1], 2), 
                             latitude = rep(coordinates(buoy)[2], 2), 
                             time = c('2002-06-01','2018-12-31'),
                             fields = 'analysed_sst')

# Extract date/time and sst values

sst.buoy <- list(erddap = sst.buoy)

sst.buoy$data <- tibble(date = as.Date(sst.buoy$erddap$data$time,
                                       origin = '1970-01-01', tz = "GMT"),
                                       sst = sst.buoy$erddap$data$analysed_sst)

# Plot time series

pdf("figures/Figure-S1a.pdf", height = 4, width = 10)
par(las=1)
plot.ts(ts(sst.buoy$data$sst, frequency = 365, start = c(2002,1)),
        xlab = NA, ylab = "SST (°C)", axes = FALSE)
time.labels <- seq(as.Date('2002-06-01'), as.Date('2018-12-31'), by = "12 months") 
axis(2)
axis(1, labels = lubridate::year(time.labels), 
     at = seq(from = 2002, by = 1, length.out = length(time.labels)))
box()
box()
dev.off() 


# Perform wavelet analysis

wavelet.timor <- WaveletComp::analyze.wavelet(my.data = sst.buoy$data,
                                              my.series = "sst",
                                              loess.span = 0,
                                              dt = 1, 
                                              lowerPeriod = 1,
                                              upperPeriod = 1500,
                                              make.pval = FALSE,
                                              n.sim = 100,
                                              date.format = "%Y-%m-%d")

# Plot wavelet power spectrum

pdf("figures/Figure-S1b.pdf", height = 5, width = 8)
WaveletComp::wt.image(wavelet.timor, 
                      exponent = 1,
                      plot.coi = TRUE,
                      color.key = "interval", 
                      color.palette = "rev(pals::parula(n.levels))",
                      n.levels = 100,
                      show.date = TRUE,
                      timelab = "",
                      spec.time.axis = list(at = time.labels, 
                                            labels = lubridate::year(time.labels)),
                      periodlab = "",
                      legend.params = list(n.ticks = 15, 
                                           label.digits = 1))
dev.off()

# Wavelet power averages across time

WaveletComp::wt.avg(wavelet.timor, siglvl = 0.001, sigcol = "blue")

# SST climatology

sst_climg <- get_climatology(dataset = 'jplMURSST41', # MUR high resolution SST
                            variable.name = 'analysed_sst',
                            time.window = 8, # 8-day average
                            no.years = 15, # Over a period of 15 years
                            region.shp = study.area,
                            rmatch = TRUE, # Resample output rasters
                            raster.dest = depth, # To the same resolution as depth raster
                            date.start = '2017-04-07', # Start date of period of interest
                            date.end = '2017-05-07') # End date of period of interest

#'---------------------------------------------
# Fronts
#'---------------------------------------------
  
# https://cran.r-project.org/web/packages/imager/vignettes/gettingstarted.html#example-2-edge-detection
# https://www.rdocumentation.org/packages/wvtool/versions/1.0/topics/edge.detect
# https://rdrr.io/github/galuardi/boaR/man/boaR-package.html
# https://github.com/LuisLauM/grec
  
  
  
  















