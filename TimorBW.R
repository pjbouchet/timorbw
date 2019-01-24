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

#' ====================================
# LIBRARIES, SETTINGS ====
#' ====================================

#'--------------------------------------------------------------------------
# Loops through the  list and loads (+ installs if needed) all packages
#'--------------------------------------------------------------------------

# Note: sf requires GDAL, which can be tricky to install on some OS
# For Mac, see: https://stackoverflow.com/questions/44973639/trouble-installing-sf-due-to-gdal


# require(devtools)
# devtools::install_github("seananderson/ggsidekick")
# devtools::install_github("dkahle/ggmap", ref = "tidyup")

pack<-c("tidyverse", # Tidyverse
        "raster", # Raster and GIS operations
        "rgdal", # Working with geospatial data
        "rgeos", # Topology operations on geometries
        "pals", # Colour palettes (parula)
        "sf", # Simple features
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
  p<-pack[i]
  if(as.character(p) %in% rownames(installed.packages())==TRUE){
    library(p, character.only = TRUE)
  }else{install.packages(p)
    library(p, character.only = TRUE)}
}

# Gets the latest install of ggmap
# if(!requireNamespace("devtools")) install.packages("devtools")
# devtools::install_github("dkahle/ggmap", ref = "tidyup", force=TRUE)

# Also, Google has changed its policies. For use after July 2018, need to
# log on to Google Cloud Platform, create project, generate API key and enable billing
# https://stackoverflow.com/questions/52565472/get-map-not-passing-the-api-key-http-status-was-403-forbidden

library(ggmap)

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

#' ====================================
# FUNCTIONS ====
#' ====================================

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
# Cleans column names
#'---------------------------------------------

bw <- bw %>% janitor::clean_names()
gps <- purrr::map(.x = gps, .f = function(x) janitor::clean_names(x))

#'---------------------------------------------
# Formats columns
#'---------------------------------------------

bw$date <- as.Date(bw$date, "%yy-%mm-%dd") # date

bw <- chr2fac(bw) # Converts all character variables to factors
gps <- purrr::map(.x = gps, .f = function(x) chr2fac(x))

#'---------------------------------------------
# Adds year column
#'---------------------------------------------

bw <- bw %>% 
  mutate(year = lubridate::year(date))

#'---------------------------------------------
# Extracts relevant data
#'---------------------------------------------

bw.pts <- bw %>% 
  dplyr::filter(species == "Blue whale") %>% 
  dplyr::filter(obs_type == "Dedicated") %>% 
  dplyr::filter(resighted == "Not a duplicate")

#' =============================
# GIS SHAPEFILES ====
#' =============================

#'---------------------------------------------
# Defines projections
#'---------------------------------------------

CRSll <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
CRSutm <- sp::CRS("+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

#'---------------------------------------------
# Downloads world polygon
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
# Defines study area extent
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
# Crops land
#'---------------------------------------------

timor <- raster::crop(land, bw.ext)
timor_utm <- sp::spTransform(timor, CRSutm)

timor_sf <- sf::st_as_sf(timor)
timor_sf_utm <- sf::st_as_sf(timor_utm)  

#'---------------------------------------------
# Creates spatial lines for each survey day
#'---------------------------------------------

# Need a character string for each date

gps <- purrr::map(.x = gps, .f = function(x) mutate(.data = x, datechr = as.character(date)))

# Generates lines

gps.trks <- purrr::map(gps, ~.x %>% 
             split(.$datechr) %>% 
             purrr::map(., ~ createLines(.)))

# Combines all tracks from 2007/2008

gps.07 <- do.call(rbind, gps.trks$`2007`)
gps.08 <- do.call(rbind, gps.trks$`2008`)

# Necessary for saving shapefiles to disk

gps.07@data$time <- as.character(gps.07@data$time)
gps.08@data$time <- as.character(gps.08@data$time)

# Exports shapefiles

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
# Downloads basemap
#'---------------------------------------------

gmap.timor<-ggmap::get_googlemap(center = c(lon=126,lat=-9.4),
                          zoom=8,
                          size = c(640, 640),
                          scale=2,
                          maptype="roadmap",
                          style = 'feature:road|element:all|visibility:simplified&style=feature:administrative.locality|element:labels|visibility:simplified')

# Quick visualisation
# ggmap(gmap.timor)

#'---------------------------------------------
# Fortifies lines for plotting
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

depth <- raster::raster(file.path("env", "gebco30sec.asc"))
raster::projection(depth) <- CRSll

depth <- raster::mask(x = depth, mask = study.area)
depth <- raster::projectRaster(depth, crs = CRSutm) # Project to UTM

depth <- raster::crop(x = depth, y = extent(study.area_utm))

# Quick plot

plot(depth, col = pals::parula(100))
plot(coast_utm, add=T, col="lightgrey")
  
#'---------------------------------------------
# Seabed slope
#'---------------------------------------------

# Returns values representing the 'rise over run'
# Convert to degrees by taking ATAN ( output ) * 57.29578
  
seabed_slope <- SDMTools::slope(mat = depth, latlon = F)
seabed_slope <- atan(seabed_slope)*180/pi

# Quick plot

plot(seabed_slope, col = pals::parula(100))
plot(coast_utm, add=T, col="lightgrey")

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
                      byid=TRUE)

# To get the nearest distance to any feature, apply min over rows
  
dist_coast <- apply(dc, 1, min)
dist_coast <- cbind(sp::coordinates(rdf), dist_coast)
dist_coast <- raster::rasterFromXYZ(xyz = dist_coast, crs = CRSutm)

plot(dist_coast, col = pals::parula(100))
plot(coast_utm, add = T, col = "lightgrey")
  
  
  
  
  
  
  















