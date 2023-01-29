##%######################################################################################%##
#                                                                                       #
####           Evidence of likely foraging by pygmy blue whales                         ####
####     in the Timor Trough during the late austral winter and early austral spring    ####
#                                                                                      #
##%######################################################################################%##

# Libraries and global settings ------------------------------------------

#'---------------------------------------------
# Load (+ install if needed) all packages
#'---------------------------------------------

# Note: sf requires GDAL, which can be tricky to install on some OS
# For Mac, see: https://stackoverflow.com/questions/44973639/trouble-installing-sf-due-to-gdal

# require(remotes)
# remotes::install_github("dkahle/ggmap", ref = "tidyup")
# remotes::install_github("beckyfisher/FSSgam_package")
# remotes::install_github("seananderson/ggsidekick")

pacman::p_load(tidyverse, # Tidyverse
               raster, # Raster and GIS operations
               geoR, # Geostatistical analysis including variogram-based, likelihood-based and Bayesian methods
               rgdal, # Working with geospatial data
               rgeos, # Topology operations on geometries
               pals, # Colour palettes (parula)
               plyr, # Split-apply-combine paradigm
               tidyr, # Tidy messy data
               sf, # Simple features
               terra, # Manipulate geographic (spatial) data in "raster" and "vector" form
               ncf, # Spatial covariance functions
               xtable, # Export tables to laTeX
               mgcv, # Generalised Additive Models (GAMs)
               WaveletComp, # Wavelet analysis
               GGally, # ggplot2 extensions
               SDMTools, # Tools for processing data in SDMs
               lubridate, # Date handling
               ggplot2, # Graphics
               ggcorrplot, # Correlation graph
               ggsidekick, # Simple theme for ggplot2, 
               meteo, # Spatio-Temporal Analysis and Mapping of Meteorological Observations 
               rworldmap, # Mapping global data
               gstat, # Spatio-Temporal Geostatistical Modelling, Prediction and Simulation
               purrrlyr, # Intersection of purrr/dplyr
               scales, # For comma format on plots
               cowplot, # For combining multiple ggplots
               rerddap, # R client for working with ERDDAP servers
               sp, # Classes and Methods for Spatial Data
               smoothr, # Smooth and Tidy Spatial Features in R
               pROC, # Display and Analyze ROC Curves
               ROCR, # Visualizing the Performance of Scoring Classifiers
               janitor, # Data cleaning 
               FSSgam) # All subsets GAMs

# Note that Google has recently changed its policies. For use after July 2018, need to
# log on to Google Cloud Platform, create a project, generate API key and enable billing
# https://stackoverflow.com/questions/52565472/get-map-not-passing-the-api-key-http-status-was-403-forbidden

# library(ggmap)

set.seed(20221202) # Set the random seed (reproducible results)

source("R/Burton_etal_Functions.R")

#'---------------------------------------------
# Set tibble options
#'---------------------------------------------

options(tibble.width = Inf) # All tibble columns shown
options(pillar.neg = FALSE) # No colouring negative numbers
options(pillar.subtle = TRUE)
options(pillar.sigfig = 5)

#'---------------------------------------------
# Set time zone
#'---------------------------------------------

Sys.setenv(TZ = "Australia/West")

# Data preparation -------------------------------------------------------------

# Define coordinate projection systems
CRSll <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
CRSutm <- sp::CRS("+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

#'---------------------------------------------
# GPS tracks
#'---------------------------------------------

gps <- list('2007' = NULL, '2008' = NULL)

gps$`2007` <- readr::read_csv("data/gps_albacora07_corrected.csv") %>%
  dplyr::mutate(Date = lubridate::dmy(Date))

gps$`2008` <- readr::read_tsv("data/gps_trks_bicuda_08.txt") %>%
  dplyr::mutate(Date = lubridate::mdy(Date))

gps <- purrr::map(.x = gps, .f = function(x) janitor::clean_names(x))
gps <- purrr::map(.x = gps, .f = function(x) chr2fac(x))

#'---------------------------------------------
# Filter data by month
#'---------------------------------------------

# Focus on the austral winter and early austral spring
gps.dates <- purrr::map(.x = gps, 
                        .f = ~.x %>% 
                          dplyr::filter(lubridate::month(date) <= 9) %>% 
                          dplyr::pull(date) %>% unique(.)) %>% 
  Reduce(c, .) %>% sort(.)

gps <- purrr::map(.x = gps, .f = ~dplyr::filter(.x, date %in% gps.dates))

#'---------------------------------------------
# Import sightings data
#'---------------------------------------------

timor.dat <- readr::read_csv("data/timorbw_data.csv") %>%
  janitor::clean_names() %>%
  dplyr::mutate(date = as.Date(date, "%yy-%mm-%dd"))

# Convert all character variables to factors
timor.dat <- chr2fac(timor.dat) 

# Add a year column
timor.dat <- timor.dat %>% 
  dplyr::mutate(year = lubridate::year(date))

#'---------------------------------------------
# Seasons
#'---------------------------------------------

timor.dat <- timor.dat %>% 
  dplyr::mutate(month = lubridate::month(date)) 

# Number of sightings (all species) per month
timor.dat %>% dplyr::pull(month) %>% table()

timor.dat <- timor.dat %>% 
  dplyr::mutate(season = as.factor(ifelse(month %in% 6:8, "Winter", 
                                          ifelse(month %in% 9:11, "Spring", "Summer"))))

#'---------------------------------------------
# Extract and filter blue whale data
#'---------------------------------------------

# Note: 1 sighting made on September 1st, 2008 when end of ops and transit back to Dili does not have an associated GPS track.
bw <- timor.dat %>% 
  dplyr::filter(species == "Blue whale", obs_type == "Dedicated", resighted == "Not a duplicate") %>%
  dplyr::filter(date %in% gps.dates) %>%
  dplyr::select(-species, -obs_type, -resighted, -block, -group) %>%
  add_xy()

bw %>% dplyr::pull(month) %>% table(.)

# GIS layers ---------------------------------------------------------------------

# Import shapefiles

# Extracted from GSHHS - http://www.soest.hawaii.edu/wessel/gshhg/
coast <- raster::shapefile(file.path("gis", "timor_leste.shp"))
coast_utm <- sp::spTransform(coast, CRSutm)

# Define regional extent
bw.ext <- as(extent(c(121, 132, -16, -7.5)), "SpatialPolygons")
sp::proj4string(bw.ext) <- CRSll

# Survey area
study.area <- raster::shapefile(file.path("gis", "study_area.shp"))
sp::proj4string(study.area) <- CRSll
extent(study.area)
study.area_utm <- sp::spTransform(study.area, CRSutm)

# Crop land
timor <- raster::crop(coast, bw.ext)
timor_utm <- sp::spTransform(timor, CRSutm)

#'---------------------------------------------
# Submarine canyons
#'---------------------------------------------

# From Harris and Whiteway (2011). Marine Geology, 285(1-4): 69-86.
# Canyon types: 1. shelf-incising and river associated; 2. shelf-incising; and 3. blind.

canyons <- raster::shapefile("gis/global_canyons.shp")
canyons <- raster::crop(x = canyons, y = rgeos::gBuffer(spgeom = study.area, width = 0.2))
canyons_utm <- sp::spTransform(canyons, CRSutm)

#'---------------------------------------------
# Identify effort segments for each survey day
#'---------------------------------------------

gps$`2007`$timediff <- gps$`2008`$timediff <- NA
gps$`2007`$line_id <- gps$`2008`$line_id <-NA

# Calculate time difference between sequential points
# Effort segments are identified as sets of GPS positions less than 5 minutes apart

gps$`2007`$timediff[1] <- gps$`2008`$timediff[1] <- 0

for(i in 2:nrow(gps$`2007`)){
  gps$`2007`$timediff[i] <- as.numeric(lubridate::as.duration(lubridate::interval(gps$`2007`$time[i-1], gps$`2007`$time[i])), "minutes")
}

for(i in 2:nrow(gps$`2008`)){
  gps$`2008`$timediff[i] <- as.numeric(lubridate::as.duration(lubridate::interval(gps$`2008`$time[i-1], gps$`2008`$time[i])), "minutes")
}

# Assign identifier to each segment of effort

# 2007

out <- split(gps$`2007`, f = gps$`2007`$date)

for(j in 1:length(out)){
  
  counter <- 1
  x <- out[[j]]
  # x$id <- 1:nrow(x)
  x <- x %>% dplyr::arrange(date, time)
  for(i in 1:nrow(x)){
    if(i==1){ x$line_id[i] <- paste0(as.character(x$date[i]), "-L", 1)
    } else {
      if(abs(x$timediff[i])>5){
        counter <- counter + 1
        x$line_id[i] <- paste0(as.character(x$date[i]), "-L", counter)
      }else{
        x$line_id[i] <- x$line_id[i-1]
      }
    }
  }
  x$line_id <- as.factor(x$line_id)
  out[[j]] <- x  
} # End for loop

gps$`2007` <- do.call(rbind, out)

# 2008

out <- split(gps$`2008`, f = gps$`2008`$date)

for(j in 1:length(out)){
  
  counter <- 1
  x <- out[[j]]
  
  for(i in 1:nrow(x)){
    if(i==1){ x$line_id[i] <- paste0(as.character(x$date[i]), "-L", 1)
    }else{
      if(x$timediff[i]>5){
        counter <- counter + 1
        x$line_id[i] <- paste0(as.character(x$date[i]), "-L", counter)
      }else{
        x$line_id[i] <- x$line_id[i-1]
      }
    }
  }
  x$line_id <- as.factor(x$line_id)  
  out[[j]] <- x  
} # End for loop

gps$`2008` <- do.call(rbind, out)

# Need a character string for each date

gps <- purrr::map(.x = gps, .f = function(x) mutate(.data = x, line_id = as.character(line_id)))
gps.all <- do.call(rbind, gps)

#'---------------------------------------------
# Create line objects for each on-effort segment
#'---------------------------------------------

oneffort_sp <- gps.all %>% 
  split(x = ., f = .$date) %>% 
  purrr::map(., ~.x %>% 
               split(x = ., f = .$line_id) %>% 
               purrr::map(., ~ createLines(.)))

# Project to UTM
oneffort_sp <- purrr::map_depth(
  .x = oneffort_sp, .depth = 2,
  .f = ~ sp::spTransform(.x, CRSobj = CRSutm)
) %>% purrr::map(.x = ., .f = ~ do.call(rbind, .))

# Discard short segments
oneffort_sp <- purrr::map(.x = oneffort_sp, .f = ~smoothr::drop_crumbs(x =.x, threshold = 250))

gps.lines <- list(`2007` = do.call(rbind, oneffort_sp[which(grepl("2007", names(oneffort_sp)))]),
                  `2008` = do.call(rbind, oneffort_sp[which(grepl("2008", names(oneffort_sp)))]))

#'---------------------------------------------
# Create 3 km buffers around tracks
#'---------------------------------------------

bufferonland <- purrr::map(.x = oneffort_sp, 
                           .f = ~ rgeos::gBuffer(spgeom = .x, 
                                                 byid = TRUE,
                                                 width = 3000, 
                                                 capStyle = "ROUND", 
                                                 joinStyle = "ROUND"))

# Clip buffers overlapping land
buffer_pol <- purrr::map(.x = bufferonland, .f = ~{
  p1 <- raster::erase(.x, timor_utm)
  raster::intersect(p1, study.area_utm)})

# plot(study.area_utm, axes = TRUE)
# plot(timor_utm, col = "grey", add = TRUE)
# purrr::walk(.x = bufferonland, .f = ~plot(.x, col = "yellow", add = TRUE))
# purrr::walk(.x = buffer_pol, .f = ~plot(.x, col = "green", add = TRUE))
# 
# plot(study.area_utm, axes = TRUE, xlim = c(930000, 950000), ylim = c(8900000, 8920000))
# purrr::walk(.x = bufferonland, .f = ~plot(.x, col = "yellow", add = TRUE))
# plot(study.area_utm, axes = TRUE, xlim = c(920000, 950000), ylim = c(8900000, 8950000), add = T)
# purrr::walk(.x = buffer_pol, .f = ~plot(.x, col = "green", add = TRUE))

# Environmental data ------------------------------------------------------

#'---------------------------------------------
# Bathymetry ====
#'---------------------------------------------

# From GEBCO 15s gridded dataset (see terms of use for acknowledgements)
# https://www.gebco.net/data_and_products/gridded_bathymetry_data/
bathy <- raster::raster(file.path("env", "gebco15sec.tif"))
proj4string(bathy) <- CRSll
bathy <- raster::projectRaster(bathy, crs = CRSutm) # Project to UTM
depth <- raster::mask(x = bathy, mask = study.area_utm)
depth <- raster::crop(x = depth, y = extent(study.area_utm))

#'---------------------------------------------
# Seabed slope ====
#'---------------------------------------------

# Values represent the 'rise over run'
# Convert to degrees by taking ATAN ( output ) * 57.29578

seabed_slope <- SDMTools::slope(mat = bathy, latlon = F) # Must be in proj units, not latlon
seabed_slope <- atan(seabed_slope)*180/pi
seabed_slope <- raster::mask(x = seabed_slope, mask = study.area_utm)
seabed_slope <- raster::crop(x = seabed_slope, y = extent(study.area_utm))

#'---------------------------------------------
# Distance to nearest canyon ====
#'---------------------------------------------

rdf <- raster::as.data.frame(depth, xy = T)
rdf <- rdf[complete.cases(rdf),] # Removes NAs
rdf <- sp::SpatialPoints(rdf[,1:2], CRSutm) # Converts to SpatialPts

# All canyons included

dcanyons <-  rgeos::gDistance(spgeom1 = canyons_utm, 
                              spgeom2 = rdf, 
                              byid = TRUE)

# To get the nearest distance to any feature, apply min over rows

dist_canyons <- apply(dcanyons, 1, min)
dist_canyons <- cbind(sp::coordinates(rdf), dist_canyons)
dist_canyons <- raster::rasterFromXYZ(xyz = dist_canyons, crs = CRSutm)
dist_canyons <- dist_canyons/1000 # Convert to km

# Only shelf-incising canyon

dcanyons.incising <- rgeos::gDistance(spgeom1 = canyons_utm[canyons_utm@data$class == 2,], spgeom2 = rdf, byid = TRUE)

# To get the nearest distance to any feature, apply min over rows

dist_incising <- apply(dcanyons.incising, 1, min)
dist_incising <- cbind(sp::coordinates(rdf), dist_incising)
dist_incising <- raster::rasterFromXYZ(xyz = dist_incising, crs = CRSutm)
dist_incising <- dist_incising/1000 # Convert to km

rm(rdf)

#'---------------------------------------------
# Sea Surface Temperature (SST) ====
#'---------------------------------------------

# https://cran.r-project.org/web/packages/rerddap/vignettes/Using_rerddap.html

# First determine the temporal variability of the system.
# MUR (Multi-scale Ultra-high Resolution) is an analyzed SST product 
# at 0.01-degree resolution going back to 2002.

sstInfo <- rerddap::info('jplMURSST41')

# Define a location at which to assess temperature time series (virtual buoy)

buoy <- sp::SpatialPoints(coords = cbind(126.25, -9.55))
plot(study.area); points(buoy) # Quick plot

# attach('someFile.RData'); someObj <- someObj; detach('file:someFile.RData')

# Retrieve daily sst between 2002 and 2020

sst.buoy <- rerddap::griddap(datasetx = sstInfo, 
                             longitude = rep(coordinates(buoy)[1], 2), 
                             latitude = rep(coordinates(buoy)[2], 2), 
                             time = c('2002-06-02','2020-12-31'),
                             fields = 'analysed_sst')

# Extract date/time and sst values

sst.buoy <- list(erddap = sst.buoy)

sst.buoy$data <- tibble(date = as.Date(sst.buoy$erddap$data$time, origin = '1970-01-01', tz = "GMT"),
                        sst = sst.buoy$erddap$data$analysed_sst)

# Perform wavelet analysis

wavelet.timor <- WaveletComp::analyze.wavelet(my.data = sst.buoy$data,
                                              my.series = "sst",
                                              loess.span = 0,
                                              dt = 1, 
                                              method = "Fourier.rand",
                                              lowerPeriod = 1,
                                              upperPeriod = 1500,
                                              make.pval = FALSE,
                                              n.sim = 100,
                                              date.format = "%Y-%m-%d")

# Wavelet power averages across time

WaveletComp::wt.avg(wavelet.timor)

# Generate SST climatology

sst_climg <- get_wkclimatology(remote.dataset = 'jplMURSST41', # MUR high resolution SST
                               variable.name = 'analysed_sst',
                               t.year = 2007,
                               start.month = 1, 
                               end.month = 12,
                               no.years = 10, # Over a period of 10 years
                               climg.direction = 'centered', 
                               region.shp = study.area,
                               rmatch = TRUE, # Resample output rasters
                               raster.dest = depth, # To the same resolution as depth raster
                               error400 = FALSE) 

# Check that values exist

sst_climg[[1]]@data@values %>% summary()

# Quick plot

sst_climg %>% raster::stack(.) %>% 
  raster::calc(x = ., mean) %>% 
  make_map(input.raster = ., plot.title = 'SST')

#'---------------------------------------------
# Chlorophyll-a ====
#'---------------------------------------------

# Aqua MODIS, NPP, L3SMI, Global, 4km, Science Quality, 2003-present (1 Day Composite)

chlaInfo <- rerddap::info('erdMH1chlamday')

# Retrieve daily Chl-a

chla.buoy <- rerddap::griddap(datasetx = chlaInfo, 
                              longitude = rep(coordinates(buoy)[1], 2), 
                              latitude = rep(coordinates(buoy)[2], 2), 
                              time = c('2003-01-01','2020-12-31'),
                              fields = 'chlorophyll')

# Extract date/time and chla values

chla.buoy <- list(erddap = chla.buoy)

chla.buoy$data <- tibble(date = as.Date(chla.buoy$erddap$data$time, origin = '1970-01-01', tz = "GMT"),
                         chla = chla.buoy$erddap$data$chlorophyll)

# Create climatology

# Mean chlorophyll 

chla_climg <- get_wkclimatology(remote.dataset = 'erdMH1chla1day', 
                                variable.name = 'chlorophyll',
                                t.year = 2007,
                                start.month = 1, 
                                end.month = 12,
                                no.years = 10,
                                summary.method =  "mean",
                                climg.direction = 'centered', 
                                region.shp = study.area,
                                rmatch = TRUE,
                                raster.dest = depth,
                                error400 = TRUE)

# Max chlorophyll

chlamax_climg <- get_wkclimatology(remote.dataset = 'erdMH1chla1day', 
                                   variable.name = 'chlorophyll',
                                   t.year = 2007,
                                   start.month = 1, 
                                   end.month = 12,
                                   no.years = 10,
                                   summary.method = "max",
                                   climg.direction = 'centered', 
                                   region.shp = study.area,
                                   rmatch = TRUE,
                                   raster.dest = depth,
                                   error400 = TRUE)

# Check that values exist

chla_climg[[1]]@data@values %>% summary()
chlamax_climg[[1]]@data@values %>% summary()

# Take the log10 of chla values

chla_climg <- purrr::map(.x = chla_climg, .f = ~log10(.x))
chlamax_climg <- purrr::map(.x = chlamax_climg, .f = ~log10(.x))

# Fill gaps in remote-sensed products

chla_climg_filled <- purrr::map(.x = chla_climg, .f = ~meteo::rfillspgaps(.x, study.area_utm))
chlamax_climg_filled <- purrr::map(.x = chlamax_climg, .f = ~meteo::rfillspgaps(.x, study.area_utm))

#'---------------------------------------------
# Frequency of Chlorophyll Peak Index (FCPI) ====
#'---------------------------------------------

fcpi <- calc_fcpi(chla.climg = chla_climg_filled, log.10 = FALSE)
fcpi$fcpi <- raster::resample(fcpi$fcpi, depth)

# Pseudo-absences ------------------------------------------------

# Summary of survey effort
survey_time <- gps.all %>% 
  dplyr::group_by(date, line_id) %>%
  dplyr::summarise(time_on_line = difftime(dplyr::last(time),dplyr::first(time), units = "mins")) %>%
  dplyr::ungroup()

survey_time <- survey_time %>% 
  dplyr::group_by(date) %>%
  dplyr::summarise(time_on_effort = sum(time_on_line)) %>%
  dplyr::ungroup()

hist(as.numeric(survey_time$time_on_effort)/60, xlab="time ON-effort (hours)", main = "")
summary(as.numeric(survey_time$time_on_effort)/60)

Total.on_effort <- sum(as.numeric(survey_time$time_on_effort))/60

# Length of effort chunks
L <- purrr::map(.x = oneffort_sp, .f = ~rgeos::gLength(.x, byid = TRUE))

# Number of pseudo-absence points 
N.pseudo <- 1000

PA <- lapply(seq_along(buffer_pol), function(x){ 
  
  day <- names(buffer_pol)[x]

  # Filter sightings by day
  bw.sights <- bw %>% dplyr::filter(date == day) 
  N <- nrow(bw.sights)
  
  # If any sightings were made, create 3 km buffer around them
  
  if(N > 0){
    bw.sights.utm <- sp::SpatialPoints(coords = bw.sights[, c("x", "y")], proj4string = CRSutm)
    exclusion.zones <- rgeos::gBuffer(spgeom = bw.sights.utm, byid = TRUE, width = 3000)
    dat.utm <- purrr::map(.x = 1:N, .f = ~sp::spsample(x = exclusion.zones[.x,], n = 1, type = "random")) %>%
      do.call(rbind, .)
  } else {
    dat.utm <-  NULL
    exclusion.zones <- NULL
  }
  
  # Time on effort that day
  t.on <- survey_time$time_on_effort[x]
  
  # Generate a regular grid of points separated by 3 km
  spts <- sp::spsample(buffer_pol[[x]], type = "regular", iter = 1000, cellsize = 3000)
  spts <- sp::spTransform(spts, CRSutm)
  
  # Prevent background points from being sampled within exclusion zones around presence points
  if(N > 0) spts <- rgeos::gDifference(spts, exclusion.zones)
  
  # Randomly sample points
  spts.out <- sample(spts, round((N.pseudo/Total.on_effort)*t.on/60,0), replace = FALSE) 
  
  return(list(absence = spts.out, presence = dat.utm,
              exclusion = exclusion.zones, 
              buffer = buffer_pol[[x]],
              segment = oneffort_sp[[x]]))
})

# Plot a random day to check
check_points(PA, 20)

# Format absence points
pseudo_sp <- purrr::map(PA, "absence") %>% purrr::set_names(nm = names(buffer_pol))

bw.abs <- pseudo_sp %>% 
  do.call(rbind,.) %>% 
  as_tibble(.)

bw.abs$date <- purrr::map_dbl(.x = pseudo_sp, .f = ~length(.x)) %>% 
  purrr::map2(.x = names(.), .y = ., .f = ~rep(.x, .y)) %>% do.call(c, .)

names(bw.abs) <- c("x","y", "date")
bw.abs$date <- lubridate::as_date(bw.abs$date)

# Find week corresponding to each date -------------------------------------------------------

bw.pres <- bw %>% 
  dplyr::mutate(climgweek = purrr::map(.x = bw$date, .f = ~find_week(.x)) %>% do.call("c", .),
  climgweek.lag1 = purrr::map(.x = bw$date, .f = ~find_week(.x, week.lag = 1)) %>% do.call("c", .),
  climgweek.lag2 = purrr::map(.x = bw$date, .f = ~find_week(.x, week.lag = 2)) %>% do.call("c", .),
  climgweek.lag3 = purrr::map(.x = bw$date, .f = ~find_week(.x, week.lag = 3)) %>% do.call("c", .),
  climgweek.lag4 = purrr::map(.x = bw$date, .f = ~find_week(.x, week.lag = 4)) %>% do.call("c", .))

bw.abs <- bw.abs %>% 
  dplyr::mutate(climgweek = purrr::map(.x = bw.abs$date, .f = ~find_week(.x)) %>% do.call("c", .),
  climgweek.lag1 = purrr::map(.x = bw.abs$date, .f = ~find_week(.x, week.lag = 1)) %>% do.call("c", .),
  climgweek.lag2 = purrr::map(.x = bw.abs$date, .f = ~find_week(.x, week.lag = 2)) %>% do.call("c", .),
  climgweek.lag3 = purrr::map(.x = bw.abs$date, .f = ~find_week(.x, week.lag = 3)) %>% do.call("c", .),
  climgweek.lag4 = purrr::map(.x = bw.abs$date, .f = ~find_week(.x, week.lag = 4)) %>% do.call("c", .))

# Extract covariate values at background points ------------------------------------------------

static.env <- raster::stack(depth, seabed_slope, dist_canyons, dist_incising, fcpi$fcpi)
names(static.env) <- c('depth', 'slope', 'dcanyons', 'dcanyonsInc', 'fcpi')

bw.pres <- dplyr::bind_cols(bw.pres, tibble::as_tibble(raster::extract(x = static.env, y = bw.pres[, c('x', 'y')])))
bw.abs <- dplyr::bind_cols(bw.abs, tibble::as_tibble(raster::extract(x = static.env, y = bw.abs[, c('x', 'y')])))

# Retrieve values of dynamic covariates
bw.pres <- bw.pres %>% 
  dplyr::mutate(sst = extract_climg(dat = ., climg = sst_climg, var.name = "climgweek"),
                chla = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek"),
                chla.max = extract_climg(dat = ., climg = chlamax_climg_filled, var.name = "climgweek")) %>% 
  dplyr::mutate(chla.lag1 = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek.lag1"),
                chla.lag2 = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek.lag2"),
                chla.lag3 = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek.lag3"),
                chla.lag4 = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek.lag4")) %>% 
  dplyr::mutate(chla.max.lag1 = extract_climg(dat = ., climg = chlamax_climg_filled, var.name = "climgweek.lag1"),
                chla.max.lag2 = extract_climg(dat = ., climg = chlamax_climg_filled, var.name = "climgweek.lag2"),
                chla.max.lag3 = extract_climg(dat = ., climg = chlamax_climg_filled, var.name = "climgweek.lag3"),
                chla.max.lag4 = extract_climg(dat = ., climg = chlamax_climg_filled, var.name = "climgweek.lag4"))

bw.abs <- bw.abs %>% 
  dplyr::mutate(sst = extract_climg(dat = ., climg = sst_climg, var.name = "climgweek"),
                chla = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek"),
                chla.max = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek")) %>% 
  dplyr::mutate(chla.lag1 = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek.lag1"),
                chla.lag2 = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek.lag2"),
                chla.lag3 = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek.lag3"),
                chla.lag4 = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek.lag4")) %>% 
  dplyr::mutate(chla.max.lag1 = extract_climg(dat = ., climg = chlamax_climg_filled, var.name = "climgweek.lag1"),
                chla.max.lag2 = extract_climg(dat = ., climg = chlamax_climg_filled, var.name = "climgweek.lag2"),
                chla.max.lag3 = extract_climg(dat = ., climg = chlamax_climg_filled, var.name = "climgweek.lag3"),
                chla.max.lag4 = extract_climg(dat = ., climg = chlamax_climg_filled, var.name = "climgweek.lag4"))

# Create dataset ------------------------------------------------

covariate.names <- c("sst", "chla", "chla.lag2", "chla.lag4", "chla.max", "depth", "slope", "dcanyons", "dcanyonsInc", "fcpi")

bw.pres <- bw.pres %>% 
  dplyr::select(date, x, y, all_of(covariate.names)) %>%
  dplyr::mutate(presence = 1) %>%
  dplyr::mutate(ID = 1:nrow(.))

bw.abs <- bw.abs %>% 
  dplyr::select(date, x, y, all_of(covariate.names)) %>%
  dplyr::mutate(presence = 0) %>%
  dplyr::mutate(ID = 1:nrow(.))

bw.df <- rbind(bw.pres, bw.abs) %>%
  dplyr::distinct(.)

# Add weights
bw.df <- bw.df %>% dplyr::mutate(weights = automatic_weights_creation(presence))

# Retrieve times
gps.pts <- do.call(rbind, gps) %>% add_xy()
gps.pts <- sp::SpatialPointsDataFrame(coords = gps.pts[, c("x", "y")], data = gps.pts, proj4string = CRSutm)
gps.pts <- sf::st_as_sf(gps.pts)

bw.pts <- sp::SpatialPointsDataFrame(coords = bw.df[, c("x", "y")], data = bw.df, proj4string = CRSutm)
bw.pts <- sf::st_as_sf(bw.pts)

pb <- dplyr::progress_estimated(nrow(bw.pts))
closest.times <- purrr::map(.x = seq_len(nrow(bw.pts)),
                              .f = ~{
                                pb$tick()$print()
                                day <- bw.pts[.x, ]$date
                                gps.subset <- subset(gps.pts, date == day)
                                id <- sf::st_nearest_feature(bw.pts[.x,], gps.subset)
                                list(day = gps.subset[id,]$date, time = gps.subset[id,]$time)
                                })

bw.df$time <- do.call(c, purrr::map(closest.times, "time"))
bw.df$datetime <- lubridate::as_datetime(paste0(bw.df$date, " ", bw.df$time))
bw.df <- dplyr::arrange(bw.df, datetime)

# Model fitting ------------------------------------------------------------------

corr <- cor(bw.df[, covariate.names], use = "complete.obs")
ggcorrplot::ggcorrplot(corr, hc.order = FALSE, type = "lower", lab = TRUE)

# Generate model list using a full subsets approach
bw.modelset <- na.omit(bw.df) %>%
  list.model.sets(response.var = "presence", 
                               dat = ., 
                               covariates = covariate.names,
                               corr.cutoff = 0.25,
                               max.predictors = 3,
                               k = 4,
                               use.weights = TRUE,
                               verbose = TRUE)

# Fit candidate models
bw.gamlist <- fit.model.set(model.set.list = bw.modelset, max.models = 1000)

# Extract model selection results
bw.gams <- FSS_results(bw.gamlist, plot.models = FALSE)
bestmodel <- bw.gams$best.model$obj
summary(bestmodel)
plot(bestmodel, scale = 0, pages = 1, shade = TRUE)
purrr::walk(.x = 1:3, .f = ~{
  pdf(paste0("plot", .x, ".pdf"), height = 5, width = 5)
  plot(bestmodel, scale = 0, shade = TRUE, select = .x)
  abline(h = 0, lty = 2)
  dev.off()
})

# Model selection summary
all.gams <- bw.gams$model.table[1:10, c(revised.covariates, "AICc", "wi.AICc", "delta.AICc")] |> 
  dplyr::select(depth, fcpi, slope, dcanyons, dcanyonsInc, chla, chla.lag2, chla.lag4, chla.max, AICc, wi.AICc, delta.AICc)
row.names(all.gams) <- NULL
xtable::xtable(all.gams)

# Indidividual terms
fit_and_plot("chla.max", k = 4)
fit_and_plot("dcanyonsInc", k = 4)
fit_and_plot("slope", k = 4)

# Covariate contributions
# --> fit full and reduced models
# See: https://stat.ethz.ch/pipermail/r-help/2011-November/295324.html

# Model without Chl-a max
b1 <- mgcv::gam(presence ~ s(dcanyonsInc, k = 4, bs = "tp") + s(slope, k = 4, bs = "tp"), 
                sp = bestmodel$sp[2:3],
                data = bw.df,
                weights = bw.df$weights,
                family = binomial(link = "logit"),
                method = "REML")

# Model without distInc
b2 <- mgcv::gam(presence ~ s(chla.max, k = 4, bs = "tp") + s(slope, k = 4, bs = "tp"), 
                sp = bestmodel$sp[c(1,3)],
                data = bw.df,
                weights = bw.df$weights,
                family = binomial(link = "logit"),
                method = "REML")

# Model without slope
b3 <- mgcv::gam(presence ~ s(chla.max, k = 4, bs = "tp") + s(dcanyonsInc, k = 4, bs = "tp"), 
                sp = bestmodel$sp[1:2],
                data = bw.df,
                weights = bw.df$weights,
                family = binomial(link = "logit"),
                method = "REML")

b4 <- mgcv::gam(presence ~ s(chla.max, k = 4, bs = "tp"), 
                sp = bestmodel$sp[1],
                data = bw.df,
                weights = bw.df$weights,
                family = binomial(link = "logit"),
                method = "REML")

b5 <- mgcv::gam(presence ~ s(dcanyonsInc, k = 4, bs = "tp"), 
                sp = bestmodel$sp[2],
                data = bw.df,
                weights = bw.df$weights,
                family = binomial(link = "logit"),
                method = "REML")

b6 <- mgcv::gam(presence ~ s(slope, k = 4, bs = "tp"), 
                sp = bestmodel$sp[3],
                data = bw.df,
                weights = bw.df$weights,
                family = binomial(link = "logit"),
                method = "REML")

# Null model
b0 <- mgcv::gam(presence ~ 1, 
                data = bw.df,
                weights = bw.df$weights,
                family = binomial(link = "logit"),
                method = "REML")

b <- bestmodel
dev.1 <- (deviance(b1)-deviance(b))/deviance(b0) ## prop explained by chla
dev.2 <- (deviance(b2)-deviance(b))/deviance(b0) ## prop explained by distInc
dev.3 <- (deviance(b3)-deviance(b))/deviance(b0) ## prop explained by slope

dev.1[2]<- (deviance(b0)-deviance(b4))/deviance(b0) 
dev.2[2]<- (deviance(b0)-deviance(b5))/deviance(b0)
dev.3[2] <- (deviance(b0)-deviance(b6))/deviance(b0)

dev.1 <- mean(dev.1)
dev.2 <- mean(dev.2)
dev.3 <- mean(dev.3)

dev.1 + dev.2 + dev.3
# summary(bestmodel)$dev.expl

# Model validation -------------------------------------------------------------

# Residual checks
par(mfrow = c(3,2))
gam.check(bestmodel)
qq.gam(bestmodel, rep = 1000, pch = 16)
acf(residuals(bestmodel, "pearson"))
par(mfrow = c(1,1))

# Spatial autocorrelation 
Vario.D <- 25000
Vario.bw <- geoR::variog(coords = na.omit(bw.df)[, c("x", "y")],
                         data = residuals(bestmodel, type = "pearson"),
                         max.dist = Vario.D,
                         breaks = seq(0, Vario.D, by = 500))

# Model performance -------------------------------------------------------------

# The final model is used to predict the data on the response scale (i.e. a value between 0 and 1)
pr <- predict(bestmodel, na.omit(bw.df), type = "response")                          

# Specify the vector of predictions (pr) and the vector of labels (i.e. the observed values "Pres")
pred <- ROCR::prediction(predictions = as.numeric(pr), labels = na.omit(bw.df)$presence)          

# Model performance in the form of the true positive rate and the false positive rate
perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr") 

plot(perf, colorize = TRUE) # to plot the ROC curve
abline(0, 1, lty = 2)

# Metric 1 (AUC)
pROC::roc(na.omit(bw.df)$presence, pr)$auc
# Same as ROCR::performance(pred, "auc")@y.values[[1]]

# Build the confusion matrix
conf.mat <- SDMTools::confusion.matrix(obs = na.omit(bw.df)$presence, pred = pr, 
                                       threshold = SDMTools::optim.thresh(obs = na.omit(bw.df)$presence, pred = pr)$`sensitivity=specificity`[1])

# Metric 2 (TP)
TP <- conf.mat[2,2] 
# / sum(bw.df$presence)
FN <- conf.mat[1,2]
FP <- conf.mat[2,1]
TN <- conf.mat[1,1]

# Metric 3 (Sørensen)
Sorensen <- 2 * TP/(FN + 2 *TP + FP)

# Metric 4 (TSS)
tss <- SDMTools::sensitivity(conf.mat) + SDMTools::specificity(conf.mat) - 1
# Same as: TP/sum(bw.df$presence == 1) + TN/sum(bw.df$presence == 0) - 1

#'---------------------------------------------
# Cross-validation settings
#'---------------------------------------------

nCV <- 100
cv.ratio <- 0.25

#'---------------------------------------------
# Create training and testing datasets
#'---------------------------------------------

training.df.day <- list() 
evaluation.df.day <- list() 

bw.df = na.omit(bw.df)
for (i in 1:nCV){
  
  print(i)
  
  dsummary <- bw.df %>% dplyr::group_by(date) %>% dplyr::summarise(n = sum(presence))
  
  int <- 0
  while(int==0){
    days_out <- unique(bw.df$date)[sample(c(1:length(unique(bw.df$date))), 
                                          round(cv.ratio*length(unique(bw.df$date))), 
                                          replace = FALSE)] 
    if(sum(dsummary[dsummary$date%in%days_out,]$n)>0) int <- 1
  }
  
  bw.training <- bw.df[!(bw.df$date %in% days_out),] %>% 
    dplyr::mutate(response = presence)
  
  bw.evaluation <- bw.df[bw.df$date %in% days_out,] %>% 
    dplyr::mutate(response = presence)
  
  if(sum(bw.evaluation$response)==0)
    
  # Add weights 
  bw.training <- bw.training %>% dplyr::mutate(weights = automatic_weights_creation(presence))
  bw.evaluation <- bw.evaluation %>% dplyr::mutate(weights = automatic_weights_creation(presence))
  
  # Save in list
  training.df.day <- c(training.df.day, list(bw.training))
  evaluation.df.day <- c(evaluation.df.day, list(bw.evaluation))
}

rm(dsummary, days_out)

#'---------------------------------------------
# Set up tibble object to hold results
#'---------------------------------------------

mod_cv <- tibble::tibble(cv.run = NA, 
                         AUC = NA, 
                         TP = NA,
                         Sørensen = NA,
                         TSS = NA)
models_cv <- list()

#'---------------------------------------------
# Perform cross-validation
#'---------------------------------------------

for(i in 1:nCV){
  
  print(i)
  
  m <- mgcv::gam(bestmodel$formula, 
                 weights = training.df.day[[i]]$weights, 
                 family = binomial(link = "logit"), 
                 method = "REML", 
                 data = training.df.day[[i]]) 
  
  models_cv[[i]] <- m # Save models in list
  
  # Predictions to evaluation and training
  pred.ext <- predict(m, evaluation.df.day[[i]], type = "response") 
  # pred.int <- predict(m, training.df.day[[i]], type = "response")
  
  # AUC
  # suppressMessages(int.AUC <- round(pROC::roc(training.df.day[[i]]$presence, pred.int)$auc, 3))
  suppressMessages(ext.AUC <- round(pROC::roc(evaluation.df.day[[i]]$presence, pred.ext)$auc, 3))
  # diff.AUC <- round(int.AUC - ext.AUC, 3)
  
  # Confusion matrix
  conf.mat.ext <- SDMTools::confusion.matrix(obs = evaluation.df.day[[i]]$presence, pred = pred.ext,
                  threshold = SDMTools::optim.thresh(obs = evaluation.df.day[[i]]$presence, 
                                                     pred = pred.ext)$`sensitivity=specificity`[1])
  
  # Metric 2 (TP)
  TP.ext <- conf.mat.ext[2,2] 
  FN.ext <- conf.mat.ext[1,2]
  FP.ext <- conf.mat.ext[2,1]
  
  # Metric 3 (Sørensen)
  Sorensen.ext <- 2 * TP.ext/(FN.ext + 2 *TP.ext + FP.ext)
  
  # Metric 4 (TSS)
  tss.ext <- SDMTools::sensitivity(conf.mat.ext) + SDMTools::specificity(conf.mat.ext) - 1
  
  mod_cv[i, ] <- list(i, ext.AUC, TP.ext / sum(evaluation.df.day[[i]]$presence == 1),
                      Sorensen.ext, tss.ext)
}

rm(m, pred.ext, conf.mat.ext, tss.ext, Sorensen.ext, TP.ext, FN.ext, FP.ext)

mod_cv <- mod_cv %>% tidyr::pivot_longer(cols = c(2:5), names_to = "metric")
mod_cv %>% dplyr::group_by(metric) %>% dplyr::summarise(min = min(value), max = max(value), mean = mean(value), sd = sd(value))

# Fluke up dives ----------------------------------------------------------

## Records of Fluke up dives

fud <- readr::read_csv("data/All data Eni WWR to S&P 22 FUDs.csv", show_col_types = FALSE) |> 
  janitor::clean_names() |> 
  dplyr::mutate(date = as.Date(date, format = "%d/%m/%Y", tz = "Australia/West", usetz = TRUE)) |> 
  dplyr::select(-x17, -id) |> 
  dplyr::slice(1:17) |> 
  dplyr::select(-block, -obs_type, -species, -group, -resighted, -depth_cat, -season, -sight_no, -sight_time) |> 
  dplyr::rename(N_adults = adults, N_calves = calves, depth = depth_m, N_FUDs = tu_dives) |> 
  dplyr::mutate(N_calves = ifelse(is.na(N_calves), 0, N_calves)) |> 
  dplyr::mutate(N_ind = N_adults + N_calves) |> 
  dplyr::select(-adcalves) |> 
  dplyr::relocate(N_ind, .before = N_adults)

fud <- add_xy(fud) %>%
  dplyr::select(-depth) %>%
  dplyr::mutate(climgweek = purrr::map(.x = fud$date, .f = ~find_week(.x)) %>% do.call("c", .))

fud <- fud %>% dplyr::mutate(
  sst = extract_climg(dat = ., climg = sst_climg, var.name = "climgweek"),
  chla = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek"),
  chla.max = extract_climg(dat = ., climg = chlamax_climg_filled, var.name = "climgweek")
)

fud <- dplyr::bind_cols(fud, tibble::as_tibble(raster::extract(x = static.env, y = fud[, c("x", "y")])))
dplyr::mutate(depth = -1 * raster::extract(depth.r, cbind(x, y))) |>
  dplyr::mutate(distInc = raster::extract(dist_canyons, cbind(x, y))) |>
  dplyr::select(-x, -y)

bw.fud <- bw %>%
  dplyr::mutate(climgweek = purrr::map(.x = bw$date, .f = ~ find_week(.x)) %>% do.call("c", .)) %>%
  dplyr::mutate(
    sst = extract_climg(dat = ., climg = sst_climg, var.name = "climgweek"),
    chla = extract_climg(dat = ., climg = chla_climg_filled, var.name = "climgweek"),
    chla.max = extract_climg(dat = ., climg = chlamax_climg_filled, var.name = "climgweek")
  )
bw.fud <- dplyr::bind_cols(bw.fud, tibble::as_tibble(raster::extract(x = static.env, y = bw.fud[, c("x", "y")])))

bw.fud <- bw.fud %>%
  dplyr::left_join(., fud, by = c("longitude", "latitude")) |>
  dplyr::mutate(FUD = ifelse(is.na(N_FUDs), 0, 1)) |>
  dplyr::mutate(FUD = as.factor(FUD))

ggplot2::ggplot(bw.fud, aes(x = FUD, y = dcanyonsInc, fill = FUD)) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter(
    shape = 15,
    color = "steelblue",
    position = position_jitter(0.21)
  ) +
  ggplot2::theme_classic()

kruskal.test(depth.x ~ FUD, data = bw.fud)
kruskal.test(dcanyonsInc ~ FUD, data = bw.fud)