##%######################################################%##
#                                                          #
####      Evidence of pygmy blue whale foraging in      ####
####     the Timor Trough during the austral winter     ####
#                                                          #
##%######################################################%##

# Libraries and global settings ------------------------------------------

#'---------------------------------------------
# Load (+ install if needed) all packages
#'---------------------------------------------

# Note: sf requires GDAL, which can be tricky to install on some OS
# For Mac, see: https://stackoverflow.com/questions/44973639/trouble-installing-sf-due-to-gdal

require(remotes)
remotes::install_github("seananderson/ggsidekick")
remotes::install_github("dkahle/ggmap", ref = "tidyup")
remotes::install_github("beckyfisher/FSSgam_package")

pacman::p_load(tidyverse, # Tidyverse
               raster, # Raster and GIS operations
               rgdal, # Working with geospatial data
               rgeos, # Topology operations on geometries
               pals, # Colour palettes (parula)
               plyr, # Split-apply-combine paradigm
               sf, # Simple features
               ncf, # Spatial covariance functions
               mgcv, # Generalised Additive Models (GAMs)
               WaveletComp, # Wavelet analysis
               GGally, # ggplot2 extensions
               SDMTools, # Tools for processing data in SDMs
               lubridate, # Date handling
               ggplot2, # Graphics
               rworldmap, # Mapping global data
               gstats,# Spatial and Spatio-Temporal Geostatistical Modelling, Prediction and Simulation
               purrrlyr, # Intersection of purrr/dplyr
               scales, # For comma format on plots
               cowplot, # For combining multiple ggplots
               ggsidekick, # Nice ggplot theme by Sean Anderson 
               rerddap, # R client for working with ERDDAP servers
               smoothr, # Smooth and Tidy Spatial Features in R. 
               janitor, # Data cleaning 
               FSSgam) # All subsets GAMs

# Note that Google has recently changed its policies. For use after July 2018, need to
# log on to Google Cloud Platform, create a project, generate API key and enable billing
# https://stackoverflow.com/questions/52565472/get-map-not-passing-the-api-key-http-status-was-403-forbidden

library(ggmap)

set.seed(87) # Set the random seed (reproducible results)

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

#'---------------------------------------------
# Import sightings data
#'---------------------------------------------

timor.dat <- readr::read_csv("data/timorbw_data.csv")

#'---------------------------------------------
# GPS tracks
#'---------------------------------------------

gps <- list('2007' = NULL, '2008' = NULL)

gps$`2007` <- readr::read_csv("data/gps_albacora07.csv")
gps$`2008` <- readr::read_tsv("data/gps_bicuda08.txt")

#'---------------------------------------------
# Clean up column names
#'---------------------------------------------

timor.dat <- timor.dat %>% janitor::clean_names()
gps <- purrr::map(.x = gps, .f = function(x) janitor::clean_names(x))

#'---------------------------------------------
# Filter data by month
#'---------------------------------------------

gps.dates <- purrr::map(.x = gps, 
                        .f = ~.x %>% 
                          dplyr::filter(lubridate::month(date)<=9) %>% 
                          dplyr::pull(date) %>% unique(.)) %>% 
  Reduce(c, .) %>% sort(.)

#'---------------------------------------------
# Format columns
#'---------------------------------------------

timor.dat$date <- as.Date(timor.dat$date, "%yy-%mm-%dd") # date
timor.dat <- chr2fac(timor.dat) # Converts all character variables to factors
gps <- purrr::map(.x = gps, .f = function(x) chr2fac(x))

#'---------------------------------------------
# Add a year column
#'---------------------------------------------

timor.dat <- timor.dat %>% dplyr::mutate(year = lubridate::year(date))

#'---------------------------------------------
# Filter out sightings with no effort data
#'---------------------------------------------

timor.dat <- timor.dat %>% 
  dplyr::filter(date%in%gps.dates, obs_type == "Dedicated")

#'---------------------------------------------
# Add seasons
#'---------------------------------------------

# Winter defined as June to September
# Inclusion of September justified from the tagging data described in Double et al (2014).
# whereby one individual continued transmitting after August 1st and was shown to remain within Banda
# Sea until September before starting migration.

timor.dat <- timor.dat %>% 
  dplyr::mutate(month = lubridate::month(date)) 

# Number of sightings (all species) per month

timor.dat %>% dplyr::pull(month) %>% table()

timor.dat <- timor.dat %>% 
  dplyr::mutate(season = as.factor(ifelse(month%in%6:9, "Winter", ifelse(month%in%10:12, "Spring", "Summer"))))

#'---------------------------------------------
# Extract blue whale data
#'---------------------------------------------

bw <- timor.dat %>% 
  dplyr::filter(species == "Blue whale") %>% 
  dplyr::filter(resighted == "Not a duplicate")

bw %>% dplyr::pull(month) %>% table(.)

#'---------------------------------------------
# Timeline of surveys and sightings
#'---------------------------------------------

data.frame(date = c(gps$`2007`$date, gps$`2008`$date)) %>% 
  dplyr::mutate(top = 1) %>% 
  dplyr::mutate(season = as.factor(ifelse(lubridate::month(date)%in%6:9, "Winter", ifelse(lubridate::month(date)%in%10:12, "Spring", "Summer")))) %>% 
  ggplot(data = ., aes(x = date, y = top)) + 
  geom_rug(sides = "b")+
  geom_point(data = bw %>% dplyr::mutate(top = 1), color = "blue")

# GIS layers ---------------------------------------------------------------------

#'---------------------------------------------
# Define coordinate projection systems
#'---------------------------------------------

CRSll <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
CRSutm <- sp::CRS("+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

#'---------------------------------------------
# Shapefile for Timor leste
#'---------------------------------------------

# Extracted from GSHHS - http://www.soest.hawaii.edu/wessel/gshhg/

coast <- raster::shapefile(file.path("gis", "timor_leste.shp"))
coast_utm <- sp::spTransform(coast, CRSutm)

#'---------------------------------------------
# Define regional extent
#'---------------------------------------------

bw.ext <- as(extent(c(121, 132, -16, -7.5)), "SpatialPolygons")
sp::proj4string(bw.ext) <- CRSll

#'---------------------------------------------
# Survey area
#'---------------------------------------------

study.area <- raster::shapefile(file.path("gis", "study_area.shp"))
sp::proj4string(study.area) <- CRSll
extent(study.area)
study.area_utm <- sp::spTransform(study.area, CRSutm)

#'---------------------------------------------
# Crop land
#'---------------------------------------------

timor <- raster::crop(coast, bw.ext)
timor_utm <- sp::spTransform(timor, CRSutm)

#'---------------------------------------------
# Submarine canyons
#'---------------------------------------------

# From Harris and Whiteway (2011). Marine Geology, 285(1-4): 69-86.
# Canyon types: 1. shelf-incising and river associated; 2. shelf-incising; and 3. blind.

canyons <- raster::shapefile("gis/global_canyons.shp")
canyons <- raster::crop(x = canyons, 
                        y = rgeos::gBuffer(spgeom = study.area, width = 0.2))
canyons_utm <- sp::spTransform(canyons, CRSutm)

#'---------------------------------------------
# Identify effort segments for each survey day
#'---------------------------------------------

gps$`2007`$timediff <- gps$`2008`$timediff <- NA
gps$`2007`$line_id <- gps$`2008`$line_id <-NA

# Calculate time difference between sequential points
# Effort segments are identified as sets of GPS positions less than 5 minutes apart

for(i in 1:nrow(gps$`2007`)){
  if(i==1){ gps$`2007`$timediff[i] <- 0
  }else{ gps$`2007`$timediff[i] <- as.numeric(lubridate::as.duration(lubridate::interval(gps$`2007`$time[i-1], gps$`2007`$time[i])), "minutes")
  }   
}

for(i in 1:nrow(gps$`2008`)){
  if(i==1){ gps$`2008`$timediff[i] <- 0
  }else{ gps$`2008`$timediff[i] <- as.numeric(lubridate::as.duration(lubridate::interval(gps$`2008`$time[i-1], gps$`2008`$time[i])), "minutes")
  }   
}

# Assign identifier to each segment of effort

# 2007

out <- split(gps$`2007`, f = gps$`2007`$date)

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
gps.all <- do.call(rbind, gps) %>% dplyr::filter(date%in%gps.dates)

#'---------------------------------------------
# Create line objects for each on-effort segment
#'---------------------------------------------

oneffort_sp <- gps.all %>% 
  split(x = ., f = .$date) %>% 
  purrr::map(., ~.x %>% 
               split(x = ., f = .$line_id) %>% 
               purrr::map(., ~ createLines(.)))

#'---------------------------------------------
# Project to UTM
#'---------------------------------------------

oneffort_sp <- purrr::map_depth(.x = oneffort_sp, .depth = 2, 
                                .f = ~sp::spTransform(.x, CRSobj = CRSutm)) %>% 
  purrr::map(.x = ., .f = ~do.call(rbind, .))

#'---------------------------------------------
# Discard short segments
#'---------------------------------------------

oneffort_sp <- purrr::map(.x = oneffort_sp, .f = ~smoothr::drop_crumbs(x =.x, threshold = 250))

#'---------------------------------------------
# Create 6 km buffers around tracks
#'---------------------------------------------

bufferonland <- purrr::map(.x = oneffort_sp, 
                           .f = ~ rgeos::gBuffer(spgeom = .x, 
                                                 byid = FALSE,
                                                 width = 6000, 
                                                 capStyle = "ROUND", 
                                                 joinStyle = "ROUND"))

plot(timor_utm, col = "grey")
purrr::walk(.x = bufferonland, .f = ~plot(.x, col = "yellow", add = TRUE))

#'---------------------------------------------
# Clip buffers overlapping land
#'---------------------------------------------

buffer_pol <- plyr::llply(bufferonland, function(x){Fun_robustBufferDiff(x, timor_utm)})

#'---------------------------------------------
# Check the clipping
#'---------------------------------------------

plot(timor_utm, col = "transparent")
purrr::walk(.x = buffer_pol,
            .f = ~plot(.x, add = TRUE, col = "green"))


# Environmental data ------------------------------------------------------

#'---------------------------------------------
# Bathymetry ====
#'---------------------------------------------

# From GEBCO 30s gridded dataset (see terms of use for acknowledgements)

bathy <- raster::raster(file.path("env", "gebco30sec.asc"))
proj4string(bathy) <- CRSll

bathy <- raster::projectRaster(bathy, crs = CRSutm) # Project to UTM

depth <- raster::mask(x = bathy, mask = study.area_utm)
depth <- raster::crop(x = depth, y = extent(study.area_utm))

# Quick plot

make_map(input.raster = depth, show.tracks = TRUE, plot.title = "Depth (m)", legend.title = FALSE)

#'---------------------------------------------
# Seabed slope ====
#'---------------------------------------------

# Values represent the 'rise over run'
# Convert to degrees by taking ATAN ( output ) * 57.29578

seabed_slope <- SDMTools::slope(mat = bathy, latlon = F) # Must be in proj units, not latlon
seabed_slope <- atan(seabed_slope)*180/pi

seabed_slope <- raster::mask(x = seabed_slope, mask = study.area_utm)
seabed_slope <- raster::crop(x = seabed_slope, y = extent(study.area_utm))

# Quick plot

make_map(input.raster = seabed_slope, show.tracks = TRUE, plot.title = "Slope (deg)", legend.title = FALSE)

#'---------------------------------------------
# Distance to coast ====
#'---------------------------------------------

# Use rgeos to compute point-polygon distances
# This outputs a matrix with a column for each feature in spgeom1.

rdf <- raster::as.data.frame(depth, xy = T)
rdf <- rdf[complete.cases(rdf),] # Removes NAs
rdf <- sp::SpatialPoints(rdf[,1:2], CRSutm) # Converts to SpatialPts

# Cartesian distances

dc <-  rgeos::gDistance(spgeom1 = coast_utm, 
                        spgeom2 = rdf, 
                        byid = TRUE)

# To get the nearest distance to any feature, apply min over rows

dist_coast <- apply(dc, 1, min)
dist_coast <- cbind(sp::coordinates(rdf), dist_coast)
dist_coast <- raster::rasterFromXYZ(xyz = dist_coast, crs = CRSutm)
dist_coast <- dist_coast/1000 # Convert to km

# Quick plot

make_map(input.raster = dist_coast, show.tracks = TRUE, 
         plot.title = "Distance to coast (km)", legend.title = FALSE)
rm(dc)

#'---------------------------------------------
# Distance to nearest canyon ====
#'---------------------------------------------

# All canyons included

dcanyons <-  rgeos::gDistance(spgeom1 = canyons_utm, 
                              spgeom2 = rdf, 
                              byid = TRUE)

# To get the nearest distance to any feature, apply min over rows

dist_canyons <- apply(dcanyons, 1, min)
dist_canyons <- cbind(sp::coordinates(rdf), dist_canyons)
dist_canyons <- raster::rasterFromXYZ(xyz = dist_canyons, crs = CRSutm)
dist_canyons <- dist_canyons/1000 # Convert to km

# Quick plot

make_map(input.raster = dist_canyons, show.tracks = TRUE, plot.title = "Distance to nearest canyon (km)", legend.title = FALSE)

# Only shelf-incising canyon

dcanyons.incising <- rgeos::gDistance(spgeom1 = canyons_utm[canyons_utm@data$class == 2,], spgeom2 = rdf, byid = TRUE)

# To get the nearest distance to any feature, apply min over rows

dist_incising <- apply(dcanyons.incising, 1, min)
dist_incising <- cbind(sp::coordinates(rdf), dist_incising)
dist_incising <- raster::rasterFromXYZ(xyz = dist_incising, crs = CRSutm)
dist_incising <- dist_incising/1000 # Convert to km

# Quick plot

make_map(input.raster = dist_incising, show.tracks = TRUE, plot.title = "Distance to nearest shelf-incising canyon (km)", legend.title = FALSE)

#'---------------------------------------------
# Sea Surface Temperature (SST) ====
#'---------------------------------------------

# https://cran.r-project.org/web/packages/rerddap/vignettes/Using_rerddap.html

# First determine the temporal variability of the system.
# MUR (Multi-scale Ultra-high Resolution) is an analyzed SST product at 0.01-degree resolution going back to 2002.

sstInfo <- rerddap::info('jplMURSST41')

# Define a location at which to assess temperature time series (virtual buoy)

buoy <- sp::SpatialPoints(coords = cbind(126.25, -9.55))
plot(study.area); points(buoy) # Quick plot

# Retrieve daily sst between 2002 and 2018

sst.buoy <- rerddap::griddap(x = sstInfo, 
                             longitude = rep(coordinates(buoy)[1], 2), 
                             latitude = rep(coordinates(buoy)[2], 2), 
                             time = c('2002-06-02','2018-12-31'),
                             fields = 'analysed_sst')

# Extract date/time and sst values

sst.buoy <- list(erddap = sst.buoy)

sst.buoy$data <- tibble(date = as.Date(sst.buoy$erddap$data$time,
                                       origin = '1970-01-01', tz = "GMT"),
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

# Retrieve daily sst

chla.buoy <- rerddap::griddap(x = chlaInfo, 
                              longitude = rep(coordinates(buoy)[1], 2), 
                              latitude = rep(coordinates(buoy)[2], 2), 
                              time = c('2003-01-01','2018-12-31'),
                              fields = 'chlorophyll')

# Extract date/time and chla values

chla.buoy <- list(erddap = chla.buoy)

chla.buoy$data <- tibble(date = as.Date(chla.buoy$erddap$data$time,
                                        origin = '1970-01-01', tz = "GMT"),
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

#'---------------------------------------------
# Frequency of Chlorophyll Peak Index (FCPI) ====
#'---------------------------------------------

fcpi <- calc_fcpi(chla.climg = chla_climg, log.10 = FALSE)

# Background sampling ------------------------------------------------

#'---------------------------------------------
# Summarise survey effort
#'---------------------------------------------

survey_time <- as.data.frame(table(gps.all$date))
names(survey_time) <- c("date", "time_on_effort")

hist(survey_time$time_on_effort, xlab="time ON-effort (min)", main = "")
mean(survey_time$time_on_effort)

Total.on_effort <- sum(survey_time$time_on_effort)/60

#'---------------------------------------------
# Generate pseudo-absence (background) points
#'---------------------------------------------

N.background <- 1000 # Number of points to generate

i <- 0

background_sp <- lapply(buffer_pol,function(x){ 
  
  i <<- i + 1 # Increase counter
  
  day <- as.character(names(buffer_pol)[i]) # Retrieve day
  bw.sights <- bw %>% dplyr::filter(date==day) # Filter sightings by day
  
  # If any sigtings were made, create 4 km buffer around them
  
  if(nrow(bw.sights)>0){
    bw.sights.utm <- sp::spTransform(sp::SpatialPoints(coords = bw.sights[, c("longitude", "latitude")], 
                                                       proj4string = CRSll), CRSobj = CRSutm)
    exclusion.zones <- rgeos::gBuffer(spgeom = bw.sights.utm, byid = TRUE, width = 4000)
  }
  
  # Time on effort that day
  t.on <- survey_time$time_on_effort[i]
  
  # Generate a regular grid of points separated by 4 km
  spts <- sp::spsample(x, type = "regular", iter = 1000, cellsize = 4000)
  spts <- sp::spTransform(spts, CRSutm)
  
  # Prevent background points from being sampled within exclusion zones around presence points
  if(nrow(bw.sights)>0) spts <- rgeos::gDifference(spts, exclusion.zones)
  
  # Randomly sample points
  spts.2 <- sample(spts, round((N.background/Total.on_effort)*t.on/60,0), replace = FALSE) 
  
  return(spts.2)
})

#'---------------------------------------------
# Plot a random day to check
#'---------------------------------------------

check_bg()

#'---------------------------------------------
# Convert bg points to a single dataframe object
#'---------------------------------------------

background_df <- background_sp %>% do.call(rbind,.) %>% 
  sp::spTransform(., "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% 
  as_tibble(.) %>% 
  dplyr::mutate(presence = 0)

#'---------------------------------------------
# Add date and clean up
#'---------------------------------------------

background_df$date <- purrr::map_dbl(.x = background_sp, .f = ~length(.x)) %>% 
  purrr::map2(.x = names(.), 
              .y = .,
              .f = ~rep(.x, .y)) %>% do.call(c, .)

names(background_df) <- c("longitude","latitude","presence", "date")
background_df$date <- as.Date(background_df$date, format = "%Y-%m-%d",
                              tz = "Australia/West", usetz = TRUE)

# Extract covariate values at presence points -------------------------------------------------------

#'---------------------------------------------
# Find week corresponding to each date
#'---------------------------------------------

bw <- bw %>% dplyr::mutate(climgweek = purrr::map(.x = bw$date, .f = ~find_week(.x)) %>% do.call("c", .))

#'---------------------------------------------
# Retrieve values of static covariates
#'---------------------------------------------

static.env <- raster::stack(depth, seabed_slope, dist_coast, dist_canyons, dist_incising, fcpi$fcpi)
names(static.env) <- c('depth', 'slope', 'dcoast', 'dcanyons', 'dcanyonsInc', 'fcpi')
static.env <- raster::projectRaster(from = static.env, crs = CRSll)

bw.env <- tibble::as_tibble(raster::extract(x = static.env, y = bw[, c('longitude', 'latitude')]))

#'---------------------------------------------
# Retrieve values of dynamic covariates
#'---------------------------------------------

bw <- bw %>% dplyr::mutate(sst = extract_climg(dat = ., climg = sst_climg, var.name = "climgweek"),
                           chla = extract_climg(dat = ., climg = chla_climg, var.name = "climgweek"),
                           chla.max = extract_climg(dat = ., climg = chlamax_climg, var.name = "climgweek"))

bw <- as_tibble(cbind(bw, bw.env))


# Extract covariate values at background points ------------------------------------------------

#'---------------------------------------------
# Find week corresponding to each date
#'---------------------------------------------

bw.abs <- background_df %>% 
  dplyr::mutate(climgweek = purrr::map(.x = background_df$date, .f = ~find_week(.x)) %>% do.call("c", .))

#'---------------------------------------------
# Retrieve values of static covariates
#'---------------------------------------------

bw.abs.env <- tibble::as_tibble(raster::extract(x = static.env, y = bw.abs[, c('longitude', 'latitude')]))

#'---------------------------------------------
# Retrieve values of dynamic covariates
#'---------------------------------------------

bw.abs <- bw.abs %>% dplyr::mutate(sst = extract_climg(dat = ., climg = sst_climg, var.name = "climgweek"),
                                   chla = extract_climg(dat = ., climg = chla_climg, var.name = "climgweek"),
                                   chla.max = extract_climg(dat = ., climg = chla_climg, var.name = "climgweek"))

bw.abs <- tibble::as_tibble(cbind(bw.abs, bw.abs.env))

# Create binomial dataset ------------------------------------------------

#'---------------------------------------------
# Prepare presence points
#'---------------------------------------------

bw.pres <- bw %>% 
  dplyr::select(date, longitude, latitude, 
                sst, chla, chla.max, depth, slope, dcoast, dcanyons, dcanyonsInc, fcpi) %>%
  dplyr::mutate(presence = 1)

head(bw.pres)

#'---------------------------------------------
# Prepare background points
#'---------------------------------------------

bw.abs <- bw.abs %>% 
  dplyr::select(date, longitude, latitude, 
                sst, chla, chla.max, depth, slope, dcoast, dcanyons, dcanyonsInc, fcpi) %>% 
  dplyr::mutate(presence = 0)

#'---------------------------------------------
# Combine presence and background points
#'---------------------------------------------

bw.pb <- rbind(bw.pres, bw.abs)

# Add weights -------------------------------------------------------------

bw.pb <- na.omit(bw.pb) %>% 
  dplyr::distinct(.)

bw.pb <- bw.pb %>% dplyr::mutate(weights = automatic_weights_creation(presence))


# Model selection ------------------------------------------------------------------

#'---------------------------------------------
# Retrieve UTM coordinates
#'---------------------------------------------

bw.xy <- sp::SpatialPoints(coords = bw.pb[, c("longitude", "latitude")], proj4string = CRSll) %>% 
  sp::spTransform(., CRSobj = CRSutm) %>% 
  as_tibble(.) %>% dplyr::rename(x = longitude, y = latitude)

bw.pb <- cbind(bw.pb, bw.xy)

bw.pb <- bw.pb %>% dplyr::mutate(depth = abs(depth))

#'---------------------------------------------
# Max chlorophyll-a only
#'---------------------------------------------

mod0 <- mgcv::gam(presence ~ s(chla.max, k = 4, bs = "tp"),
                  data = bw.pb,
                  weights = bw.pb$weights,
                  family = binomial(link = "logit"),
                  method = "REML")

summary(mod0)
plot(mod0, shade = TRUE, scale = 0); abline (h = 0, lty = 2)

#'---------------------------------------------
# Generate model list using a full subsets approach
#'---------------------------------------------

bw.modelset <- list.model.sets(response.var = "presence", 
                               dat = bw.pb, 
                               covariates = c("sst", "chla.max", "chla", "depth",
                                              "slope", "dcanyonsInc", "fcpi", "dcoast"),
                               corr.cutoff = 0.5,
                               max.predictors = 4,
                               k = 5,
                               use.weights = TRUE,
                               verbose = TRUE)

#'---------------------------------------------
# Fit candidate models
#'---------------------------------------------

bw.gamlist <- fit.model.set(model.set.list = bw.modelset, max.models = 1000)

#'---------------------------------------------
# Extract model selection results
#'---------------------------------------------

bw.gams <- FSS_results(bw.gamlist, plot.models = FALSE)
bestmodel <- bw.gams$best.model$obj

#'---------------------------------------------
# Perform model checks
#'---------------------------------------------

summary(bestmodel)

par(mfrow = c(2,2))
gam.check(bestmodel)
qq.gam(bestmod, rep = 1000, pch = 16)

testmodel <- gam(presence ~ s(fcpi, k = 3) + s(chla.max, k = 3) + s(slope, k = 3) + s(dcanyonsInc, k = 3, bs = "tp"),
                 weights = bw.pb$weights,
                 REML = TRUE,
                 family = binomial(link = "logit"),
                 data = bw.pb)

summary(testmodel)
acf(residuals(testmodel, "pearson"))
pacf(residuals(testmodel, "pearson"))
plot(bestmodel, scale = 0, shade = TRUE, pages = 1)
qq.gam(testmodel, rep = 1000, pch = 16)
par(mfrow = c(2,2))
gam.check(testmodel)
par(mfrow = c(1,1))

# Spatial autocorrelation -------------------------------------------------------------

correlog.bw <- plot_correlogram(model = bestmodel, data = bw.pb, maxD = 50000)
plot_variogram(model = bestmodel, dat = bw.pb)

Vario1 <- geoR::variog(coords = bw.pb[, c("x", "y")],
                       data = residuals(bestmodel, type = "pearson"),
                       max.dist = 25000, breaks = seq(0, 25000, by = 1000))
plot(Vario1)

# generate standadized residuals

resids <- residuals(bestmodel, type = "pearson")
var.dat_resid <- gstat::variogram(resids~1, loc= ~x+y, data=bw.pb)
variog.fit <- gstat::fit.variogram(var.dat_resid, gstat::vgm(c("Sph", "Exp")))
plot(var.dat_resid, variog.fit)  


inv_dist = with(bw.pb, 1/dist(cbind(x, y), diag = T, upper = T))
inv_dist = as.matrix(inv_dist)
ape::Moran.I(bw.pb$presence, weight = inv_dist, scaled=T)

plot(nlme::Variogram(bestmodel, form = ~ x + y, data = bw.pb))

resids <- as.numeric(residuals(bestmodel, type = "pearson"))
resids.df <- data.frame(bw.pb, resids = resids)
ggplot(resids.df[resids.df$resids > -1,], aes(x, y, colour = resids)) + 
  geom_point(size = 5, alpha = 0.5) + scale_color_gradient2()

lr.fit <- mgcViz::getViz(testmodel)
mgcViz::qq(lr.fit, rep = 100, level = .9, CI = "quantile")

correlog1.1 <- ncf::correlog(bw.pb$x, bw.pb$y, residuals(bestmodel, type = "pearson"),
                        na.rm=T, increment=1, resamp=0)
plot(correlog1.1$correlation[1:50], type="b", pch=16, cex=1.5, lwd=1.5,
     xlab="distance", ylab="Moran's I", cex.lab=2, cex.axis=1.5, ylim = c(-1, 1)); abline(h=0)
plot(bw.pb$x, bw.pb$y, col=c("blue", "red")[sign(resid(bestmodel))/2+1.5], pch=19,
     cex=abs(resid(bestmodel))/max(resid(bestmodel))*2, xlab="geographical xcoordinates", ylab="geographical y-coordinates")


gamm_spat = gamm(presence ~ s(chla.max, k = 5, bs = "tp") + 
                   s(fcpi, k = 5, bs = "tp"), 
                 family = binomial(link = "logit"),
                 data=bw.pb, 
                 correlation = corGaus(1, form = ~ x + y))
plot(gamm_spat)

acf(residuals(bestmodel))
pacf(residuals(bestmodel))

bw.pb$fda <- factor(bw.pb$date)


# Model performance -------------------------------------------------------------

# The final model is used to predict the data on the response scale (i.e. a value between 0 and 1)
pr <- predict(bestmodel, bw.pb, type = "response")                          

# Specify the vector of predictions (pr) and the vector of labels (i.e. the observed values "Pres")
pred <- ROCR::prediction(predictions = as.numeric(pr), labels = bw.pb$presence)                           
# to assess model performance in the form of the true positive rate and the false positive rate
perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr") 

plot(perf, colorize=TRUE) # to plot the ROC curve
abline(0,1, lty = 2)

auc <- ROCR::performance(pred, "auc")
auc <- auc@y.values[[1]]

myroc <- pROC::roc(bw.pb$presence, pr)
plot(myroc)
ci(myroc)

# Choice of the best cut-off probability

y<-as.data.frame(perf@y.values)
x<-as.data.frame(perf@x.values)
fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45° line and the line joining the origin with the point (x;y) on the ROC curve
L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
d <- L*sin(fi)                                                     # to calculate the distance between the 45° line and the ROC curve
d <- tibble::tibble(d)
names(d) <- "value"
# The table should then be opened in Microsoft Excel to find the maximum distance with the command "Sort", and the relative position (i.e. the number of the corresponding record)
# MAX d= 0.277124605038953 --> position 2098

which(d == max(d, na.rm = TRUE))
alpha<-as.data.frame(perf@alpha.values)    # the alpha values represent the corresponding cut-offs
alpha[which(d == max(d, na.rm = TRUE)),]      


DATA<-matrix(0,nrow(bw.pb),3)     # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 4973 - the number of rows can be checked with dim(dat)) 
DATA<-as.data.frame(DATA)
names(DATA)<-c("plotID","Observed","Predicted")
DATA$plotID<-1:nrow(bw.pb)                                           # the first column is filled with an ID value that is unique for each row
DATA$Observed<-bw.pb$presence                                            # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(bestmodel,bw.pb,type="response")                 # the third column reports the predictions
confusmat <- PresenceAbsence::cmx(DATA, threshold = alpha[which(d == max(d, na.rm = TRUE)),])       

PresenceAbsence::presence.absence.accuracy(DATA, threshold = alpha[which(d == max(d, na.rm = TRUE)),] )
PresenceAbsence::auc.roc.plot(DATA, color = TRUE, legend.cex = 1.4, main = "")


sorensen <- 2*confusmat[1,1]/(confusmat[2,1] + 2 * confusmat[1,1] + confusmat[1,2])

modEvA::HLfit(bestmodel, bin.method = "n.bins", n.bins = 10, main = "Hosmer-Lemeshow GOF, N bins")
Dsquared(bestmodel)
AUC(bestmodel)

#'---------------------------------------------
# Cross-validation settings
#'---------------------------------------------

nCV <- 100
cv.ratio <- 0.3

#'---------------------------------------------
# Create training and testing datasets
#'---------------------------------------------

training.df.day <- list() 
evaluation.df.day <- list() 

for (i in 1:nCV){
  
  print(i)
  
  dsummary <- bw.pb %>% dplyr::group_by(date) %>% dplyr::summarise(n = sum(presence))
  
  int <- 0
  while(int==0){
    days_out <- unique(bw.pb$date)[sample(c(1:length(unique(bw.pb$date))), 
                                          round(cv.ratio*length(unique(bw.pb$date))), 
                                          replace = FALSE)] 
    if(sum(dsummary[dsummary$date%in%days_out,]$n)>0) int <- 1
  }
  
  bw.training <- bw.pb[!(bw.pb$date %in% days_out),] %>% 
    dplyr::mutate(response = presence)
  
  bw.evaluation <- bw.pb[bw.pb$date %in% days_out,] %>% 
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
                         int.AUC = NA, 
                         ext.AUC = NA, 
                         diff.AUC = NA, 
                         threshold = NA, 
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
  pred.int <- predict(m, training.df.day[[i]], type = "response")
  
  # AUC = Threshold-independent metric of performance
  suppressMessages(int.AUC <- round(pROC::roc(training.df.day[[i]]$presence, pred.int)$auc, 3))
  suppressMessages(ext.AUC <- round(pROC::roc(evaluation.df.day[[i]]$presence, pred.ext)$auc, 3))
  diff.AUC <- round(int.AUC - ext.AUC, 3)
  
  # TSS = threshold-dependent metric calculated on confusion matrix from pred.ext
  thresh <- SDMTools::optim.thresh(obs = evaluation.df.day[[i]]$presence,
                                   pred = pred.ext)$`max.sensitivity+specificity`[1]
  
  conf.mat <- SDMTools::confusion.matrix(obs = evaluation.df.day[[i]]$presence, 
                                         pred = pred.ext, threshold = thresh)
  
  tss <- SDMTools::sensitivity(conf.mat) + SDMTools::specificity(conf.mat) - 1
  
  mod_cv[i, ] <- c(i, int.AUC, ext.AUC, diff.AUC, thresh, tss)
}

rm(m, pred.ext, pred.int, diff.AUC, thresh, conf.mat, tss)

summary(mod_cv)

# Boxplot

mod_cv %>% dplyr::select(-diff.AUC, -cv.run, -threshold) %>% 
  tidyr::pivot_longer(data = ., cols = c("int.AUC", "ext.AUC", "TSS"), names_to = "metric") %>% 
  dplyr::mutate(metric = factor(metric)) %>% 
  ggplot + aes(x = metric, y = value) + geom_boxplot()

# Model predictions -------------------------------------------------------------

#'---------------------------------------------
# Project static layers to UTM
#'---------------------------------------------

static.env.utm <- raster::projectRaster(from = static.env, crs = CRSutm)

#'---------------------------------------------
# Ensure all layers have same resolution and extent
#'---------------------------------------------

chla.rs <- purrr::map(.x = chlamax_climg, .f = ~ raster::resample(x = .x, y = static.env.utm$depth))

#'---------------------------------------------
# Get model predictions
#'---------------------------------------------

winter.weeks <- 22:35 

weekly.predictions <- get_predictions(model = bestmodel, time.span = winter.weeks) %>% 
  purrr::set_names(x = ., nm = names(chla.rs)[winter.weeks])

#'---------------------------------------------
# Plot average predictions
#'---------------------------------------------

bw.preds <- raster::calc(x = weekly.predictions, mean) 
plot(bw.preds, col = pals::viridis(n = 100))
points(bw.pb[bw.pb$presence==1,]$x, bw.pb[bw.pb$presence==1,]$y, pch = 16)
points(bw.pb[bw.pb$presence==0,]$x, bw.pb[bw.pb$presence==0,]$y, col = "grey", pch = '.')

#'---------------------------------------------
# Point biserial correlation
#'---------------------------------------------

pbis.df <- bw.pb %>%
  dplyr::mutate(HS = raster::extract(bw.preds, .[, c("x", "y")]))

pbis.df1 <- pbis.df %>% dplyr::filter(presence == 1)
pbis.df0 <- pbis.df %>% dplyr::filter(presence == 0)

dismo::evaluate(p = pbis.df1, a = pbis.df0, model = bestmodel)


ltm::biserial.cor(x = pbis.df$HS, y = pbis.df$presence, use = c("complete"), level = 2)
cor.test(pbis.df$HS, pbis.df$presence)
