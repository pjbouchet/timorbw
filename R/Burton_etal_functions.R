#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------
#' 
#'  Winter distribution and habitat use of blue whales (Balaenoptera musculus sp.) 
#'  in the Timor Trough, south of Timor-Leste during 2007-08
#'
#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------

#'---------------------------------------------
# Change all character columns to factors and vice versa
#'---------------------------------------------

chr2fac <- function(y){y %>% dplyr::mutate_if(is.character, as.factor)}
fac2chr <- function(y){y %>% dplyr::mutate_if(is.factor, as.character)}

#'---------------------------------------------
# Clip buffers 
#'---------------------------------------------

Fun_robustBufferDiff = function(spgeom1, spgeom2) {
  tryCatch(expr = rgeos::gDifference(spgeom1, spgeom2, byid = F),
           error = function(e) {w <- gSimplify(spgeom1, tol = 1000, topologyPreserve = T);
           w <- gBuffer(w, byid = TRUE, width = 100);
           ww <- rgeos::gDifference(w, spgeom2, byid = F);
           return(ww)
           }
  )
}

#'---------------------------------------------
# Quick convenience function to plot points and buffers
#'---------------------------------------------

check_bg <- function(){
  b.ind <- sample(x = seq_along(buffer_pol), size = 1)
  plot(buffer_pol[[b.ind]]); plot(background_sp[[b.ind]], add = TRUE, pch = 16)
  if(nrow(bw[bw$date==names(buffer_pol[b.ind]),])>0) points(sp::spTransform(sp::SpatialPoints(coords = bw[bw$date==names(buffer_pol[b.ind]), c("longitude", "latitude")], proj4string = CRSll), CRSobj = CRSutm),
                                                            col = "blue", pch = 16)}

#'---------------------------------------------
# Create a SpatialLinesDataFrame from a data.frame
#'---------------------------------------------

createLines <- function(df){
  
  # The df object must contain two columns named "longitude" and "latitude".
  
  #'--------------------------------------------------------------------
  # PARAMETERS
  #'--------------------------------------------------------------------
  #' @param df Input data.frame object.
  #'--------------------------------------------------------------------
  
  df.L <- df[, c("longitude", "latitude")] %>%  
    as.matrix(.) %>% 
    raster::spLines(., crs = CRSll)
  
  Sldf <- sp::SpatialLinesDataFrame(df.L, data = df)
  return(Sldf)}


#'---------------------------------------------
# Generate weekly climatologies from ERRDAPP data
#'---------------------------------------------

get_wkclimatology <- function(remote.dataset,
                              variable.name,
                              t.year = 2017,
                              start.month = 4, 
                              end.month = 4,
                              no.years = 15,
                              climg.direction = 'backwards',
                              summary.method = "mean",
                              region.shp = NULL,
                              rmatch = TRUE,
                              raster.dest = NULL,
                              error400 = TRUE){
  
  
  #'-----------------------------------------------------------
  # PARAMETERS
  #'-----------------------------------------------------------
  #' @param remote.dataset Name of the ERDDAP dataset of interest.
  #' @param variable.name Name of the variable of interest.
  #' @param t.year Focal year.
  #' @param start.month Start month for the climatology.
  #' @param end.month End month for the climatology.
  #' @param no.years Length of the climatology (in years).
  #' @param climg.direction Year span for the climatology, relative to the focal year t.year. One of 'backwards', 'forwards', or 'centered'.
  #' @param summary.method Desired summary statistic.
  #' @param region.shp Shapefile of the study area, used to define the spatial extent of the ERDDAP data query. 
  #' @param rmatch Whether to resample the remote sensing data to the same resolution as an existing raster of choice.
  #' @param raster.dest Raster to match.
  #' @param error400 Set to TRUE if error 400 occurs.
  #'-----------------------------------------------------------
  
  #'---------------------------------------------
  # Weekly climatologies
  #'---------------------------------------------
  
  time.window <-  7 # in days
  
  #'---------------------------------------------
  # Perform function checks
  #'---------------------------------------------
  
  if(!class(remote.dataset)=="character") stop('Dataset name must be a character vector')
  if(!climg.direction%in%c('backwards', 'forwards', 'centered')) stop('Unrecognised direction')
  if(!is.numeric(t.year) | !is.numeric(no.years) | !is.numeric(start.month) | !is.numeric(end.month)) stop('Non-numeric inputs to time.window')
  if(!class(raster.dest)=='RasterLayer') stop('Input raster must be of class RasterLayer')
  if(!class(region.shp)%in%c('SpatialPolygonsDataFrame', 'SpatialPolygons')) stop('Input raster must be of class RasterLayer')
  
  message('Retrieving dataset details ...')
  
  #'---------------------------------------------
  # Get information on desired ERDDAP dataset.
  #'---------------------------------------------
  
  datInfo <- rerddap::info(remote.dataset)
  
  message('Generating dates ...')
  
  #'---------------------------------------------
  # Generate weekly date windows
  #'---------------------------------------------
  
  # Start January 1 and end December 24
  
  t.start <- lubridate::ymd(paste0(t.year, '-', 1, '-', 1))
  t.end <- lubridate::ymd(paste0(t.year, '-', 12, '-', 24))
  
  start.dates <- seq(t.start, t.end, by = paste0(time.window, ' days'))
  
  # Adjust for end month in subsequent year
  
  if(end.month < start.month) start.dates <- c(start.dates, start.dates + lubridate::years(1))
  
  # Generate start dates with one week leeway on either side
  
  start.dates <- start.dates[(min(which(lubridate::month(start.dates)==start.month))-1):(max(which(lubridate::month(start.dates)==end.month))+1)]
  
  rnames <- gsub(pattern = paste0(t.year, '-'), replacement = '', x = start.dates)
  rnames <- na.omit(rnames)
  
  # Generate corresponding end dates
  
  end.dates <- start.dates + lubridate::days(time.window-1)
  
  start.dates <- na.omit(start.dates)
  end.dates <- na.omit(end.dates)
  # if(end.month < start.month) t.end <- lubridate::ymd(paste0(t.year + 1, '-', 12, '-', 24))
  # if(end.month >= start.month) t.end <- lubridate::ymd(paste0(t.year, '-', 12, '-', 24))
  # 
  # start.dates <- seq(t.start, t.end, by = paste0(time.window, ' days'))
  # 
  # t.start <- lubridate::ymd(paste0(t.year, '-', start.month, '-', 1))
  # 
  # 
  # if(end.month%in%c(1, 3, 5, 7, 8, 10, 12)){
  #   
  #   if(end.month < start.month) t.end <- lubridate::ymd(paste0(t.year + 1, '-', end.month, '-', 31))
  #   if(end.month >= start.month) t.end <- lubridate::ymd(paste0(t.year, '-', end.month, '-', 31))
  #   
  # }else{
  #   
  #   if(end.month < start.month) t.end <- lubridate::ymd(paste0(t.year + 1, '-', end.month, '-', 30))
  #     if(end.month >= start.month) t.end <- lubridate::ymd(paste0(t.year, '-', end.month, '-', 30))
  # }
  # 
  # start.dates <- rnames <- seq(t.start, t.end, by = paste0(time.window, ' days'))
  
  # wkstart <- lubridate::floor_date(as.Date(date.start), unit = "week") + 1
  # wkend <- lubridate::floor_date(as.Date(date.end), unit = "week") + 1
  
  # start.dates <- rnames <- seq(wkstart, wkend, by = paste0(time.window, ' days'))
  
  # if(lubridate::ymd(date.end)-lubridate::ymd(date.start)<time.window){
  #   
  #   end.dates <- lubridate::ymd(wkstart) + lubridate::days(time.window-1)
  #   
  # }else{
  #   
  #   end.dates <- seq(lubridate::ymd(wkstart) + lubridate::days(time.window-1), 
  #                    lubridate::ymd(wkend), 
  #                    by = paste0(time.window, ' days'))
  # }
  
  
  #'---------------------------------------------
  # Adjust end dates if necessary
  #'---------------------------------------------
  
  # if(length(end.dates)<length(start.dates)) end.dates <- c(end.dates, start.dates[length(start.dates)]+ days(time.window-1))
  
  #'---------------------------------------------
  # Retain only dates corresponding to data
  #'---------------------------------------------
  
  # if(!is.null(match.df)){
  #   
  #   start.dates <- start.dates[start.dates%in%match.df$wkstart]
  #   end.dates <- end.dates[end.dates%in%match.df$wkstart]
  # }
  
  #'---------------------------------------------
  # Expand date list according to direction of climatology
  #'---------------------------------------------
  
  if(climg.direction == 'backwards'){
    
    start.dates <- purrr::map(.x = start.dates,
                              .f = ~rev(seq(.x, length = no.years, by = "-1 year"))) 
    
    end.dates <- purrr::map(.x = end.dates,
                            .f = ~rev(seq(.x, length = no.years, by = "-1 year"))) 
    
  }else if(climg.direction == 'forwards'){
    
    start.dates <- purrr::map(.x = start.dates,
                              .f = ~seq(.x, length = no.years, by = "1 year"))
    
    end.dates <- purrr::map(.x = end.dates,
                            .f = ~seq(.x, length = no.years, by = "1 year"))
    
  }else if(climg.direction == 'centered'){
    
    start.dates <- purrr::map(.x = start.dates,
                              .f = ~c(rev(seq(.x, length = floor(no.years/2), by = "-1 year")),
                                      seq(.x + years(1), length = ceiling(no.years/2), by = "1 year")))
    
    end.dates <- purrr::map(.x = end.dates,
                            .f = ~c(rev(seq(.x, length = floor(no.years/2), by = "-1 year")),
                                    seq(.x + years(1), length = ceiling(no.years/2), by = "1 year")))
    
  }
  
  #'---------------------------------------------
  # Compile dates
  #'---------------------------------------------
  
  dates.list <- purrr::map2(.x = start.dates,
                            .y = end.dates,
                            .f = ~tibble(.x, .y))
  
  #'---------------------------------------------
  # Correct errors in latitude queries
  #'---------------------------------------------
  
  if(error400){
    
    # Solution found on https://github.com/ropensci/rerddap/issues/68
    
    spacing_string <- unlist(strsplit(datInfo$alldata$latitude$value[1], ","))
    spacing <- unlist(strsplit(spacing_string[3], "="))
    spacing <- as.numeric(spacing[2])
    if (spacing < 0)  latSouth <- FALSE
    
    latVal <- datInfo$alldata$latitude[datInfo$alldata$latitude$attribute_name == "actual_range", "value"]
    latVal2 <- as.numeric(strtrim(strsplit(latVal, ",")[[1]], width = 100))
    tempLat <- paste0(latVal2[2], ',', latVal2[1])
    datInfo$alldata$latitude[datInfo$alldata$latitude$attribute_name == "actual_range", "value"] <- tempLat
  }
  
  #'---------------------------------------------
  # Download the ERDDAP data
  #'---------------------------------------------
  
  message('Downloading data ...')
  
  listlg <- purrr::map(.x = dates.list, .f = ~split(.x, seq(nrow(.x)))) %>% flatten(.) %>% length(.)
  
  pb <- dplyr::progress_estimated(listlg) # Set up progress bar
  
  # Wrapper around rerddap::griddap that prints a progress bar when called with purrr
  
  satellite.dat <- function(data.info, lon.extent, lat.extent, timespan, field.name){
    
    pb$tick()$print()
    
    rerddap::griddap(x = data.info,
                     longitude = lon.extent,
                     latitude = lat.extent,
                     time = timespan, 
                     fields = field.name)
    
  }
  
  dat <- purrr::map(.x = dates.list,
                    .f = ~split(.x, seq(nrow(.x))) %>% 
                      purrr::map(.x = ., .f = ~ as.data.frame(.x)) %>%
                      purrr::map(.x = ., .f = ~c(.$.x, .$.y)) %>% 
                      purrr::map(.x = ., 
                                 .f = ~satellite.dat(data.info = datInfo,
                                                     lon.extent = c(extent(region.shp)[1], 
                                                                    extent(region.shp)[2]),
                                                     lat.extent = c(extent(region.shp)[3], 
                                                                    extent(region.shp)[4]),
                                                     timespan = .x, 
                                                     field.name = variable.name),
                                 .pb = pb))
  
  
  dat <- purrr::map(.x = dat,
                    .f = ~purrr::map(.x = ., .f = ~.x$data) %>% 
                      bind_rows(.) %>% as_tibble(.))
  
  #'---------------------------------------------
  # Calculate climatology
  #'---------------------------------------------
  
  message('\nCreating climatology ...')
  
  # Function to summarise the data into climatology.
  # Prints a progress bar when called with purrr
  
  build.climg <- function(dat, fieldname){
    
    pb$tick()$print()
    
    dat %>% 
      dplyr::group_by(lon, lat) %>% 
      dplyr::summarize_at(.vars = fieldname, .funs = summary.method, na.rm = TRUE) %>% 
      dplyr::ungroup()
    
  }
  
  pb <- dplyr::progress_estimated(length(dat)) # Set up progress bar
  
  climg <- purrr::map(.x = dat, .f = ~build.climg(dat = .x, fieldname = variable.name), .pb = pb)
  
  #'---------------------------------------------
  # Set up empty raster
  #'---------------------------------------------
  
  # Via an extent object derived from the data
  
  e <- climg[[1]] %>% 
    dplyr::select(lon, lat) %>% 
    dplyr::rename(x = lon, y = lat) %>% 
    raster::extent()
  
  r <- raster::raster(e, nrows = sqrt(nrow(climg[[1]])), ncols = sqrt(nrow(climg[[1]])))
  proj4string(r) <- CRSll
  
  message('\nProducing rasters ...')
  
  #'---------------------------------------------
  # Generate rasters
  #'---------------------------------------------
  
  rasOut <- purrr::map(.x = climg,
                       .f = ~ raster::rasterize(x = .x[, c("lon", "lat")],
                                                y = r,
                                                field = .[,c(variable.name)],
                                                fun = mean) %>%
                         raster::projectRaster(from = ., crs = CRSutm) %>%
                         raster::crop(., sp::spTransform(region.shp, CRSutm)) %>%
                         raster::mask(., sp::spTransform(region.shp, CRSutm)))
  
  #'---------------------------------------------
  # Resample rasters
  #'---------------------------------------------
  
  if(rmatch)  rasOut <- purrr::map(.x = rasOut, .f = ~ raster::resample(x = .x, y = raster.dest))
  
  
  # rasOut <- purrr::map(.x = rasOut, .f = ~raster::stack(unlist(.x)))
  
  names(rasOut) <- rnames
  
  # seq(lubridate::ymd(date.start), 
  #                    lubridate::ymd(date.end), 
  #                    by = paste0(time.window, ' days'))
  
  return(rasOut)
  
}


#'---------------------------------------------
# Find calendar week corresponding to a date
#'---------------------------------------------

find_week <- function(input.date, week.lag = NULL){
  
  #'-----------------------------------------------------------
  # PARAMETERS
  #'-----------------------------------------------------------
  #' @param input.date Target date.
  #' @param week.lag Time lag in weeks.
  #'-----------------------------------------------------------
  
  #'---------------------------------------------
  # Perform function checks
  #'---------------------------------------------
  
  if(!class(input.date)=="Date") stop("Input not of class Date")
  
  time.window <- 7
  
  #'---------------------------------------------
  # Time lag
  #'---------------------------------------------
  
  if(!is.null(week.lag)) input.date <- input.date - 7*week.lag
  
  #'---------------------------------------------
  # Find the associate year
  #'---------------------------------------------
  
  date.yr <- lubridate::year(input.date)
  
  # Deal with leap years
  
  if(lubridate::leap_year(input.date)) {
    date.yr <- date.yr - 1
    input.date <- as.Date(paste0(date.yr, "-", lubridate::month(input.date), "-", ifelse(lubridate::day(input.date)==29, 28, lubridate::day(input.date))))
    # input.date <- (input.date - lubridate::days(1))-lubridate::years(1)
  }
  
  #'---------------------------------------------
  # Generate week sequence
  #'---------------------------------------------
  
  t.start <- lubridate::ymd(paste0(date.yr, '-', 1, '-', 1))
  t.end <- lubridate::ymd(paste0(date.yr, '-', 12, '-', 24))
  
  start.dates <- seq(t.start, t.end, by = paste0(time.window, ' days'))
  
  #'---------------------------------------------
  # Find matching week
  #'---------------------------------------------
  
  if(input.date == start.dates[1]){
    gsub(pattern = paste0(date.yr, '-'), replacement = '', x = input.date)
  }else{
    max(start.dates[start.dates<input.date]) %>% 
      gsub(pattern = paste0(date.yr, '-'), replacement = '', x = .)}
}

#'---------------------------------------------
# Extract values from climatology rasters
#'---------------------------------------------

extract_climg <- function(dat,
                          climg,
                          var.name){
  
  #'-----------------------------------------------------------
  # PARAMETERS
  #'-----------------------------------------------------------
  #' @param dat Input data. Must contain a column with the week date as returned by find_week().
  #' @param climg Climatology. Must be a list of rasters as returned by get_wkclimatology().
  #' @param var.name Name of column containing week IDs.
  #'-----------------------------------------------------------
  
  #'---------------------------------------------
  # Perform function checks
  #'---------------------------------------------
  
  if(!class(climg)%in%c("list")) stop('The climatology should be a list object')
  if(!var.name%in%names(dat)) stop('Unrecognised column')
  
  #'---------------------------------------------
  # Extract values from correct raster
  #'---------------------------------------------
  
  pb <- dplyr::progress_estimated(nrow(dat)) # Set up progress bar
  
  purrr::map_dbl(.x = 1:nrow(dat), 
                 .f = ~{
                   pb$tick()$print()
                   climg[names(climg)==dat[.x,] %>% dplyr::pull(var.name)][[1]] %>% 
                     raster::projectRaster(from = ., crs = CRSll) %>% 
                     raster::extract(x = ., y = dat[.x, c('longitude', 'latitude')])
                 }, .pb = pb)
  
} # End function

#'---------------------------------------------
# Calculate the FCPI from a climatology of chlorophyll values
#'---------------------------------------------

calc_fcpi <- function(chla.climg, log.10 = FALSE){
  
  #'--------------------------------------------------------------------
  # PARAMETERS
  #'--------------------------------------------------------------------
  #' @param chla.climg Input climatology.
  #' @param log.10 If TRUE, climatology values are log10 transformed.
  #'--------------------------------------------------------------------
  
  # Step 1: standardise chlorophyll-a values for each pixel using log10 and z-score
  
  if(log.10) chla.climg <- purrr::map(.x = chla.climg, .f = ~log10(.x)) 
  
  chla.stack <- raster::stack(chla.climg)
  chla.mean <- raster::calc(chla.stack, fun = function(x) mean(x, na.rm = TRUE))
  chla.sd <- raster::calc(chla.stack, fun = function(x) sd(x, na.rm = TRUE))
  
  chla.1 <- (chla.stack - chla.mean)/chla.sd
  
  # Step 2: spatial means for each month
  
  chla.2 <- purrr::map_dbl(.x = as.list(chla.1), .f = ~mean(raster::getValues(.x), na.rm = TRUE))
  chla.2 <- tibble::enframe(chla.2) %>% 
    dplyr::mutate(week = dplyr::row_number()) %>% 
    dplyr::mutate(x1 = sin(week*2*pi/52), 
                  x2 = cos(week*2*pi/52),
                  x3 = sin(week*2*pi/26),
                  x4 = cos(week*2*pi/26))
  
  chla.mod <- lm(value ~ x1 + x2 + x3 + x4 + week, data = chla.2) 
  chla.preds <- predict(chla.mod, newdata = chla.2, se = TRUE)
  
  par(mfrow = c(2,1))
  plot(value ~ week, data = chla.2, ylab = "Standardised Chl-a", xlab = "Week")
  lines(chla.2$week, chla.preds$fit, col = "orange")
  lines(chla.2$week, chla.preds$fit+chla.preds$se.fit, col = "darkgrey", lty = 2)
  lines(chla.2$week, chla.preds$fit-chla.preds$se.fit, col = "darkgrey", lty = 2)
  
  # Step 3: Indentify % time > 1 SD of area-wide model
  
  upper <- chla.preds$fit+1*chla.preds$se.fit
  chla.1 <- as.list(chla.1)
  
  chla.3 <- purrr::map2(.x = chla.1, .y = upper, .f = ~.x > .y)
  chla.4 <- raster::stack(chla.3)
  chla.4 <- raster::calc(x = chla.4, fun = function(x) 100*sum(x, na.rm = TRUE)/length(chla.3))
  chla.4[chla.4==0] <- NA
  
  plot(chla.4, col = viridis(100))
  
  return(list(fcpi = chla.4, means = chla.2, model = chla.mod, preds = chla.preds))
}

#'---------------------------------------------
# Assign weights to presence and pseudo-absence points
#'---------------------------------------------

# Taken from the biomod2 package

automatic_weights_creation <- function(resp, prev = 0.5, subset = NULL){
  
  #'--------------------------------------------------------------------
  # PARAMETERS
  #'--------------------------------------------------------------------
  #' @param resp Vector of presences and absences, coded as 1 and 0 respectively.
  #' @param prev Desired prevalence (defaults to 0.5).
  #' @param subset Vector indicating which elements of resp to subset for calculations.
  #'--------------------------------------------------------------------
  
  if(is.null(subset)) subset <- rep(TRUE, length(resp))
  
  nbPres <- sum(resp[subset], na.rm = TRUE)
  
  # The number of true absences + pseudo absences to maintain true value of prevalence
  nbAbsKept <- sum(subset, na.rm = TRUE) - sum(resp[subset], na.rm = TRUE) 
  Yweights <- rep(1, length(resp))
  
  if(nbAbsKept > nbPres){ # code absences as 1
    Yweights[which(resp>0)] <- (prev * nbAbsKept) / (nbPres * (1-prev))
  } else{ # code presences as 1
    Yweights[which(resp==0 | is.na(resp))] <- (nbPres * (1-prev)) / (prev * nbAbsKept)
  }
  Yweights = round(Yweights[])
  Yweights[!subset] <- 0
  
  return(Yweights)
}

#'---------------------------------------------
# Plot spline correlogram and semivariogram
#'---------------------------------------------

plot_correlogram <- function(model, data, maxD = 15000, nboot = 500){
  
  #'--------------------------------------------------------------------
  # PARAMETERS
  #'--------------------------------------------------------------------
  #' @param model Input model, from which to extract residuals.
  #' @param dat Data object containing x and y coordinates for each data point.
  #' @param maxD Maximum range of the correlogram (m).
  #' @param nboot Number of bootstrap iterations.
  #'--------------------------------------------------------------------
  
  # For the data
  
  Correlog.data <- ncf::spline.correlog(x = data$x, y = data$y, z = data$presence, xmax = maxD, resamp = nboot)
  
  # For model residuals
  
  Correlog.residuals <- ncf::spline.correlog(x = data$x, y = data$y, z = residuals(model, type = "pearson"), xmax = maxD, resamp = nboot)
  
  # Plot results
  par(mfrow = c(2,1))
  ncf::plot.spline.correlog(Correlog.data)
  ncf::plot.spline.correlog(Correlog.residuals)
  par(mfrow = c(1,1))
}

plot_variogram <- function(model, dat, maxD = 15000){
  
  #'--------------------------------------------------------------------
  # PARAMETERS
  #'--------------------------------------------------------------------
  #' @param model Input model, from which to extract residuals.
  #' @param dat Data object containing x and y coordinates for each data point.
  #' @param maxD Range over which to assess the semi-variance (m).
  #'--------------------------------------------------------------------
  
  vario.data <- data.frame(residuals(model, type = "pearson"), x = dat$x, y = dat$y)
  coordinates(vario.data) <- c("x", "y")
  bw.variogram <- gstat::variogram(residuals(model, type = "pearson") ~ 1, vario.data, 
                                   cutoff = maxD,
                                   alpha = c(0, 45, 90, 135))
  plot(bw.variogram)
  
}

#'---------------------------------------------
# Extract model selection results
#'--------------------------------------------- 

# Modified version of FSSgam::extract.mod.dat, which also includes AIC in addition to AICc

extract.model.dat <- function (mod.fit, r2.type. = r2.type) {
  
#'--------------------------------------------------------------------
# PARAMETERS
#'--------------------------------------------------------------------
#' See the FSSgam package.
#'--------------------------------------------------------------------
 
   mod.dat = list(AICc = NA, AIC = NA, BIC = NA, r2.vals = NA, r2.vals.unique = NA, 
                 edf = NA, edf.less.1 = NA)
  if (class(mod.fit)[[1]] != "try-error") {
    # mod.dat$AIC = MuMIn::QAIC(mod.fit, chat = mod.fit$deviance/mod.fit$df.residual)
    mod.dat$AIC = AIC(mod.fit)
    mod.dat$AICc = AICc(mod.fit)
    mod.dat$BIC = BIC(mod.fit)
    tempOut = NA
    if (class(mod.fit)[1] == "gam" & r2.type. == "dev") {
      tempOut = summary(mod.fit)$dev.expl
    }
    if (class(mod.fit)[1] == "gam" & r2.type. == "r2") {
      tempOut = summary(mod.fit)$r.sq
    }
    if (class(mod.fit)[1] == "gam" & r2.type. == "r2.lm.est") {
      tempOut = summary(lm(mod.fit$y ~ predict(mod.fit)))$r.sq
    }
    if (class(mod.fit)[[1]] == "gamm4" & r2.type. == "dev") {
      tempOut = summary(mod.fit$gam)$dev.expl
      if (length(tempOut) == 0) {
        tempOut = NA
      }
    }
    if (class(mod.fit)[[1]] == "gamm4" & r2.type. == "r2") {
      tempOut = summary(mod.fit$gam)$r.sq
    }
    if (class(mod.fit)[[1]] == "gamm" & r2.type. == "r2") {
      tempOut = summary(mod.fit$gam)$r.sq
    }
    if (class(mod.fit)[[1]] == "gamm4" & r2.type. == "r2.lm.est") {
      tempOut = summary(lm(attributes(mod.fit$mer)$frame$y ~ 
                             predict(mod.fit[[1]], re.form = NA, type = "response")))$r.sq
    }
    if (is.null(tempOut)) {
      tempOut = NA
    }
    mod.dat$r2.vals = round(tempOut, 5)
    if (class(mod.fit)[1] == "gam") {
      edf.m = summary(mod.fit)$edf
      p.coeff.m = summary(mod.fit)$p.coeff
    }
    else {
      edf.m = mod.fit$gam$edf
      p.coeff.m = mod.fit$gam$p.coeff
    }
    edf.m[which(edf.m < 1)] = 1
    mod.dat$edf = round(sum(c(edf.m, length(p.coeff.m))), 
                        2)
    if (class(mod.fit)[1] == "gam") {
      edf.m = summary(mod.fit)$edf
    }
    else {
      edf.m = mod.fit$gam$edf
    }
    mod.dat$edf.less.1 = length(which(edf.m < 0.25))
  }
  return(mod.dat)
}

#'---------------------------------------------
# Full subsets listing
#'--------------------------------------------- 

list.model.sets <- function(response.var, 
                            dat,
                            covariates,
                            use.weights = TRUE,
                            corr.cutoff = 0.28, 
                            max.predictors, 
                            verbose = TRUE,
                            k = NULL){
  
  #'--------------------------------------------------------------------
  # PARAMETERS
  #'--------------------------------------------------------------------
  #' @param response.var Response variable.
  #' @param dat Input data.
  #' @param covariates Input covariates.
  #' @param include.xy Whether to weigh presence and pseudo-absence points.
  #' @param corr.cutoff Maximm correlation between covariates.
  #' @param max.predictors Maximum number of predictors in a single model.
  #' @param verbose If TRUE, prints function status during execution.
  #' @param k Basis dimension for GAM smoothers. Sets the upper limit on the degrees of freedom associated with an s smooth.
  #'--------------------------------------------------------------------
  
  # Prepare the data
  
  if(verbose) message('Preparing data ...')
  
  use.dat <- dat %>% 
    dplyr::rename(response = response.var)
  
  if(use.weights) gamweights <<- use.dat %>% dplyr::pull(weights)
  
  # Initial model
  
  if(verbose) message('Building null model ...')
  
  # Build model formula
  
  if(is.null(k)){
    gterms <- paste0("s(", covariates, ") ", collapse = "+ ")
  }else{
    gterms <- paste0("s(", covariates, ", k = ", k," ) ", collapse = "+ ")}

  
  if(use.weights){
    
    null.model <- mgcv::gam(formula = as.formula(paste0("response ~ ", gterms, collapse = "+")), data = use.dat, weights = gamweights, REML = TRUE, family = binomial(link = "logit"))
    
  }else{
    
    null.model <- mgcv::gam(formula = as.formula(paste0("response ~ ", gterms, collapse = "+")), data = use.dat, REML = TRUE, family = binomial(link = "logit"))}
  
  if(verbose) message('Generating model set ...')
  
  if(is.null(k)){
    
    model.set <- FSSgam::generate.model.set(use.dat = use.dat,
                                            test.fit = null.model,
                                            cov.cutoff = corr.cutoff,
                                            non.linear.correlations = FALSE,
                                            pred.vars.cont = covariates,
                                            pred.vars.fact = NA, 
                                            max.predictors = max.predictors,
                                            bs.arg = "'tp'")
  }else{
    
    model.set <- FSSgam::generate.model.set(use.dat = use.dat,
                                            test.fit = null.model,
                                            cov.cutoff = corr.cutoff,
                                            non.linear.correlations = FALSE,
                                            pred.vars.cont = covariates,
                                            pred.vars.fact = NA, 
                                            max.predictors = max.predictors,
                                            k = k,
                                            bs.arg = "'tp'")
  }
  
  if(verbose) message('Done!')
  return(model.set)
  
}

#'---------------------------------------------
# Full subsets model fitting
#'--------------------------------------------- 

# Modified version of FSSgam::fit.mod.set

fit.model.set <- function (model.set.list, max.models = 200, save.model.fits = T, parallel = F, n.cores = 4, r2.type = "r2.lm.est", report.unique.r2 = F) 
{
  
  #'--------------------------------------------------------------------
  # PARAMETERS
  #'--------------------------------------------------------------------
  #' See FSSgam for details.
  #'--------------------------------------------------------------------
  
  use.datModSet <- model.set.list$used.data
  n.mods = length(model.set.list$mod.formula)
  mod.formula = model.set.list$mod.formula
  test.fit = model.set.list$test.fit
  included.vars = model.set.list$included.vars
  if (n.mods > max.models) {
    warning(paste("You have ", n.mods, " models. Individual models fits will not be saved.\n        If you want to save all the model fits all of these you need to\n        increase 'max.models' from ", 
                  max.models, ".", sep = ""))
    save.model.fits = F
  }
  require(MuMIn)
  if (save.model.fits == T) {
    pb <- txtProgressBar(max = length(mod.formula), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    if (parallel == T) {
      require(doSNOW)
      cl = makeCluster(n.cores)
      registerDoSNOW(cl)
      opts <- list(progress = progress)
      out.dat <- foreach(l = 1:length(mod.formula), .packages = c("mgcv", 
                                                                  "gamm4", "MuMIn", "FSSgam"), .errorhandling = "pass", 
                         .options.snow = opts) %dopar% {
                           fit.mod.l(mod.formula[[l]], test.fit. = test.fit, 
                                     use.dat = use.datModSet)
                         }
      close(pb)
      stopCluster(cl)
      registerDoSEQ()
    }
    else {
      out.dat <- vector("list", length(mod.formula))
      for (l in 1:length(mod.formula)) {
        mod.l = fit.mod.l(mod.formula[[l]], test.fit. = test.fit, 
                          use.dat = use.datModSet)
        out.dat[[l]] = mod.l
        setTxtProgressBar(pb, l)
      }
    }
    close(pb)
    names(out.dat) = names(mod.formula)[1:n.mods]
    model.success = lapply(lapply(out.dat, FUN = class), 
                           FUN = function(x) {
                             length(grep("gam", x)) > 0
                           })
    failed.models = mod.formula[which(model.success == F)]
    success.models = out.dat[which(model.success == T)]
    if (length(success.models) == 0) {
      stop("None of your models fitted successfully. Please check your input objects.")
    }
    var.inclusions = build.inclusion.mat(included.vars = included.vars, 
                                         formula.list = success.models)
    mod.data.out = data.frame(modname = names(success.models))
    mod.data.out$formula = unlist(lapply(success.models, 
                                         FUN = function(x) {
                                           as.character(formula(x)[3])
                                         }))
    mod.data.out = cbind(mod.data.out, do.call("rbind", lapply(success.models, 
                                                               FUN = function(x) {
                                                                 unlist(extract.model.dat(x, r2.type. = r2.type))
                                                               })))
  }
  else {
    var.inclusions = build.inclusion.mat(included.vars = included.vars, 
                                         formula.list = mod.formula)
    mod.data.out = data.frame(modname = names(mod.formula))
    mod.data.out$formula = unlist(lapply(mod.formula, FUN = function(x) {
      as.character(formula(x))[2]
    }))
    pb <- txtProgressBar(max = length(mod.formula), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    if (parallel == T) {
      require(doSNOW)
      cl = makeCluster(n.cores)
      registerDoSNOW(cl)
      opts <- list(progress = progress)
      mod.dat <<- foreach(l = 1:length(mod.formula), .packages = c("mgcv", 
                                                                   "gamm4", "MuMIn", "FSSgam"), .options.snow = opts) %dopar% 
        {
          unlist(extract.model.dat(fit.mod.l(mod.formula[[l]], 
                                             test.fit. = test.fit, use.dat = use.datModSet), 
                                   r2.type. = r2.type))
        }
      close(pb)
      stopCluster(cl)
      registerDoSEQ()
    }
    else {
      mod.dat = vector("list", length(mod.formula))
      for (l in 1:length(mod.formula)) {
        mod.l = fit.mod.l(mod.formula[[l]], test.fit. = test.fit, 
                          use.dat = use.datModSet)
        out = unlist(extract.model.dat(mod.l, r2.type. = r2.type))
        mod.dat[[l]] = out
        setTxtProgressBar(pb, l)
      }
    }
    close(pb)
    names(mod.dat) = names(mod.formula[1:n.mods])
    mod.data.out = cbind(mod.data.out, do.call("rbind", mod.dat))
    failed.models = mod.formula[which(is.na(mod.data.out$AICc) == 
                                        T)]
    success.models = mod.formula[which(is.na(mod.data.out$AICc) == 
                                         F)]
  }
  mod.data.out$delta.AIC = round(mod.data.out$AIC - min(mod.data.out$AIC, 
                                                        na.rm = T), 3)
  
  mod.data.out$delta.AICc = round(mod.data.out$AICc - min(mod.data.out$AICc, 
                                                          na.rm = T), 3)
  mod.data.out$delta.BIC = round(mod.data.out$BIC - min(mod.data.out$BIC, 
                                                        na.rm = T), 3)
  mod.data.out$wi.AIC = round(wi(mod.data.out$AIC), 3)
  mod.data.out$wi.AICc = round(wi(mod.data.out$AICc), 3)
  mod.data.out$wi.BIC = round(wi(mod.data.out$BIC), 3)
  if (report.unique.r2 == T) {
    null.r2 = mod.data.out$r2.vals[which(mod.data.out$modname == 
                                           "null")]
    mod.data.out$r2.vals.unique = mod.data.out$r2.vals - 
      null.r2
  }
  mod.data.out = cbind(mod.data.out, var.inclusions)
  min.mods = min(colSums(mod.data.out[, included.vars]))
  var.weights = unlist(lapply(included.vars, FUN = function(x) {
    sum(sort(mod.data.out$wi.AICc[which(mod.data.out[, x] == 
                                          1)], decreasing = T)[1:min.mods])
  }))
  var.weights = unlist(lapply(included.vars, FUN = function(x) {
    sum(sort(mod.data.out$wi.AIC[which(mod.data.out[, x] == 
                                         1)], decreasing = T)[1:min.mods])
  }))
  names(var.weights) = included.vars
  variable.weights.raw = var.weights
  aic.var.weights = list(variable.weights.raw = variable.weights.raw)
  var.weights = unlist(lapply(included.vars, FUN = function(x) {
    sum(sort(mod.data.out$wi.BIC[which(mod.data.out[, x] == 
                                         1)], decreasing = T)[1:min.mods])
  }))
  names(var.weights) = included.vars
  variable.weights.raw = var.weights
  bic.var.weights = list(variable.weights.raw = variable.weights.raw)
  return(list(mod.data.out = mod.data.out, failed.models = failed.models, 
              success.models = success.models, variable.importance = list(aic = aic.var.weights, 
                                                                          bic = bic.var.weights)))
}

#'---------------------------------------------
# Compile results from full subsets 
#'--------------------------------------------- 

FSS_results <- function(FSS.list, plot.models, AIC = FALSE, verbose = TRUE){
  
  #'--------------------------------------------------------------------
  # PARAMETERS
  #'--------------------------------------------------------------------
  #' @param FSS.list Output from fit.model.set.
  #' @param plot.models Whether to plot smooths for the top-scoring models.
  #' @param AIC If TRUE, model ranking/selection is based on the AIC, rather than the AICc.
  #' @param verbose Whether to print function status during execution.
  #'--------------------------------------------------------------------
  
  fss.obj <- list()
  
  # Retrieve failed and successful models
  
  fss.obj$failed.models <- FSS.list$failed.models
  fss.obj$success.models <- FSS.list$success.models
  
  if(verbose) message(paste0('Number of failed models: ', length(fss.obj$failed.models), ' (', 100 * round(length(fss.obj$failed.models)/(length(fss.obj$failed.models)+length(fss.obj$success.models))), '%)'))
  
  if(verbose) message(paste0('Number of successful models: ', length(fss.obj$success.models), ' (', 100 * round(length(fss.obj$success.models)/(length(fss.obj$failed.models)+length(fss.obj$success.models))), '%)'))
  
  # Draw comparative table
  
  fss.obj$model.table <- FSS.list$mod.data.out
  
  if(AIC)  fss.obj$model.table <-   fss.obj$model.table[order(fss.obj$model.table$AIC),]
  if(!AIC)  fss.obj$model.table <-   fss.obj$model.table[order(fss.obj$model.table$AICc),]
  
  # Variable importance
  
  fss.obj$var.importance <- FSS.list$variable.importance$aic$variable.weights.raw
  fss.obj$var.importance <- sort(fss.obj$var.importance, decreasing = TRUE)
  
  # All models within 2 AICc units
  
  if(AIC) fss.obj$top.models <- fss.obj$model.table[which(fss.obj$model.table$delta.AIC<2),]
  if(!AIC) fss.obj$top.models <- fss.obj$model.table[which(fss.obj$model.table$delta.AICc<2),]
  
  # Model with minimum AICc
  
  if(AIC)   fss.obj$best.model <- list(summary = fss.obj$model.table[which(fss.obj$model.table$AIC==min(fss.obj$model.table$AIC, na.rm = TRUE)),])
  if(!AIC)   fss.obj$best.model <- list(summary = fss.obj$model.table[which(fss.obj$model.table$AICc==min(fss.obj$model.table$AICc, na.rm = TRUE)),])
  
  fss.obj$best.model$obj <- fss.obj$success.models[[as.character(fss.obj$best.model$summary$modname)]]
  
  if(verbose) message(paste0('Top-ranking model: ', fss.obj$best.model$summary$modname))
  if(verbose) message(paste0('Most important variable: ', names(fss.obj$var.importance)[1]))
  
  if(plot.models){
    
    par(oma = c(1,1,4,1))
    
    for(r in 1:nrow(fss.obj$top.models)){
      
      best.model.name <- as.character(fss.obj$top.models$modname[r])
      best.model <- FSS.list$success.models[[best.model.name]]
      
      plot(best.model, all.terms = TRUE, pages = 1, residuals = FALSE, scale = 0, shade = TRUE)
      mtext(side = 3, text = resp.vars[i], outer = TRUE)
    }
  }
  
  return(fss.obj)
  
}

#'---------------------------------------------
# Make predictions from a fitted model
#'---------------------------------------------

get_predictions <- function(model, time.span){
  
  #'-----------------------------------------------------------
  # PARAMETERS
  #'-----------------------------------------------------------
  #' @param model Fitted model object.
  #' @param time.span Time period for which predictions are sought (week numbers).
  #'-----------------------------------------------------------
  
  pb <- dplyr::progress_estimated(length(time.span)) # Set up progress bar
  r.out <- purrr::map(.x = time.span,
                      .f = ~{
                        pb$tick()$print()
                        # preds.stack <- raster::stack(static.env.utm, r.x, r.y, chla.rs[[.x]])
                        # names(preds.stack) <- c(names(static.env.utm), "x", "y", "chla.0")
                        preds.stack <- raster::stack(static.env.utm, chla.rs[[.x]])
                        names(preds.stack) <- c(names(static.env.utm), "chla.max")
                        predict(object = preds.stack, model = model, type = "response")},
                      .pb = pb)
  raster::stack(r.out)
}

#'---------------------------------------------
# Make a map of an input raster over the study area
#'---------------------------------------------

make_map <- function(input.raster, 
                     col.ramp = pals::parula(100),
                     show.canyons = FALSE,
                     show.tracks = FALSE,
                     show.sightings = FALSE,
                     sighting.size = 1,
                     plot.title,
                     legend.title = TRUE,
                     legend.size = 0.5,
                     legend.x = 0.1,
                     legend.y = 0.2,
                     axis.txt.x = 12,
                     axis.txt.y = 12){
  
  #'--------------------------------------------------------------------
  # PARAMETERS
  #'--------------------------------------------------------------------
  #' @param input.raster Input raster.
  #' @param col.ramp Colour scheme for the raster.
  #' @param show.canyons If TRUE, overlay canyon locations.
  #' @param show.tracks If TRUE, overlay survey tracklines.
  #' @param show.sightings If TRUE, overlay survey sightings.
  #' @param plot.title Plot title.
  #' @param legend.title Title of legend.
  #' @param legend.size Size of legend text area.
  #' @param legend.x Horizontal placement of legend.
  #' @param legend.y Vertical placement of legend.
  #' @param axis.txt.x Font size for x-axis.
  #' @param axis.txt.y Font size for y-axis.
  #'--------------------------------------------------------------------
  
  #'---------------------------------------------
  # Name of raster as legend title
  #'---------------------------------------------
  
  raster.name <- deparse(substitute(input.raster))
  raster.name <- paste(toupper(substr(raster.name, 1, 1)), 
                       substr(raster.name, 2, nchar(raster.name)), sep="")
  raster.name <- gsub(pattern = "_", replacement = " ", raster.name)
  
  #'---------------------------------------------
  # Ensure that raster is in lat/lon
  #'---------------------------------------------
  
  if(!proj4string(input.raster)==as.character(CRSll)){
    input.raster <- raster::projectRaster(from = input.raster, crs = CRSll)
  }
  
  #'---------------------------------------------
  # Stops execution if basemap does not exist
  #'---------------------------------------------
  
  if(exists("gmap.timor")==FALSE) stop("Missing basemap")
  
  if(show.canyons) canyons.plot <- fortify(canyons)
  if(show.tracks){
    gps.07f <- fortify(gps.07)
    gps.08f <- fortify(gps.08)}
  
  
  #'---------------------------------------------
  # Convert to data.frame for plotting
  #'---------------------------------------------
  
  input.raster <- raster::as.data.frame(input.raster, xy = TRUE)
  names(input.raster) <- c("lon", "lat", "z")
  
  suppressMessages(ggmap(gmap.timor)+ # basemap
                     
                     coord_equal() + # Needs to be ### to produce map in right dimension pair
                     
                     geom_tile(data = input.raster, aes(x = lon, y = lat, fill = z))+ # raster
                     
                     {if(show.canyons) geom_path(data = canyons.plot, aes(long, lat, group=group))}+
                     
                     # GPS tracks
                     
                     {if(show.tracks) geom_path(data = gps.07f, aes(long, lat, group = group), alpha = 0.25)}+
                     {if(show.tracks) geom_path(data = gps.08f, aes(long, lat, group = group), alpha = 0.25)}+
                     
                     # Sightings
                     
                     {if(show.sightings) geom_point(data = bw, 
                                                    aes(longitude, latitude), 
                                                    pch = 21, fill = "black", size = sighting.size, alpha = 1)}+
                     
                     scale_fill_gradientn(colors = col.ramp,
                                          na.value = 'transparent',
                                          limits = range(input.raster$z),
                                          breaks = pretty(input.raster$z),
                                          labels = pretty(input.raster$z))+
                     
                     {if(legend.title) labs(fill = raster.name)} +
                     {if(!legend.title) labs(fill = "")} +
                     
                     ggtitle(plot.title) +
                     
                     xlab("")+
                     ylab("")+
                     
                     scale_x_continuous(limits = range(xval), 
                                        breaks= xval,
                                        labels = lab.x)+
                     
                     scale_y_continuous(limits = range(-yval), 
                                        breaks= rev(-yval),
                                        labels = rev(lab.y), expand = c(0,0))+
                     
                     theme_sleek() + # ggsidekick magic happens here
                     
                     theme(axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black", size = axis.txt.y),
                           axis.text.x = element_text(size = axis.txt.x, colour = "black"),
                           plot.title = element_text(colour = "black"),
                           panel.border = element_rect(colour = "black", fill = NA, size=0.8),
                           legend.key = element_rect(fill = "transparent"),
                           legend.position = c(legend.x, legend.y),
                           legend.key.width = unit(legend.size,"cm"),
                           legend.background = element_rect(fill = "transparent", size = 2),
                           legend.text = element_text(size = 12, colour = "black"),
                           legend.key.size = unit(0.75,"cm")) +
                     
                     {if(legend.title) theme(legend.title = element_text(size = 13, face = "bold", colour = "black"))})
  
}