#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------
#' 
#'  Winter distribution and habitat use of blue whales (Balaenoptera musculus sp.) 
#'  in the Timor Trough, south of Timor-Leste during 2007-08
#'
#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------

# Register Google API Key
# ggmap::register_google(key = "xxxxxx")

# Table 1 ---------------------------------------------------------------

# Number of sightings
sightings.per.month <- bw %>% 
  dplyr::group_by(year, month) %>% 
  tally(species == "Blue whale") %>% dplyr::ungroup()

# Number of calves
bw %>% 
  dplyr::group_by(year, month) %>% 
  tally(calves) %>% dplyr::ungroup()

# Number of individuals
bw %>% 
  dplyr::group_by(year, month) %>% 
  tally(adcalves) %>% dplyr::ungroup()

# Number of effort days

gps$`2007` %>% 
  dplyr::mutate(month = factor(lubridate::month(date))) %>%
  dplyr::group_by(year, month) %>% 
  dplyr::summarise(effort_days = dplyr::n_distinct(date), min_date = min(date), max_date = max(date)) %>% dplyr::ungroup()

gps$`2008` %>% 
  dplyr::mutate(month = factor(lubridate::month(date))) %>%
  dplyr::group_by(year, month) %>% 
  dplyr::summarise(effort_days = dplyr::n_distinct(date), min_date = min(date), max_date = max(date)) %>% dplyr::ungroup()

gps.lines <- rbind(gps$`2007`, gps$`2008`) %>% 
  split(x = ., f = .$date) %>% 
  purrr::map(., ~.x %>% split(x = ., f = .$line_id) %>% purrr::map(., ~ createLines(.))) %>% 
  purrr::map_depth(.x = ., .depth = 2, .f = ~sp::spTransform(.x, CRSobj = CRSutm)) %>% 
  purrr::map(.x = ., .f = ~do.call(rbind, .)) %>% 
  purrr::map(.x = ., .f = ~smoothr::drop_crumbs(x =.x, threshold = 250))

depth.along.lines <- 
 purrr::map(.x = gps.lines, 
            .f = ~{
              tmp <- sp::spsample(.x, n = 1000, type = "regular") %>%
                raster::extract(x = depth, y = .)
              range(tmp)
                }) %>% tibble::enframe() %>% 
  dplyr::mutate(minDepth = purrr::map_dbl(.x = value, .f = ~.x[[1]]),
                maxDepth = purrr::map_dbl(.x = value, .f = ~.x[[2]])) %>%
  dplyr::mutate(date = lubridate::as_date(name)) %>%
  dplyr::select(-value, name) %>%
  dplyr::mutate(year = lubridate::year(date), month = lubridate::month(date))
  
depth.along.lines %>% 
  dplyr::group_by(year, month) %>% 
  dplyr::summarise(minDepth = min(minDepth), maxdepth = max(maxDepth)) %>%
  dplyr::ungroup()

daily.effort <- purrr::map_df(.x = gps.lines, .f = ~rgeos::gLength(.x)) %>% 
  t() %>% as_tibble(rownames = "date") %>% 
  dplyr::rename(length = V1) %>% 
  dplyr::mutate(km = length / 1000,
                month = lubridate::month(date),
                year = lubridate::year(date))

monthly.tally <- daily.effort %>% 
  dplyr::group_by(year, month) %>% 
  dplyr::summarise(total = sum(km)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(x = ., y = sightings.per.month, by = c("year", "month")) %>% 
  dplyr::mutate(n = ifelse(is.na(n), 0, n))


# Figure 1 (R version) ---------------------------------------------------------------

# Code taken from https://egallic.fr/en/maps-with-r/

# World map
worldMap <- rworldmap::getMap(resolution = "high")
world.points <- ggplot2::fortify(worldMap)
world.points$region <- world.points$id
world.df <- world.points[,c("long","lat","group", "region")]

w <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group), 
               fill = "grey20", color = "grey20") +
  coord_map("ortho", orientation = c(12, 125, -2)) + 
  theme_void() + xlab("") + ylab("") +
  theme(panel.background = element_rect(fill = "grey95", colour = "grey95"))

ggsave("/Users/philippebouchet/Google Drive/Documents/git/timorbw/fig/globe.pdf", plot = w)

#'---------------------------------------------
# Basemap
#'---------------------------------------------

gmap.timor <- readRDS("/Users/philippebouchet/Google Drive/Documents/git/timorbw/gis/timor_gmap.rds")

ggmap(gmap.timor) # Quick visualisation

#'---------------------------------------------
# Axis labels
#'---------------------------------------------

xval <- seq(124.75,127.25,0.5)
yval <- seq(8,10.5,0.5)
lab.x <- c(paste(xval, "°E",sep=""))
lab.y <- c(paste(yval, "°S",sep=""))

#'---------------------------------------------
# ggplot map
#'---------------------------------------------

fig.1 <- ggmap::ggmap(gmap.timor) + # basemap
  
  # GPS tracks
  geom_path(data = fortify(gps.07), aes(long, lat, group = group), colour = "#d8b365", alpha = 0.85)+
  geom_path(data = fortify(gps.08), aes(long, lat, group = group), colour = "#5ab4ac", alpha = 0.85)+
  
  coord_equal() + # Needs to be ### to produce map in right dimension pair
  
  geom_point(data = bw, aes(longitude, latitude, fill = as.factor(year)), pch = 21, size = 2, alpha = 0.85)+
  
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"), name = NULL)+
  
  xlab("")+
  ylab("")+
  
  scale_x_continuous(limits = range(xval), 
                     breaks= xval,
                     labels = lab.x)+
  
  scale_y_continuous(limits = range(-yval), 
                     breaks= rev(-yval),
                     labels = rev(lab.y), expand = c(0,0))+
  
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
  
  ggplot2::guides(fill = guide_legend(override.aes = list(size=6))) +
  
  
  #'---------------------------------------------
  # North arrow and scale bar
  #'---------------------------------------------
  
  # Data argument HAS to match first data set of ggplot/ggmap
  # Also requires coord_equal otherwise returns an error about need for Cartesian coordinates
  
  ggsn::north(data = fortify(gps.07),
              location = "bottomleft",
              symbol = 12, 
              scale = 0.15,
              anchor = c(x = 126.95, y = -10.2)) +
  
  ggsn::scalebar(data = fortify(gps.07),
                 x.min = 125,
                 x.max = 127.5,
                 y.min = -10.5,
                 y.max = -9,
                 dist = 25, 
                 transform = TRUE,
                 dist_unit = "km",
                 height = 0.03,
                 model = "WGS84",
                 anchor = c(x = 127.25, y = -10.3),
                 st.dist = 0.05, border.size = 0.2, nudge_y = -0.025) +
  
  annotate("text", 
           label = "N", 
           x = 127.025, 
           y = -10, 
           size = 5) + 
  
  theme(axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))


ggsave(filename = file.path("fig/Figure1.png"), width = 25, height = 20, units = "cm")


# Figure 2 ----------------------------------------------------------------

png("fig/Figure2.png", res = 450, width = 2500, height = 2750)
par(mfrow = c(2,2))
purrr::walk(.x = 1:4, .f = ~{
  plot(bestmodel, shade = TRUE, scale = 0, select = .x, xlab = c("Max chlorophyll-a concentration", "Distance to incising canyon", "FCPI", "Slope")[.x]); abline(h = 0, lty = 2)})
dev.off()

# Figure S1 ---------------------------------------------------------------

#'---------------------------------------------
# (A) Time series
#'---------------------------------------------

png("fig/Figure-S1a.png", res = 450, height = 1800, width = 3500)
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

#'---------------------------------------------
# (B) Wavelet power spectrum
#'---------------------------------------------

png("fig/Figure-S1b.png", res = 450, height = 2500, width = 3500)
WaveletComp::wt.image(wavelet.timor, 
                      plot.coi = TRUE,
                      plot.contour = TRUE,
                      exponent = 1,
                      siglvl = 0.1,
                      color.key = "interval", 
                      color.palette = "rev(pals::parula(n.levels))",
                      n.levels = 100,
                      show.date = TRUE,
                      timelab = "",
                      spec.time.axis = list(at = time.labels, labels = lubridate::year(time.labels)),
                      periodlab = "",
                      legend.params = list(n.ticks = 15, label.digits = 1))
dev.off()


# gmap.timor <- ggmap::get_googlemap(center = c(lon=126,lat=-9.4),
#                           zoom=8,
#                           size = c(640, 640),
#                           scale=2,
#                           maptype="roadmap",
#                           style = 'feature:road|element:all|visibility:simplified&style=feature:administrative.locality|element:labels|visibility:simplified')


# Figure S2 ---------------------------------------------------------------

# Covariate maps

covariate.maps <- purrr::pmap(list(
  
  x = list(depth = depth, 
           slope = seabed_slope, 
           dcoast = dist_coast, 
           dcanyons = dist_canyons, 
           dincising = dist_incising,
           sst = raster::stack(sst_climg) %>% mean(., na.rm = TRUE),
           chla = raster::stack(chla_climg) %>% mean(., na.rm = TRUE)), 
  
  # Colour ramps
  y = list(depth = rev(pals::brewer.blues(100)),
           slope = pals::viridis(100),
           dcoast = pals::cividis(100),
           dcanyons = pals::parula(100),
           dincising = pals::parula(100),
           sst = pals::coolwarm(100),
           chla = pals::ocean.speed(100)),
  
  z = as.list(c('Depth (m)', 'Seabed slope (°)', 'Distance to coast (km)',
                'Distance to canyon (km)', 'Distance to incising canyon (km)', 'Sea surface temperature (°C)',
                'log10(chla) (mg/m3)'))),
  
  .f = function(x, y, z) { make_map(input.raster = x, 
                                    col.ramp = y, 
                                    show.canyons = FALSE, 
                                    show.tracks = TRUE, 
                                    show.sightings = TRUE,
                                    sighting.size = 1,
                                    plot.title = z,
                                    legend.title = FALSE,
                                    legend.size = 0.5,
                                    legend.x = 0.1,
                                    legend.y = 0.2,
                                    axis.txt.x = 12,
                                    axis.txt.y = 12)})
require(patchwork)
fig.s2 <- cowplot::ggdraw() +
  cowplot::draw_plot(covariate.maps$depth, x = 0, y = 0.6, width = 0.5, height = 0.5) +
  cowplot::draw_plot(covariate.maps$slope, x = 0.5, y = 0.6, width = 0.5, height = 0.5) +
  cowplot::draw_plot(covariate.maps$dcoast, x = 0, y = 0.4, width = 0.5, height = 0.5) +
  cowplot::draw_plot(covariate.maps$dcanyons, x = 0.5, y = 0.4, width = 0.5, height = 0.5) +
  cowplot::draw_plot(covariate.maps$dincising, x = 0, y = 0.2, width = 0.5, height = 0.5) +
  cowplot::draw_plot(covariate.maps$sst, x = 0.5, y = 0.2, width = 0.5, height = 0.5) +
  cowplot::draw_plot(covariate.maps$chla, x = 0, y = -0.2, width = 0.5, height = 0.5)



fig.s2 <- (covariate.maps$depth | covariate.maps$slope) /
  (covariate.maps$dcoast | covariate.maps$dcanyons) /
    (covariate.maps$sst | covariate.maps$chla)

ggsave(plot = fig.s2, filename = file.path("fig", "FigureS2.png"), width = 30, height = 30, units = "cm")

# Correlation -------------------------------------------------------------

# Variance inflation factor functions
# https://highstat.com/Books/Book2/HighstatLibV10.R

corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

# Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}



bw.pb %>% dplyr::select(depth, slope, dcoast, dcanyons, dcanyonsInc, sst, chla) %>% 
  cor(., use = 'complete.obs', method = 'spearman')

bw.pb %>% dplyr::select(depth, slope, dcoast, dcanyons, dcanyonsInc, sst, chla) %>% 
  GGally::ggcorr(data = ., method = c('complete.obs', 'pearson'), label = TRUE)

bw.pb %>% dplyr::select(depth, slope, dcoast, dcanyons, dcanyonsInc, sst, chla.0, chla.4, chla.8, chla.12, chla.15) %>% 
  GGally::ggcorr(., geom = "blank", label = TRUE, hjust = 0.75, method = c('complete.obs', 'pearson')) +
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient) > 0.6)) +
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
  guides(color = FALSE, alpha = FALSE)

bw.pb %>% dplyr::select(depth, slope, dcoast, dcanyons, dcanyonsInc, sst, chla) %>% 
  GGally::ggcorr(data = ., method = c('complete.obs', 'spearman'), label = TRUE)

bw.pb %>% dplyr::select(depth, slope, dcoast, dcanyons, dcanyonsInc, sst, chla) %>% 
  GGally::ggcorr(., geom = "blank", label = TRUE, hjust = 0.75, method = c('complete.obs', 'spearman')) +
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient) > 0.6)) +
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
  guides(color = FALSE, alpha = FALSE)

bw.pb %>% dplyr::select(depth, slope, dcanyons, sst, chla) %>% 
  corvif(.)
