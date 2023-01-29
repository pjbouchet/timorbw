##%######################################################################################%##
#                                                                                       #
####           Evidence of likely foraging by pygmy blue whales                         ####
####     in the Timor Trough during the late austral winter and early austral spring    ####
#                                                                                      #
##%######################################################################################%##


# Register Google API Key
# ggmap::register_google(key = "xxxxxx")

# Summary of results ---------------------------------------------------------------

# Number of sightings
sightings.per.month <- bw %>% 
  dplyr::group_by(year, month) %>% 
  tally(adcalves>0) %>% dplyr::ungroup()

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

gps_lines <- rbind(gps$`2007`, gps$`2008`) %>% 
  split(x = ., f = .$date) %>% 
  purrr::map(., ~.x %>% split(x = ., f = .$line_id) %>% purrr::map(., ~ createLines(.))) %>% 
  purrr::map_depth(.x = ., .depth = 2, .f = ~sp::spTransform(.x, CRSobj = CRSutm)) %>% 
  purrr::map(.x = ., .f = ~do.call(rbind, .)) %>% 
  purrr::map(.x = ., .f = ~smoothr::drop_crumbs(x =.x, threshold = 250))

depth.along.lines <- 
 purrr::map(.x = gps_lines, 
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

daily.effort <- purrr::map_df(.x = gps_lines, .f = ~rgeos::gLength(.x)) %>% 
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

# Figure 3 ----------------------------------------------------------------

pdf("Figure3.pdf", width = 8, height = 5)
par(mfrow = c(1,2))
purrr::walk(.x = 1:4, .f = ~{
  plot(bestmodel, shade = TRUE, scale = 0, select = .x, xlab = c("Max chlorophyll-a concentration (mg/m3)", "Distance to incising canyon (km)", "Seabed slope (degrees)")[.x]); abline(h = 0, lty = 2)})
dev.off()

# Figure S1 ---------------------------------------------------------------

make_map(input.raster = depth, show.canyons = TRUE, plot.title = "Depth (m)", legend.title = FALSE)
ggsave("canyon_map.pdf")

# Figure S2 ---------------------------------------------------------------

check_points(PA, 20)

# Figure S3 ---------------------------------------------------------------

quickplot(depth, save = TRUE, sightings = TRUE)
quickplot(seabed_slope, breaks = c(0:15, 20,25,30,40), save = TRUE, sightings = TRUE)
quickplot(fcpi$fcpi, save = TRUE, sightings = TRUE)
quickplot(dist_canyons, save = TRUE, sightings = TRUE)
quickplot(dist_incising, save = TRUE, sightings = TRUE)
quickplot(raster::calc(raster::stack(chlamax_climg_filled[27:39]),median), save = TRUE, sightings = TRUE)
quickplot(raster::calc(raster::stack(chla_climg_filled[27:39]),median), save = TRUE, sightings = TRUE)

# Figure S4 ---------------------------------------------------------------

# (A) Time series

png("fig/Figure-S1a.png", res = 450, height = 1800, width = 3500)
par(las = 1)
plot.ts(ts(sst.buoy$data$sst, frequency = 365, start = c(2002,1)),
        xlab = NA, ylab = "SST (°C)", axes = FALSE)
time.labels <- seq(as.Date('2002-06-01'), as.Date('2020-12-31'), by = "12 months") 
axis(2)
axis(1, labels = lubridate::year(time.labels), 
     at = seq(from = 2002, by = 1, length.out = length(time.labels)))
box()
box()
dev.off() 

# (B) Wavelet power spectrum

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

# Figure S5 ---------------------------------------------------------------

ggcorrplot::ggcorrplot(corr, hc.order = FALSE, type = "lower", lab = TRUE)

# Figure S6 ---------------------------------------------------------------

mod_cv <- mod_cv %>% dplyr::mutate(metric = ifelse(metric == "AUC", "ROC-AUC", metric))
ggplot(mod_cv, aes(metric, value)) + geom_boxplot(fill = "grey") +
  ylab("") + xlab("") +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        axis.ticks.x = element_blank())

# Figure S7 ---------------------------------------------------------------

plot_variogram(model = bestmodel, dat = na.omit(bw.df), maxD = 25000, pch = 16, col = "black", ylim = c(-1, 100))

ggplot(data.frame(x = Vario.bw$u, y = Vario.bw$v), aes(x, y)) +
  geom_point() +
  geom_smooth(span = 1) +
  scale_x_continuous(breaks = seq(0, Vario.D, 5000), labels = format(seq(0, Vario.D, 5000), big.mark = ",")) +
  scale_y_continuous(breaks = seq(0, 3500, 500)) +
  ylab("Semivariance") + xlab("Distance (m)") + 
  ggplot2::theme(text = element_text(size = 12),
                 axis.text.y = element_text(angle = 90, hjust = 0.5),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
                 axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
                 axis.text = element_text(color = 'black', size = 11),
                 plot.margin = margin(1, 1, 1, 1, "cm"))
