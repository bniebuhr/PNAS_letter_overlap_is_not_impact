# Re-analysis of wind power development overlap/impacts with important conservation areas
# Using maps of priority areas for expansion of RE and PA from Dunnet et al 2022 
#
# Here we propose a measure of impact overlap ratio (impact OR), to be compared 
# to a simple area overlap ratio (OR) proposed by Dunnet et al

# packages
library(dplyr)
library(purrr)
library(terra)
library(sf)

library(tmap)
library(grid)
library(ggplot2)

library(oneimpact)

#-----
# Functions
binarize <- function(x, 
                     type = c("threshold", "higher", "lower", "quantile")[1],
                     # breaks = 0.7, 
                     breaks = list(c(0.7), c(0.99)),
                     na.rm = TRUE) {
  
  # v <- terra::values()
  # if(is.list(breaks)) {
    # check if the number of elements is the same number of layers
    qg <- purrr::map2(as.list(x), breaks, function(a, b) 
      terra::global(a, quantile, prob = b, na.rm = na.rm))
  # } else {
  #   qg <- terra::global(as.list(x), quantile, prob = breaks, na.rm = na.rm)
  # }
  
  x_select <- sapply(seq_along(qg), 
                     function(i) {
                       rr <- sapply(qg[[i]][1,], function(y) x[[i]] > y)
                       rast(rr)
                     })
                       
  # terra::freq(x_select) %>% 
  #   tibble::as_tibble() %>% 
  #   dplyr::mutate(prop = count/sum(count))
  x_select <- rast(x_select)
  names(x_select) <- names(x)
  
  x_select
}

# codes
# 0: no priority
# 1: wind/infrastructure expansion priority
# 2: pa expansion priority
# 3: overlap 1 and 2
# 5: influence area of predicted infrastructure with no overlap
# 7: influence area of predicted infrastructure with overlap with pa expansion areas
overlap_pa_infra <- function(x, pa_layer = 1, infra_layer = 2, infra_zoi = NULL) {
  over_pa_infra <- x[[infra_layer]] + (2 * x[[pa_layer]])
  if(!is.null(infra_zoi)) {
    r <- 5 * (x[[infra_zoi]] - x[[infra_layer]])
    over_pa_infra <- over_pa_infra + r
  }
  over_pa_infra
}


overlap_ratio <- function(x) {
  df <- terra::freq(x) %>%
    tibble::as_tibble() %>% 
    dplyr::mutate(count_pa = ifelse(value %in% c(2,3,7), count, 0),
                  freq_pa = count_pa/sum(count_pa),
                  count_infra = ifelse(value %in% c(1,3), count, 0),
                  perc_land_infra = sum(count_infra)/sum(count),
                  overlap_ratio = freq_pa/perc_land_infra)
  list(overlap_ratio_orig = df$overlap_ratio[4], 
       overlap_ratio_zoi = df$overlap_ratio[6],
       overlap_ratio_orig_zoi = df$overlap_ratio[4] + df$overlap_ratio[6],
       df = df)
}

get_overlap_ratio <- function(x) x[1:3]

zoi_around <- function(x, zoi, layer = 2) {
  
  r <- x[[layer]]
  
  m <- focalMat(r, zoi, "circle")
  m <- m/max(m)
  
  ff <- terra::focal(r, m, fun = "max", na.rm = T)
  names(ff) <- paste0(names(r), "_buff", zoi)
  
  out <- x
  add(out) <- ff
  out
}

#----- 
# example theoretical - zoi buffer
r <- rast(ncols=36, nrows=18, crs = "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
")
v <- rep(0, ncell(r))
v[c(250,500)] <- 1
values(r) <- v
plot(r)

z <- zoi_around(r, zoi = 50, layer = 1)
plot(z)

#-----
# example for 1 region and 1 threshold value of % PA and RE priority area

# read example rasters
m_n_eu <- terra::rast("data/global_wind_solar_conservation_priorities_2022/1km/n.europe_mask_1km.tif")
w_n_eu <- terra::rast("data/global_wind_solar_conservation_priorities_2022/1km/n.europe_wind_1km.tif")
# plot(w_n_eu)
w_glob <- terra::rast("data/global_wind_solar_conservation_priorities_2022/1km/global_wind_1km.tif")
# plot(w_glob)

# pa
pa <- terra::rast("data/global_wind_solar_conservation_priorities_2022/1km/pa_priority_exp_1km.tif")
# plot(pa)
pa_n_eu <- oneimpact::rescale_vals(terra::crop(pa, m_n_eu) * m_n_eu)
# plot(pa_n_eu)

# input <- c(pa_n_eu, w_n_eu)
# input <- c(pa_sas, w_sas)
# 
# plot(input)
# max_30 <- binarize(input, breaks = list(c(0.7), c(0.7)))
# plot(max_30)
# b30 <- zoi_around(max_30, zoi = 1000)
# # over <- overlap_pa_infra(max_30)
# over <- overlap_pa_infra(b30, infra_zoi = 3)
# plot(over, col = c("grey", "yellow", "lightblue", "red"))
# overlap_ratio_all <- overlap_ratio(over)
# over_ratio <- get_overlap_ratio(overlap_ratio_all)

# in a single go
input <- c(pa_n_eu, w_n_eu)

# no buffer/zoi
over_ratio <- input %>% 
  binarize(breaks = list(0.7, 0.7)) %>% 
  overlap_pa_infra(infra_zoi = NULL) %>% 
  overlap_ratio() %>% 
  get_overlap_ratio()
over_ratio

# with buffer
over_ratio <- input %>% 
  binarize(breaks = 0.7) %>% 
  zoi_around(zoi = 1000) %>% 
  overlap_pa_infra(infra_zoi = 3) %>% 
  overlap_ratio() %>% 
  get_overlap_ratio()
over_ratio

#-----------
# create table and figure equal the one in the paper
# regions
regions <- c("c_europe", "far_east", "middle_east", 
             "n_america", "n_europe", "rus_balkans", 
             "s_america", "s_asia", "s_europe")

# masks
masks <- list.files("data/global_wind_solar_conservation_priorities_2022/1km/", 
                    pattern = "mask", full.names = T) %>% 
  sapply(terra::rast)
names(masks) <- regions

# wind
wind_index <- list.files("data/global_wind_solar_conservation_priorities_2022/1km/",
                         pattern = "wind") %>% 
  grep(pattern = "global|priority|bivmap|windmat", invert = TRUE)

wind <- list.files("data/global_wind_solar_conservation_priorities_2022/1km/", 
                   pattern = "wind", full.names = T)[wind_index] %>% 
  sapply(terra::rast)
names(wind) <- regions

# pa
pa_global <- terra::rast("data/global_wind_solar_conservation_priorities_2022/1km/pa_priority_exp_1km.tif")
pa <- sapply(masks, function(x) terra::crop(pa_global, x) * x)
# plot(pa[[2]])

#-----
# run original analysis

# parameters
parms <- expand.grid(x = regions,
                     y = c(0.99, 0.9, 0.7),
                     z = seq(0.99, 0.7, -0.01))

over_ratio <- parms %>% 
  purrr::pmap_dbl(
  function(x, y, z, ...) {
    print(paste(x, y, z))
    input <- c(pa[[x]], wind[[x]])
    or <- input %>% 
      binarize(breaks = list(y, z)) %>% 
      overlap_pa_infra() %>% 
      overlap_ratio() %>% 
      get_overlap_ratio()
    or
  })

# data.frame
regions_labels <- c("Central Europe", "East Asia", "Middle East",
                    "North America", "Northern Europe", "Russia and the Balkans",
                    "South America", "South Asia", "Southern Europe")
cols <- c("#5F4690", "#1D6996", "#7DC4C4",
          "#4AA37E", "#73AF48", "#EEB41F",
          "#E79535", "#CE5948", "#6D6D6D")

df_overlap <- parms %>%
  tibble::as_tibble() %>% 
  dplyr::rename(region = 1, pa_threshold = 2, wind_threshold = 3) %>% 
  dplyr::mutate(region_lab = factor(region, levels = regions, labels = regions_labels),
                pa_th = 100*(1-pa_threshold), 
                pa_th = factor(pa_th, labels = c("1% PA", "10% PA", "30% PA")),
                wind_th = 1-wind_threshold,
                or_original = over_ratio)

# save(df_overlap, file = "output/df_overlap_original.rda")
# write.csv(df_overlap, file = "output/df_overlap_original.csv", row.names = F)

# plot
df_overlap %>% 
  ggplot(aes(wind_th, or_original, color = region_lab)) +
  geom_line(size = 1.5) +
  facet_wrap(~pa_th) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Percentage of land devoted to renewables",
       y = "Overlap ratio", 
       color = "") +
  scale_color_manual(values = cols) +
  theme_bw()

#-----
# run analysis zoi 1 km
  
# parameters
zoi <- 1000 # in meters
parms <- expand.grid(x = regions,
                     y = c(0.99, 0.9, 0.7),
                     z = seq(0.99, 0.7, -0.01))

over_ratio_1km <- parms %>% 
  purrr::pmap(
    function(x, y, z, ...) {
      print(paste(x, y, z))
      input <- c(pa[[x]], wind[[x]])
      or <- input %>% 
        binarize(breaks = list(y, z)) %>% 
        zoi_around(zoi = zoi) %>% 
        overlap_pa_infra(infra_zoi = 3) %>% 
        overlap_ratio() %>% 
        get_overlap_ratio()
      or
    })

over1km_orig <- purrr::map_dbl(over_ratio_1km, ~ .[[1]])
over1km_zoi <- purrr::map_dbl(over_ratio_1km, ~ .[[2]])
over1km_orig_zoi <- purrr::map_dbl(over_ratio_1km, ~ .[[3]])
  
# read original data.frame
load("output/df_overlap_original.rda")
df_overlap

# check
all(df_overlap$or_original == over1km_orig) # TRUE, ok!

# data.frame
df_overlap_1km <- df_overlap %>%
  dplyr::mutate(or_zoi_1km = over1km_zoi,
                or_original_zoi_1km = over1km_orig_zoi)

save(df_overlap_1km, file = "output/df_overlap_original_zoi1km.rda")
write.csv(df_overlap_1km, file = "output/df_overlap_original_zoi1km.csv", row.names = F)

# plot
df_overlap_1km %>% 
  ggplot(aes(wind_th, or_original_zoi_1km, color = region_lab)) +
  geom_line(size = 1.2) +
  facet_wrap(~pa_th) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Percentage of land devoted to renewables",
       y = "Overlap ratio with 1 km ZoI", 
       color = "") +
  scale_color_manual(values = cols) +
  theme_bw()

df_overlap_1km %>% 
  ggplot(aes(wind_th, or_original_zoi_1km/or_original, color = region_lab)) +
  # ggplot(aes(wind_th, or_zoi_1km, color = region_lab)) +
  geom_line(size = 1.2) +
  facet_wrap(~pa_th) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Percentage of land devoted to renewables",
       y = "Overlap ratio increase from 0 to 1km ZoI", 
       color = "") +
  scale_color_manual(values = cols) +
  theme_bw()
  
# plot both original and buffer
df_overlap_1km %>% 
  tidyr::pivot_longer(cols = c(or_original, or_original_zoi_1km), 
                      names_to = "or_type", values_to = "or") %>% 
  dplyr::mutate(or_type = factor(or_type, levels = c("or_original", "or_original_zoi_1km"),
                                 labels = c("ZoI = 0", "ZoI = 1 km"))) %>% 
  ggplot(aes(wind_th, or, color = region_lab)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, color = "grey", linetype = 2) + 
  facet_grid(or_type~pa_th) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Percentage of land devoted to renewables",
       y = "Overlap ratio", 
       color = "") +
  scale_color_manual(values = cols) +
  theme_bw()

#-----
# run analysis zoi 5 km

# parameters
zoi <- 5000 # in meters
parms <- expand.grid(x = regions,
                     y = c(0.99, 0.9, 0.7),
                     z = seq(0.99, 0.7, -0.01))

over_ratio_5km <- parms %>% 
  purrr::pmap(
    function(x, y, z, ...) {
      print(paste(x, y, z))
      input <- c(pa[[x]], wind[[x]])
      or <- input %>% 
        binarize(breaks = list(y, z)) %>% 
        zoi_around(zoi = zoi) %>% 
        overlap_pa_infra(infra_zoi = 3) %>% 
        overlap_ratio() %>% 
        get_overlap_ratio()
      or
    })

over5km_orig <- purrr::map_dbl(over_ratio_5km, ~ .[[1]])
over5km_zoi <- purrr::map_dbl(over_ratio_5km, ~ .[[2]])
over5km_orig_zoi <- purrr::map_dbl(over_ratio_5km, ~ .[[3]])

# read original data.frame
load("output/df_overlap_original_zoi1km.rda")
df_overlap_1km

# check
all(df_overlap_1km$or_original == over5km_orig) # TRUE, ok!

# data.frame
df_overlap_5km <- df_overlap_1km %>%
  dplyr::mutate(or_zoi_5km = over5km_zoi,
                or_original_zoi_5km = over5km_orig_zoi)

# save(df_overlap_5km, file = "output/df_overlap_original_zoi5km.rda")
# write.csv(df_overlap_5km, file = "output/df_overlap_original_zoi5km.csv", row.names = F)
load("output/df_overlap_original_zoi5km.rda")

# plot
df_overlap_5km %>% 
  ggplot(aes(wind_th, or_original_zoi_5km, color = region_lab)) +
  geom_line(size = 1.2) +
  facet_wrap(~pa_th) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Percentage of land devoted to renewables",
       y = "Overlap ratio with 5 km ZoI", 
       color = "") +
  scale_color_manual(values = cols) +
  theme_bw()

df_overlap_5km %>% 
  ggplot(aes(wind_th, or_original_zoi_5km/or_original, color = region_lab)) +
  # ggplot(aes(wind_th, or_zoi_1km, color = region_lab)) +
  geom_line(size = 1.2) +
  facet_wrap(~pa_th) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Percentage of land devoted to renewables",
       y = "Overlap ratio increase from 0 to 5km ZoI", 
       color = "") +
  scale_color_manual(values = cols) +
  theme_bw()

# plot both original and buffer
df_overlap_5km %>% 
  tidyr::pivot_longer(cols = c(or_original, or_original_zoi_1km, or_original_zoi_5km), 
                      names_to = "or_type", values_to = "or") %>% 
  dplyr::mutate(or_type = factor(or_type, levels = c("or_original", "or_original_zoi_1km",
                                                     "or_original_zoi_5km"),
                                 labels = c("ZoI = 0", "ZoI = 1 km", "ZoI = 5 km"))) %>% 
  ggplot(aes(wind_th, or, color = region_lab)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, color = "grey", linetype = 2) + 
  facet_grid(or_type~pa_th, scales = "free_y") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Percentage of land devoted to renewables",
       y = "Overlap ratio", 
       color = "") +
  scale_color_manual(values = cols) +
  theme_bw()

cols <- c("slateblue", "skyblue4", "deepskyblue3",
          "springgreen4", "#33CC00", "goldenrod1",
          "darkorange1", "brown2", "darkgrey")
cols <- c("#5F4690", "#1D6996", "#7DC4C4",
          "#4AA37E", "#73AF48", "#EEB41F",
          "#E79535", "#CE5948", "#6D6D6D")
wind_zoi_1_5 <- df_overlap_5km %>% 
  tidyr::pivot_longer(cols = c(or_original_zoi_1km, or_original_zoi_5km), 
                      names_to = "or_type", values_to = "or") %>% 
  dplyr::mutate(or_type = factor(or_type, levels = c("or_original_zoi_1km",
                                                     "or_original_zoi_5km"),
                                 labels = c("ZoI = 1 km", "ZoI = 5 km"))) %>%
  ggplot(aes(wind_th, or, color = region_lab)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 1, color = "darkred", linetype = 2) +
  facet_grid(or_type~pa_th, scales = "free_y") +
  scale_x_continuous(breaks = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3), 
                     labels = scales::percent_format(accuracy = 1L)) +
  scale_color_manual(values = cols) +
  labs(x = "Percentage of land devoted to renewables",
       y = "Impact overlap ratio", 
       color = "") +
  theme_light() +
  theme(legend.position = "top") +
  guides(color = guide_legend(ncol = 3, keyheight = 0.1))
wind_zoi_1_5

ggsave("overlap_ratio_zoi_1_5.png", plot = wind_zoi_1_5, path = "output",
       width = 20, height = 15, units = "cm", dpi = 300)
# ggsave("overlap_ratio_zoi_1_5.tiff", plot = wind_zoi_1_5, path = "output",
#        width = 20, height = 15, units = "cm", dpi = 300)
ggsave("overlap_ratio_zoi_1_5.pdf", plot = wind_zoi_1_5, path = "output",
       width = 20, height = 15, units = "cm", dpi = 300)

##########################
# Global map with ZoI

# masks
masks <- list.files("data/global_wind_solar_conservation_priorities_2022/1km/", 
                    pattern = "mask", full.names = T) %>% 
  sapply(terra::rast)
names(masks) <- regions

# wind
wind_index <- list.files("data/global_wind_solar_conservation_priorities_2022/1km/",
                         pattern = "wind") %>% 
  grep(pattern = "global|priority|bivmap|windmat", invert = TRUE)

wind <- list.files("data/global_wind_solar_conservation_priorities_2022/1km/", 
                   pattern = "wind", full.names = T)[wind_index] %>% 
  sapply(terra::rast)
names(wind) <- regions

# pa
pa_global <- terra::rast("data/global_wind_solar_conservation_priorities_2022/1km/pa_priority_exp_1km.tif")
pa <- sapply(masks, function(x) terra::crop(pa_global, x) * x)
# plot(pa[[2]])

# buffer 1km
zoi <- 1000 # in meters
map_1km <- purrr::map(seq_along(wind),
                      function(i) {
                        print(paste(i))
                        input <- c(pa[[i]], wind[[i]])
                        or <- input %>% 
                          binarize(breaks = list(0.7, 0.7)) %>% 
                          zoi_around(zoi = zoi) %>% 
                          overlap_pa_infra(infra_zoi = 3) 
                        or
                      })

map_1km_merge <- terra::merge(map_1km[[1]], map_1km[[2]], map_1km[[3]], 
                              map_1km[[4]], map_1km[[5]], map_1km[[6]], 
                              map_1km[[7]], map_1km[[8]], map_1km[[9]])
plot(map_1km_merge)
map_1km_merge46 <- terra::classify(map_1km_merge, rbind(c(5, 4), c(7, 6)))
freq(map_1km_merge46) %>% tibble::as_tibble() %>% mutate(prop = count/sum(count))

# buffer
zoi <- 5000 # in meters
map_5km <- purrr::map(seq_along(wind),
    function(i) {
      print(paste(i))
      input <- c(pa[[i]], wind[[i]])
      or <- input %>% 
        binarize(breaks = list(0.7, 0.7)) %>% 
        zoi_around(zoi = zoi) %>% 
        overlap_pa_infra(infra_zoi = 3) 
      or
    })

map_5km_merge <- terra::merge(map_5km[[1]], map_5km[[2]], map_5km[[3]], 
                              map_5km[[4]], map_5km[[5]], map_5km[[6]], 
                              map_5km[[7]], map_5km[[8]], map_5km[[9]])
names(map_5km_merge) <- "wind_pa_overlap"
plot(map_5km_merge)
freq(map_5km_merge) %>% tibble::as_tibble() %>% mutate(prop = count/sum(count))
# map_5km_moll = terra::project(map_5km_merge, "+proj=moll")
# plot(map_5km_moll)

# export
terra::writeRaster(map_1km_merge, filename = "output/map_1km_4overout_6overinpa.tif")
terra::writeRaster(map_5km_merge, filename = "output/map_5km_5overout_7overinpa.tif", overwrite = TRUE)

# read
map_5km_merge <- terra::rast("output/map_5km_5overout_7overinpa.tif")*1
zoom_areas <- sf::st_read("data/zoom_areas.shp") %>% 
  sf::st_transform(sf::st_crs(pa_global)) %>% 
  dplyr::arrange(area)

# plot
cols5 <- c("grey", "yellow", "dodgerblue1", "violetred4", NA, "seagreen2", NA, "orange")
plot(map_5km_merge, col = cols5)

cols1 <- c("grey", "yellow", "dodgerblue1", "violetred4", "seagreen2", NA, "orange")
plot(map_1km_merge, col = cols1)

# plot with tmap
map_5km_moll <- terra::project(map_5km_merge, "+proj=moll", method = "near")
plot(map_5km_moll, col = cols5)
map_5km_eck4 <- terra::project(map_5km_merge, terra::crs(pa_global), method = "near")
plot(map_5km_eck4, col = cols5)

# zoom maps
m1 <- terra::crop(map_5km_eck4, zoom_areas[1,])
plot(m1, col = cols5)
m2 <- terra::crop(map_5km_eck4, zoom_areas[2,])
plot(m2, col = cols5)
m3 <- terra::crop(map_5km_eck4, zoom_areas[3,])
plot(m3, col = cols5)
m4 <- terra::crop(map_5km_eck4, zoom_areas[4,])
plot(m4, col = cols5)
m5 <- terra::crop(map_5km_eck4, zoom_areas[5,])
plot(m5, col = cols5)

library(spData)
library(rnaturalearthdata)
world_moll <- st_transform(world, crs = "+proj=moll")
world_eck4 <- st_transform(world, crs = terra::crs(pa_global))

# colors and legend
brk <- c(0, 0.5, 1.5, 2.5, 3.5, 5.5, 7.1)
labs <- c("No priority", "RE priority", "PA priority", 
         "Overlap", "ZoI 5km", 
         "Overlap - ZoI 5km")

# stretch in x direction
bb <- st_bbox(world_eck4)
xrang <- diff(bb[c(1,3)])*0.07
# diff(bb[c(1,2)])*0.25
add <- c(-xrang, 0, +xrang, 0)

map5km_exp <- tm_shape(world_eck4, bbox = st_bbox(world_eck4) + add)  +
  tm_fill(col = grey(0.9)) +
  tm_shape(map_5km_eck4) +
  tm_raster(title = "", style = "fixed", breaks = brk, 
            palette = cols5[!is.na(cols5)],
            labels = labs) +
  tm_shape(zoom_areas) +
  tm_borders(col = "black") +
  tm_scale_bar(breaks = c(0, 2000, 4000), 
               position = c("center", "BOTTOM")) +
  tm_compass(position = c(0.6, 0.05)) +
  tm_layout(frame = FALSE, legend.bg.color = "white")
map5km_exp

tmap_save(map5km_exp, filename = "output/map_5km_tmap.png", 
          width = 20, height = 12, units = "cm", dpi = 300)

world_detailed <- rnaturalearthdata::countries50 %>% 
  sf::st_as_sf() %>% 
  sf::st_transform(crs = st_crs(pa_global))

map5km_zoom1 <- tm_shape(world_eck4, bbox = zoom_areas[1,]) +
  tm_polygons(col = grey(0.9)) +
  # tm_shape(map_5km_eck4) +
  tm_shape(m1) +
  tm_raster(title = "", style = "fixed", breaks = brk, 
            palette = cols5[!is.na(cols5)], 
            legend.show = FALSE)
map5km_zoom1

map5km_zoom2 <- tm_shape(world_detailed, bbox = zoom_areas[2,]) +
  tm_polygons(col = grey(0.9)) +
  # tm_shape(map_5km_eck4) +
  tm_shape(m2, ) +
  tm_raster(title = "", style = "fixed", breaks = brk, 
            palette = cols5[!is.na(cols5)], 
            legend.show = FALSE)
map5km_zoom2

map5km_zoom3 <- tm_shape(world_detailed, bbox = zoom_areas[3,]) +
  tm_polygons(col = grey(0.9)) +
  # tm_shape(map_5km_eck4) +
  tm_shape(m3) +
  tm_raster(title = "", style = "fixed", breaks = brk, 
            palette = cols5[!is.na(cols5)], 
            legend.show = FALSE)
map5km_zoom3

map5km_zoom4 <- tm_shape(world_detailed, bbox = zoom_areas[4,]) +
  tm_polygons(col = grey(0.9)) +
  # tm_shape(map_5km_eck4) +
  tm_shape(m4) +
  tm_raster(title = "", style = "fixed", breaks = brk, 
            palette = cols5[!is.na(cols5)], 
            legend.show = FALSE)
map5km_zoom4

map5km_zoom5 <- tm_shape(world_detailed, bbox = zoom_areas[5,]) +
  tm_polygons(col = grey(0.9)) +
  # tm_shape(map_5km_eck4) +
  tm_shape(m5) +
  tm_raster(title = "", style = "fixed", breaks = brk, 
            palette = cols5[!is.na(cols5)], 
            legend.show = FALSE)
map5km_zoom5


xy <- st_bbox(zoom_areas[1,])
asp2 <- (xy$xmax - xy$xmin)/(xy$ymax - xy$ymin)
w <- 0.35
h <- w/asp2
vp1 <- viewport(x=-0.05, y=0.85, width = w, height=h, just=c("left", "top"))

xy <- st_bbox(zoom_areas[2,])
asp2 <- (xy$xmax - xy$xmin)/(xy$ymax - xy$ymin)
w <- 0.35
h <- w/asp2
vp2 <- viewport(x=-0.05, y=0.6, width = w, height=h, just=c("left", "top"))

xy <- st_bbox(zoom_areas[3,])
asp2 <- (xy$xmax - xy$xmin)/(xy$ymax - xy$ymin)
w <- 0.35
h <- w/asp2
vp3 <- viewport(x=1.05, y=0.85, width = w, height=h, just=c("right", "top"))

xy <- st_bbox(zoom_areas[4,])
asp2 <- (xy$xmax - xy$xmin)/(xy$ymax - xy$ymin)
w <- 0.35
h <- w/asp2
vp4 <- viewport(x=1.05, y=0.6, width = w, height=h, just=c("right", "top"))

xy <- st_bbox(zoom_areas[5,])
asp2 <- (xy$xmax - xy$xmin)/(xy$ymax - xy$ymin)
w <- 0.35
h <- w/asp2
vp5 <- viewport(x=1.05, y=0.3, width = w, height=h, just=c("right", "top"))


tmap_save(map5km_exp, filename = "output/map_5km_tmap.png", 
          width = 25, height = 12, units = "cm", dpi = 300,
          insets_tm = list(map5km_zoom1, map5km_zoom2, map5km_zoom3, map5km_zoom4, map5km_zoom5), 
          insets_vp = list(vp1, vp2, vp3, vp4, vp5))
# tmap_save(map5km_exp, filename = "output/map_5km_tmap.tiff", 
#           width = 25, height = 12, units = "cm", dpi = 300,
#           insets_tm = list(map5km_zoom1, map5km_zoom2, map5km_zoom3, map5km_zoom4, map5km_zoom5), 
#           insets_vp = list(vp1, vp2, vp3, vp4, vp5))
tmap_save(map5km_exp, filename = "output/map_5km_tmap.pdf", 
          width = 25, height = 12, units = "cm", dpi = 300,
          insets_tm = list(map5km_zoom1, map5km_zoom2, map5km_zoom3, map5km_zoom4, map5km_zoom5), 
          insets_vp = list(vp1, vp2, vp3, vp4, vp5))

