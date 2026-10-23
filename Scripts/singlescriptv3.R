# =============================================================================
# REPRODUCE PAPER: PM2.5-ATTRIBUTABLE MORTALITY IN YEREVAN, 2023-2026
# =============================================================================


## ---- 0. SETUP --------------------------------------------------------------
library(readr); library(dplyr); library(sf); library(gstat); library(spacetime)
library(raster); library(sp); library(parallel); library(pbapply)
library(terra); library(ggplot2); library(lubridate); library(zoo)
library(viridis); library(tidyr); library(patchwork); library(stringr)
library(scales); library(ggspatial); library(readxl); library(grid)
library(leaflet); library(leafem); library(htmlwidgets)
library(maptiles); library(tidyterra); library(cowplot); library(maps)

# >>> EDIT THESE PATHS ONCE <<<
PROJ_DIR  <- "C:/Users/mprhc/OneDrive/Documents/GitHub/ArmeniaPollutionAnalysis"
DATA_DIR  <- "C:/Users/mprhc/Downloads"
OUT_DIR   <- file.path(PROJ_DIR, "Interpolation outputs_ST")
PRED_DIR  <- file.path(OUT_DIR, "pred_layers")
VAR_DIR   <- file.path(OUT_DIR, "var_layers")
# ---- Versioned run folders ----
# Every distinct GEMM coefficient file gets its own run folder (keyed by file hash),
# so changing the coefficients never collides with a previous run's cached outputs.
# Kriging products (PRED_DIR/VAR_DIR under OUT_DIR) are input-independent and stay shared.
GEMM_XLSX   <- file.path(PROJ_DIR, "GEMMcoefficients_NCDLRI_GHE2023.xlsx")
run_id <- if (file.exists(GEMM_XLSX)) substr(unname(tools::md5sum(GEMM_XLSX)), 1, 8) else "MISSING_COEFFS"
RUN_DIR <- file.path(PROJ_DIR, "Runs", paste0("run_", run_id))
CACHE_DIR <- file.path(RUN_DIR, "Processed_Cache")
QUAD_DIR  <- file.path(RUN_DIR, "Quadrant_Maps")
YEAR_PLOTS_DIR <- file.path(RUN_DIR, "Yearly_Plots")
dir.create(PRED_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(VAR_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QUAD_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(YEAR_PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
message("Run folder: ", RUN_DIR)

READINGS_CSV   <- file.path(DATA_DIR, "readings_clean_2023_2026.csv")
SENSORS_CSV    <- file.path(DATA_DIR, "sensors_filtered_2023_2026.csv")
BOUNDARY_LVL4  <- file.path(DATA_DIR, "boundary.gpkg")
BOUNDARY_LVL5  <- file.path(PROJ_DIR, "Yerevan Boundary/boundary-polygon-lvl5.gpkg")
POP_TIF        <- file.path(PROJ_DIR, "WorldPopData/arm_age_data/pop.tif")

PERIOD_LABELS <- c("June 2023 - June 2024",
                   "June 2024 - June 2025",
                   "June 2025 - June 2026")

# Common helpers --------------------------------------------------------------
translate_districts <- function(x) {
  str_replace_all(x, c(
    "Կենտրոն"="Kentron","Արաբկիր"="Arabkir","Քանաքեռ"="Kanaker","Կանաքեռ"="Kanaker",
    "Քանաքեր"="Kanaker","Կանաքեր"="Kanaker","Զեյթուն"="Zeytun","Նորք"="Nork",
    "Նոր"="Nor","Մարաշ"="Marash","Շենգավիթ"="Shengavit","Աջափնյակ"="Ajapnyak",
    "Դավթաշեն"="Davtashen","Էրեբունի"="Erebuni","Մալաթիա"="Malatia",
    "Սեբաստիա"="Sebastia","Ավան"="Avan","Նուբարաշեն"="Nubarashen"))
}

load_districts <- function(path = BOUNDARY_LVL5) {
  d <- st_read(path, layer = "boundary-polygon-lvl5", quiet = TRUE) |>
    st_transform(4326)
  nm <- grep("name|district|nom", colnames(d), ignore.case = TRUE, value = TRUE)[1]
  d[[nm]] <- translate_districts(d[[nm]])
  d[[nm]] <- str_replace(d[[nm]], "Kanaker[- ]Zeytun", "Kanaker-Zeytun")
  d[[nm]] <- str_replace(d[[nm]], "Nork[- ]Marash",    "Nork-\nMarash")
  d[[nm]] <- str_replace(d[[nm]], "Nor[- ]?Nork",      "Nor-\nNork")
  list(sf = d, name_col = nm)
}

# Resolve namespace conflicts (terra over raster)
extract  <- terra::extract; mask <- terra::mask; values <- terra::values
resample <- terra::resample; project <- terra::project; crop <- terra::crop

## ============================================================================
## 1. FIGURE 1 - Temporal Coverage of Sensor Network
## ============================================================================
figure_1 <- function() {
  message("=== Figure 1: Sensor network temporal coverage ===")
  out_file <- file.path(OUT_DIR, "Figure1_sensor_temporal_coverage.png")
  if (file.exists(out_file)) {
    message("Output already exists. Skipping.")
    return(invisible(NULL))
  }
  
  sensors <- read_csv(SENSORS_CSV, show_col_types = FALSE)
  sensors$first <- as.POSIXct(sensors$first_measurement_time)
  sensors$last  <- as.POSIXct(sensors$last_measurement_time)
  
  sensors_time_clean <- sensors %>%
    mutate(
      provider_clean = case_when(
        provider == "clarity"    ~ "Clarity",
        provider == "purple_air" ~ "PurpleAir",
        TRUE                     ~ stringr::str_to_title(gsub("[_-]", " ", provider))
      )
    ) %>%
    dplyr::filter(!is.na(first) & !is.na(last))
  
  days <- seq(min(sensors_time_clean$first, na.rm = TRUE),
              max(sensors_time_clean$last, na.rm = TRUE),
              by = "day")
  
  unique_providers <- unique(sensors_time_clean$provider_clean)
  
  df_list <- lapply(unique_providers, function(p) {
    sub_sensors <- sensors_time_clean[sensors_time_clean$provider_clean == p, ]
    active_count <- sapply(days, function(d) {
      sum(sub_sensors$first <= d & sub_sensors$last >= d, na.rm = TRUE)
    })
    data.frame(date = days, provider = p, active_sensors = active_count)
  })
  df_time <- do.call(rbind, df_list)
  
  provider_order <- df_time %>%
    group_by(provider) %>%
    summarise(max_active = max(active_sensors, na.rm = TRUE)) %>%
    arrange(desc(max_active)) %>%
    pull(provider)
  
  df_time$provider <- factor(df_time$provider, levels = provider_order)
  df_time <- df_time %>% dplyr::filter(date >= as.POSIXct("2022-01-01"))
  
  p1 <- ggplot(df_time, aes(x = date, y = active_sensors, fill = provider)) +
    geom_area(position = "stack", alpha = 0.85, color = "white", linewidth = 0.1) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    scale_fill_viridis_d(option = "plasma", name = "Sensor Provider",
                         guide = guide_legend(reverse = FALSE)) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey92"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey40", size = 8),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    ) +
    labs(
      title = "Temporal Coverage of the Yerevan Air Quality Sensor Network",
      subtitle = "The total height shows overall active sensors; colored bands represent individual provider contributions from 2022 onward.",
      x = "Timeline Year",
      y = "Total Active Sensors (Cumulative)"
    )
  
  ggsave(out_file, plot = p1, width = 11, height = 6, dpi = 300)
}

## ============================================================================
## 2. FIGURE 2 - Spatiotemporal Expansion of Sensor Network (Publication Map)
## ============================================================================
figure_2 <- function() {
  message("=== Figure 2: Sensor network spatiotemporal expansion ===")
  out_file <- file.path(OUT_DIR, "Figure2_sensor_spatial_expansion.png")
  if (file.exists(out_file)) {
    message("Output already exists. Skipping.")
    return(invisible(NULL))
  }
  
  options(HTTPUserAgent = "Mozilla/5.0")
  options(timeout = 120)
  
  # Viewport bounds for a balanced, squarish quadrant view
  trim_xmin <- 44.43; trim_xmax <- 44.59
  trim_ymin <- 40.13; trim_ymax <- 40.24
  
  # --- 1. Load sensor data ---
  sensors <- read.csv(SENSORS_CSV)
  sensors_all <- sensors %>%
    dplyr::filter(provider %in% c("clarity", "purple_air", "yerevan_gis")) %>%
    dplyr::filter(longitude >= trim_xmin, longitude <= trim_xmax,
                  latitude  >= trim_ymin, latitude  <= trim_ymax) %>%
    mutate(
      install_time = as.POSIXct(first_measurement_time, format = "%Y-%m-%d %H:%M:%S"),
      install_time_num = as.numeric(install_time),
      type = case_when(
        provider == "clarity"      ~ "Residential Clarity Sensors",
        provider == "purple_air"   ~ "PurpleAir",
        provider == "yerevan_gis"  ~ "Municipality Construction Clarity Sensors"
      ),
      type = factor(type, levels = c("Residential Clarity Sensors",
                                     "Municipality Construction Clarity Sensors",
                                     "PurpleAir")),
      district_en = case_when(
        district_slug == "ajapnyak"          ~ "Ajapnyak",
        district_slug == "arabkir"           ~ "Arabkir",
        district_slug == "avan"              ~ "Avan",
        district_slug == "kentron"           ~ "Kentron",
        district_slug == "davtashen"         ~ "Davtashen",
        district_slug == "erebuni"           ~ "Erebuni",
        district_slug == "malatia-sebastia"  ~ "Malatia-Sebastia",
        district_slug == "nor-nork"          ~ "Nor Nork",
        district_slug == "nork-marash"        ~ "Nork-Marash",
        district_slug == "nubarashen"        ~ "Nubarashen",
        district_slug == "shengavit"         ~ "Shengavit",
        district_slug == "kanaker-zeytun"    ~ "Kanaker-Zeytun",
        TRUE ~ tools::toTitleCase(gsub("-", " ", district_slug))
      )
    ) %>%
    dplyr::filter(install_time >= as.POSIXct("2023-01-01"))
  
  sensors_sf <- st_as_sf(sensors_all, coords = c("longitude","latitude"), crs = 4326)
  bbox_vals <- c(xmin=trim_xmin, xmax=trim_xmax, ymin=trim_ymin, ymax=trim_ymax)
  yerevan_bbox <- st_as_sfc(st_bbox(bbox_vals, crs = 4326))
  sensors_sf <- st_intersection(sensors_sf, yerevan_bbox)
  
  # --- 2. Load district boundaries ---
  districts_sf <- st_read(BOUNDARY_LVL5, layer = "boundary-polygon-lvl5", quiet = TRUE) %>%
    st_transform(4326)
  
  # --- 3. Basemap tiles ---
  tiles <- get_tiles(yerevan_bbox, provider = "OpenTopoMap", zoom = 12, crop = TRUE)
  
  # --- 4. Main map ---
  p_main <- ggplot() +
    geom_spatraster_rgb(data = tiles, alpha = 0.88) +
    geom_sf(data = districts_sf, fill = alpha("#2c3e50", 0.04),
            color = "#2c3e50", linewidth = 0.4, linetype = "solid") +
    geom_sf(data = sensors_sf, aes(fill = install_time_num, shape = type),
            color = "#111111", stroke = 0.8, size = 2.1, alpha = 0.95) +
    scale_shape_manual(
      name = "Monitoring Network",
      values = c("Residential Clarity Sensors" = 23,
                 "Municipality Construction Clarity Sensors" = 21,
                 "PurpleAir" = 24)
    ) +
    scale_fill_gradientn(
      name = "First Measurement Year",
      colors = c("#03045E", "#6A0DAD", "#D81B60", "#F4A261", "#FFEA00"),
      limits = c(as.numeric(as.POSIXct("2023-01-01")),
                 as.numeric(as.POSIXct("2026-01-01"))),
      oob = scales::squish,
      breaks = as.numeric(as.POSIXct(c("2023-01-01","2024-01-01","2025-01-01","2026-01-01"))),
      labels = c("2023","2024","2025","2026"),
      guide = guide_colorbar(direction = "horizontal", title.position = "top",
                             barwidth = unit(6, "cm"), barheight = unit(0.4, "cm"), order = 1)
    ) +
    annotation_scale(location = "br", width_hint = 0.25, text_size = 7,
                     pad_x = unit(0.4, "in")) +
    annotation_north_arrow(location = "br", pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
                           height = unit(0.8, "cm"), width = unit(0.8, "cm"),
                           style = north_arrow_orienteering(fill = c("grey10","white"),
                                                            text_col = "grey10", text_size = 6)) +
    scale_x_continuous(name = "Longitude",
                       breaks = seq(trim_xmin + 0.02, trim_xmax - 0.02, length.out = 4),
                       labels = label_dollar(prefix = "", suffix = "°E")) +
    scale_y_continuous(name = "Latitude",
                       breaks = seq(trim_ymin + 0.01, trim_ymax - 0.01, length.out = 4),
                       labels = label_dollar(prefix = "", suffix = "°N")) +
    coord_sf(xlim = c(trim_xmin, trim_xmax), ylim = c(trim_ymin, trim_ymax), expand = FALSE) +
    labs(title = "Spatiotemporal Expansion of the Air Quality Sensor Network and Districts in Yerevan, Armenia") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_line(color = alpha("grey20", 0.12), linewidth = 0.4, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      axis.text   = element_text(size = 7.5, color = "#444444"),
      axis.title  = element_text(size = 8.5, color = "#222222", face = "bold"),
      axis.ticks  = element_line(color = "#888888", linewidth = 0.3),
      plot.title  = element_text(face = "bold", size = 11, hjust = 0),
      plot.title.position = "plot",
      plot.margin = margin(t = 12, r = 10, b = 10, l = 12),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "top",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9)
    )
  
  # --- 5. Inset maps ---
  world_maps <- sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
  armenia_sf <- world_maps %>% dplyr::filter(ID == "Armenia")
  yerevan_point <- st_sfc(st_point(c(44.5152, 40.1792)), crs = 4326)
  yerevan_shadow_sf <- st_convex_hull(st_union(sensors_sf))
  
  p_inset_world <- ggplot() +
    geom_sf(data = world_maps, fill = "#ECECEC", color = "#FFFFFF", size = 0.1) +
    geom_sf(data = armenia_sf, fill = "#1E8449", color = NA) +
    coord_sf(xlim = c(28, 62), ylim = c(30, 50), expand = FALSE) +
    theme_void()
  
  p_inset_regional <- ggplot() +
    geom_sf(data = world_maps, fill = "#E5E8E8", color = "#FFFFFF", size = 0.1) +
    geom_sf(data = armenia_sf, fill = "#F9F9F9", color = "#7F8C8D", linewidth = 0.3) +
    geom_sf(data = yerevan_shadow_sf, fill = "#1E8449", alpha = 0.25, color = NA) +
    geom_sf(data = yerevan_point, color = "#1E8449", size = 2.2, shape = 16) +
    coord_sf(xlim = c(36, 52), ylim = c(36, 44), expand = FALSE) +
    theme_void()
  
  # --- 6. Combine ---
  final_plot <- ggdraw() +
    draw_plot(p_main) +
    draw_plot(p_inset_world,    x = 0.61, y = 0.1, width = 0.17, height = 0.17) +
    draw_plot(p_inset_regional, x = 0.79, y = 0.1, width = 0.17, height = 0.17)
  
  ggsave(out_file, plot = final_plot, width = 9, height = 7, dpi = 600)
  message("Saved: ", out_file)
}

## ============================================================================
## 3. SPATIO-TEMPORAL KRIGING  
## ============================================================================
run_kriging <- function() {
  message("=== Section 3: Spatio-temporal kriging ===")
  nc_file <- file.path(OUT_DIR, "pm25_yerevan_ST_stack.nc")
  if (file.exists(nc_file) && length(list.files(PRED_DIR, "*.tif$")) > 0) {
    message("Kriging outputs already exist. Skipping.")
    return(invisible(NULL))
  }
  
  data <- read_csv(READINGS_CSV, show_col_types = FALSE)
  data$timestamp <- as.POSIXct(gsub("\u00A0", " ", data$timestamp), format = "%Y-%m-%d %H:%M:%S")
  sensors <- read_csv(SENSORS_CSV, show_col_types = FALSE)
  
  sensors.c <- dplyr::filter(sensors, provider %in% c("clarity","purple_air","yerevan_gis") & !is_suspicious)
  lon_col <- intersect(c("longitude","long","lon"), names(sensors.c))[1]
  lat_col <- intersect(c("latitude","lat"), names(sensors.c))[1]
  sensors_coords <- dplyr::select(sensors.c, id, !!sym(lon_col), !!sym(lat_col))
  sensors_coords <- dplyr::rename(sensors_coords, longitude = !!sym(lon_col), latitude = !!sym(lat_col))
  
  data.v <- data |>
    dplyr::filter(sensor_id %in% sensors.c$id) |>
    left_join(sensors_coords, by = c("sensor_id"="id")) |>
    dplyr::filter(!sensor_id %in% c(288,79,92,45)) |>
    dplyr::select(timestamp, sensor_id, pm2.5, longitude, latitude) |>
    dplyr::filter(!is.na(timestamp), !is.na(pm2.5), !is.na(longitude), !is.na(latitude))
  
  yerevan_boundary <- st_read(BOUNDARY_LVL4, layer = "boundary-polygon-lvl4", quiet = TRUE)
  yerevan_proj <- st_transform(yerevan_boundary, 32638)
  grid_centers <- st_make_grid(yerevan_proj, cellsize = 100, what = "centers")
  grid_inside  <- st_intersection(grid_centers, yerevan_proj)
  if (length(grid_inside) == 0) grid_inside <- grid_centers
  grid_sp <- as(grid_inside, "Spatial"); proj4string(grid_sp) <- CRS("EPSG:32638")
  yerevan_proj_sp <- as(yerevan_proj, "Spatial")
  
  sensor_locs <- dplyr::select(data.v, sensor_id, longitude, latitude) |> distinct()
  sf_sensors_utm <- st_transform(st_as_sf(sensor_locs, coords = c("longitude","latitude"), crs = 4326), 32638)
  utm_coords <- st_coordinates(sf_sensors_utm)
  sensor_locs$x_utm <- utm_coords[, "X"]; sensor_locs$y_utm <- utm_coords[, "Y"]
  
  set.seed(42)
  dvar <- data.v[sample(min(nrow(data.v), 1000)), ]
  dvar <- dvar[complete.cases(dvar[, "pm2.5"]), ]
  sill <- var(dvar$pm2.5); if (!is.finite(sill) || sill == 0) sill <- 100
  vgm_fit <- vgmST("sumMetric",
                   space = vgm(sill*0.8, "Sph", 2000, sill*0.2),
                   time  = vgm(sill*0.8, "Exp",   12, sill*0.2),
                   joint = vgm(sill*0.3, "Sph", 1000, 0),
                   stAni = 100)
  
  WINDOW_HOURS <- 24; MAX_OBS <- 200
  
  krigeST_one <- function(target_ts) {
    ts_str <- format(target_ts, "%Y%m%d_%H%M%S")
    res <- list(timestamp = ts_str, success = FALSE, pred_file = NA, var_file = NA, reason = "unknown")
    tryCatch({
      hw <- WINDOW_HOURS * 3600 / 2
      dw <- data.v[data.v$timestamp >= target_ts - hw & data.v$timestamp <= target_ts + hw, ]
      if (nrow(dw) < 3) { res$reason <- "too_few_obs"; return(res) }
      if (nrow(dw) > MAX_OBS) { set.seed(42); dw <- dw[sample(nrow(dw), MAX_OBS), ] }
      dw <- merge(dw, sensor_locs, by = c("sensor_id","longitude","latitude"))
      dw <- dw[complete.cases(dw[, c("x_utm","y_utm","pm2.5","timestamp")]), ]
      dw <- dw[!duplicated(dw[, c("x_utm","y_utm","timestamp")]), ]
      if (nrow(dw) < 3) { res$reason <- "too_few_after_clean"; return(res) }
      
      sp_win <- SpatialPoints(dw[, c("x_utm","y_utm")], proj4string = CRS("EPSG:32638"))
      stData  <- STIDF(sp_win, dw$timestamp, data.frame(pm2.5 = dw$pm2.5))
      stGrid  <- STF(sp = as(grid_sp, "SpatialPoints"), time = c(target_ts, target_ts + 3600))
      st_pred <- tryCatch(krigeST(pm2.5 ~ 1, stData, stGrid, model = vgm_fit, computeVar = TRUE, debug.level = 0), error = function(e) NULL)
      if (is.null(st_pred) || nrow(st_pred@data) == 0) { res$reason <- "krigeST_failed"; return(res) }
      
      n_grid <- length(grid_sp)
      spdf <- SpatialPixelsDataFrame(
        points = SpatialPoints(coordinates(grid_sp), proj4string = CRS("EPSG:32638")),
        data = data.frame(var1.pred = st_pred@data$var1.pred[1:n_grid], var1.var = st_pred@data$var1.var[1:n_grid]))
      pr <- tryCatch(mask(crop(raster(spdf["var1.pred"]), yerevan_proj_sp), yerevan_proj_sp), error = function(e) NULL)
      vr <- tryCatch(mask(crop(raster(spdf["var1.var"]),  yerevan_proj_sp), yerevan_proj_sp), error = function(e) NULL)
      if (is.null(pr) || is.null(vr)) { res$reason <- "crop_mask_failed"; return(res) }
      
      pf <- file.path(PRED_DIR, paste0("pm25_", ts_str, ".tif"))
      vf <- file.path(VAR_DIR,  paste0("var_",  ts_str, ".tif"))
      writeRaster(pr, pf, overwrite = TRUE); writeRaster(vr, vf, overwrite = TRUE)
      res$pred_file <- pf; res$var_file <- vf; res$success <- TRUE; res$reason <- "success"; res
    }, error = function(e) { res$reason <- conditionMessage(e); res })
  }
  
  ts_all <- sort(unique(data.v$timestamp))
  message(sprintf("Processing %d timestamps...", length(ts_all)))
  unlink(file.path(PRED_DIR, "*.tif")); unlink(file.path(VAR_DIR, "*.tif"))
  
  n_cores <- max(1, detectCores() - 1)
  cl <- makeCluster(n_cores, type = ifelse(.Platform$OS.type == "windows", "PSOCK","FORK"))
  clusterExport(cl, c("data.v","sensor_locs","grid_sp","yerevan_proj_sp","vgm_fit","krigeST_one","PRED_DIR","VAR_DIR","WINDOW_HOURS"))
  clusterEvalQ(cl, { library(raster); library(sp); library(gstat); library(spacetime) })
  res <- pblapply(ts_all, krigeST_one, cl = cl)
  stopCluster(cl)
  
  pf <- sort(list.files(PRED_DIR, "*.tif$", full.names = TRUE))
  vf <- sort(list.files(VAR_DIR,  "*.tif$", full.names = TRUE))
  if (length(pf)) {
    writeCDF(rast(pf), nc_file, overwrite = TRUE, compression = 5)
    writeCDF(rast(vf), file.path(OUT_DIR, "pm25_variance_ST_stack.nc"), overwrite = TRUE, compression = 5)
  }
}

## ============================================================================
## 4. GEMM MORTALITY COMPUTATION (Annual / June-June Periods)
## ============================================================================
run_GEMM <- function() {
  message("=== Section 4: GEMM Mortality Computation ===")
  deaths_cache <- file.path(CACHE_DIR, "deaths_rasters.tif")
  se_cache <- file.path(CACHE_DIR, "se_rasters.tif")
  eff_cache <- file.path(CACHE_DIR, "eff_rasters.tif")
  summary_csv <- file.path(RUN_DIR, "yerevan_mortality_all_periods.csv")
  
  if (file.exists(deaths_cache) && file.exists(se_cache) && file.exists(eff_cache)) {
    message("GEMM outputs already exist. Skipping.")
    return(invisible(NULL))
  }
  
  pop_raster_raw <- rast(POP_TIF); crs(pop_raster_raw) <- "EPSG:4326"
  if (!file.exists(GEMM_XLSX)) {
    alt <- file.path(dirname(PROJ_DIR), basename(GEMM_XLSX))
    if (file.exists(alt)) GEMM_XLSX <- alt
  }
  if (!file.exists(GEMM_XLSX)) {
    stop("Coefficient file not found: ", GEMM_XLSX,
         "\nPlace GEMMcoefficients_NCDLRI_GHE2023.xlsx in PROJ_DIR (", PROJ_DIR, ").")
  }
  coeff_all <- read_excel(GEMM_XLSX)
  colnames(coeff_all) <- c("Cause","AgeGroup","Theta","Theta_SD","Alpha","Mu","Nu","Deaths23","Pop23","DeathRate")
  coeff_all <- coeff_all[complete.cases(coeff_all), ]
  coeff_df  <- coeff_all %>% dplyr::filter(Cause == "NCD+LRI")
  
  # Burnett GEMM NCD+LRI applies to adults 25+ ONLY. Children contribute zero
  # attributable deaths by design; they map to NA and are skipped below. The old
  # pooled "0-25" row was a mislabel of the pooled >25 coefficient and is removed
  # from the coefficient file.
  age_map <- c("0-9"=NA,"1-10"=NA,"5-14"=NA,"10-19"=NA,"15-24"=NA,"20-29"="25-30","25-34"="30-35",
               "30-39"="35-40","35-44"="40-45","40-49"="45-50","45-54"="50-55","50-59"="55-60","55-64"="60-65","60-69"="65-70",
               "60-70"="65-70","70-79"="75-80","75-84"="80+","80-89"="80+","85-94"="80+","90-99"="80+")
  age_order <- c("25-30","30-35","35-40","40-45","45-50","50-55","55-60","60-65","65-70","70-75","75-80","80+")
  causes_list <- unique(coeff_df$Cause)
  coef_ages   <- unique(coeff_df$AgeGroup)
  coef_matrix <- array(NA, dim = c(length(causes_list), length(coef_ages), 5),
                       dimnames = list(causes_list, coef_ages, c("theta","theta_sd","alpha","mu","nu")))
  for (i in seq_along(causes_list)) {
    for (j in seq_along(coef_ages)) {
      row <- coeff_df[coeff_df$Cause == causes_list[i] & coeff_df$AgeGroup == coef_ages[j], ]
      if (nrow(row) > 0) coef_matrix[i, j, ] <- c(row$Theta[1], row$Theta_SD[1], row$Alpha[1], row$Mu[1], row$Nu[1])
    }
  }
  
  calculate_gemm_rr <- function(pm25, theta, alpha, mu, nu, cf = 2.4) {
    z <- pm25 - cf; z[z < 0] <- 0
    omega <- 1 / (1 + exp(-(z - mu) / nu))
    exp(theta * log((z / alpha) + 1) * omega)
  }
  
  pop_raster_raw <- project(pop_raster_raw, "EPSG:32638")
  first_r <- rast(list.files(PRED_DIR, pattern = "*.tif$", full.names = TRUE)[1])
  pop_raster <- crop(pop_raster_raw, ext(first_r))
  
  pred_files <- list.files(PRED_DIR, pattern = "*.tif$", full.names = TRUE)
  pred_dates <- as.POSIXct(gsub("pm25_|.tif","",basename(pred_files)), format="%Y%m%d_%H%M%S", tz="UTC")
  
  periods <- list(
    list(start = as.POSIXct("2023-06-01 00:00:00", tz="UTC"), end = as.POSIXct("2024-06-01 00:00:00", tz="UTC")),
    list(start = as.POSIXct("2024-06-01 00:00:00", tz="UTC"), end = as.POSIXct("2025-06-01 00:00:00", tz="UTC")),
    list(start = as.POSIXct("2025-06-01 00:00:00", tz="UTC"), end = as.POSIXct("2026-06-01 00:00:00", tz="UTC"))
  )
  
  deaths_rasters_list <- se_rasters_list <- eff_rasters_list <- list()
  all_periods_summary <- data.frame()
  
  for (p_idx in seq_along(periods)) {
    p <- periods[[p_idx]]
    p_label <- PERIOD_LABELS[p_idx]
    p_files <- pred_files[pred_dates >= p$start & pred_dates < p$end]
    if(length(p_files) == 0) next
    message(sprintf("Processing %s (%d files)...\n", p_label, length(p_files)))
    
    pm25_annual_avg <- mean(rast(p_files), na.rm = TRUE)
    pm25_resampled  <- resample(pm25_annual_avg, pop_raster[[1]], method = "bilinear")
    
    death_rates <- coeff_df[, c("AgeGroup","DeathRate")]
    n_pop_cells <- ncell(pop_raster[[1]])
    deaths_by_pop_cell <- rep(0, n_pop_cells)
    var_by_pop_cell <- rep(0, n_pop_cells)
    pm25_pop_values <- values(pm25_resampled)
    
    age_deaths_summary <- data.frame(RasterAgeGroup=character(), CoefAge=character(), Deaths=numeric(), BaselineDeaths=numeric(), MeanRR=numeric(), MeanPAF=numeric(), MeanAbsoluteRisk=numeric(), stringsAsFactors = FALSE)
    
    for (ag in names(pop_raster)) {
      if (!ag %in% names(age_map)) next
      pop_values <- values(pop_raster[[ag]])
      coef_age <- unname(age_map[ag])
      if (is.na(coef_age)) next  # under-25: outside GEMM NCD+LRI scope
      age_dr <- death_rates$DeathRate[death_rates$AgeGroup == coef_age]
      if (length(age_dr) == 0 || is.na(age_dr[1])) next
      age_dr <- age_dr[1]
      idx <- which(coef_ages == coef_age)
      if (length(idx) == 0 || is.na(coef_matrix[1, idx, 1])) next
      
      theta <- coef_matrix[1, idx, 1]; theta_sd <- coef_matrix[1, idx, 2]
      alpha <- coef_matrix[1, idx, 3]; mu <- coef_matrix[1, idx, 4]; nu <- coef_matrix[1, idx, 5]
      
      rr_values <- calculate_gemm_rr(pm25_pop_values, theta, alpha, mu, nu)
      paf_values <- (rr_values - 1) / rr_values; paf_values[is.na(paf_values)] <- 0
      baseline_death_values <- pop_values * age_dr
      death_values <- paf_values * baseline_death_values
      # Delta-method: Var(D) = (B * T(z)/RR)^2 * Var(theta), where
      # T(z) = log1p(z/alpha) * omega(z); previously the T(z)/RR factor was omitted,
      # massively overstating uncertainty for low PM2.5.
      z_var <- pm25_pop_values - 2.4; z_var[is.na(z_var) | z_var < 0] <- 0
      omega_var <- 1 / (1 + exp(-(z_var - mu) / nu))
      tz_values <- log1p(z_var / alpha) * omega_var
      var_values <- (baseline_death_values * tz_values / rr_values)^2 * theta_sd^2
      var_values[is.na(var_values)] <- 0; death_values[is.na(death_values)] <- 0
      
      deaths_by_pop_cell <- deaths_by_pop_cell + death_values
      var_by_pop_cell <- var_by_pop_cell + var_values
      age_deaths_summary <- rbind(age_deaths_summary, data.frame(RasterAgeGroup=ag, CoefAge=coef_age, Deaths=sum(death_values, na.rm=TRUE), BaselineDeaths=sum(baseline_death_values, na.rm=TRUE), MeanRR=mean(rr_values, na.rm=TRUE), MeanPAF=mean(paf_values, na.rm=TRUE), MeanAbsoluteRisk=mean(paf_values * age_dr, na.rm=TRUE)))
    }
    
    deaths_raster <- pop_raster[[1]]; values(deaths_raster) <- deaths_by_pop_cell
    se_raster <- pop_raster[[1]]; values(se_raster) <- sqrt(var_by_pop_cell)
    pop_total <- sum(pop_raster, na.rm = TRUE)
    pop_total_values <- values(pop_total)
    
    deaths_per_cap_pm25_values <- ifelse(!is.na(pm25_pop_values) & pm25_pop_values > 2.4 & !is.na(pop_total_values) & pop_total_values > 0 & deaths_by_pop_cell > 0, (deaths_by_pop_cell / pop_total_values) / pm25_pop_values, NA)
    eff_raster <- deaths_raster; values(eff_raster) <- deaths_per_cap_pm25_values * 100000
    values(eff_raster)[!is.finite(values(eff_raster))] <- NA
    
    deaths_rasters_list[[p_idx]] <- deaths_raster
    se_rasters_list[[p_idx]] <- se_raster
    eff_rasters_list[[p_idx]] <- eff_raster
    
    total_deaths <- sum(deaths_by_pop_cell, na.rm = TRUE)
    total_se <- sqrt(sum(var_by_pop_cell, na.rm = TRUE))
    all_periods_summary <- rbind(all_periods_summary, data.frame(Period=p_label, TotalDeaths=total_deaths, SE=total_se, CI_Lower=total_deaths - 1.96*total_se, CI_Upper=total_deaths + 1.96*total_se))
    
    # Generate period-specific plots
    df1 <- rast_to_wgs84_df(pm25_resampled, "PM25")
    p1 <- safe_spatial_plot(df1, "PM25", "PM2.5\n(µg/m³)", paste0("Annual Average PM2.5 (", p_label, ")"), palette = "YlOrRd")
    if (!is.null(p1)) ggsave(file.path(YEAR_PLOTS_DIR, paste0("p", p_idx, "_plot1_pm25_map.png")), p1, width=9, height=7, dpi=300)
  }
  
  writeRaster(rast(deaths_rasters_list), deaths_cache, overwrite=TRUE)
  writeRaster(rast(se_rasters_list), se_cache, overwrite=TRUE)
  writeRaster(rast(eff_rasters_list), eff_cache, overwrite=TRUE)
  write.csv(all_periods_summary, summary_csv, row.names=FALSE)
  cat("\n================ PM2.5-Attributable Deaths per Period ================\n")
  print(all_periods_summary)
}

# Helper for GEMM plotting
rast_to_wgs84_df <- function(r, value_name = "value") {
  r[!is.finite(r)] <- NA
  r_wgs84 <- tryCatch(project(r, "EPSG:4326", method = "bilinear"), error = function(e) NULL)
  if (is.null(r_wgs84)) return(data.frame())
  df <- as.data.frame(r_wgs84, xy = TRUE, na.rm = FALSE)
  if (ncol(df) < 3) return(data.frame())
  colnames(df)[3] <- value_name
  df[is.finite(df$x) & is.finite(df$y) & is.finite(df[[value_name]]), ]
}
safe_spatial_plot <- function(df, fill_col, fill_name, title, subtitle = NULL, palette = "YlOrRd", direction = 1, is_viridis = FALSE) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  p_base <- ggplot(df, aes(x = x, y = y, fill = .data[[fill_col]])) + geom_raster(alpha = 0.75) + coord_sf(crs = 4326, expand = FALSE) + labs(title = title, subtitle = subtitle) + theme_minimal() + theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_text(size = 8, color = "gray40"), legend.position = "right", plot.title = element_text(face = "bold", size = 12))
  if (is_viridis) p_base <- p_base + scale_fill_viridis_c(na.value = NA, name = fill_name)
  else p_base <- p_base + scale_fill_distiller(palette = palette, direction = direction, na.value = NA, name = fill_name)
  p_osm <- tryCatch({ p_base + annotation_map_tile(type = "osm", zoom = 12, progress = "none", quiet = TRUE) }, error = function(e) NULL)
  if (!is.null(p_osm)) if (tryCatch({ ggplotGrob(p_osm); TRUE }, error = function(e) FALSE)) return(p_osm)
  return(p_base)
}

## ============================================================================
## 5. FIGURE 3 - Spatial PM2.5 heatmaps (January & July, by year)
## ============================================================================
figure_3 <- function() {
  message("=== Figure 3: PM2.5 spatial heatmaps (Jan/July) ===")
  out_file <- file.path(OUT_DIR, "Figure3_seasonal_pm25.png")
  if (file.exists(out_file)) {
    message("Output already exists. Skipping.")
    return(invisible(NULL))
  }
  
  pred_files <- list.files(PRED_DIR, "*.tif$", full.names = TRUE)
  t_strs <- gsub("pm25_|.tif", "", basename(pred_files))
  t_info <- as.POSIXct(t_strs, format = "%Y%m%d_%H%M%S", tz = "UTC")
  df_files <- data.frame(file = pred_files, time = t_info)
  df_files$year  <- year(df_files$time); df_files$month <- month(df_files$time)
  
  idx <- which(df_files$year %in% 2024:2026 & df_files$month %in% c(1, 7))
  r_sub  <- rast(df_files$file[idx])
  groups <- paste(year(df_files$time[idx]), month(df_files$time[idx]), sep = "_")
  r_monthly <- tapp(r_sub, groups, fun = "mean", na.rm = TRUE)
  
  df_long <- as.data.frame(r_monthly, xy = TRUE, na.rm = TRUE) |>
    pivot_longer(-c(x, y), names_to = "Yr_Mo", values_to = "pm25") |>
    separate(Yr_Mo, c("Year","Month"), sep = "_") |>
    mutate(Year = as.integer(Year), Month = factor(ifelse(Month == "1", "January", "July"), levels = c("January","July")))
  
  p <- ggplot(df_long, aes(x, y, fill = pm25)) + geom_raster() +
    scale_fill_viridis_c(option = "inferno", direction = -1) +
    facet_grid(Year ~ Month) + coord_equal() + theme_minimal() +
    labs(title = "Yerevan Seasonal PM2.5 Comparison", x = "", y = "") +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), strip.text = element_text(size = 12, face = "bold"))
  ggsave(out_file, p, width = 8, height = 12, dpi = 300)
}

## ============================================================================
## 6. FIGURE 4 - Citywide pop-weighted time series
## ============================================================================
figure_4 <- function() {
  message("=== Figure 4: Citywide pop-weighted time series ===")
  out_file <- file.path(OUT_DIR, "Figure4_citywide_timeseries.png")
  ts_cache <- file.path(CACHE_DIR, "citywide_pop_weighted_means.csv")
  
  if (file.exists(out_file)) {
    message("Output already exists. Skipping.")
    return(invisible(NULL))
  }
  
  if (file.exists(ts_cache)) {
    df <- read.csv(ts_cache)
    df$time <- as.POSIXct(df$time, tz="UTC")
  } else {
    pred_files <- list.files(PRED_DIR, "*.tif$", full.names = TRUE)
    t_info <- as.POSIXct(gsub("pm25_|.tif","",basename(pred_files)), format = "%Y%m%d_%H%M%S", tz = "UTC")
    
    pop_r <- rast(POP_TIF)
    if (nlyr(pop_r) > 1) pop_r <- sum(pop_r, na.rm = TRUE)
    sample_r <- rast(pred_files[1])
    pop_aligned <- if (same.crs(pop_r, sample_r)) resample(pop_r, sample_r, method = "bilinear") else project(pop_r, sample_r)
    pop_aligned <- mask(pop_aligned, sample_r)
    pop_total <- global(pop_aligned, "sum", na.rm = TRUE)[1, 1]
    
    batch <- 200; n <- length(pred_files); means <- numeric(n)
    pb <- txtProgressBar(0, n, style = 3)
    for (i in seq(1, n, by = batch)) {
      ei <- min(i + batch - 1, n); rb <- rast(pred_files[i:ei])
      means[i:ei] <- global(rb * pop_aligned, "sum", na.rm = TRUE)[, 1] / pop_total
      setTxtProgressBar(pb, ei)
    }; close(pb)
    
    df <- data.frame(time = t_info, pm25 = means) |> dplyr::filter(!is.na(time), !is.na(pm25))
    write.csv(df, ts_cache, row.names = FALSE)
  }
  
  df$pm25_7d <- rollmean(df$pm25, k = 168, fill = NA, align = "center")
  
  p <- ggplot(df, aes(time, pm25)) +
    geom_point(aes(color = "Pop-Weighted Hourly Exposure"), alpha = 0.3, size = 0.5) +
    geom_line(aes(y = pm25_7d, color = "7-Day Moving Average"), linewidth = 1) +
    geom_hline(aes(yintercept = 15, color = "WHO Guideline: Maximum Daily Mean (15 µg/m³)"), linetype = "dashed", linewidth = 0.6) +
    scale_color_manual(values = c("Pop-Weighted Hourly Exposure" = "gray60", "7-Day Moving Average" = "firebrick", "WHO Guideline: Maximum Daily Mean (15 µg/m³)" = "black")) +
    scale_x_datetime(date_breaks = "3 months", date_labels = "%b %Y") + theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(title = "Average PM2.5 Levels in Yerevan, 2023-2026 (Constant Population-Weighted Exposure)", x = "Date", y = "PM2.5 (µg/m³)")
  ggsave(out_file, p, width = 12, height = 5, dpi = 300)
}

## ============================================================================
## 7. FIGURE 5 - Excess mortality quadrant maps (3 periods)
## ============================================================================
figure_5 <- function() {
  message("=== Figure 5: Excess mortality maps (3 periods) ===")
  out_file <- file.path(QUAD_DIR, "Figure5_mortality_quadrant.png")
  if (file.exists(out_file)) {
    message("Output already exists. Skipping.")
    return(invisible(NULL))
  }
  
  district_info <- load_districts()
  yerevan_districts <- district_info$sf; name_col <- district_info$name_col
  districts_sv <- vect(st_transform(yerevan_districts, 32638))
  
  deaths_cache <- file.path(CACHE_DIR, "deaths_rasters.tif")
  if (!file.exists(deaths_cache)) stop("Run GEMM section first - ", deaths_cache, " missing.")
  deaths_stack <- rast(deaths_cache); crs(deaths_stack) <- "EPSG:32638"
  if (crs(deaths_stack) != crs(districts_sv)) deaths_stack <- project(deaths_stack, districts_sv)
  
  pop_raw <- rast(POP_TIF); crs(pop_raw) <- "EPSG:4326"; pop_raw <- project(pop_raw, "EPSG:32638")
  first_r <- rast(list.files(PRED_DIR, "*.tif$", full.names = TRUE)[1])
  pop_r <- crop(pop_raw, ext(first_r)); pop_total <- pop_r[[1]]
  values(pop_total) <- rowSums(values(pop_r), na.rm = TRUE)
  
  deaths_rasters <- lapply(seq_len(nlyr(deaths_stack)), function(i) deaths_stack[[i]])
  deaths_pc_rasters <- lapply(deaths_rasters, function(r) {
    pop_al <- resample(pop_total, r, method = "near")
    rv <- values(r); pv <- values(pop_al)
    values(r) <- (rv / pv) * 1000; values(r)[!is.finite(values(r))] <- NA; r
  })
  
  make_map <- function(r, lims, title) {
    r <- mask(r, districts_sv); df <- rast_to_wgs84_df(r, "val")
    ggplot(df, aes(x, y, fill = val)) + geom_raster() +
      scale_fill_gradientn(colors = colorRampPalette(c("#2ECC71","#F1C40F","#E67E22","#E74C3C","#9B59B6"))(100), limits = lims, oob = squish, guide = "none") +
      geom_sf(data = yerevan_districts, fill = NA, color = "black", linewidth = 0.4, inherit.aes = FALSE) +
      geom_sf_text(data = yerevan_districts, aes(label = .data[[name_col]]), color = "black", size = 2.5, fontface = "bold", inherit.aes = FALSE, check_overlap = FALSE, lineheight = 0.8) +
      annotation_north_arrow(location = "tr", which_north = "true",
                             style = north_arrow_fancy_orienteering(),
                             height = unit(1.2,"cm"), width = unit(1.2,"cm")) +
      annotation_scale(location = "bl", width_hint = 0.3, unit_category = "metric") +
      coord_sf(crs = 4326, expand = FALSE) +
      labs(title = title) + theme_minimal() +
      theme(panel.grid = element_blank(), axis.title = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 14, face = "bold"),
            plot.margin = margin(0,0,0,0))
  }
  
  all_vals <- unlist(lapply(deaths_pc_rasters, function(r) values(r)[!is.na(values(r))]))
  lims <- c(0, max(all_vals, na.rm = TRUE))
  maps <- lapply(seq_along(deaths_pc_rasters), function(i) make_map(deaths_pc_rasters[[i]], lims, PERIOD_LABELS[i]))
  
  dummy <- ggplot(data.frame(val = lims), aes(1, 1, fill = val)) + geom_tile(alpha = 0) +
    scale_fill_gradientn(colors = colorRampPalette(c("#2ECC71","#F1C40F","#E67E22","#E74C3C","#9B59B6"))(100), limits = lims,
                         guide = guide_colorbar(barwidth = unit(2.25,"cm"), barheight = unit(9.4,"cm"),
                                                title.position = "top", title.hjust = 0.5,
                                                frame.colour = "black", ticks.colour = "black"),
                         name = "Annual Excess\nDeaths per\n1,000\nResidents", oob = squish) +
    theme_void() +
    theme(legend.position = c(0.5, 0.5),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12))
  
  p <- (maps[[1]] | maps[[2]]) / (maps[[3]] | dummy)
  ggsave(out_file, p, width = 12, height = 12, dpi = 300)
}
## ============================================================================
## 8. FIGURE 6 - District-level monthly PM2.5 (pop-weighted)
## ============================================================================
figure_6 <- function() {
  message("=== Figure 6: District monthly pop-weighted PM2.5 ===")
  out_file <- file.path(OUT_DIR, "Figure6_district_pm25.png")
  # PM2.5-only product: shared location, not run-specific
  pm25_cache <- file.path(OUT_DIR, "pm25_popweighted_by_district_month.csv")
  
  # Skip only if BOTH the PNG and the CSV cache exist (figure_7 needs the CSV)
  if (file.exists(out_file) && file.exists(pm25_cache)) {
    message("Output already exists. Skipping.")
    return(invisible(NULL))
  }
  
  district_info <- load_districts()
  yerevan_districts <- district_info$sf; district_names <- yerevan_districts[[district_info$name_col]]
  districts_sv <- vect(st_transform(yerevan_districts, 32638))
  
  if (file.exists(pm25_cache)) {
    pm25_df <- read.csv(pm25_cache); pm25_df$month_date <- as.Date(pm25_df$month_date)
  } else {
    pop_raw <- rast(POP_TIF); crs(pop_raw) <- "EPSG:4326"; pop_raw <- project(pop_raw, "EPSG:32638")
    first_r <- rast(list.files(PRED_DIR, "*.tif$", full.names = TRUE)[1])
    pop_r <- resample(pop_raw, first_r, method = "bilinear")
    pop_total <- sum(pop_r, na.rm = TRUE) |> mask(first_r)
    
    pred_files <- list.files(PRED_DIR, "*.tif$", full.names = TRUE)
    file_months <- format(as.POSIXct(gsub("pm25_|.tif","",basename(pred_files)), format = "%Y%m%d_%H%M%S", tz = "UTC"), "%Y-%m")
    unique_months <- sort(unique(file_months)); pm25_list <- list()
    
    for (m in unique_months) {
      m_files <- pred_files[file_months == m]; pm25_avg <- mean(rast(m_files), na.rm = TRUE)
      pop_eff <- mask(pop_total, pm25_avg); pxp <- pm25_avg * pop_eff
      z_num <- extract(pxp, districts_sv, fun = "sum", na.rm = TRUE)
      z_den <- extract(pop_eff, districts_sv, fun = "sum", na.rm = TRUE)
      pm25_list[[m]] <- data.frame(district = district_names, month = m, month_date = as.Date(paste0(m, "-01")), pm25 = z_num[[2]] / z_den[[2]])
    }
    pm25_df <- do.call(rbind, pm25_list); write.csv(pm25_df, pm25_cache, row.names = FALSE)
  }
  
  district_order <- pm25_df |> group_by(district) |> summarise(avg = mean(pm25, na.rm = TRUE)) |> arrange(desc(avg)) |> pull(district)
  pm25_df$district <- factor(pm25_df$district, levels = district_order)
  colors <- rev(colorRampPalette(c("#2ECC71","#F1C40F","#E67E22","#E74C3C","#9B59B6"))(length(district_order)))
  names(colors) <- district_order
  
  p <- ggplot(pm25_df, aes(month_date, pm25, color = district, group = district)) + geom_line(linewidth = 0.7) +
    scale_color_manual(values = colors, name = "District") + scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
    coord_cartesian(xlim = as.Date(c("2024-06-01", NA))) +
    labs(title = "Population-Weighted PM2.5 Exposure by District Over Time", subtitle = "Monthly population-weighted average PM2.5 (µg/m³)", x = NULL, y = "PM2.5 (µg/m³)") +
    theme_minimal() + theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 10, color = "gray40"), axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "right")
  ggsave(out_file, p, width = 12, height = 6, dpi = 300)
}

## ============================================================================
## 9. FIGURE 7 - District-level monthly mortality heatmap
## ============================================================================
figure_7 <- function() {
  message("=== Figure 7: District monthly mortality heatmap ===")
  out_file <- file.path(RUN_DIR, "Figure7_mortality_heatmap.png")
  if (file.exists(out_file)) {
    message("Output already exists. Skipping.")
    return(invisible(NULL))
  }
  
  pm25_cache <- file.path(OUT_DIR, "pm25_popweighted_by_district_month.csv")
  deaths_cache <- file.path(CACHE_DIR, "deaths_rasters.tif")
  if (!file.exists(pm25_cache) || !file.exists(deaths_cache)) stop("Run Sections 4 & 8 first.")
  
  district_info <- load_districts()
  yerevan_districts <- district_info$sf; district_names <- yerevan_districts[[district_info$name_col]]
  districts_sv <- vect(st_transform(yerevan_districts, 32638))
  
  pop_raw <- rast(POP_TIF); crs(pop_raw) <- "EPSG:4326"; pop_raw <- project(pop_raw, "EPSG:32638")
  first_r <- rast(list.files(PRED_DIR, "*.tif$", full.names = TRUE)[1])
  pop_r <- crop(pop_raw, ext(first_r)); pop_total <- pop_r[[1]]
  values(pop_total) <- rowSums(values(pop_r), na.rm = TRUE)
  z_pop <- extract(pop_total, districts_sv, fun = "sum", na.rm = TRUE)
  district_pop <- z_pop[[2]]
  
  deaths_stack <- rast(deaths_cache); crs(deaths_stack) <- "EPSG:32638"
  if (crs(deaths_stack) != crs(districts_sv)) deaths_stack <- project(deaths_stack, districts_sv)
  
  mortality_df <- data.frame()
  for (i in seq_len(min(nlyr(deaths_stack), length(PERIOD_LABELS)))) {
    r <- deaths_stack[[i]]; z <- extract(r, districts_sv, fun = "sum", na.rm = TRUE)
    mortality_df <- rbind(mortality_df, data.frame(district = district_names, period = PERIOD_LABELS[i], period_idx = i, deaths_per_1k = (z[[2]] / district_pop) * 1000))
  }
  
  pm25_df <- read.csv(pm25_cache)
  district_order <- pm25_df |> group_by(district) |> summarise(avg = mean(pm25, na.rm = TRUE)) |> arrange(desc(avg)) |> pull(district)
  mortality_df$district <- factor(mortality_df$district, levels = rev(district_order))
  mortality_df$period <- factor(mortality_df$period, levels = PERIOD_LABELS)
  val_mid <- median(mortality_df$deaths_per_1k, na.rm = TRUE)
  mortality_df$text_color <- ifelse(mortality_df$deaths_per_1k > val_mid, "white", "black")
  
  p <- ggplot(mortality_df, aes(period, district, fill = deaths_per_1k)) + geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = sprintf("%.2f", deaths_per_1k), color = text_color), size = 3.5, fontface = "bold", show.legend = FALSE) +
    scale_fill_gradientn(colors = colorRampPalette(c("#2ECC71","#F1C40F","#E67E22","#E74C3C","#9B59B6"))(100), name = "Excess Deaths\nper 1,000") +
    scale_color_manual(values = c(black = "black", white = "white"), guide = "none") +
    scale_x_discrete(position = "top", labels = c("June 2023-\nJune 2024", "June 2024-\nJune 2025", "June 2025-\nJune 2026")) +
    labs(title = "Excess Mortality by District Over Time", x = NULL, y = NULL) + theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 10), legend.position = "right", panel.grid = element_blank())
  ggsave(out_file, p, width = 10, height = 8, dpi = 300)
}

## ============================================================================
## 10. SUPPLEMENTARY: Percent Elderly (60+) by District
## ============================================================================
figure_supp_elderly <- function() {
  message("=== Supplementary: Elderly % Heatmap ===")
  out_file <- file.path(PROJ_DIR, "yerevan_pct_elderly_by_district.png")
  if (file.exists(out_file)) {
    message("Output already exists. Skipping.")
    return(invisible(NULL))
  }
  
  district_info <- load_districts()
  yerevan_districts <- district_info$sf; name_col <- district_info$name_col
  
  pop_raster <- rast(POP_TIF); crs(pop_raster) <- "EPSG:4326"
  pop_raster <- crop(pop_raster, ext(yerevan_districts))
  band_names <- names(pop_raster)
  
  get_lower <- function(nm) as.numeric(sub("-.*", "", nm))
  lower_bounds <- sapply(band_names, get_lower)
  elderly_bands <- band_names[!is.na(lower_bounds) & lower_bounds >= 60]
  if (length(elderly_bands) == 0) elderly_bands <- tail(band_names, 4)
  
  pop_total <- sum(pop_raster, na.rm = TRUE); pop_elderly <- sum(pop_raster[[elderly_bands]], na.rm = TRUE)
  districts_sv <- vect(yerevan_districts)
  
  z_total <- extract(pop_total, districts_sv, fun = "sum", na.rm = TRUE)
  z_elderly <- extract(pop_elderly, districts_sv, fun = "sum", na.rm = TRUE)
  yerevan_districts$pct_elderly <- 100 * z_elderly[[2]] / z_total[[2]]
  
  p <- ggplot(yerevan_districts) + geom_sf(aes(fill = pct_elderly), color = "grey20", linewidth = 0.4) +
    geom_sf_text(aes(label = .data[[name_col]]), color = "black", size = 3, fontface = "bold", check_overlap = TRUE, lineheight = 0.8) +
    scale_fill_gradient(low = "red", high = "green", name = "% aged\n60+", na.value = "grey80") +
    coord_sf(crs = 4326, expand = FALSE) + labs(title = "Age Structure of Yerevan by District", subtitle = "Percentage of population aged 60+ (WorldPop age-structured raster, aggregated to districts)") +
    theme_minimal() + theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.key.height = unit(1.6, "cm"), plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 10, color = "gray40"))
  ggsave(out_file, p, width = 9, height = 8, dpi = 300)
}
## ============================================================================
## 11. GRAPHICAL ABSTRACT - Figure 5 maps + workflow pipeline (Elsevier: 531 x 1328 px)
## ============================================================================
graphical_abstract <- function() {
  message("=== Graphical Abstract ===")
  out_file <- file.path(RUN_DIR, "GraphicalAbstract.png")
  out_tif  <- file.path(RUN_DIR, "GraphicalAbstract.tif")
  
  # ---- Data prep (mirrors Figure 5) ----
  district_info <- load_districts()
  yerevan_districts <- district_info$sf; name_col <- district_info$name_col
  districts_sv <- vect(st_transform(yerevan_districts, 32638))
  
  deaths_cache <- file.path(CACHE_DIR, "deaths_rasters.tif")
  if (!file.exists(deaths_cache)) stop("Run GEMM section first - ", deaths_cache, " missing.")
  deaths_stack <- rast(deaths_cache); crs(deaths_stack) <- "EPSG:32638"
  if (crs(deaths_stack) != crs(districts_sv)) deaths_stack <- project(deaths_stack, districts_sv)
  
  pop_raw <- rast(POP_TIF); crs(pop_raw) <- "EPSG:4326"; pop_raw <- project(pop_raw, "EPSG:32638")
  first_r <- rast(list.files(PRED_DIR, "*.tif$", full.names = TRUE)[1])
  pop_r <- crop(pop_raw, ext(first_r)); pop_total <- pop_r[[1]]
  values(pop_total) <- rowSums(values(pop_r), na.rm = TRUE)
  
  deaths_pc <- lapply(seq_len(nlyr(deaths_stack)), function(i) {
    r <- deaths_stack[[i]]
    pop_al <- resample(pop_total, r, method = "near")
    rv <- values(r); pv <- values(pop_al)
    values(r) <- (rv / pv) * 1000
    values(r)[!is.finite(values(r))] <- NA
    mask(r, districts_sv)
  })
  
  # Same palette, shared limits and district labeling as Figure 5
  pal <- colorRampPalette(c("#2ECC71","#F1C40F","#E67E22","#E74C3C","#9B59B6"))(100)
  all_vals <- unlist(lapply(deaths_pc, function(r) values(r)[!is.na(values(r))]))
  lims <- c(0, max(all_vals, na.rm = TRUE))
  
  ga_map <- function(r, title) {
    df <- rast_to_wgs84_df(r, "val")
    ggplot(df, aes(x, y, fill = val)) + geom_raster() +
      scale_fill_gradientn(colors = pal, limits = lims, oob = squish, guide = "none") +
      geom_sf(data = yerevan_districts, fill = NA, color = "black", linewidth = 0.25,
              inherit.aes = FALSE) +
      geom_sf_text(data = yerevan_districts, aes(label = .data[[name_col]]),
                   color = "black", size = 2.5, fontface = "bold",
                   inherit.aes = FALSE, check_overlap = FALSE, lineheight = 0.8) +
      coord_sf(crs = 4326, expand = FALSE) +
      labs(title = title) +
      theme_void(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5,
                                      margin = margin(b = 2)),
            plot.margin = margin(1, 2, 1, 2))
  }
  
  # ---- Side legend: identical colorbar to Figure 5 ----
  legend_plot <- ggplot(data.frame(val = lims), aes(1, 1, fill = val)) + geom_tile(alpha = 0) +
    scale_fill_gradientn(colors = pal, limits = lims, oob = squish,
                         name = "Annual Excess\nDeaths per\n1,000\nResidents",
                         guide = guide_colorbar(barwidth = unit(1.6, "cm"),
                                                barheight = unit(6.5, "cm"),
                                                title.position = "top", title.hjust = 0.5,
                                                frame.colour = "black", ticks.colour = "black")) +
    theme_void() +
    theme(legend.position = c(0.5, 0.5),
          legend.title = element_text(face = "bold", size = 13),
          legend.text  = element_text(size = 11))
  
  # ---- Total excess deaths per period (labels placed below each map) ----
  summary_csv <- file.path(RUN_DIR, "yerevan_mortality_all_periods.csv")
  if (file.exists(summary_csv)) {
    period_totals <- round(read.csv(summary_csv)$TotalDeaths)
  } else {
    period_totals <- round(sapply(deaths_pc, function(r) sum(values(r), na.rm = TRUE)))
  }
  period_totals <- period_totals[seq_len(min(length(period_totals), length(deaths_pc)))]
  
  total_label <- function(n) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("Total Excess Deaths: ", format(n, big.mark = ",")),
               size = 4.5, fontface = "bold") +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_void()
  }
  blank_panel <- ggplot() + theme_void()
  totals_row <- NULL  # folded into the single 8-panel grid below
  
  # ---- Workflow pipeline strip (bottom) ----
  steps <- data.frame(
    x    = c(1.35, 4.05, 6.75, 9.45),
    y    = rep(1, 4),
    lab  = c("333 low-cost\nPM2.5 sensors",
             "Spatio-temporal\nkriging (100 m)",
             "GEMM health\nimpact model",
             "Excess mortality\nmaps by period"),
    fill = c("#DEEBF7", "#9ECAE1", "#4292C6", "#08519C")
  )
  workflow <- ggplot(steps) +
    geom_tile(aes(x, y, width = 2.35, height = 0.85, fill = fill)) +
    geom_text(aes(x, y, label = lab), size = 6.5, fontface = "bold", lineheight = 0.95,
              data = steps[1:3, ], color = "black") +
    geom_text(aes(x, y, label = lab), size = 6.5, fontface = "bold", lineheight = 0.95,
              data = steps[4, ], color = "white") +
    annotate("segment", x = c(2.55, 5.25, 7.95), xend = c(2.85, 5.55, 8.25),
             y = 1, yend = 1,
             arrow = arrow(length = unit(0.4, "cm"), type = "closed"), linewidth = 1.4) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0, 10.8), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.35, 1.65), expand = c(0, 0)) +
    theme_void(base_size = 1) +
    theme(plot.margin = margin(2, 0, 0, 0))
  
  # ---- Compose: single 8-panel grid (avoids deeply nested patchwork), workflow below ----
  maps_row <- ga_map(deaths_pc[[1]], PERIOD_LABELS[1]) +
    ga_map(deaths_pc[[2]], PERIOD_LABELS[2]) +
    ga_map(deaths_pc[[3]], PERIOD_LABELS[3]) +
    legend_plot +
    total_label(period_totals[1]) +
    total_label(period_totals[ifelse(length(period_totals) >= 2, 2, 1)]) +
    total_label(period_totals[ifelse(length(period_totals) >= 3, 3, 1)]) +
    blank_panel +
    plot_layout(ncol = 4, widths = c(1, 1, 1, 0.52), heights = c(3, 0.28))
  ga <- maps_row / workflow + plot_layout(heights = c(1, 0.8))
  
  ggsave(out_file, ga, width = 1600, height = 1100, units = "px", dpi = 96)
  ggsave(out_tif,  ga, width = 1600, height = 1100, units = "px", dpi = 96, compression = "lzw")
  message("Saved: ", out_file)
}
## ============================================================================
## MAIN EXECUTION ORDER
## ============================================================================
# Safe to run all at once. Functions will automatically skip if their 
# specific output files (PNGs, TIFs, or CSVs) already exist.
figure_1()              # Sensor temporal coverage
figure_2()              # Sensor spatial expansion (publication-ready map)
run_kriging()           # Kriging (Skips if pm25_yerevan_ST_stack.nc exists)
run_GEMM()              # GEMM (Skips if deaths_rasters.tif exists)
figure_3()              # Spatial heatmaps
figure_4()              # Time series chart
figure_5()              # Excess mortality quadrant maps
figure_6()              # District monthly PM2.5 line chart
figure_7()              # District mortality heatmap
figure_supp_elderly()   # Elderly demographic map
graphical_abstract()    # Graphical abstract
