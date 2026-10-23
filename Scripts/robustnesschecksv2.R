# =============================================================================
# ROBUSTNESS CHECKS FOR PM2.5-ATTRIBUTABLE MORTALITY IN YEREVAN
# =============================================================================
# Caching Logic: If an output already exists, the check is skipped.
# =============================================================================

library(readr); library(dplyr); library(sf); library(gstat); library(spacetime)
library(raster); library(sp); library(terra); library(ggplot2); library(lubridate)
library(pbapply); library(stringr)

# Ensure paths from main script are loaded, or define them:
if (!exists("PROJ_DIR")) {
  PROJ_DIR  <- "C:/Users/mprhc/OneDrive/Documents/GitHub/ArmeniaPollutionAnalysis"
  DATA_DIR  <- "C:/Users/mprhc/Downloads"
  READINGS_CSV <- file.path(DATA_DIR, "readings_clean_2023_2026.csv")
  SENSORS_CSV  <- file.path(DATA_DIR, "sensors_filtered_2023_2026.csv")
  BOUNDARY_LVL4 <- file.path(DATA_DIR, "boundary.gpkg")
  POP_TIF <- file.path(PROJ_DIR, "WorldPopData/arm_age_data/pop.tif")
  GEMM_XLSX <- file.path(PROJ_DIR, "GEMMcoefficients.xlsx")
  PRED_DIR <- file.path(PROJ_DIR, "Interpolation outputs_ST/pred_layers")
}

ROB_DIR <- file.path(PROJ_DIR, "Robustness")
dir.create(ROB_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper to run a lightweight GEMM calculation on an annual mean raster
calculate_gemm_annual <- function(pm25_annual_rast, pop_raster, coeff_df) {
  pm25_resampled <- resample(pm25_annual_rast, pop_raster[[1]], method = "bilinear")
  pm25_pop_values <- values(pm25_resampled)
  n_pop_cells <- ncell(pop_raster[[1]])
  deaths_by_pop_cell <- rep(0, n_pop_cells)
  
  age_map <- c("0-9"="0-25","1-10"="0-25","5-14"="0-25","10-19"="0-25","15-24"="0-25",
               "20-29"="25-30","25-34"="30-35","30-39"="35-40","35-44"="40-45","40-49"="45-50",
               "45-54"="50-55","50-59"="55-60","55-64"="60-65","60-69"="65-70","60-70"="65-70",
               "70-79"="75-80","75-84"="80+","80-89"="80+","85-94"="80+","90-99"="80+")
  causes_list <- unique(coeff_df$Cause)
  coef_ages <- unique(coeff_df$AgeGroup)
  coef_matrix <- array(NA, dim = c(length(causes_list), length(coef_ages), 5),
                       dimnames = list(causes_list, coef_ages, c("theta","theta_sd","alpha","mu","nu")))
  for (i in seq_along(causes_list)) {
    for (j in seq_along(coef_ages)) {
      row <- coeff_df[coeff_df$Cause == causes_list[i] & coeff_df$AgeGroup == coef_ages[j], ]
      if (nrow(row) > 0) coef_matrix[i, j, ] <- c(row$Theta[1], row$Theta_SD[1], row$Alpha[1], row$Mu[1], row$Nu[1])
    }
  }
  
  calc_rr <- function(pm25, theta, alpha, mu, nu, cf = 2.4) {
    z <- pm25 - cf; z[z < 0] <- 0
    omega <- 1 / (1 + exp(-(z - mu) / nu))
    exp(theta * log((z / alpha) + 1) * omega)
  }
  
  death_rates <- coeff_df[, c("AgeGroup","DeathRate")]
  
  for (ag in names(pop_raster)) {
    if (!ag %in% names(age_map)) next
    pop_values <- values(pop_raster[[ag]])
    coef_age <- unname(age_map[ag])
    age_dr <- death_rates$DeathRate[death_rates$AgeGroup == coef_age]
    if (length(age_dr) == 0 || is.na(age_dr[1])) next
    age_dr <- age_dr[1]
    idx <- which(coef_ages == coef_age)
    if (length(idx) == 0 || is.na(coef_matrix[1, idx, 1])) next
    
    theta <- coef_matrix[1, idx, 1]; alpha <- coef_matrix[1, idx, 3]
    mu <- coef_matrix[1, idx, 4]; nu <- coef_matrix[1, idx, 5]
    
    rr_values <- calc_rr(pm25_pop_values, theta, alpha, mu, nu)
    paf_values <- (rr_values - 1) / rr_values; paf_values[is.na(paf_values)] <- 0
    baseline_death_values <- pop_values * age_dr
    death_values <- paf_values * baseline_death_values
    death_values[is.na(death_values)] <- 0
    deaths_by_pop_cell <- deaths_by_pop_cell + death_values
  }
  return(sum(deaths_by_pop_cell, na.rm = TRUE))
}
# =============================================================================
# CHECK 1: LEAVE-ONE-PROVIDER-OUT (Drop Yerevan GIS) — OPTIMIZED v2
# =============================================================================
check_lopo <- function(grid_m = 500, n_cores = max(1, parallel::detectCores() - 1),
                       chunk_size = 30) {
  message("\n=== Robustness Check 1: Leave-Yerevan-GIS-Out (optimized v2) ===")
  
  out_tif <- file.path(ROB_DIR, "lopo_no_gis_annual_pm25.tif")
  if (file.exists(out_tif)) {
    message("Output already exists. Skipping. Delete it to re-run.")
    return(invisible(NULL))
  }
  
  data <- read_csv(READINGS_CSV, show_col_types = FALSE)
  data$timestamp <- as.POSIXct(gsub("\u00A0", " ", data$timestamp), format = "%Y-%m-%d %H:%M:%S")
  sensors <- read_csv(SENSORS_CSV, show_col_types = FALSE)
  
  sensors.c <- dplyr::filter(sensors, provider %in% c("clarity","purple_air") & !is_suspicious)
  lon_col <- intersect(c("longitude","long","lon"), names(sensors.c))[1]
  lat_col <- intersect(c("latitude","lat"), names(sensors.c))[1]
  
  sensors_coords <- dplyr::select(sensors.c, id, !!sym(lon_col), !!sym(lat_col))
  sensors_coords <- dplyr::rename(sensors_coords, longitude = !!sym(lon_col), latitude = !!sym(lat_col))
  
  data.v <- data |>
    dplyr::filter(sensor_id %in% sensors.c$id) |>
    left_join(sensors_coords, by = c("sensor_id"="id")) |>
    dplyr::filter(!is.na(timestamp), !is.na(pm2.5), !is.na(longitude), !is.na(latitude))
  
  start_date <- as.POSIXct("2024-06-01 00:00:00", tz="UTC")
  end_date   <- as.POSIXct("2025-06-01 00:00:00", tz="UTC")
  data.v <- data.v[data.v$timestamp >= start_date & data.v$timestamp < end_date, ]
  
  # --- Compute UTM coords ONCE and pre-join outside the loop ---
  sensor_locs <- dplyr::select(data.v, sensor_id, longitude, latitude) |> distinct()
  sf_sensors_utm <- st_transform(st_as_sf(sensor_locs, coords = c("longitude","latitude"), crs = 4326), 32638)
  utm_coords <- st_coordinates(sf_sensors_utm)
  sensor_locs$x_utm <- utm_coords[, "X"]; sensor_locs$y_utm <- utm_coords[, "Y"]
  
  data.v <- dplyr::left_join(data.v, sensor_locs, by = c("sensor_id","longitude","latitude"))
  data.v$timestamp <- as.POSIXct(data.v$timestamp, origin = "1970-01-01", tz = "UTC")
  data.v$date <- as.Date(data.v$timestamp)
  
  # --- Aggregate hourly -> DAILY means per sensor ---
  daily_data <- data.v |>
    group_by(sensor_id, x_utm, y_utm, date) |>
    summarise(pm2.5 = mean(pm2.5, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(!is.na(pm2.5))
  
  # FIX #1: Convert to plain data.frame so SpatialPoints() accepts column subsets
  daily_data <- as.data.frame(daily_data)
  
  cat(sprintf("Aggregated to %d daily sensor-observations.\n", nrow(daily_data)))
  
  # --- Boundary + grid ---
  yerevan_boundary <- st_read(BOUNDARY_LVL4, layer = "boundary-polygon-lvl4", quiet = TRUE)
  yerevan_proj <- st_transform(yerevan_boundary, 32638)
  grid_centers <- st_make_grid(yerevan_proj, cellsize = grid_m, what = "centers")
  grid_inside  <- st_intersection(grid_centers, yerevan_proj)
  if (length(grid_inside) == 0) grid_inside <- grid_centers
  grid_sp <- as(grid_inside, "Spatial"); proj4string(grid_sp) <- CRS("EPSG:32638")
  yerevan_proj_sp <- as(yerevan_proj, "Spatial")
  
  cat(sprintf("Grid: %d points at %dm resolution.\n", length(grid_sp), grid_m))
  
  # --- Fit a SPATIAL variogram ---
  set.seed(42)
  dvar <- daily_data[sample(min(nrow(daily_data), 1500)), ]
  sp_sample <- SpatialPointsDataFrame(
    SpatialPoints(as.matrix(dvar[, c("x_utm","y_utm")]), proj4string = CRS("EPSG:32638")),
    data.frame(pm2.5 = dvar$pm2.5)
  )
  
  sill <- var(dvar$pm2.5, na.rm=TRUE); if (!is.finite(sill) || sill == 0) sill <- 100
  vgm_start <- vgm(sill, "Sph", 2000, sill * 0.1)
  emp_vgm <- variogram(pm2.5 ~ 1, sp_sample, cutoff = 20000, width = 1000)
  vgm_fit <- tryCatch(fit.variogram(emp_vgm, vgm_start), error = function(e) vgm_start)
  cat(sprintf("Spatial variogram: model=%s, psill=%.1f, range=%.0f\n",
              vgm_fit$model[2], vgm_fit$psill[2], vgm_fit$range[2]))
  
  # --- Process days in CHUNKS so progress bar actually updates ---
  unique_days <- sort(unique(daily_data$date))
  n_days <- length(unique_days)
  chunks <- split(seq_along(unique_days), ceiling(seq_along(unique_days) / chunk_size))
  
  cat(sprintf("Processing %d days on %d cores (in %d chunks of %d)...\n",
              n_days, n_cores, length(chunks), chunk_size))
  
  library(parallel)
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("daily_data", "grid_sp", "vgm_fit"), envir = environment())
  clusterEvalQ(cl, {
    library(sp); library(gstat); library(raster)
    NULL
  })
  
  daily_means_list <- list()
  pb <- txtProgressBar(0, length(chunks), style = 3)
  
  for (ci in seq_along(chunks)) {
    idx_chunk <- chunks[[ci]]
    days_chunk <- unique_days[idx_chunk]
    
    # Returns a list of rasters (or NULL) named by date
    chunk_results <- parLapply(cl, days_chunk, function(d) {
      day_data <- daily_data[daily_data$date == d, ]
      if (nrow(day_data) < 3) return(NULL)
      
      # FIX: as.matrix() on the coords subset — tibble-safe
      sp_train <- SpatialPointsDataFrame(
        SpatialPoints(as.matrix(day_data[, c("x_utm","y_utm")]),
                      proj4string = CRS("EPSG:32638")),
        data.frame(pm2.5 = day_data$pm2.5)
      )
      
      pred <- tryCatch(
        krige(pm2.5 ~ 1, sp_train, grid_sp, model = vgm_fit, debug.level = 0),
        error = function(e) NULL
      )
      if (is.null(pred)) return(NULL)
      
      r <- raster(SpatialPixelsDataFrame(
        points = SpatialPoints(coordinates(grid_sp), proj4string = CRS("EPSG:32638")),
        data = data.frame(pred = pred@data$var1.pred)
      ))
      r
    })
    
    names(chunk_results) <- as.character(days_chunk)
    daily_means_list <- c(daily_means_list, chunk_results)
    setTxtProgressBar(pb, ci)
  }
  close(pb)
  stopCluster(cl)
  
  daily_means_list <- daily_means_list[!sapply(daily_means_list, is.null)]
  cat(sprintf("\nSuccessfully kriged %d / %d days.\n",
              length(daily_means_list), n_days))
  
  annual_mean_rast <- mean(stack(daily_means_list), na.rm = TRUE)
  annual_mean_rast <- mask(annual_mean_rast, yerevan_proj_sp)
  
  # --- GEMM mortality calculation ---
  pop_raster_raw <- rast(POP_TIF); crs(pop_raster_raw) <- "EPSG:4326"
  pop_raster_raw <- project(pop_raster_raw, "EPSG:32638")
  first_r <- rast(list.files(PRED_DIR, "*.tif$", full.names = TRUE)[1])
  pop_raster <- crop(pop_raster_raw, ext(first_r))
  
  coeff_all <- readxl::read_excel(GEMM_XLSX)
  colnames(coeff_all) <- c("Cause","AgeGroup","Theta","Theta_SD","Alpha","Mu","Nu","Deaths23","Pop23","DeathRate")
  coeff_df <- coeff_all[complete.cases(coeff_all), ] |> dplyr::filter(Cause == "NCD+LRI")
  
  deaths_lopo <- calculate_gemm_annual(rast(annual_mean_rast), pop_raster, coeff_df)
  
  baseline_csv <- file.path(PROJ_DIR, "yerevan_mortality_all_periods.csv")
  if (file.exists(baseline_csv)) {
    baseline_df <- read.csv(baseline_csv)
    baseline_deaths <- baseline_df$TotalDeaths[baseline_df$Period == "June 2024 - June 2025"]
    diff_pct <- ((deaths_lopo - baseline_deaths) / baseline_deaths) * 100
    cat(sprintf("Baseline (June 24-25): %.0f deaths | No-GIS: %.0f deaths (%.1f%% change)\n",
                baseline_deaths, deaths_lopo, diff_pct))
  } else {
    cat(sprintf("No-GIS model estimated %.0f deaths. (Run main script to get baseline comparison)\n", deaths_lopo))
  }
  
  writeRaster(rast(annual_mean_rast), out_tif, overwrite = TRUE)
  cat("Saved:", out_tif, "\n")
  
  # --- Publication-ready map PNG (the raw GeoTIFF does not render in Word) ---
  out_png <- file.path(ROB_DIR, "lopo_no_gis_annual_pm25_map.png")
  map_df <- as.data.frame(rast(annual_mean_rast), xy = TRUE, na.rm = TRUE)
  names(map_df)[3] <- "pm25"
  p_map <- ggplot(map_df, aes(x = x, y = y, fill = pm25)) +
    geom_raster() +
    geom_sf(data = yerevan_proj, fill = NA, color = "black", linewidth = 0.4,
            inherit.aes = FALSE) +
    scale_fill_viridis_c(option = "plasma", name = "PM2.5 (µg/m³)") +
    coord_sf(expand = FALSE) +
    theme_minimal(base_size = 11) +
    labs(x = NULL, y = NULL) +
    theme(axis.text = element_blank(), panel.grid = element_blank())
  ggsave(out_png, p_map, width = 7, height = 6.2, dpi = 300)
  cat("Saved map:", out_png, "\n")
}# =============================================================================
# CHECK 2: CONSTANT SENSOR LOCATIONS (Network Stability)
# =============================================================================
check_constant_sensors <- function() {
  message("\n=== Robustness Check 2: Constant Sensor Locations ===")
  
  out_file <- file.path(ROB_DIR, "check_constant_sensors_trend.png")
  out_csv <- file.path(ROB_DIR, "check_constant_sensors_data.csv")
  
  if (file.exists(out_file) && file.exists(out_csv)) {
    message("Outputs already exist. Skipping.")
    return(invisible(NULL))
  }
  
  data <- read_csv(READINGS_CSV, show_col_types = FALSE)
  data$timestamp <- as.POSIXct(gsub("\u00A0", " ", data$timestamp), format = "%Y-%m-%d %H:%M:%S")
  sensors <- read_csv(SENSORS_CSV, show_col_types = FALSE)
  
  sensors.c <- dplyr::filter(sensors, provider %in% c("clarity","purple_air","yerevan_gis") & !is_suspicious)
  data.c <- dplyr::filter(data, sensor_id %in% sensors.c$id)
  
  data.c$year <- year(data.c$timestamp)
  sensor_years <- data.c |> group_by(sensor_id) |> summarise(n_years = n_distinct(year))
  stable_sensors <- sensor_years$sensor_id[sensor_years$n_years == 3]
  
  cat(sprintf("Found %d sensors active in 2023, 2024, and 2025.\n", length(stable_sensors)))
  
  if(length(stable_sensors) < 5) {
    message("Not enough stable sensors for robust kriging. Try relaxing the criteria.")
    return(NULL)
  }
  
  sensors_sub <- dplyr::select(sensors.c, id, longitude, latitude)
  
  data.v <- data.c |> 
    dplyr::filter(sensor_id %in% stable_sensors) |>
    left_join(sensors_sub, by = c("sensor_id"="id")) |>
    dplyr::filter(!is.na(timestamp), !is.na(pm2.5), !is.na(longitude), !is.na(latitude))
  
  citywide_stable <- data.v |> group_by(timestamp) |> summarise(pm25 = mean(pm2.5, na.rm=TRUE))
  citywide_all    <- data.c |> group_by(timestamp) |> summarise(pm25_all = mean(pm2.5, na.rm=TRUE))
  
  df_compare <- left_join(citywide_stable, citywide_all, by = "timestamp")
  df_compare <- dplyr::filter(df_compare, !is.na(pm25), !is.na(pm25_all))
  df_compare$month <- floor_date(df_compare$timestamp, "month")
  
  monthly_compare <- df_compare |> group_by(month) |>
    summarise(Stable_Sensors = mean(pm25, na.rm=TRUE), All_Sensors = mean(pm25_all, na.rm=TRUE))
  
  write.csv(monthly_compare, out_csv, row.names = FALSE)
  
  monthly_compare$year <- year(monthly_compare$month)
  annual_means <- monthly_compare |> 
    group_by(year) |> 
    summarise(Stable_Annual_Mean = mean(Stable_Sensors, na.rm=TRUE),
              All_Annual_Mean = mean(All_Sensors, na.rm=TRUE))
  
  if(2024 %in% annual_means$year && 2025 %in% annual_means$year) {
    s24 <- annual_means$Stable_Annual_Mean[annual_means$year == 2024]
    s25 <- annual_means$Stable_Annual_Mean[annual_means$year == 2025]
    pct_increase <- ((s25 - s24) / s24) * 100
    cat(sprintf("Trend Check: Constant sensors increased %.1f%% from 2024 to 2025.\n", pct_increase))
  }
  
  p <- ggplot(monthly_compare, aes(x = month)) +
    geom_line(aes(y = Stable_Sensors, color = "Constant Sensors Only"), linewidth = 1) +
    geom_line(aes(y = All_Sensors,    color = "Full Network"),          linewidth = 1) +
    scale_color_manual(values = c("Constant Sensors Only" = "blue",
                                  "Full Network"          = "red")) +
    theme_minimal() + 
    labs(title = "Robustness: PM2.5 Trend (Constant vs Full Network)", 
         x = "Date", y = "PM2.5 (µg/m³)", color = "") +
    theme(legend.position = "bottom")
  
  ggsave(out_file, p, width = 10, height = 5, dpi = 300)
  cat("Saved trend comparison plot and data.\n")
}

# =============================================================================
# CHECK 3: SPATIO-TEMPORAL CROSS-VALIDATION
# =============================================================================
check_cv <- function() {
  message("\n=== Robustness Check 3: 5-Fold Spatio-Temporal CV ===")
  
  out_png <- file.path(ROB_DIR, "check_cv_scatterplot.png")
  out_csv <- file.path(ROB_DIR, "cv_predictions.csv")
  
  if (file.exists(out_png) && file.exists(out_csv)) {
    message("Outputs already exist. Skipping.")
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
    dplyr::filter(!is.na(timestamp), !is.na(pm2.5), !is.na(longitude), !is.na(latitude))
  
  sensor_locs <- dplyr::select(data.v, sensor_id, longitude, latitude) |> distinct()
  sf_sensors_utm <- st_transform(st_as_sf(sensor_locs, coords = c("longitude","latitude"), crs = 4326), 32638)
  utm_coords <- st_coordinates(sf_sensors_utm)
  sensor_locs$x_utm <- utm_coords[, "X"]; sensor_locs$y_utm <- utm_coords[, "Y"]
  
  set.seed(42)
  sample_ts <- sample(unique(data.v$timestamp), 100)
  data.v <- data.v[data.v$timestamp %in% sample_ts, ]
  data.v <- dplyr::left_join(data.v, sensor_locs, by = c("sensor_id","longitude","latitude"))
  # Defensive: ensure timestamp stays POSIXct (left_join preserves it, but be safe)
  data.v$timestamp <- as.POSIXct(data.v$timestamp, origin = "1970-01-01", tz = "UTC")
  # FIX: convert to data.frame so SpatialPoints() accepts subsets
  data.v <- as.data.frame(data.v)
  
  unique_sensors <- unique(sensor_locs$sensor_id)
  folds <- sample(rep(1:5, length.out = length(unique_sensors)))
  names(folds) <- unique_sensors
  
  data.v$fold <- folds[as.character(data.v$sensor_id)]
  
  cv_results <- data.frame()
  pb <- txtProgressBar(0, 5, style=3)
  
  for(k in 1:5) {
    train <- data.v[data.v$fold != k, ]
    test  <- data.v[data.v$fold == k, ]
    
    if(nrow(train) == 0 || nrow(test) == 0) next
    
    sill <- var(train$pm2.5, na.rm=TRUE); if (!is.finite(sill) || sill == 0) sill <- 100
    # FIX: stAni = 100 (was 0 in my earlier paste — invalid)
    vgm_fit <- vgmST("sumMetric",
                     space = vgm(sill*0.8, "Sph", 2000, sill*0.2),
                     time  = vgm(sill*0.8, "Exp", 12,   sill*0.2),
                     joint = vgm(sill*0.3, "Sph", 1000, 0),
                     stAni = 100)
    
    for(ts in unique(test$timestamp)) {
      train_ts <- train[train$timestamp == ts, ]
      test_ts  <- test[test$timestamp == ts, ]
      
      if(nrow(train_ts) < 3 || nrow(test_ts) == 0) next
      
      # Ensure ts is a time-based class for STF()
      ts <- as.POSIXct(ts, origin = "1970-01-01", tz = "UTC")
      
      # Use as.matrix to be tibble-safe (data.v is already data.frame, but defensive)
      sp_train <- SpatialPoints(as.matrix(train_ts[, c("x_utm","y_utm")]),
                                proj4string = CRS("EPSG:32638"))
      stData <- STIDF(sp_train, train_ts$timestamp, data.frame(pm2.5 = train_ts$pm2.5))
      
      sp_test <- SpatialPoints(as.matrix(test_ts[, c("x_utm","y_utm")]),
                               proj4string = CRS("EPSG:32638"))
      stGrid <- STF(sp = sp_test, time = c(ts, ts + 3600))
      
      pred <- tryCatch(krigeST(pm2.5 ~ 1, stData, stGrid, model = vgm_fit,
                               computeVar = FALSE, debug.level = 0),
                       error = function(e) NULL)
      
      if(!is.null(pred)) {
        cv_results <- rbind(cv_results, data.frame(
          observed  = test_ts$pm2.5,
          predicted = pred@data$var1.pred[1:nrow(test_ts)]
        ))
      }
    }
    setTxtProgressBar(pb, k)
  }
  close(pb)
  
  cv_results <- cv_results[complete.cases(cv_results), ]
  cv_results$residual <- cv_results$observed - cv_results$predicted
  
  rmse <- sqrt(mean(cv_results$residual^2))
  mae  <- mean(abs(cv_results$residual))
  ss_res <- sum(cv_results$residual^2)
  ss_tot <- sum((cv_results$observed - mean(cv_results$observed))^2)
  r2 <- 1 - (ss_res / ss_tot)
  
  cat(sprintf("ST-CV Results: R2 = %.3f | RMSE = %.2f µg/m³ | MAE = %.2f µg/m³\n", r2, rmse, mae))
  
  p <- ggplot(cv_results, aes(x = observed, y = predicted)) +
    geom_point(alpha = 0.2, color = "darkblue") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
    theme_minimal() + labs(title = "Robustness: Spatio-Temporal Cross-Validation (5-Fold)",
                           subtitle = sprintf("R² = %.2f | RMSE = %.2f | MAE = %.2f", r2, rmse, mae),
                           x = "Observed PM2.5", y = "Predicted PM2.5")
  ggsave(out_png, p, width = 6, height = 6, dpi = 300)
  write.csv(cv_results, out_csv, row.names = FALSE)
}

# =============================================================================
# CHECK 4: VARIOGRAM SENSITIVITY (Manual vs Auto-Fit) — FIXED
# =============================================================================
check_variogram <- function() {
  message("\n=== Robustness Check 4: Variogram Fitting Sensitivity ===")
  
  out_pdf <- file.path(ROB_DIR, "check_variogram_comparison.pdf")
  if (file.exists(out_pdf)) {
    message("Output already exists. Skipping.")
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
    dplyr::filter(!is.na(timestamp), !is.na(pm2.5), !is.na(longitude), !is.na(latitude))
  
  sensor_locs <- dplyr::select(data.v, sensor_id, longitude, latitude) |> distinct()
  sf_sensors_utm <- st_transform(st_as_sf(sensor_locs, coords = c("longitude","latitude"), crs = 4326), 32638)
  utm_coords <- st_coordinates(sf_sensors_utm)
  sensor_locs$x_utm <- utm_coords[, "X"]; sensor_locs$y_utm <- utm_coords[, "Y"]
  
  data.v <- dplyr::left_join(data.v, sensor_locs, by = c("sensor_id","longitude","latitude"))
  data.v$timestamp <- as.POSIXct(data.v$timestamp, origin = "1970-01-01", tz = "UTC")
  # FIX: convert to data.frame so SpatialPoints() accepts subsets
  data.v <- as.data.frame(data.v)
  
  # FIX #2: larger sample for a more stable empirical ST variogram
  set.seed(42)
  sample_n <- min(nrow(data.v), 10000)
  dvar <- data.v[sample(nrow(data.v), sample_n), ]
  cat(sprintf("Using %d observations for variogram estimation.\n", sample_n))
  
  sp_win <- SpatialPoints(as.matrix(dvar[, c("x_utm","y_utm")]), proj4string = CRS("EPSG:32638"))
  stData <- STIDF(sp_win, dvar$timestamp, data.frame(pm2.5 = dvar$pm2.5))
  
  sill <- var(dvar$pm2.5, na.rm=TRUE); if (!is.finite(sill) || sill == 0) sill <- 100
  
  # FIX #3: more conservative starting ranges (5000m, 24h, 2000m joint)
  vgm_manual <- vgmST("sumMetric",
                      space = vgm(sill*0.8, "Sph", 5000, sill*0.2),
                      time  = vgm(sill*0.8, "Exp", 24,  sill*0.2),
                      joint = vgm(sill*0.3, "Sph", 2000, 0),
                      stAni = 100)
  
  cat("Computing empirical ST variogram...\n")
  emp_vgm <- variogramST(pm2.5 ~ 1, stData, tunit = "hours", cutoff = 20000, width = 1000)
  
  cat("Fitting StVariogram (method=1, bounded)...\n")
  vgm_auto <- tryCatch({
    fit.StVariogram(emp_vgm, vgm_manual,
                    fit.method = 1,
                    stAni = 100,
                    lower = c(0.01, 100, 0.01, 1, 0.01, 100, 1),
                    upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf))
  }, error = function(e) {
    message("Auto-fit failed (method=1): ", conditionMessage(e))
    message("Falling back to method=7 with relaxed bounds...")
    tryCatch(
      fit.StVariogram(emp_vgm, vgm_manual, fit.method = 7, stAni = 100,
                      lower = c(0.01, 100, 0.01, 1, 0.01, 100, 1),
                      upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf)),
      error = function(e2) {
        message("Auto-fit failed entirely. Using manual fit only.")
        vgm_manual
      }
    )
  })
  
  # Compute MSEs safely (the manual fit may not have optim.output)
  mse_manual <- tryCatch(attr(vgm_manual, "optim.output")$value, error = function(e) NA)
  mse_auto   <- tryCatch(attr(vgm_auto, "MSE"), error = function(e) {
    # If vgm_auto fell back to vgm_manual, attr may be missing — compute MSE manually
    tryCatch(
      mean((emp_vgm$np - predict(vgm_auto, emp_vgm$dist, emp_vgm$timelag))^2, na.rm = TRUE),
      error = function(e2) NA
    )
  })
  
  cat(sprintf("Manual Fit MSE: %s\n", ifelse(is.na(mse_manual), "N/A", sprintf("%.2f", mse_manual))))
  cat(sprintf("Auto Fit MSE:   %s\n", ifelse(is.na(mse_auto),   "N/A", sprintf("%.2f", mse_auto))))
  
  pdf(out_pdf, width = 10, height = 5)
  plot(emp_vgm, vgm_manual, main = "Manual Variogram Fit")
  plot(emp_vgm, vgm_auto, main = "Automated Variogram Fit")
  dev.off()
  
  cat("Saved variogram comparison plots to Robustness directory.\n")
}
# =============================================================================
# EXECUTE CHECKS
# =============================================================================
check_lopo()             # Leave-One-Provider-Out (Takes ~10-15 mins)
check_constant_sensors() # Constant sensor locations (Takes ~1 min)
check_cv()               # 5-Fold CV (Takes ~5-10 mins)
check_variogram()        # Variogram sensitivity (Takes ~2 mins)
