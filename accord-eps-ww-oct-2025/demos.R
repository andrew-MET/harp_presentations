# Demos

library(harp)
library(dplyr)

### Batch running
new_point_verif_project("~/harp/demo")


### Plotting

wspd <- read_forecast(
  2025072400,
  "meps",
  "ws10m",
  lead_time     = seq(0, 24, 3),
  members       = c(0, 1, 2, 9, 12),
  file_path     = "~/data/meps/ens",
  file_template = "member_{MBR2}/{fcst_model}_sfc_{LDT2}_{YYYY}{MM}{DD}T{HH}Z.nc",
  file_format_opts = netcdf_opts("met_norway_det"),
  return_data   = TRUE
)

plot(wspd, col = meps_mbr000)

plot(
  filter(pivot_members(wspd), lead_time == 24),
  facet_col = member
)

perts <- select(
  mutate(wspd, across(contains("meps_mbr"), ~.x - meps_mbr000)),
  -meps_mbr000
)

plot(
  filter(pivot_members(perts), lead_time == 0),
  facet_col = member
) +
  scale_fill_gradient2(limits = abs_range)


# Neighbourhood probabilities
pcp3h <- read_forecast(
  2025072400,
  "meps",
  "pcp",
  lead_time     = seq(0, 24, 3),
  members       = c(0, 1, 2, 9, 12),
  file_path     = "~/data/meps/ens",
  file_template = "member_{MBR2}/{fcst_model}_sfc_{LDT2}_{YYYY}{MM}{DD}T{HH}Z.nc",
  file_format_opts = netcdf_opts("met_norway_det"),
  return_data   = TRUE
)

plot(pcp3h, col = meps_mbr000) +
  scale_fill_viridis_c(
    limits = c(0.125, NA),
    trans = "log2",
    oob = censor_low_squish_high,
    na.value = NA,
    option = "G",
    direction = -1
  )

prob8 <- ens_prob(pcp3h, 8)

plot(prob8, col = prob_ge_8) +
  scale_fill_viridis_c(option = "H", limits = c(0.1, 1), na.value = NA)

nbhd_at <- ens_stats(
  mutate(pcp3h, across(contains("meps_mbr"), ~nbhd_smooth(.x, 7, 8))),
  sd = FALSE
)

plot(nbhd_at, col = ens_mean) +
  scale_fill_viridis_c(option = "H", limits = c(0.1, 1), na.value = NA)

nbhd_within <- ens_stats(
  mutate(pcp3h, across(contains("meps_mbr"), ~nbhd_smooth(.x, 7, 8) > 0)),
  sd = FALSE
)

plot(nbhd_within, col = ens_mean) +
  scale_fill_viridis_c(option = "H", limits = c(0.1, 1), na.value = NA)


# Wind vectors
wind <- read_forecast(
  2025072400,
  "meps",
  c("u10m", "v10m"),
  lead_time     = 12,
  file_path     = "~/data/meps/ens",
  file_template = "member_00/{fcst_model}_sfc_{LDT2}_{YYYY}{MM}{DD}T{HH}Z.nc",
  file_format_opts = netcdf_opts("met_norway_det"),
  return_data   = TRUE
) |>
  pivot_parameters_wider(fcst) |>
  mutate(wind_speed = sqrt(u10m ^ 2 + v10m ^ 2))

plot(wind, col = wind_speed) +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  geom_geowindvec(
    aes(u = u10m, v = v10m),
    wind,
    skip = 25
  )

plot(wind, col = wind_speed) +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  geom_geowindvec(
    aes(u = u10m, v = v10m),
    wind,
    skip = 25,
    max_vec = 5
  ) +
  theme(panel.border = element_rect(colour = NA))

# Cross sections
oslo    <- c(10.75, 59.91)
ålesund <- c(6.15, 62.47)

my_params <- add_param_def(
  "clw",
  netcdf = new_netcdf_param("mass_fraction_of_cloud_condensed_water_in_air_ml")
)

my_params <- add_param_def(
  "cli",
  netcdf = new_netcdf_param("mass_fraction_of_cloud_ice_in_air_ml"),
  param_defs = my_params
)


os_ål <- read_forecast(
  2025072400,
  "meps",
  c("T", "Q", "psfc", "sfc_geo", "caf", "clw", "cli", "w"),
  lead_time           = 12,
  file_path           = "~/data/meps/det",
  file_template       = "{fcst_model}_det_2_5km_{YYYY}{MM}{DD}T{HH}Z.nc",
  transformation      = "xsection",
  transformation_opts = xsection_opts(oslo, ålesund),
  param_defs          = my_params,
  vertical_coordinate = "model",
  return_data         = TRUE
) |>
  pivot_parameters_wider(fcst)

ggplot(os_ål, aes(xs = caf)) +
  geom_xs()

ggplot(os_ål, aes(xs = caf)) +
  geom_xs(interpolate = TRUE) +
  scale_fill_viridis_c(option = "F") +
  scale_y_reverse() +
  coord_cartesian(expand = FALSE) +
  theme_bw()

ggplot(os_ål, aes(xs = caf, psfc = psfc)) +
  geom_xs_pressure(interpolate = TRUE) +
  scale_fill_viridis_c(option = "F") +
  scale_y_reverse() +
  coord_cartesian(expand = FALSE) +
  theme_bw()

ggplot(os_ål, aes(xs = caf, psfc = psfc)) +
  geom_xs_pressure(
    interpolate = TRUE,
    scale_pressure = 1 / 100,
    scale_distance = 1 / 1000,
    topo_colour = "grey",
    topo_linewidth = 0.5
  ) +
  scale_fill_viridis_c(option = "F") +
  scale_y_reverse() +
  coord_cartesian(expand = FALSE) +
  theme_bw() +
  labs(
    x        = "Oslo -> Ålesund [km]",
    y        = "Pressure [hPa]",
    fill     = "Cloud\nFraction",
    title    = "Cloud Fraction",
    subtitle = "for a cross section between Oslo and Ålesund"
  )

ggplot(os_ål, aes(xs = caf, psfc = psfc, temp = T, spec_hum = Q, topo = sfc_geo)) +
  geom_xs_height(
    interpolate = TRUE,
    scale_distance = 1 / 1000,
    topo_colour = "grey",
    topo_linewidth = 0.5
  ) +
  scale_fill_viridis_c(option = "F") +
  coord_cartesian(expand = FALSE) +
  theme_bw() +
  labs(
    x        = "Oslo -> Ålesund [km]",
    y        = "Height [m]",
    fill     = "Cloud\nFraction",
    title    = "Cloud Fraction",
    subtitle = "for a cross section between Oslo and Ålesund"
  )

ggplot(os_ål, aes(xs = caf, psfc = psfc, temp = T, spec_hum = Q, topo = sfc_geo)) +
  geom_xs_height(
    interpolate = TRUE,
    scale_distance = 1 / 1000,
    height_lims = c(NA, 12000),
    topo_colour = "grey",
    topo_linewidth = 0.5
  ) +
  scale_fill_viridis_c(option = "F") +
  coord_cartesian(expand = FALSE) +
  theme_bw() +
  labs(
    x        = "Oslo -> Ålesund [km]",
    y        = "Height [m]",
    fill     = "Cloud\nFraction",
    title    = "Cloud Fraction",
    subtitle = "for a cross section between Oslo and Ålesund"
  )

ggplot(os_ål, aes(xs = T)) +
  geom_xs_map() +
  coord_cartesian(expand = FALSE)

ggplot(os_ål, aes(xs = T)) +
  geom_xs_map(
    section_labels = c("Oslo", "Ålesund"),
    section_label_hjust = c(-0.2, 1.1),
    map_fill = "grey"
  ) +
  coord_equal(expand = FALSE) +
  theme_harp_map() +
  theme(panel.background = element_rect(fill = "skyblue"))




ggplot(os_ål, aes(xs = w, psfc = psfc)) +
  geom_xs_pressure(
    interpolate = TRUE,
    scale_pressure = 1 / 100,
    scale_distance = 1 / 1000,
    topo_fill = "grey",
    topo_colour = "grey20",
    topo_linewidth = 0.5
  ) +
  scale_fill_gradient2(limits = abs_range) +
  scale_y_reverse() +
  coord_cartesian(expand = FALSE) +
  theme_bw() +
  labs(
    x        = "Oslo -> Ålesund [km]",
    y        = "Pressure [hPa]",
    fill     = bquote(m.s^{-1}),
    title    = "Veritcal Velocity",
    subtitle = "for a cross section between Oslo and Ålesund"
  )
