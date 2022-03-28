# cannibalized code from non-stationary paper...
library(sdmTMB)
library(sp)
library(ggplot2)
library(ggsidekick)
library(dplyr)
species <- read.csv("survey_data/species_list.csv")
bio <- readRDS("survey_data/wcbts_bio_2019-08-01.rds")
haul <- readRDS("survey_data/wcbts_haul_2019-08-01.rds")
# catch = readRDS("survey_data/wcbts_catch_2019-08-01.rds")
# names(catch) = tolower(names(catch))
names(bio) <- tolower(names(bio))
bio$date_yyyymmdd <- as.character(bio$date_yyyymmdd)
bio$trawl_id <- as.numeric(bio$trawl_id)
bio$sampling_end_hhmmss <- as.character(bio$sampling_end_hhmmss)
bio$sampling_start_hhmmss <- as.character(bio$sampling_start_hhmmss)

dat <- dplyr::left_join(bio, haul)
dat <- dplyr::filter(dat,sex=="F")

# filter species from list
dat <- dplyr::filter(dat, scientific_name %in% species$Scientific.Name)
# convert to UTM
coordinates(dat) <- c("longitude_dd", "latitude_dd")
proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
newproj <- paste("+proj=utm +zone=10 ellps=WGS84")
dat <- spTransform(dat, CRS(newproj))
dat <- as.data.frame(dat)
dat$lon <- dat$longitude_dd/1000
dat$lat <- dat$latitude_dd/1000
#dat$year = as.numeric(substr(dat$date_yyyymmdd,1,4))

# set up df of cells surveyed. all we really care about are epsilon_st by year, so depth can be 0
unique_coords <- strsplit(unique(paste(floor(dat$lon), floor(dat$lat))), " ")
grid_dat <- expand.grid(depth_scaled = 0, depth_scaled2 = 0,
  year = unique(dat$year),
  i = seq(1, length(unique_coords)))
grid_dat$lon <- as.numeric(unlist(lapply(unique_coords, getElement, 1)))[grid_dat$i]
grid_dat$lat <- as.numeric(unlist(lapply(unique_coords, getElement, 2)))[grid_dat$i]

spp = dplyr::group_by(dat, scientific_name) %>%
  dplyr::summarise(common_name = common_name[1])
saveRDS(spp,"spp.rds")
for(i in 1:length(unique(dat$scientific_name))) {

  sub = dplyr::filter(dat, scientific_name==unique(dat$scientific_name)[i])
  spde <- make_mesh(sub, c("lon", "lat"), cutoff = 20)
  sub$depth_scaled = as.numeric(scale(sub$depth_m))
  sub$depth_scaled2 = sub$depth_scaled^2
  sub$log_length = log(sub$length_cm)
  m = sdmTMB(log_length ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = sub,
    time = "year",
    mesh = spde,
    spatial = "on",
    spatiotemporal = "AR1"
  )
  saveRDS(m, paste0("bio_runs/",i,".rds"))

}
