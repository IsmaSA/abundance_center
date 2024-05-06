
## "Abundance-center" theory test:


suppressMessages({
  library(dplyr, quiet = TRUE, warn.conflicts = FALSE)
  library(reshape, quiet = TRUE, warn.conflicts = FALSE)
  library(tidyr)  
  library(stringr)
  library(raster)
  library(rnaturalearth)
  library(xlsx)
  library(readxl)
  library(writexl)
  library(rgdal)
  library(devtools)
  Sys.setenv(LANGUAGE = "en")
})


# Global data
setwd("C:/Users/Propietario/Desktop/Escritorio/full databse/DATA")
df <- read_xlsx("Global_dataset.xlsx", sheet = "Sheet1")
df1<- df %>% filter(Alien =="Y")


# Check species, time series and occ

data_s<- df1 %>% group_by(taxon) %>% summarise(
  time_series = n_distinct(site_id),
  occ = n(), 
  Countries = n_distinct(country)) #%>% filter(time_series >4)

data_s <- df1 %>% group_by(country) %>% summarise(species = n_distinct(taxon))


# Map to see the occ

europe_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(continent == "Europe" | name == "Turkey" | 
           name %in% c("Morocco", "Algeria", "Tunisia", "Libya", "Egypt"))

europe_map$color <- ifelse(europe_map$name %in% c("Morocco", "Algeria", "Tunisia", "Libya", "Egypt", "Russia", "Turkey"),
                           "grey80", "Wheat")


ggplot() +
  geom_sf(data = europe_map, aes(fill = color))+   scale_fill_identity()+ 
  coord_sf() + theme_bw()+ xlim(-11,30) + ylim(36,60)+
  geom_point(data = df1, aes(x = Longitude_X, y = Latitude_Y), 
   size = 1, shape = 21, color = "black", fill = "red", stroke = 0.2)


# Fill the gaps ----

time_series <- unique(df1$site_id)
# Replace 'corbicula "fluminalis"' with 'corbicula fluminalis' in the dataframe
df1$taxon <- gsub('corbicula "fluminalis"', 'corbicula fluminalis2', df1$taxon)

df1<- df1 %>% filter(!taxon=="dikerogammarus haemobaphes\villosus")
df1 <- df1 %>% filter(!grepl("dikerogammarus haemobaphes|villosus", taxon))

for(i in 1:length(time_series)){
  
  subset_df <- df1 %>% filter(site_id == time_series[i]) 
  subset_df2 <- df %>% filter(site_id == time_series[i]) 
  
  species <- unique(subset_df$taxon)
  
  for(j in 1:length(species)){
    
    # Subset the data frame to include only the current species
    subset_species_df <- subset_df %>% filter(taxon == species[j])
    
    # check if the species has less than 2  values
    if(sum(!is.na(subset_species_df$abundance[subset_species_df$taxon == species[j]])) < 2) {
      next
    }
    
    
    subset_species_df <- subset_species_df %>%  group_by(site_id)  %>% ## Add missing years e.g. 2010
      complete(year = seq(min(year), max(year), period = 1))   
    
    sampling_years <- unique(subset_df2$year)
    
    subset_species_df2 <- subset_species_df %>% filter(year %in% sampling_years)
    
    subset_species_df2$taxon <- species[j]
    subset_species_df2<- subset_species_df2[,c(1:4)]
    subset_species_df2[is.na(subset_species_df2)] = 0
    
    subset_species_df3 <- subset_species_df2 %>%  group_by(site_id)  %>% ## Add missing years e.g. 2010
      complete(year = seq(min(year), max(year))) 
    
    subset_species_df3$taxon <- species[j]
    
    subset_species_df3$abundance<- as.numeric(subset_species_df3$abundance)
    Paride <- complete(mice(subset_species_df3, m = 5, maxit = 50, method = "pmm", seed = 123))
    write_xlsx(Paride,paste0(time_series[i], species[j],".xlsx"))
    cat("\n The loop for-->",time_series[i])
  }
  
}

files <- list.files()

data <- do.call(rbind, lapply(files, read_excel))
data$abundance[is.na(data$abundance)] <- 0
write_xlsx(data,"data_filled.xlsx")


# Mann-Kendall trend test modified: 

df1
names(df1)
df1.1 <- data


#MK modified
My.mmkh=function (x, ci = 0.95) 
{
  x = x
  z = NULL
  z0 = NULL
  pval = NULL
  pval0 = NULL
  S = 0
  Tau = NULL
  essf = NULL
  ci = ci
  if (is.vector(x) == FALSE) {
    stop("Input data must be a vector")
  }
  if (any(is.finite(x) == FALSE)) {
    x <- x[-c(which(is.finite(x) == FALSE))]
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }
  n <- length(x)
  V <- rep(NA, n * (n - 1)/2)
  k = 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      k = k + 1
      V[k] = (x[j] - x[i])/(j - i)
    }
  }
  slp <- median(V, na.rm = TRUE)
  t = 1:length(x)
  xn <- (x[1:n]) - ((slp) * (t))
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      S = S + sign(x[j] - x[i])
    }
  }
  ro <- acf(rank(xn), lag.max = (n - 1), plot = FALSE)$acf[-1]
  sig <- qnorm((1 + ci)/2)/sqrt(n)
  rof <- rep(NA, length(ro))
  for (i in 1:(length(ro))) {
    if (ro[i] > sig || ro[i] < -sig) {
      rof[i] <- ro[i]
    }
    else {
      rof[i] = 0
    }
  }
  cte <- 2/(n * (n - 1) * (n - 2))
  ess = 0
  for (i in 1:(n - 1)) {
    ess = ess + (n - i) * (n - i - 1) * (n - i - 2) * rof[i]
  }
  essf = 1 + ess * cte
  var.S = n * (n - 1) * (2 * n + 5) * (1/18)
  if (length(unique(x)) < n) {
    aux <- unique(x)
    for (i in 1:length(aux)) {
      tie <- length(which(x == aux[i]))
      if (tie > 1) {
        var.S = var.S - tie * (tie - 1) * (2 * tie + 
                                             5) * (1/18)
      }
    }
  }
  VS = var.S * essf
  if (S == 0) {
    z = 0
    z0 = 0
  }
  if (S > 0) {
    z = (S - 1)/sqrt(VS)
    z0 = (S - 1)/sqrt(var.S)
  }
  else {
    z = (S + 1)/sqrt(VS)
    z0 = (S + 1)/sqrt(var.S)
  }
  pval = 2 * pnorm(-abs(z))
  pval0 = 2 * pnorm(-abs(z0))
  Tau = S/(0.5 * n * (n - 1))
  return(c("Corrected Zc" = z, "new P-value" = pval, "N/N*" = essf, 
           "Original Z" = z0, "old P.value" = pval0, "Tau" = Tau, 
           "Sen s slope" = slp, "old.variance" = var.S, "new.variance" = VS, "S statistic" = S, "n" = n))
}

res <- data.frame()
for (sp in unique(df1$taxon)) {
  df2 <- df1.1 %>% filter(taxon ==sp)
  
  time_series <- unique(df2$site_id)
  
  for (ts in time_series) {
    df3 <- df2 %>% filter(site_id ==ts)
    
    if(nrow(df3) < 3) {
      cat("Records", sp,":", ts, "/n")
    } else {
      
      xy.list <- split(df3$abundance, df3$site_id) 
      xy.list <- xy.list[lengths(xy.list) >= 3]  

      if(var(xy.list[[1]], na.rm = TRUE) == 0) {
        next
      }
      MK <-as.data.frame(do.call(rbind,lapply(xy.list[1:1],function(x)unlist(My.mmkh(x))))) #96 are time series,
      
        if (!is.null(MK) && nrow(MK) > 0) {
          MK$site_id <- ts
          MK$taxon <- sp
         S <- MK$`S statistic`
         Var <- MK$new.variance
         P <- MK$`new P-value`
         Tau <- MK$Tau
          row <- data.frame("Species" = sp, "site_id" = ts, "Trend" = S, 
                            "Variance" = Var, "P-value" = P, 
                            "N" = MK$n, "Tau" = Tau)
          res <- rbind(res, row)
        }
      }
    }
    cat("Processed species: ", sp, "\n")
  }
}

write_xlsx(res,"S_trends.xlsx")


#  EMILI approach------

spn <- unique(res$Species)
n<- spn[1]

keys <- tibble(Species = character(), key = integer())
for (n in spn) {
  tryCatch({
    if( n=="corbicula fluminalis2"){
      next
    }
    if( n=="rhyacophila dorsalis ssp."){
      next
    }  
    sp <- occ_data(scientificName = n, hasCoordinate = TRUE, occurrenceStatus = "PRESENT",
                   limit = 5)
    sp <- sp[["data"]]
    
    if (nrow(sp) > 0) {
      key <- sp$taxonKey[1]
    } else {
      key <- NA
    }
    
    keys <- rbind(keys, tibble(Species = n, key = key))
    cat(n, "Extracted\n")
  })
  
}

unique(keys$Species)

sp <- unique(keys$key)
s<- sp[1]
dois <- tibble(Species = character(), DOI = character())
for(s in sp){
  if(s == "corbicula fluminalis2" || s == "rhyacophila dorsalis ssp.") {
    next
  }
  x = occ_download(
    pred_in("basisOfRecord", c("OCCURRENCE", "HUMAN_OBSERVATION", "OBSERVATION")),
    pred_in("taxonKey", s),
    pred_and(pred_gte("year", 1970), pred_lte("year", 2020)),
    pred("hasCoordinate", TRUE),
    pred_lte("coordinateUncertaintyInMeters", 9999), 
    user = "ismaelsoto",
    pwd = "Ismaputas123.",
    email = "isma-sa@hotmail.com"
  )
  status <- occ_download_meta(x)$status
  
  while(status != "SUCCEEDED") {
    Sys.sleep(15) 
    status <- occ_download_meta(x)$status
  }
  
  # Once the download is complete, get the DOI
  z <- occ_download_meta(x)
  doi2 <- z$doi 
  dois <- rbind(dois, tibble(Species = s, DOI = doi2))
  cat(s, "downloaded\n")
}

write_xlsx(dois, "dois.xlsx")

# Check data ---

files <- list.files(pattern = ".zip")
data_gbif <- data.frame()
for (file in files) {
  dir_name <- gsub("\\.zip$", "", file)
  if (dir.exists(dir_name)) {
    unlink(dir_name, recursive = TRUE)
  }
  dir.create(dir_name)
  unzip(file, exdir = dir_name)
  unzipped_files <- list.files(dir_name, full.names = TRUE)
  target_file <- file.path(dir_name, "occurrence.txt")
  if (file.exists(target_file)) {
    data <- read_delim(target_file, delim = "\t")
    
    data1 <- data[,c("species","taxonKey","decimalLatitude", "decimalLongitude", "occurrenceStatus","countryCode")]
    data_gbif <- rbind(data_gbif, data1)
  } else {
    warning("File  with occurrence OUT: ", file)
  }
}

nrow(data_gbif) #502,551

names(data_gbif)
# records from our data set
names(df1)
df11<- df1 %>% filter(!taxon %in% c("dikerogammarus haemobaphes/villosus","rhyacophila dorsalis ssp.",
                             "corbicula \"fluminalis\"","corophium sowinskyi","physella heterostropha","dikerogammarus bispinosus")) 

gbif<- unique(data_gbif$species) %>% tolower()
gbif[grepl("cavimana", gbif, ignore.case = TRUE)]

df11$taxon[df11$taxon=="dugesia tigrina"] <- "girardia tigrina"
df11$taxon[df11$taxon=="corophium curvispinum"] <- "chelicorophium curvispinum"
df11$taxon[df11$taxon=="corophium robustum"] <- "chelicorophium robustum"
df11$taxon[df11$taxon=="echinogammarus trichiatus"] <- "chaetogammarus trichiatus"
df11$taxon[df11$taxon=="balanus improvisus"] <- "amphibalanus improvisus"
df11$taxon[df11$taxon=="orconectes limosus"] <- "faxonius limosus"
df11$taxon[df11$taxon=="echinogammarus ischnus"] <- "chaetogammarus ischnus"
df11$taxon[df11$taxon=="orchestia cavimana"] <- "cryptorchestia cavimana"
df11$taxon[df11$taxon=="ferrissia clessiniana"] <- "pettancylus clessinianus"
df11$taxon[df11$taxon=="psychomyia fragilis"] <- "metalype fragilis"


setdiff(df11$taxon, gbif)


# Clean GBIF data ----

data <- data_gbif
data$countryCode <- as.character(data$countryCode)
data$countryCode <- countrycode(data$countryCode[!is.na(data$countryCode)], origin = 'iso2c', destination = 'iso3c')

data <- data %>% dplyr::select(decimalLongitude, decimalLatitude, countryCode, species)
data$decimalLatitude <- as.numeric(data$decimalLatitude)
data$decimalLongitude <- as.numeric(data$decimalLongitude)

clean_data <- data.frame()
for(dat in unique(data$species)){
  da <- data[data$species==dat,]
  flags<- CoordinateCleaner::clean_coordinates(x = da,
                                               lon = "decimalLongitude",
                                               lat = "decimalLatitude",
                                               countries = "countryCode",
                                               species = "species",
                                               tests = c("capitals", "centroids", "equal","gbif",
                                                         "zeros")) 
  dat_cl <- da[flags$.summary,]
  
  dat_cl<- dat_cl %>% as.data.frame() 
  clean_data <- rbind(clean_data, dat_cl) 
}

nrow(clean_data) # 410, 890

# Now I will combine ts and gbif data

df11 <- df11[!duplicated(df11[c("site_id", "taxon")]), ]

df12<- df11[,c(3,8,9)]
colnames(df12) <- c("species","decimalLatitude","decimalLongitude")

clean_data1 <- clean_data[,c(4,2,1)]

names(df12)
names(clean_data1)  

df12$source <- "time_series"
clean_data1$source <- "gbif"

data <- rbind(df12, clean_data1)

nrow(data) # 414219



# Place GBIF records in GRID ----
bio <- raster::getData('worldclim', var='bio', res=5)

grid_extent <- extent(-180, 180, -60, 90)  

#Grid
grid <- raster(grid_extent, res = 0.08333333)
grid_count <- raster(grid)
res<- data.frame()
for (z in unique(data$species)) {
  data0 <- data[data$species==z,]
  
  occ_points <- sp::SpatialPoints(data.frame(lon = data0$decimalLongitude, 
                                             lat = data0$decimalLatitude))
  grid_count <- rasterize(occ_points, grid_count, fun = 'count')
  freq_table <- freq(grid_count, useNA='no') %>%
    as.data.frame() %>%
    filter(value > 0) 
  
  # Remove overlap
  valid_cells <- sum(!is.na(grid_count[]))
  ids <- rep(NA, ncell(grid_count))
  ids[!is.na(grid_count[])] <- 1:valid_cells
  grid_count[] <- ids
  raster_cells<- grid_count@data@values %>% as.data.frame()
  cell_ids <- extract(grid_count, occ_points)
  point_counts_per_cell <- as.data.frame(table(cell_ids))
  colnames(point_counts_per_cell) <- c("Cell_ID", "Count")
  points_df <- as.data.frame(occ_points)
  points_df$Cell_ID <- extract(grid_count, occ_points, cellnumbers=TRUE)[,1]

  points_df$Cell_ID <- as.character(points_df$Cell_ID)
  point_counts_per_cell$Cell_ID <- as.character(point_counts_per_cell$Cell_ID)
  points_with_cell_ids <- points_df %>%
    left_join(point_counts_per_cell, by = "Cell_ID")
  
  #select just one point:
  unique_points_per_cell <- points_with_cell_ids %>%
    group_by(Cell_ID) %>%
    dplyr::slice(1) %>%
    ungroup()

  unique_points_per_cell$pecies <- z
  res <- rbind(res, unique_points_per_cell)
}
write_xlsx(res, "cell_ID_presence.xlsx")

# Remove the species that are in < 10 cells

cells <- read_xlsx("cell_ID_presence.xlsx")

cells2<- cells %>% group_by(pecies) %>% summarise(count=n()) %>% filter(count <10)
remove_cells <- unique(cells2$pecies)
#"Chelicorophium robustum"   "Dreissena rostriformis"   
# "amphibalanus improvisus"   "baetis liebenauae"        
# "barbronia weberi"          "chaetogammarus trichiatus"
# "corbicula fluminalis"      "dreissena rostriformis"   
# "echinogammarus berilloni"  "gammarus roeselii"        
# "gyraulus parvus"           "hemimysis anomala"        
# "katamysis warpachowskyi"   "menetus dilatatus"        
# "metalype fragilis"         "musculium transversum"    
# "obesogammarus obesus"      "pectinatella magnifica"   
# "pettancylus clessinianus"  "physa fontinalis"         
# "physella gyrina"           "piscicola haranti"        
# "pontogammarus robustoides" "potamothrix bavaricus"    
# "potamothrix heuscheri"     "potamothrix vejdovskyi"   
# "psammoryctides moravicus"  "pseudosuccinea columella" 
# "psychomyia pusilla"        "rhithropanopeus harrisii" 
# "rhyacophila munda"  

presence <- cells %>% filter(!pecies %in% remove_cells) 


# Extra Predictors ----
memory.limit(size = 900000)
gc()
setwd("C:/Users/Propietario/Desktop/abundance-center/Extra predictors")
bio
# res = 0.08333333

# Elevation
elevation <- raster("elevation_1KMmd_GMTEDmd.tif")
elevation <- projectRaster(elevation, bio, method = "bilinear") #re-proyect
elevation <- setNames(elevation, "Elevation")
elevation

# Slope
Slope <- raster("slope_1KMmd_GMTEDmd.tif")
Slope <- projectRaster(Slope, bio, method = "bilinear") #re-proyect
Slope <- setNames(Slope, "Slope")
Slope

# bio

bioc <- bio[[c(1,4,10,11,15,16,17)]]
var_stack <- terra::rast(bioc)


predictors1<- stack(var_stack, elevation,Slope,coastline, flow_accum, strahler)
predictors1<- stack(var_stack, elevation,Slope)


predictors1_spat <- rast(predictors1)
coords_vect <- vect(coordinates, geom=c("decimalLongitude", "decimalLatitude"), crs=crs(predictors1_spat))
val <- extract(predictors1_spat, coords_vect)



# Flow acc --> https://www.earthenv.org/streams (Stream length and flow accumulation)
library(ncdf4)

flow_accum <- terra::rast("C:/Users/Propietario/Desktop/abundance-center/Extra predictors/flow_acc.nc")
flow_accum <- resample(flow_accum, rast_template, method="near")

coordinates <- presence
coordinates <- read_xlsx("points_presence.xlsx")
accum <- terra::extract(flow_accum, coordinates[,c(2,3)])
accum


extract_with_buffer <- function(raster, points, max_radius = 500, step = 100) {
  # Ensure points are in a SpatVector
  if (!inherits(points, "SpatVector")) {
    points <- vect(data.frame(lon = points$decimalLongitude, lat = points$decimalLatitude), crs = crs(raster))
  }
  
  # Initialize the extraction
  extracted_values <- extract(raster, points, buffer = 0, fun = mean, na.rm = TRUE)
  
  # Initialize buffer radius
  radius <- step
  
  while (any(is.na(extracted_values)) && radius <= max_radius) {
    # Extract with current radius
    buffered_values <- extract(raster, points, buffer = radius, fun = mean, na.rm = TRUE)
    
    # Indices where original extractions are NA
    na_indices <- which(is.na(extracted_values))
    
    # Check structure of buffered_values to ensure it is a vector
    if (!is.vector(buffered_values)) {
      stop("Buffered values extraction did not return a vector as expected.")
    }
    
    # Replace NA values in the original extraction
    extracted_values[na_indices] <- buffered_values[na_indices]
    
    # Increase the radius
    radius <- radius + step
  }
  
  return(extracted_values)
}
coordinates_vect <- vect(coordinates, geom=c("decimalLongitude", "decimalLatitude"), crs=crs(flow_accum))
extracted_values <- extract_with_buffer(flow_accum, coordinates[,c(2,3)])
extracted_values <- extract_with_buffer(flow_accum, coordinates_vect)


# Distance to the Sea (https://invasiber.org/GarciaBerthou/docs/papers/Cano%E2%80%90Barbacil_GEB_2022_withSI.pdf)
coastline <- st_read("ne_50m_coastline.shp")
coordinates
coordinates_sf <- st_as_sf(coordinates, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, agr = "constant")
coordinates_sf <- st_as_sf(coordinates, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

utm_crs <- "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs"  
coordinates_utm <- st_transform(coordinates_sf, crs = utm_crs)
coastline_utm <- st_transform(coastline, crs = utm_crs)

nearest_index <- st_nearest_feature(coordinates_utm, coastline_utm)
distances <- st_distance(coordinates_utm, coastline_utm[nearest_index, ], by_element = TRUE)
distances_km <- units::set_units(distances, "km")

coordinates$distanceSea_km = distances_km

# Strahler order

strahler <- st_read("HydroRIVERS_v10.shp")

nearest_indices <- st_nearest_feature(coordinates_sf, strahler)
nearest_rivers <- strahler[nearest_indices,]

coordinates$Strahler = nearest_rivers$ORD_STRA


coordinates1 = cbind(coordinates,val)
write_xlsx(coordinates1, "presence_predictors.xlsx")



# SDMs ----
library(terra)
library(data.table)
library(sf)
library(sp)

setwd("C:/Users/Propietario/Desktop/abundance-center")
cells <- read_xlsx("cell_ID_presence.xlsx")

cells2<- cells %>% group_by(pecies) %>% summarise(count=n()) %>% filter(count <10)
remove_cells <- unique(cells2$pecies)

presence <- cells %>% filter(!pecies %in% remove_cells) 

spn <- unique(presence$pecies)
sp <- spn[1]
sp <- "caspiobdella fadejewi"

buffer_width <- 200000  # in meters


#care with this 
plot(var_stack)
ext(var_stack)
elevation_spat <- rast(elevation)
slope_spat <- rast(Slope)
var_stack <- rast("bio_ele_slope.tif")
rast_template <- var_stack[[2]]

pseudo_absences <- data.frame()

for(spn in sp){
  presence_data <- presence[presence$pecies==sp, ]
  
# pseudo-absence
  presences_sv <- vect(presence_data, geom = c("lon", "lat"), crs = "epsg:4326")
  #buff <- st_union(st_buffer(presences_sf, dist = buffer_width))
  buff <- terra::aggregate(terra::buffer(presences_sv, width = buffer_width))
  #maps::map("world", col = "grey")  
  #plot(buff, border = "blue", add = TRUE)
  #plot(presences_sv, pch = ".", col = "salmon", add = TRUE)
  
  rast_buff <- terra::mask(rast_template, buff)
  rast_buff_nopres <- rast(raster::mask(raster(rast_buff), as(presences_sv, "Spatial"), inverse = TRUE))
  
  n_pseudoabs <- nrow(presence_data)
  pseudoabsences <- vect(sampleRandom(raster(rast_buff_nopres), size = n_pseudoabs, cells = TRUE, xy = TRUE, sp = TRUE))
  
  cat("pseudoabsences for", sp,":", n_pseudoabs  ,"\n")
  
  pseudoabsences_ext <- terra::extract(var_stack, pseudoabsences)
  pseudoabsences <- data.frame(pseudoabsences[ , 1:3], pseudoabsences_ext[ , -grep("ID", names(pseudoabsences_ext))])
  pseudoabsences <- setDT(data.frame(presence = rep(0, nrow(pseudoabsences)), pseudoabsences[ , c(2, 3, 1, 4:ncol(pseudoabsences))]))
  names(pseudoabsences)[grep("cell", names(pseudoabsences))] <- "cells"  
  
  #Now with the other predictors
  abs <- pseudoabsences[,c(2,3,4)]
  
  #DistanceSea
  abs_sf <- st_as_sf(abs, coords = c("x", "y"), crs = 4326)
  coastline_utm <- st_transform(coastline, crs = utm_crs)
  nearest_index <- st_nearest_feature(abs_sf, coastline_utm)
  distances <- st_distance(abs_sf, coastline_utm[nearest_index, ], by_element = TRUE)
  distances_km <- units::set_units(distances, "km")
  
  pseudoabsences$distanceSea_km = distances_km
  
  #accumulation (2 layer)
  #accum <- terra::extract(flow_accum, abs[,c(1,2)], method="simple")
  #accum
  #pseudoabsences$Flow_1 = accum$`flow_accumulation_variable=_1`
  #pseudoabsences$Flow_2 = accum$`flow_accumulation_variable=_2`
  

  #StrahlerOrder --> need to do it separately
  #nearest_indices <- st_nearest_feature(abs_sf, strahler)
  #nearest_rivers <- strahler[nearest_indices,]
  #pseudoabsences$Strahler = nearest_rivers$ORD_STRA
  
  pseudoabsences$species <-sp
  
  pseudo_absences <- rbind(pseudo_absences,pseudoabsences )
  cat("Completed for:", sp, "\n")
  
    }


writeRaster(var_stack, "bio_ele_slope.tif")




strahler <- st_read("HydroRIVERS_v10.shp")

nearest_indices <- st_nearest_feature(abs_sf, strahler)
nearest_rivers <- strahler[nearest_indices,]

coordinates$Strahler = nearest_rivers$ORD_STRA




# distance to water ----

river= st_read("HydroRIVERS_v10.shp")
lake = st_read("HydroLAKES_polys_v10.shp")

river_clean <- st_make_valid(river)
lake_clean <- st_make_valid(lake)


if (any(!river_valid)) {
  print("There are invalid river geometries")
}
if (any(!lake_valid)) {
  print("There are invalid lake geometries")
}

comb <- st_union(st_geometry(river), st_geometry(lake))
combined_sf <- st_sf(geometry = comb)

st_write(combined_sf, "path/to/your/combined_shapefile.shp")



river_simplified <- st_simplify(river_clean, preserveTopology = TRUE)
lake_simplified <- st_simplify(lake_clean, preserveTopology = TRUE)

# Try union again
comb <- st_union(st_geometry(river_simplified), st_geometry(lake_simplified))





# SDM -----
presence <- cells %>% filter(!pecies %in% remove_cells) 
var_stack <- rast("bio_ele_slope.tif")
var_stack <- brick(var_stack)

spn <- unique(presence$pecies)
sp <- spn[1]

#Folders directory
folder_a <- "Raster_models"  # Folder to save raster files
folder_b <- "SDM_models"  # Folder to save plot PDFs

# Create the folders if they don't exist
dir.create(folder_a, showWarnings = FALSE)
dir.create(folder_b, showWarnings = FALSE)

for (sp in spn) {
  presence1<- presence %>% filter (pecies== sp)
  presence2<-presence1[,c(1,2)]
  
  presence2$species<- 1
  coordinates(presence2) <- c('lon','lat')
  
  cat("SDM running =",sp, "\n")
  sdm_data <- sdm::sdmData(species~.,presence2,
                           predictors= var_stack,
                           bg=bg_data) 

  
  m <- sdm(species~., sdm_data, methods=c('brt','gam','rf', "glm"), 
           replication='cv', cv.folds=10, n=10,
           parallelSetting=list(ncore=10,method='parallel'), 
           modelSettings = list(
             gam = list(family = binomial, method = "REML", select= TRUE, k = 5),  #
             rf = list(ntree = 500, mtry = 3, importance = TRUE),                #
             gbm = list(distribution = "bernoulli", n.trees = 1000, interaction.depth = 3, shrinkage = 0.01,
                        n.minobsinnode = 10)) )
  
  SDM_mod <- file.path(folder_b, paste0("Model_", sp, "sdm")) 
  writeRaster(m, filename = SDM_mod, format = "sdm", overwrite=TRUE)
  
  write.sdm(m,paste0("Model_",sp,".sdm"))
  cat("Model done -->",sp,"\n")
  
  en1 <- ensemble(m, newdata = var_stack, filename='en.img', overwrite =TRUE,
                  setting=list(method='weighted',stat='AUC',opt=2))
  
  prediction_file <- file.path(folder_a, paste0("Model_", sp, ".tif")) 
  writeRaster(en1, filename = prediction_file, format = "GTiff", overwrite=TRUE)
  
  cat("Model finished-->",sp,"\n")
  
  }





# bg points
n_points <- 500
random_points <- data.frame(
  longitude = runif(n_points, -180, 180),
  latitude = runif(n_points, -90, 90)
)
random_points_sf <- st_as_sf(random_points, coords = c("longitude", "latitude"), crs = 4326)


bg_data = extract( var_stack, random_points_sf)
bg_data = bg_data %>% data.frame() %>% drop_na()
bg_data$species<- 0
bg_data$lon<- runif(nrow(bg_data), -180, 180)
bg_data$lat<-runif(nrow(bg_data), -180, 180)

nrow(bg_data)
nrow(presence)

# note, all the names must match also species and coodinates in bg_data
sdm_data <- sdm::sdmData(species~.,presence,
                         predictors= var_stack,
                         bg=bg_data)

sdm_data

getmethodNames()
m <- sdm(species~., sdm_data, methods=c('brt','gam','rf', "glm"), 
         replication='cv', cv.folds=2,
         n=2,
         parallelSetting=list(ncore=10,method='parallel'), 
         modelSettings = list(
           gam = list(family = binomial, method = "REML", select= TRUE, k = 5),  #
           rf = list(ntree = 500, mtry = 3, importance = TRUE),                #
           gbm = list(distribution = "bernoulli", n.trees = 1000, interaction.depth = 3, shrinkage = 0.01,
               n.minobsinnode = 10))
         )


en1 <- ensemble(m, newdata = var_stack, filename='en.img',
       setting=list(method='weighted',stat='AUC',opt=2))


plot(en1)
points(presence, add=T, col="red")
bg_data_sf <- st_as_sf(bg_data, coords = c("lon", "lat"), crs = 4326)
points(bg_data_sf, add=T, col="purple")












e<-extent(c(-12,32,30,72)) #look at the map and check the coordinates
bioc <- crop(bio,e) # xmin, xmax, ymin, ymax
plot(bioc[[1]])

sp<- unique_points_per_cell %>% dplyr::select(lon,lat)
sp$species <- 1 
head(sp)
coordinates(sp) <- c('lon','lat')

# VIF
ex <- raster::extract(bioc,sp)
v <- vifstep(ex)
bioc <- exclude(bioc, v)


# SDM ------
library(sdm)
installAll()

d <- sdmData(species~., sp, predictors= bioc, bg = list(method='gRandom',n=500))
d

getmethodNames()

m <- sdm(species~., d, methods=c('glm','gam','rf','maxent'), replication=c('boot'),
  n=2, parallelSetting=list(ncore=8,method='parallel'))

m

gui(m)


en1 <- ensemble(m, bioc, filename='en.img', setting=list(method='weighted',stat='auc',opt=2))

plot(en1)

# representation  -----

# centroid of the distribution of GBIF data
centroid <- unique_points_per_cell %>%
  summarise(lon = mean(lon, na.rm = TRUE), lat = mean(lat, na.rm = TRUE))%>% as.data.frame()
coordinates(centroid) <- c('lon','lat')



df_ex<- df1 %>% filter(taxon=="dugesia tigrina")
df_ex = df_ex[!duplicated(df_ex$site_id), ]
points = df_ex[,c(1,2,8,9)]
coordinates(points) <- c('Longitude_X','Latitude_Y')

coordinates(points_with_cell_ids) <- c('lon','lat')

plot(en1)
points(centroid, pch=21, bg="purple", cex= 1.5)
points(points, pch=21, bg="red", cex= 0.8)
points(points_with_cell_ids, pch=21, bg="blue", cex= 0.5)


time_series <- extract(en1, points)
time_series= time_series %>% as.data.frame()
colnames(time_series) = "Suitability"

points$Suitability = time_series$Suitability


# distance form point to centroid:
crs(centroid) <- CRS("+init=epsg:4326")

coordinates(df_ex) <- ~Longitude_X+Latitude_Y
proj4string(df_ex) <- CRS("+init=epsg:4326")

distances <- spDistsN1(as.matrix(coordinates(df_ex)), coordinates(centroid), longlat = TRUE)
df_ex$distance_to_centroid <- distances


points = df_ex
points$trend <- runif(nrow(points), min=-5, max=5)
points<- points %>% as.data.frame()
points$Suitability = time_series$Suitability

hist(points$Suitability)
glm_mk <- glm(trend ~ Suitability + distance_to_centroid, family=gaussian(), data = points)
summary(glm_mk)

par(mfrow = c(2, 2))
plot(glm_mk)

points$predicted_trend <- predict(glm_mk, type = "response")

ggplot(points, aes(x = predicted_trend, y = trend)) +
  geom_point(aes(color = "Observed"), alpha = 0.6) +
  geom_smooth(method = lm, se = FALSE, color = "blue", linetype = "dashed") +
  labs(x = "Predicted Trend", y = "Observed Trend", title = "Observed vs. Predicted Trend") +
  theme_minimal() +
  scale_color_manual("", labels = c("Observed"), values = c("red")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "darkgreen")


points$residuals <- points$trend - points$predicted_trend

ggplot(points, aes(x = predicted_trend, y = residuals)) +
  geom_point(aes(color = "Residuals"), alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "darkgreen") +
  labs(x = "Predicted Trend", y = "Residuals", title = "Residuals of the Model") +
  theme_minimal() +
  scale_color_manual("", labels = c("Residuals"), values = c("orange"))

install.packages("sjPlot")
library(sjPlot)

p1<- plot_model(glm_mk, type = "eff", terms = "Suitability")
p2<- plot_model(glm_mk, type = "eff", terms = "distance_to_centroid")

p1+p2

