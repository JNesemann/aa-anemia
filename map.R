# ALTO AMAZONAS ANEMIA MAPS
# JOHN NESEMANN

# packages
library(osmdata)
library(sf)
library(ggmap)
library(tidyverse)
library(maptools)
library(ggspatial)

#### registering google key ####
register_google(key = "AIzaSyDIUlzMxLwTzuYkbqkUNeMuPu8QOam8DbU")

#### importing data ###
setwd("~/Desktop/DATA/aa-anemia")

# anemia
data <- haven::read_dta("anemia_alto_amazonas_abril.dta") %>%
  select(unique_id:hh_numb, pplinhh:gps_alt, dbs_hg, dbs_hg_na, anemia, anemia_level_sa, anemia_orig, gps_lat, gps_long) %>%
  dplyr::rename(anemia.f=anemia_orig) %>%
  group_by(community) %>%
  summarise(n_tested=sum(!is.na(age)),
            n_any_anemia=sum(anemia.f %in% c("mild","moderate","severe")),
            n_mild=sum(anemia.f=="mild"),
            # grouping moderate and severe together due to the low number of severe cases
            n_modsev=sum(anemia.f %in% c("moderate","severe")),
            p_any=n_any_anemia/n_tested,
            p_mild=n_mild/n_tested,
            p_modsev=n_modsev/n_tested,
            # going for medians to get village midpoint as often there were a few houses located very far away
            lat=median(gps_lat, na.rm = T),
            lon=median(gps_long, na.rm = T)) %>%
  # converting to sf
  st_as_sf(., coords = c("lon","lat"), crs = 4326)

# Health posts identity and location
hp <- read_csv("Data_EESS_Nacional_28-02-2020.csv") %>%
  rename(sector=Institución, uniqueid=`Código Único`, name=`Nombre del establecimiento`, classification=`Clasificación`,
         type=Tipo, department=Departamento, province=Provincia, district=Distrito, disa_code=`Código DISA`, red_code=`Código Red`,
         microred_code=`Código Microrred`, disa=DISA, red=Red, microrred=Microrred, eu_code=`Código UE`, executive_unit=`Unidad Ejecutora`,
         hp_level=Categoria, op_state=Estado, condition=Condición, long=NORTE, lat=ESTE, beds=CAMAS) %>%
  filter(department=="LORETO" & province=="ALTO AMAZONAS" & op_state=="ACTIVADO") %>%
  select(sector:district, disa_code:hp_level, condition, lat, long, beds) %>%
  # selecting only active medical facilities in alto amazonas
  filter(!is.na(lat)) %>% filter(!is.na(long)) %>%
  select(-sector, -type, -(disa_code:executive_unit), -condition, -beds) %>%
  # filtering to just get health posts (not dental, specialist, or ophthalmic/optometrist centers)
  filter(!(classification %in% c("CENTRO ODONTOLOGICO","CENTROS MEDICOS ESPECIALIZADOS","CENTROS OPTICOS",
                                 "ATENCION PRE HOSPITALARIA,SERVICIO DE TRASLADO DE PACIENTES,ATENCION DOMICILIARIA",
                                 "CONSULTORIOS MEDICOS Y DE OTROS PROFESIONALES DE LA SALUD"))) %>%
  # xtabs(data=hp,~classification, addNA=T)
  # changing to SF object
  st_as_sf(., coords = c("long","lat"), crs = 4326)

#### quering OSM for map data ###
# building a query for trunks, i.e., the road leading into Yurimaguas
q_hwy <- getbb("alto amazonas, peru") %>%
  opq() %>%
  add_osm_feature("highway","trunk")
# sending the query to the server
hwy <- osmdata_sf(q_hwy)
# smaller roads and footpaths
q_primary <- getbb("alto amazonas, peru") %>%
  opq() %>%
  add_osm_feature("highway",c("primary","secondary","tertiary","unclassified","residential","road","path"))
primary <- osmdata_sf(q_primary)
# waterways
q_waterway <- getbb("alto amazonas, peru") %>%
  opq() %>%
  add_osm_feature("waterway", available_tags("waterway"))
waterway <- osmdata_sf(q_waterway)
# natural bodies of water
q_natwater <- getbb("alto amazonas, peru") %>%
  opq() %>%
  add_osm_feature("natural", "water")
natwater <- osmdata_sf(q_natwater)
# admin boundaries
q_admin <- getbb("alto amazonas, peru") %>%
  opq() %>%
  add_osm_feature("boundary", "administrative")
admin <- osmdata_sf(q_admin)
# all villages 
q_villages <- getbb("alto amazonas, peru") %>%
  opq() %>%
  add_osm_feature("place", "village")
villages <- osmdata_sf(q_villages)
# get towns and hamlets too as they are larger than villages but smaller than cities
q_town <- getbb("alto amazonas, peru") %>%
  opq() %>%
  add_osm_feature("place", "town")
town <- osmdata_sf(q_town)
# hamlets
q_hamlet <- getbb("alto amazonas, peru") %>%
  opq() %>%
  add_osm_feature("place", "hamlet")
hamlet <- osmdata_sf(q_hamlet)

# getting basemap
aa <- get_map(location = c(lon = -76.2, lat = -5.3), zoom = 9, maptype = "terrain-background")
ggmap(aa)

#### mapping ####

# mild or worse anemia
map_mild <- ggmap(aa) +
  # roads and highways
  geom_sf(data = primary$osm_lines, inherit.aes = F, color = "#ffffd4", linetype = 3, size = 0.75, show.legend = F) + # alpha = 0.6
  geom_sf(data = hwy$osm_lines, inherit.aes = F, color = "yellow", linetype = 1, size = 1.25, show.legend = F) +
  # waterways and natural water
  geom_sf(data = waterway$osm_lines, inherit.aes = F, color = "#9ecae1", fill = "#9ecae1") +
  geom_sf(data = waterway$osm_polygons, inherit.aes = F, color = "#9ecae1", fill = "#9ecae1") +
  geom_sf(data = natwater$osm_lines, inherit.aes = F, color = "#9ecae1", fill = "#9ecae1") +
  geom_sf(data = natwater$osm_polygons, inherit.aes = F, color = "#9ecae1", fill = "#9ecae1") +
  geom_sf(data = natwater$osm_multipolygons, inherit.aes = F, color = "#9ecae1", fill = "#9ecae1") +
  # admin boundaries
  geom_sf(data = admin$osm_lines, inherit.aes = F, color = "#636363", linetype=5, show.legend = F) + #alpha = 0.7
  # other villages, towns, and hamlets
  geom_sf(data=town$osm_points, aes(color = "black", alpha = 0.5), inherit.aes = F) + 
  geom_sf(data=villages$osm_points, aes(color = "black", alpha = 0.5), inherit.aes = F) +
  geom_sf(data=hamlet$osm_points, aes(color = "black", alpha = 0.5), inherit.aes = F) +
  # health posts
  # geom_sf(data=hp, aes(color = "red", alpha = 0.5), inherit.aes = F) + 
  # visited villages with prevalence of anemia
  geom_sf(data=data, aes(fill = p_any),
          shape = 21, size = 3,
          inherit.aes = F) +
  # changing labels and making it look nice
  xlab("Longitude") + ylab("Latitude") +
  scale_fill_gradient(low = "white", 
                      high = "red",
                      labels = scales::percent,
                      limits = c(0, 0.75)) +
  labs(color = "Legend", fill = "Prevalence", alpha = NULL, title = "A) Mild-or-worse anemia") +
  scale_color_manual(name = "Legend",
                     labels = c("Other villages"), # ,"Health posts"
                     values = c("black")) + # ,"red"
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  annotation_scale(location = "tr")
map_mild

# moderate to severe anemia
map_mod <- ggmap(aa) +
  # roads and highways
  geom_sf(data = primary$osm_lines, inherit.aes = F, color = "#ffffd4", linetype = 3, size = 0.75, show.legend = F) + # alpha = 0.6
  geom_sf(data = hwy$osm_lines, inherit.aes = F, color = "yellow", linetype = 1, size = 1.25, show.legend = F) +
  # waterways and natural water
  geom_sf(data = waterway$osm_lines, inherit.aes = F, color = "#9ecae1", fill = "#9ecae1") +
  geom_sf(data = waterway$osm_polygons, inherit.aes = F, color = "#9ecae1", fill = "#9ecae1") +
  geom_sf(data = natwater$osm_lines, inherit.aes = F, color = "#9ecae1", fill = "#9ecae1") +
  geom_sf(data = natwater$osm_polygons, inherit.aes = F, color = "#9ecae1", fill = "#9ecae1") +
  geom_sf(data = natwater$osm_multipolygons, inherit.aes = F, color = "#9ecae1", fill = "#9ecae1") +
  # admin boundaries
  geom_sf(data = admin$osm_lines, inherit.aes = F, color = "#636363", linetype=5, show.legend = F) + #alpha = 0.7
  # other villages, towns, and hamlets
  geom_sf(data=town$osm_points, aes(color = "black", alpha = 0.5), inherit.aes = F) + 
  geom_sf(data=villages$osm_points, aes(color = "black", alpha = 0.5), inherit.aes = F) +
  geom_sf(data=hamlet$osm_points, aes(color = "black", alpha = 0.5), inherit.aes = F) +
  # health posts
  # geom_sf(data=hp, aes(color = "darkorchid", alpha = 0.5), inherit.aes = F) + 
  # visited villages with prevalence of anemia
  geom_sf(data=data, aes(fill = p_modsev),
          shape = 21, size = 3.5,
          inherit.aes = F) +
  # changing labels and making it look nice
  xlab("Longitude") + ylab(NULL) +
  scale_fill_gradient(low = "white", 
                      high = "red",
                      labels = scales::percent,
                      limits = c(0, 0.75)) +
  labs(color = "Legend", fill = "Prevalence", alpha = NULL, title = "B) Moderate-or-worse anemia") +
  scale_color_manual(name = "Legend",
                     labels = c("Other villages"), # "Health posts"
                     values = c("black")) + # ,"darkorchid"
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  annotation_scale(location = "tr")
map_mod

# putting it all together

map <- ggpubr::ggarrange(map_mild, map_mod,
                         nrow = 1,
                         common.legend = T,
                         legend = "right",
                         heights = c(1.1, 1),
                         widths = c(1.1, 1))
# map

# saving it
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/R projects/aa-anemia")
ggsave("figures/map_mildmod.eps", map, 
       device = cairo_ps, # this allows me to save alpha (semi-transparancy) on my mac
       width = 12, height = 7)
























