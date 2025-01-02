library(sf)
library(tidyverse)
library(ggsflabel)
library(ggspatial)
library(ggOceanMaps)
library(patchwork)

setwd("~/Dropbox/COMSA - Serosurveillance/Data/Saki/")

## Read in COMSA cluster spatial data
dat_COMSA <- st_read("comsa.zambezia_with_clusters/comsa.zambezia_with_clusters.shp")

## Read in urban-rural designation
dat_urban_rural <- read.csv("residence.csv") %>% rename(Residence = Resid)

## Pull in shapefile for Mozambique
## From: https://data.humdata.org/dataset/cod-ab-moz
Mozambique <- st_read("moz_adm_20190607b_shp/moz_admbnda_adm1_ine_20190607.shp")

## Extract Zambezia province only
Mozambique %>%
  filter(ADM1_PT == "Zambezia") -> Zambezia

## Extract Zambezia's districts
st_read("moz_adm_20190607b_shp/moz_admbnda_adm2_ine_20190607.shp") %>%
  filter(ADM1_PT=="Zambezia") -> Zambezia_districts

## Fancier map, zoomed out
basemap(data=Mozambique, land.col="white", bathymetry=FALSE) +
  geom_sf(data=Mozambique) +
  geom_sf(data=Zambezia, size=4, fill="gold") +
  #geom_sf(data=dat_COMSA %>% filter(sero==1), fill="black") +
  theme_minimal() +
 # annotation_scale(location = "br", 
                   #width_hint = 0.25, 
                   #text_cex = 1) +
  #annotation_north_arrow(location = "br", 
                         #which_north = "true", 
                         #pad_x = unit(0.15, "in"),
                         #pad_y = unit(0.3, "in"),
                         #style = north_arrow_fancy_orienteering) +
  xlab("") + ylab("") +
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.ontop = FALSE,
        panel.grid = element_blank()) -> plot_map_Moz

## Plot just the 30 DBS clusters
dat_COMSA %>%
  filter(sero==1) %>%
  group_by(cluster) %>%
  summarize(geometry = st_union(geometry)) %>% 
  ungroup() %>%
  ggplot() +
  geom_sf(data=Zambezia, size=0.05, fill="white", colour="black") +
  geom_sf(data=Zambezia_districts, size=0.05, fill=NA, colour="black") +
  geom_sf(size=0.05, fill="grey", colour=NA, alpha=0.5) +
  ## Add cluster label
  geom_sf_text_repel(data=
                        dat_COMSA %>% filter(sero==1) %>% group_by(cluster) %>% slice(1) %>% ungroup() %>% left_join(rank_sums_df, by=c("cluster"="Cluster_clean"))
                      , aes(label=Vuln_Rank, colour=Cluster_Type), force=40, size=2, segment.size=0.1) +
  theme_classic() +
  ggtitle("30 DBS clusters") +scale_color_brewer(palette = "Dark2")+
  xlab("") + ylab("")+theme(legend.title = element_text("Cluster Type")) -> plot_map_Zambezia


library(grid)
library(gridExtra)
grid.newpage()
mainmap <- viewport(width = 1, height = 1, x = 0.5, y = 0.5) # main map
insetmap <- viewport(width = 0.5, height = 0.4, x = 0.8, y = 0.7)
print(plot_map_Zambezia, vp = mainmap)
print(plot_map_Moz, vp = insetmap)

pdf(file="~/Dropbox/COMSA - Serosurveillance/Data/Saki/Final/Map_v2.pdf", width=10)
plot_map_Moz + plot_map_Zambezia + plot_annotation(tag_levels="A")
dev.off()
