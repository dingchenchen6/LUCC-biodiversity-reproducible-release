##%######################################################%##
#                                                          #
####         Organise PREDICTS data for insects         ####
#                                                          #
##%######################################################%##

# This script takes the complete PREDICTS database, selects those entries for
# insects, and organises the data for analysis.

# directories
getwd()
dataDir <- "0_data/"
outDir <- "1_PreparePREDICTSData/"

if(!dir.exists(outDir)) dir.create(outDir)

sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)

# load required libraries
library(predictsFunctions)
library(ggplot2)

sessionInfo()
# First we need to download and load the main 2016 release, as well as the 
# supplementary data released in 2022
url("https://www.dropbox.com/scl/fi/0mqoacviqiiurrtplby2m/predicts.rds?rlkey=m5ijge7w1dthkvm5gbp652uel&dl=1") %>% readRDS() %>% droplevels() -> predicts_2016

url("https://www.dropbox.com/scl/fi/3dqv4ptkqkuswxx350xbx/predicts_2022.rds?rlkey=20g6v43l0vjdh3k8e7p533hri&dl=1") %>% readRDS() %>% droplevels() -> predicts_2022

full_join(predicts_2016,predicts_2022) -> predicts_combined

# Set the path to your local copy of the database
predicts.path <- paste0(dataDir,"database.rds")
predicts_add.path <- paste0(dataDir,"database_add.rds")
# Read in the PREDICTS data
#Predicts <- ReadPREDICTS(predicts.path)
predicts <- ReadPREDICTS(predicts.path)
predicts_add <- ReadPREDICTS(predicts_add.path)
# Select only data for birds
predicts <- predicts[(predicts$Class=="Aves"),]
predicts_add <- predicts_add[(predicts_add$Class=="Aves"),]
predicts_combined <- predicts_combined[(predicts_combined$Class=="Aves"),]
# Combined the old predicts database with added new data
library(dplyr)
Predicts <- bind_rows(
  predicts %>% filter(Class == "Aves"),
  predicts_add %>% filter(Class == "Aves")
)

# Correct effort-sensitive abundance measures (assumes linear relationship between effort and recorded abundance)
predicts <- CorrectSamplingEffort(diversity = predicts)
predicts_add <- CorrectSamplingEffort(diversity = predicts_add)
Predicts <- CorrectSamplingEffort(diversity = Predicts)

predicts_combined <- CorrectSamplingEffort(diversity = predicts_combined)

# Correcting 366812 values for sensitivity to sampling effort (most of these are 0s, 19184 non-zero)

table(Predicts$Diversity_metric)


# Merge sites that have the same coordinates (e.g. multiple traps on a single transect)
Predicts <- MergeSites(diversity = Predicts)

predicts_combined <- MergeSites(diversity = predicts_combined)
# remove rows where land use or use intensity info is missing
predicts.complete <- droplevels(predicts_combined[(
  predicts_combined$Predominant_land_use!="Cannot decide"),])
predicts.complete <- droplevels(predicts.complete[(
  predicts.complete$Use_intensity!="Cannot decide"),])
saveRDS(Predicts, "Predicts.rds")
nrow(predicts.complete)
# 550,242 records

# get counts of n species in major groups
species <- unique(predicts.complete[,c('Order','Taxon_name_entered')])
species <- unique(predicts.complete[,c('Order','Taxon_name_entered','Best_guess_binomial')])
order.counts <- tapply(X = species$Taxon_name_entered,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))


# Calculate site metrics of diversity
Sites <- SiteMetrics(diversity = Predicts,
                     extra.cols = c("Predominant_land_use","Longtitude","Latitude",
                                    "SSB","SSBS", "Biome","Sampling_efforts","Rescaled_sampling_efforts"))
Sites <- SiteMetrics(diversity = predicts_combined,
                     extra.cols = c("Predominant_land_use","Longtitude","Latitude",
                                    "SSB","SSBS", "Biome","Sampling_efforts","Rescaled_sampling_efforts"))
# First, we will rearrange the land-use classification a bit
Sites$LandUse <- paste(Sites$Predominant_land_use)

# Drop classification where land use could not be identified
Sites$LandUse[(Sites$LandUse=="Cannot decide")] <- NA

# Now make the variable a factor, and set the reference level to primary vegetation
Sites$LandUse <- factor(Sites$LandUse)
Sites$LandUse <- relevel(Sites$LandUse,ref="Primary vegetation")

Sites$Use_intensity[Sites$Use_intensity=="Cannot decide"] <- NA

# combine LU and UI 
Sites$UI <- paste0(Sites$LandUse,'_',Sites$Use_intensity)
Sites$UI[grep("NA",Sites$UI)] <- NA

# recode according to land use and use intensity combinations
Sites$UI2 <- dplyr::recode(Sites$UI,
                           'Primary vegetation_Minimal use' = 'Primary vegetation',
                           'Cropland_Light use' = 'Agriculture_High',
                           'Secondary vegetation (indeterminate age)_Minimal use' = 'Secondary vegetation',
                           'Urban_Light use' = 'Urban',
                           'Secondary vegetation (indeterminate age)_Light use' = 'Secondary vegetation',
                           'Cropland_Intense use' = 'Agriculture_High',
                           'Cropland_Minimal use' = 'Agriculture_Low',
                           'Pasture_Light use' = 'Agriculture_Low',
                           'Pasture_Minimal use' = 'Agriculture_Low',
                           'Intermediate secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Secondary vegetation (indeterminate age)_Intense use' = 'Secondary vegetation',
                           'Pasture_Intense use' = 'Agriculture_High',
                           'Urban_Minimal use' = 'Urban',
                           'Primary vegetation_Light use' = 'Primary vegetation',
                           'Young secondary vegetation_Light use' = 'Secondary vegetation',
                           'Urban_Intense use' = 'Urban',
                           'Primary vegetation_Intense use' = 'Primary vegetation',
                           'Young secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Minimal use' = 'Agriculture_Low',
                           'Plantation forest_Intense use' = 'Agriculture_High',
                           'Young secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Light use' = 'Agriculture_High',
                           'Mature secondary vegetation_Light use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Light use' = 'Secondary vegetation')

Sites$Use_intensity[((Sites$LandUse=="Mature secondary vegetation") & 
                       (Sites$Use_intensity=="Intense use"))] <- "Light use"
Sites$Use_intensity[((Sites$LandUse=="Intermediate secondary vegetation") & 
                       (Sites$Use_intensity=="Intense use"))] <- "Light use"
Sites$Use_intensity[((Sites$LandUse=="Young secondary vegetation") & 
                       (Sites$Use_intensity=="Intense use"))] <- "Light use"

# remove the urban sites and sites that are NA in UI2
#sites <- sites[!sites$UI2 == "Urban", ]
Sites <- Sites[!is.na(Sites$UI2), ]


Sites <- droplevels(Sites)

# transform abundance values 
Sites$LogAbund <- log(Sites$Total_abundance+1)


# Remove sites without coordinates
Sites <- Sites[!is.na(Sites$Latitude), ]


# save the prepared dataset
saveRDS(object = sites,file = paste0(outDir,"PREDICTSSiteData.rds"))
saveRDS(object = Sites,file = paste0(outDir,"PREDICTSSiteData1.rds"))
saveRDS(object = Sites,file = paste0(outDir,"PREDICTSSiteData2.rds"))

saveRDS(object = Predicts,file = paste0(outDir,"PREDICTSDatabase1.rds"))
saveRDS(object = predicts_combined,file = paste0(outDir,"PREDICTSDatabase2.rds"))


##%######################################################%##
#                                                          #
#### redo dataset summaries after removing other sites  ####
#                                                          #
##%######################################################%##


Predicts2 <- predicts_combined[predicts_combined$SSBS %in% Sites$SSBS, ]


table(Predicts2$Diversity_metric)


nrow(Predicts2)
# 550,242 records

# get counts of n species in major groups
Species <- unique(Predicts2[,c('Order','Best_guess_binomial')])
Species <- unique(Predicts2[,c('Order','Taxon_name_entered')])
Order.counts <- tapply(X = Species$Taxon_name_entered,
                       INDEX = Species$Order,
                       FUN = function(sp) length(unique(sp)))

sum(Order.counts)
# 6,049

##%######################################################%##
#                                                          #
####            basic map of PREDICTS sites             ####
#                                                          #
##%######################################################%##


# plot the raster in ggplot
map.world <- map_data('world')

# map of sites
P1 <-ggplot() +
         geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
         geom_point(data = Sites, aes(x = Longitude, y = Latitude), col = c("#1E90FF"), fill = c("#104E8B"), shape = 21) +
         theme(axis.title = element_blank(), 
               plot.background = element_blank(), 
               panel.background = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank()) +
         ggtitle("a.")
   
P1

if (!require(devtools)) install.packages("devtools")
devtools::install_github("valentinitnelav/plotbiomes")
library(plotbiomes)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(plotbiomes)
library(ggplot2)

# 直接下载 + 直接提取 + 直接加入 Sites
library(terra)
library(geodata)


# 一次性：下载 → 建点 → 提取 → 赋值
Sites <- within(Sites, {
  bio <- rast(worldclim_global("bio", res = 10, path = "data_worldclim"))
  v   <- extract(bio[[c(1, 12)]],
                 vect(Sites, c("Longitude", "Latitude"), "EPSG:4326"))
  Temp_mean <- v[, 2] / 10   # BIO1: 0.1°C → °C
  Precip    <- v[, 3]        # BIO12: mm
})

# 直接提取并加入 Sites（0.5 分辨率=30 arc-sec）
Sites <- Sites |>
  cbind(
    Temp_mean = extract(
      geodata::worldclim_global(var = "bio", res = 0.5, path = tempdir())[[1]],
      vect(Sites, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
    )[ ,2] / 10,
    
    Precip = extract(
      geodata::worldclim_global(var = "bio", res = 0.5, path = tempdir())[[12]],
      vect(Sites, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
    )[ ,2]
  )


plot_final <- ggplot() +
  # 1. Whittaker biome polygons
  geom_polygon(
    data = Whittaker_biomes,
    aes(x = temp_c , y = precp_cm, fill = biome),
    colour = "gray98", size = 1
  ) +
  
  # 2. Biome colors
  scale_fill_manual(
    name   = "Whittaker biomes",
    breaks = names(Ricklefs_colors),
    labels = names(Ricklefs_colors),
    values = Ricklefs_colors
  ) +
  
  # 3. Sampling points
  geom_point(
    data = Sites,
    aes(x = Temp_mean *10, y = Precip/10),
    size = 1.5, shape = 21,
    colour = "white", 
    fill = "black",
    stroke = 1, alpha = 0.5
  ) +
  
  # 4. X/Y axis labels + custom tick marks
  scale_x_continuous(
    name = "Mean annual temperature (°C)",    # X轴名称
    breaks = seq(-10, 30, by = 10),            # X轴刻度
    limits = c(-10, 30)                       # 可根据需要修改
  ) +
  scale_y_continuous(
    name = "Mean annual precipitation (cm)",       # Y轴名称
    breaks = seq(0, 400, by = 100),            # Y轴刻度
    limits = c(0, 400)                        # 根据数据调整
  ) +
  
  theme_bw()

plot_final

# save plot
ggsave(filename = paste0(outDir, "/Whittaker biomes_sampling points.pdf"), height = 4, width = 6)


##==============================================================##
##  WWF Biomes + Sites 全球分布图（Nature 风格）
##  WWF 陆地生态区做背景多边形 + Sites 采样点叠加
##  WWF terrestrial ecoregions as background + Sites as points
##==============================================================##

library(sf)
library(dplyr)
library(ggplot2)

##--------------------------------------------------------------##
## 1. 读入 WWF shapefile 和 Sites 数据
##    Read WWF shapefile and Sites data
##--------------------------------------------------------------##

# 请改成你自己的 wwf_terr_ecos.shp 路径
# Change to your own path of wwf_terr_ecos.shp
wwf <- st_read("0_data/wwf_terr_ecos.shp")

# 确保是经纬度坐标系 WGS84
# Make sure CRS is WGS84 lon/lat
wwf <- st_transform(wwf, crs = 4326)

# 假定 Sites 已经在环境中，是一个 data.frame，
# 且包含列：Longitude, Latitude, BIOME (1–14)
# Assume Sites is already loaded with columns:
# Longitude, Latitude, BIOME (1–14)


##--------------------------------------------------------------##
## 2. 定义 WWF 14 个 Biome 的官方名称与统一颜色
##    Define official names & color palette for 14 WWF biomes
##--------------------------------------------------------------##

biome_lookup <- tibble::tibble(
  BIOME = 1:14,
  biome_name = c(
    "Tropical & Subtropical Moist Broadleaf Forests",        # 1
    "Tropical & Subtropical Dry Broadleaf Forests",          # 2
    "Tropical & Subtropical Coniferous Forests",             # 3
    "Temperate Broadleaf & Mixed Forests",                   # 4
    "Temperate Conifer Forests",                             # 5
    "Boreal Forests/Taiga",                                  # 6
    "Tropical & Subtropical Grasslands, Savannas & Shrublands", # 7
    "Temperate Grasslands, Savannas & Shrublands",           # 8
    "Flooded Grasslands & Savannas",                         # 9
    "Montane Grasslands & Shrublands",                       #10
    "Tundra",                                                #11
    "Mediterranean Forests, Woodlands & Scrub",              #12
    "Deserts & Xeric Shrublands",                            #13
    "Mangroves"                                              #14
  ),
  # 颜色参考常见 WWF/Olson 2001 地图配色，保证绿色系对应森林，
  # 黄/棕对应草地和荒漠，蓝色对应湿地/苔原/红树林
  # Colors roughly follow WWF/Olson 2001 style:
  biome_col = c(
    "#009A44",  # 1 深绿 Tropical moist broadleaf forest
    "#66A61E",  # 2 亮绿 Tropical dry broadleaf forest
    "#4DAC26",  # 3 绿 Tropical coniferous forest
    "#B3DE69",  # 4 淡绿 Temperate broadleaf & mixed
    "#1B9E77",  # 5 深蓝绿 Temperate conifer forest
    "#5E81AC",  # 6 蓝灰 Boreal forest/taiga
    "#FEE08B",  # 7 浅黄 Tropical grasslands
    "#FDC863",  # 8 金黄 Temperate grasslands
    "#F6E8C3",  # 9 米色 Flooded grasslands & savannas
    "#A6CEE3",  #10 淡蓝 Montane grasslands
    "#E0ECF4",  #11 冷淡蓝 Tundra
    "#F4A582",  #12 砖红 Mediterranean
    "#FEE090",  #13 浅沙色 Deserts & xeric shrublands
    "#2C7BB6"   #14 深蓝 Mangroves
  )
)

# 为后面 scale_fill_manual 准备一个 named vector
# Prepare a named vector for scale_fill_manual
biome_fill_vec <- setNames(biome_lookup$biome_col,
                           biome_lookup$biome_name)


##--------------------------------------------------------------##
## 3. 统一 WWF shapefile 与 Sites 的 Biome 名称
##    Harmonize biome names between WWF and Sites
##--------------------------------------------------------------##

# WWF shapefile 里通常 BIOME 列也是 1–14 的整数
# 一些数据集字段名可能是 "BIOME" 或 "BIOME_ID"，这里假定叫 BIOME
# If your field name is different, change "BIOME" below.

wwf <- wwf %>%
  left_join(biome_lookup, by = "BIOME")

# Sites 里 BIOME 列也是 1–14
Sites <- Sites %>%
  left_join(biome_lookup, by = "BIOME")

# 检查是否有 NA（说明有 Biome 代码不在 1:14 中）
# Check for NA (means some BIOME codes not in 1:14)
table(is.na(wwf$biome_name))
table(is.na(Sites$biome_name))


##--------------------------------------------------------------##
## 4. Nature 风格地图：Biome 做背景，多边形填色 + Sites 黑点
##    Nature-style global map: biomes as polygons + Sites points
##--------------------------------------------------------------##

# 将 Sites 转为空间对象以便 coord_sf 更好控制（可选）
# Convert Sites to sf object (optional but nice with coord_sf)
sites_sf <- st_as_sf(
  Sites,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

p_biome_sites <- ggplot() +
  # (1) Biome 多边形背景 / Biome polygons as background
  geom_sf(
    data = wwf,
    aes(fill = biome_name),
    color = NA,        # 不画边框，避免太乱；如需边界可改成 "grey70"
    size  = 0
  ) +
  
  # (2) Sites 采样点 / Sites sampling points
  geom_sf(
    data = sites_sf,
    colour = "black",
    fill   = "black",
    shape  = 21,
    size   = 1.5,
    alpha  = 0.5
  ) +
  
  # (3) 统一 Biome 填充颜色 / Unified biome fill colors
  scale_fill_manual(
    name   = "Biome",
    values = biome_fill_vec,
    breaks = biome_lookup$biome_name,
    labels = biome_lookup$biome_name,
    guide  = guide_legend(
      title.position = "top",
      ncol = 2
    )
  ) +
  
  # (4) 坐标范围和投影 / Coordinate range & projection
  coord_sf(
    xlim   = c(-180, 180),
    ylim   = c(-60, 85),
    expand = FALSE
  ) +
  
  # (5) 坐标轴标签 / Axis labels
  labs(
    x = "Longitude (°)",
    y = "Latitude (°)"
  ) +
  
  # (6) Nature 风格主题（接近 Nature 文章中的地图样式）
  #     Nature-style theme: clean background, thin axes, small text
  theme_bw(base_size = 10) +
  theme(
    panel.background   = element_rect(fill = "white", colour = NA),
    panel.grid.major   = element_line(colour = "grey90", size = 0.2),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 8),
    axis.ticks         = element_line(colour = "grey40", size = 0.3),
    #legend.position    = "bottom",
    legend.title       = element_text(size = 9, face = "bold"),
    legend.text        = element_text(size = 7),
    legend.key.width   = unit(0.7, "cm"),
    legend.key.height  = unit(0.4, "cm"),
    legend.background  = element_rect(fill = "white", colour = "grey80"),
    legend.box.margin  = margin(t = 2, r = 2, b = 2, l = 2),
    plot.margin        = margin(t = 5, r = 5, b = 5, l = 5),
    legend.position = "bottom") +
  ggtitle("a.")
  
p_biome_sites

library(cowplot)

# 1. 先创建原始图（带图例）
p_full <- p_biome_sites  # 你的原图对象名称

# 2. 抽取图例
legend <- cowplot::get_legend(
  p_full +
    guides(fill = guide_legend(ncol = 3)) +  # 图例放 4 列（可调）
    theme(
      legend.position = "bottom"
      #legend.box = "horizontal"
    )
)

# 3. 去掉原图中的图例
p_nolegend <- p_full + theme(legend.position = "none")

# 4. 拼接图 + 宽度一致的图例
final_plot <- cowplot::plot_grid(
  p_nolegend,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.35)   # 下方图例高度占比（可调整）
)

final_plot

ggsave(filename = paste0(outDir, "/p_biome_sites2.pdf"))

P1 <- ggplot() +
  # 世界灰色底图
  geom_map(data = map.world, map = map.world,
           aes(x = long, y = lat, group = group, map_id = region),
           fill = "grey90", colour = "grey60", size = 0.2) +
  
  # 采样点按 Biome 着色
  geom_point(data = Sites,
             aes(x = Longitude, y = Latitude, fill = Biome),
             shape = 21, color = "black", size = 2) +
  
  scale_fill_viridis_d(option = "turbo", name = "Biome") +
  
  theme(axis.title = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right") +
  ggtitle("a.")

P1

# 获取国家边界
world <- ne_countries(scale = "medium", returnclass = "sf")

# 获取 Whittaker biome 底图
p_whittaker <- whittaker_base_plot()  # 这是 ggplot 对象

# 在 whittaker 底图基础上添加国家边界与采样点
P1 <- p_whittaker +
  
  # 国家边界（黑色线）
  geom_sf(data = world, fill = NA, color = "black", size = 0.2) +
  
  # 采样点（黑色）
  geom_point(data = Sites,
             aes(x = Longitude, y = Latitude),
             color = "black", size = 1.5, alpha = 0.8) +
  
  ggtitle("a.") +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

P1


# save plot
ggsave(filename = paste0(outDir, "/PREDICTS_points_map1.pdf"), height = 4, width = 8)



### Basic summaries ###

# nstudies/ nsites - all
length(unique(Sites$SS)) # 138
length(unique(Sites$SSBS)) # 6462

# nstudies/nsites - abun
length(unique(Sites[!is.na(Sites$LogAbund) , 'SS'])) # 121
length(unique(Sites[!is.na(Sites$LogAbund) , 'SSBS'])) # 5817

# reviewer request, UI2 by Biome

table(Sites$Biome, Sites$UI2)

#                                                         Agriculture_High Agriculture_Low Primary vegetation Secondary vegetation
# Boreal Forests/Taiga                                                    2               6                172                   13
# Temperate Conifer Forests                                               4             103                 10                   88
# Temperate Broadleaf & Mixed Forests                                  1308             671                315                  787
# Montane Grasslands & Shrublands                                         2             200                247                   33
# Temperate Grasslands, Savannas & Shrublands                            15              11                 15                   27
# Mediterranean Forests, Woodlands & Scrub                               21              32                 96                   58
# Deserts & Xeric Shrublands                                              0              30                 16                   16
# Tropical & Subtropical Grasslands, Savannas & Shrublands               56              47                175                   78
# Tropical & Subtropical Coniferous Forests                               2              26                 32                   43
# Flooded Grasslands & Savannas                                           6               6                  0                    0
# Tropical & Subtropical Dry Broadleaf Forests                           62              66                 13                   50
# Tropical & Subtropical Moist Broadleaf Forests                        293             119                420                  283
# Mangroves                                                               8               0                  5                    7

table(Sites[!is.na(Sites$LogAbund), 'Biome'], Sites[!is.na(Sites$LogAbund), 'UI2'])

#                                                          Agriculture_High Agriculture_Low Primary vegetation Secondary vegetation
# Boreal Forests/Taiga                                                    2               6                172                   13
# Temperate Conifer Forests                                               4             103                 10                   88
# Temperate Broadleaf & Mixed Forests                                  1280             669                289                  729
# Montane Grasslands & Shrublands                                         2             200                247                   33
# Temperate Grasslands, Savannas & Shrublands                            15              11                 15                   27
# Mediterranean Forests, Woodlands & Scrub                               21              26                 92                   58
# Deserts & Xeric Shrublands                                              0              30                 16                   16
# Tropical & Subtropical Grasslands, Savannas & Shrublands               56              47                165                   70
# Tropical & Subtropical Coniferous Forests                               2              26                 32                   43
# Flooded Grasslands & Savannas                                           0               0                  0                    0
# Tropical & Subtropical Dry Broadleaf Forests                           62              66                 13                   50
# Tropical & Subtropical Moist Broadleaf Forests                        265             110                354                  204
# Mangroves                                                               8               0                  5                    7

t.end <- Sys.time()

print(round(t.end - t.start,0))

predicts_effort<-predicts[,c("SSBS","Measurement","Diversity_metric_type","Sampling_effort","Rescaled_sampling_effort","Taxon_name_entered")]

predicts_effort <- predicts_effort %>%
  dplyr::filter(Diversity_metric_type == "Abundance")

predictsSites_joined<-left_join(model_data,predicts_effort,by="SSBS")
predictSites_filtered <- predictsSites_joined %>%
  filter(Rescaled_sampling_effort == 1)
predictSites_filtered$Abund<-predictSites_filtered$Measurement


