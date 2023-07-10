source("R-code/00_preamble.R")


# aim ---------------------------------------------------------------------

# this script calculates per botanical country:
# 1) the number of accepted hybrid species
# 2) the number of neophyte species
# 3) the number of native species
# 3) the size of the country
# 4) the number of authors per country (weighted by species range size)
# 5) gdp

# the script explores the question whether there is a correlation between
# the number of neophytes and hybrids.



# hybrids -----------------------------------------------------------------
# read species list
# can be downloaded from: https://storage.googleapis.com/kew-dev-backup/world_checklist_names_and_distribution_feb_21.zip
kew_sp <- fread("Data/kew/checklist_names.txt", encoding = "UTF-8")

# kew species data
kew_hyb <- kew_sp %>% filter(species_hybrid == "×" | species_hybrid == "+")

# how many of these hybrids are accepted
kew_hyb %>% count(taxon_status)

# we want only the accepted ones, and importantly also not the artificial hybrids
kew_hyb <- kew_hyb %>% filter(taxon_status == "Accepted")

# kew data on distribution
kew_dis <- fread("Data/kew/checklist_distribution.txt")

# join
d_hyb <- left_join(kew_hyb, kew_dis)

# sum hybrid occurrence per area code l3
d_hyb <- d_hyb %>% group_by(area_code_l3) %>% summarise(num_hyb = n_distinct(plant_name_id)) 

# get layer for botanical countries
sf::sf_use_s2(FALSE)

# read shapefile
l3 <- sf::read_sf("Data/bot_countries/wgsrpd-master/level3/level3.shp")
l3 <- sf::st_crop(l3, xmin = -180, ymin = -90, xmax = 180, ymax = 90)

l3 <- sf::st_transform(l3, crs = "+proj=moll")
l3 <- st_simplify(l3, dTolerance = 1000)

# add number of hybrids
l3_merged <- left_join(l3, d_hyb, by = c("LEVEL3_COD" = "area_code_l3") )


# tmap alternative
(m1 <- tmap::tm_shape(l3_merged) +
    tm_polygons("num_hyb",
                palette=c("white","#E1BE6A", "#40B0A6"),
                #palette=c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac"),
                style= "quantile",
                title="Hybrids",
                legend.is.portrait = T) +
    tm_layout(fontfamily = "Arial Narrow", 
              outer.margins = c(0, 0, 0, 0),
              legend.outside.position = "right",
              legend.outside.size = 0.2,
              legend.outside = T,
              frame = FALSE))


tmap_save(m1, "Figures/worlmap_hybrids.png",   
          width = 6.5,
          height = 2.95,
          dpi = 600)


# ggplot alternative
# specify the number of desired classes
n <- 5
# create quantile breaks
l3_merged <- l3_merged %>% 
  mutate(quantile_num_hyb = cut(num_hyb, 
                                breaks = c(quantile(num_hyb, probs = seq(0, 1, by = 1/n), na.rm = T)),
                                include.lowest = T, dig.lab=10)) %>% 
  mutate(pretty_labs = pretty_labels(quantile_num_hyb))

(cc <- paletteer_d("dichromat::LightBluetoDarkBlue_7")[-c(4,6)])
ggplot() +
  geom_sf(data = l3_merged, aes(fill = quantile_num_hyb), size =0.1) +
  scale_fill_manual(values=cc, na.value="grey70", name = "Hybrids",
                    labels = levels(l3_merged$pretty_labs)) +
  theme_ipsum(
    axis_title_size = 12,
    plot_title_size = 12,
    axis_title_just = "mm",
    grid = "XY",
    base_family = "Helvetica"
  ) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) -> map_hybrids



# neophytes ---------------------------------------------------------
# species data
glonaf_sp <- read_delim("Data/GLONAF/Taxon_x_List_GloNAF_vanKleunenetal2018Ecology.csv", 
                        delim = "\t", escape_double = FALSE, 
                        locale = locale(encoding = "utf-16"), 
                        trim_ws = TRUE)

# select relevant columns and subsets of data from glonaf data
glonaf_sp <- glonaf_sp %>% 
  filter(status == "naturalized") %>% 
  filter(name_status == "accepted") %>% 
  select(standardized_name, hybrid, status, region_id) %>% distinct

# how many species in the database
glonaf_sp %>% select(standardized_name) %>% distinct %>% nrow

# how many hybrid species are in glonaf?
glonaf_sp %>% filter(hybrid == 1) %>% select(standardized_name) %>% distinct

# should these be removed? Probably because otherwise we correlate hybrid vs 
# hybrid
glonaf_sp <- glonaf_sp %>% filter(hybrid != 1)
glonaf_sp %>% select(standardized_name) %>% distinct %>% nrow


# join with distribution data
# read region ids for l3 areas
glonaf_reg <- read_csv("Data/GLONAF/Region_GloNAF_vanKleunenetal2018Ecology.csv")


# join glonaf species with their regions to glonaf_reg to have this on tdwg3 resolution
d_neo <- full_join(glonaf_sp, glonaf_reg)

# summarize per tdwg3 the number of naturalized species
d_neo <- d_neo %>% 
  group_by(tdwg3) %>% 
  summarise(num_neo = n_distinct(standardized_name))
d_neo %>% arrange(desc(num_neo))

# tdwg3 with numerical, won't match with powo, recode
d_neo <- d_neo %>% filter(str_detect(tdwg3,"^\\s*[0-9]*\\s*$"))

# add number of hybrids
l3_glonaf <- left_join(l3, d_neo, by = c("LEVEL3_COD" = "tdwg3") )

# tmap alternative
(m2 <- tm_shape(l3_glonaf) +
    tm_polygons("num_neo", 
                palette=c("white","#E1BE6A", "#40B0A6"),
                style= "quantile",
                title="Neophytes",
                legend.is.portrait = T) +
    tm_layout(fontfamily = "Arial Narrow", 
              outer.margins = c(0, 0, 0, 0),
              legend.outside.position = "right",
              legend.outside.size = 0.2,
              legend.outside = T,
              frame = FALSE))

tmap_save(m, "Figures/worlmap_nonnatives.png",   
          width = 6.5,
          height = 2.95,
          dpi = 600)


# ggplot alternative
# specify the number of desired classes
n <- 5
# create quantile breaks
l3_glonaf <- l3_glonaf %>% 
  mutate(quantile_num_neo = cut(num_neo, 
                                breaks = c(quantile(num_neo, probs = seq(0, 1, by = 1/n), na.rm = T)),
                                include.lowest = T, dig.lab=10)) %>% 
  # create prettier labels for legend
  mutate(pretty_labs = pretty_labels(quantile_num_neo))

(cc <- paletteer_d("dichromat::LightBluetoDarkBlue_7")[-c(4,6)])
ggplot() +
  geom_sf(data = l3_glonaf, aes(fill = quantile_num_neo), size =0.1) +
  scale_fill_manual(values=cc, na.value="grey70", name = "Neophytes",
                    labels = levels(l3_glonaf$pretty_labs)) +
  theme_ipsum(
    axis_title_size = 12,
    plot_title_size = 12,
    axis_title_just = "mm",
    grid = "XY",
    base_family = "Helvetica"
  ) +
  theme(plot.margin = margin(0, 0, 0,0, "cm")) -> map_neophytes


plot_grid(map_hybrids + theme(legend.justification = c(0,1)), 
          map_neophytes + theme(legend.justification = c(0,1)), 
          ncol = 1, labels = c("a", "b"), align = "v")


# join neophyte and hybrid data
d <- full_join(d_hyb %>% rename(tdwg3 = area_code_l3), d_neo)



# calculate total number species from Powo --------------------------------
acc_sp <- kew_sp %>% 
  filter(taxon_rank == "Species") %>% 
  filter(taxon_status == "Accepted") %>% 
  filter(species_hybrid == "") %>% 
  filter(infraspecific_rank == "") %>% 
  filter(infraspecies == "") %>% 
  select(accepted_plant_name_id) %>% 
  distinct

kew_dis_native <- kew_dis %>% 
  filter(location_doubtful == 0) %>% 
  filter(extinct == 0) %>% 
  filter(introduced == 0)

kew_tot <- left_join(acc_sp, kew_dis_native, 
                     by = c("accepted_plant_name_id" = "plant_name_id"))

d_nat <- kew_tot %>% 
  group_by(area_code_l3) %>% 
  summarise(num_nat = n_distinct(accepted_plant_name_id))

# join native spp. number to d
d <- full_join(d, d_nat %>% rename(tdwg3 = area_code_l3))



# calculate tdwg3 country size as a confound ----------------------------
# calculate polygon area in l3
l3$area <- st_area(l3)

# convert into km2
l3$area <- l3$area / 1e+6

# join area to d
d <- left_join(d, l3 %>% st_drop_geometry() %>% 
                 select(area, LEVEL3_COD) %>% 
                 mutate(area = drop_units(area)), 
               by = c("tdwg3" = "LEVEL3_COD"))



# quantify sampling effort: number of authors -----------------------------
# approach according to: https://academic.oup.com/bioscience/article/60/10/798/231551
# but we attempt to adjust for spill over counts, taking into account species range
# size

# in how many botanical countries does a species occur?
kew_sp_range <- kew_tot %>% select(accepted_plant_name_id, area_code_l3) %>% distinct %>% 
  group_by(accepted_plant_name_id) %>% 
  summarize(n_countries = n_distinct(area_code_l3))

# attach author name
kew_sp_range <- left_join(kew_sp_range , kew_sp %>% select(plant_name_id, primary_author),
                          by = c("accepted_plant_name_id" = "plant_name_id"))

kew_sp_range <- kew_sp_range %>% filter(!is.na(primary_author))

# calculate the number of authors per country weighted by the average range size of the species they described
# and that occur in that country
kew_author <- left_join(kew_sp_range, kew_tot %>% 
                          select(accepted_plant_name_id, area_code_l3)) %>% 
  filter(area_code_l3 != "")

kew_author <- kew_author %>% 
  group_by(area_code_l3, primary_author) %>% 
  summarize(mean_range = mean(n_countries))

kew_author <- kew_author %>% mutate(author_endemism = 1/mean_range)
kew_author <- kew_author %>% 
  ungroup %>% 
  group_by(area_code_l3) %>% 
  summarise(author_weighted = sum(author_endemism))

# divide by number of species
d_author <- left_join(kew_author, d_nat) %>% 
  mutate(ratio = author_weighted/num_nat)

# display on map
l3_author <- left_join(l3, d_author, by = c("LEVEL3_COD" = "area_code_l3") )

(m3 <- tmap::tm_shape(l3_author) +
    tm_polygons("ratio",
                style= "quantile",
                title="Authors per species"))

# join authors to d
d <- left_join(d,
                d_author %>% select(area_code_l3, ratio, author_weighted) %>% 
                rename(tdwg3 = area_code_l3, authors_per_spp = ratio))




# quantify sampling effort: gdp -------------------------------------------
# get gdp data
gdp <- WDI(
  country = "all",
  indicator = "NY.GDP.PCAP.KD",
  start = 1960,
  end = 2020,
  extra = FALSE,
  cache = NULL,
  latest = NULL,
  language = "en"
)

# calculate average gdp over the years
gdp <- gdp %>% rename(gdp = 3) %>% group_by(iso2c, country) %>% summarize(avg_gdp = mean(gdp, na.rm = T))


# match to botanical countries (tdwg3)
data("tdwg_level4")
gdp <- left_join(tdwg_level4 %>% select(ISO_Code, Level3_cod) %>% st_drop_geometry() %>% distinct,
                 gdp %>% rename(ISO_Code = iso2c))

# if a botanical country has more than 2 value?
gdp <- gdp %>% group_by(Level3_cod) %>% summarize(avg_gdp = mean(avg_gdp, na.rm = T))

# join to d
d <- left_join(d , gdp %>% rename(tdwg3 = Level3_cod) %>% select(tdwg3, avg_gdp))



# visualization -----------------------------------------------------------
d <- d %>% 
  filter(tdwg3 != "") %>% 
  filter(!is.na(num_hyb)) %>% 
  filter(!is.na(num_neo))


# calculate a couple of ratios
d_ratios <- d %>% 
  mutate(hyb_per_nat = num_hyb/num_nat, neo_per_nat = num_neo/num_nat) %>% 
  mutate(hyb_per_km2 = num_hyb/area, neo_per_km2 = num_neo/area) %>% 
  mutate(hyb_per_author = num_hyb/author_weighted, neo_per_author = num_neo/author_weighted) %>% 
  mutate(hyb_per_gdp = num_hyb/avg_gdp, neo_per_gdp = num_neo/ avg_gdp)


# visualization
# 1) number of neophytes against number of hybrids
(ggplot(d, aes(x = num_neo, y = num_hyb)) +
    geom_point(alpha = 0.2, size = 2.4, pch = 21, fill = "#00009a") +
    geom_smooth(col = "black", fill = "#9a9a00", alpha = 0.2, method = "loess") +
    theme_ipsum(
      axis_title_size = 12,
      plot_title_size = 12,
      axis_title_just = "mm"
    ) +
    scale_x_log10()+
    scale_y_log10()+
    theme(
      text = element_text(size = 14),
      plot.margin=unit(rep(0.2,4),"cm"),
      legend.position = "none",
      aspect.ratio = 1,
      legend.title = element_blank() ,
      axis.line = element_line() ) +
    labs(y = "Number of hybrids per country", 
         x = "Number of neophytes per country") -> p1)

# 2) control for number of natives
(ggplot(d_ratios, aes(x = neo_per_nat, y = hyb_per_nat)) +
    geom_point(alpha = 0.2, size = 2.4, pch = 21, fill = "#00009a") +
    geom_smooth(col = "black", fill = "#9a9a00", alpha = 0.2, method = "loess") +
  theme_ipsum(
    axis_title_size = 12,
    plot_title_size = 12,
    axis_title_just = "mm"
  ) +
    scale_x_log10()+
    scale_y_log10()+
  theme(
    text = element_text(size = 14),
    plot.margin=unit(rep(0.2,4),"cm"),
    legend.position = "none",
    aspect.ratio = 1,
    legend.title = element_blank() ,
    axis.line = element_line() ) +
  labs(y = "Number of hybrids / number of natives", 
       x = "Number of neophytes / number of natives") -> p2)


# 3) control for area
(ggplot(d_ratios, aes(x = neo_per_km2, y = hyb_per_km2)) +
    geom_point(alpha = 0.2, size = 2.4, pch = 21, fill = "#00009a") +
    geom_smooth(col = "black", fill = "#9a9a00", alpha = 0.2, method = "loess") +
    theme_ipsum(
      axis_title_size = 12,
      plot_title_size = 12,
      axis_title_just = "mm"
    ) +
    scale_x_log10()+
    scale_y_log10()+
    theme(
      text = element_text(size = 14),
      plot.margin=unit(rep(0.2,4),"cm"),
      legend.position = "none",
      aspect.ratio = 1,
      legend.title = element_blank() ,
      axis.line = element_line() ) +
    labs(y = "Number of hybrids / km2", 
         x = "Number of neophytes / km2") -> p3)

# 4) per author
(ggplot(d_ratios, aes(x = neo_per_author, y = hyb_per_author)) +
    geom_point(alpha = 0.2, size = 2.4, pch = 21, fill = "#00009a") +
    geom_smooth(col = "black", fill = "#9a9a00", alpha = 0.2, method = "loess") +
    theme_ipsum(
      axis_title_size = 12,
      plot_title_size = 12,
      axis_title_just = "mm"
    ) +
    scale_x_log10()+
    scale_y_log10()+
    theme(
      text = element_text(size = 14),
      plot.margin=unit(rep(0.2,4),"cm"),
      legend.position = "none",
      aspect.ratio = 1,
      legend.title = element_blank() ,
      axis.line = element_line() ) +
    labs(y = "Number of hybrids / number of authors", 
         x = "Number of neophytes / number of authors") -> p4)


# per gdp, that does not make much sense, but still
(ggplot(d_ratios, aes(x = neo_per_gdp, y = hyb_per_gdp)) +
    geom_point(alpha = 0.2, size = 2.4, pch = 21, fill = "#00009a") +
    geom_smooth(col = "black", fill = "#9a9a00", alpha = 0.2, method = "loess") +
    theme_ipsum(
      axis_title_size = 12,
      plot_title_size = 12,
      axis_title_just = "mm"
    ) +
    scale_x_log10()+
    scale_y_log10()+
    theme(
      text = element_text(size = 14),
      plot.margin=unit(rep(0.2,4),"cm"),
      legend.position = "none",
      aspect.ratio = 1,
      legend.title = element_blank() ,
      axis.line = element_line() ) +
    labs(y = "Number of hybrids / GDP", 
         x = "Number of neopyhtes / GDP") -> p5)

# this all seems pretty robust
(p1|p3) /
(p4|p5)

ggsave(dpi = 600, height = 5.00, width = 7.18, bg = "white",
       "Figures/hybrid_neo_geography_supp.png")


# visualization and statistical models -----------------------------------
# number of hybrids versus number of neophytes, country area and taxonomic effort

# prepare data:
d_model <- d %>%
  select(tdwg3, num_hyb, num_neo, area, avg_gdp, author_weighted) %>% 
  na.omit %>% 
  mutate(num_hyb_log = log10(num_hyb),
         num_neo_log = log10(num_neo),
         area_log = log10(area),
         avg_gdp_log = log10(avg_gdp),
         author_log = log10(author_weighted))

# run model
b_mod <- brm(data = d_model, 
             family = gaussian,
             num_hyb_log ~ num_neo_log + area_log + author_log + avg_gdp_log,
             iter = 2000, warmup = 500, chains = 4, cores = 4, file = "Models/geography_numbers.rds")

b_mod <- readRDS("Models/geography_numbers.rds")

# interpret the coefficient: 
# Both dependent/response variable and independent/predictor variable(s) are 
# log-transformed. Interpret the coefficient as the percent increase in the 
# dependent variable for every 1% increase in the independent variable. 
# Example: the coefficient is 0.198. For every 1% increase in the independent 
# variable, our dependent variable increases by about 0.20%. For x percent 
# increase, calculate 1.x to the power of the coefficient, subtract 1, and 
# multiply by 100. Example: For every 20% increase in the independent variable, 
# our dependent variable increases by about (1.20 0.198 – 1) * 100 = 3.7 percent.
(1.2^0.26 -1) * 100

nd = tibble(avg_gdp_log = mean(d_model$avg_gdp_log),
            area_log = mean(d_model$area_log),
            author_log = mean(d_model$author_log),
            num_neo_log = seq(min(d_model$num_neo_log), 
                              max(d_model$num_neo_log), 
                              by = 0.1))

fit <- fitted(b_mod, newdata = nd) %>% 
    data.frame() %>% select(-Est.Error) %>% 
    bind_cols(nd)

draws <- d_model %>%
  modelr::data_grid(num_neo_log = seq_range(num_neo_log, n = 101),
                    avg_gdp_log = mean(d_model$avg_gdp_log),
                    area_log = mean(d_model$area_log),
                    author_log = mean(d_model$author_log)) %>%
  add_epred_draws(b_mod, ndraws = 100)


(ggplot(d_model, aes(x = num_neo, y = num_hyb)) +
  geom_point( alpha = 0.5, size = 2.5, pch = 21, fill = "grey70") +
  scale_x_log10() +
  scale_y_log10() +
  geom_line(inherit.aes = F,
            data = draws , 
            aes(y = 10^.epred, 
                x = 10^num_neo_log, 
                group = .draw), 
            col = "#33E4FFFF",
            alpha = .1) +
  geom_smooth(inherit.aes = F, 
              data = fit,
              aes(x = 10^num_neo_log, 
                  y = 10^Estimate, 
                  ymin = 10^Q2.5,  
                  ymax = 10^Q97.5), stat = "identity", se = FALSE, col = "#007A99FF") +
    theme_ipsum(
      axis_title_size = 12,
      plot_title_size = 12,
      axis_title_just = "mm",
      grid = "XY",
      base_family = "Helvetica"
    ) +
  theme(
    plot.margin = margin(0, 1,0,0, "cm"),
    legend.position = "none",
    aspect.ratio = 1,
    legend.title = element_blank() ,
    axis.line = element_line() ) +
  labs(y = "Number of hybrids", 
       x = "Number of neophytes") ->fig1)


ggsave(dpi = 600, height = 3.68, width = 3.68, bg = "white",
       "Figures/hybrid_neo_geography.png")




# percentage of hybrids versus percentage of neophytes, country area and taxonomic effort

# prepare data:
d_model_ratio <- d_ratios %>%
  select(tdwg3, hyb_per_nat, neo_per_nat, area, avg_gdp, author_weighted) %>% 
  na.omit %>% 
  mutate(hyb_per_nat_log = log10(hyb_per_nat),
         neo_per_nat_log = log10(neo_per_nat),
         area_log = log10(area),
         avg_gdp_log = log10(avg_gdp),
         author_log = log10(author_weighted))

# run model
b_mod_ratio <- brm(data = d_model_ratio, 
             family = gaussian,
             hyb_per_nat_log ~ neo_per_nat_log + area_log + author_log + avg_gdp_log,
             iter = 2000, warmup = 500, chains = 4, cores = 4, 
             file = "Models/geography_ratios.rds")

b_mod_ratio <- readRDS("Models/geography_ratios.rds")

# interpret result
(1.2^0.21 -1) * 100


nd_ratio = tibble(avg_gdp_log = mean(d_model_ratio$avg_gdp_log),
            area_log = mean(d_model_ratio$area_log),
            author_log = mean(d_model_ratio$author_log),
            neo_per_nat_log = seq(min(d_model_ratio$neo_per_nat_log), 
                              max(d_model_ratio$neo_per_nat_log), 
                              by = 0.1))

fit_ratio <- fitted(b_mod_ratio, newdata = nd_ratio) %>% 
  data.frame() %>% select(-Est.Error) %>% 
  bind_cols(nd_ratio)

draws_ratio <- d_model_ratio %>%
  modelr::data_grid(neo_per_nat_log = seq_range(neo_per_nat_log, n = 101),
                    avg_gdp_log = mean(d_model$avg_gdp_log),
                    area_log = mean(d_model$area_log),
                    author_log = mean(d_model$author_log)) %>%
  add_epred_draws(b_mod_ratio, ndraws = 100)


(ggplot(d_model_ratio, aes(x = neo_per_nat, y = hyb_per_nat)) +
  geom_point( alpha = 0.5, size = 2.5, pch = 21, fill = "grey70") +
  scale_x_log10() +
  scale_y_log10() +
  geom_line(inherit.aes = F,
            data = draws_ratio , 
            aes(y = 10^.epred, 
                x = 10^neo_per_nat_log, 
                group = .draw), 
            col = "#33E4FFFF",
            alpha = .1) +
  geom_smooth(inherit.aes = F, 
              data = fit_ratio,
              aes(x = 10^neo_per_nat_log, 
                  y = 10^Estimate, 
                  ymin = 10^Q2.5,  
                  ymax = 10^Q97.5), stat = "identity", se = FALSE, col = "#007A99FF") +
    theme_ipsum(
      axis_title_size = 12,
      plot_title_size = 12,
      axis_title_just = "mm",
      grid = "XY",
      base_family = "Helvetica"
    ) +
    theme(
      plot.margin = margin(0, 1, 0,0, "cm"),
      legend.position = "none",
      aspect.ratio = 1,
      legend.title = element_blank() ,
    axis.line = element_line() ) +
  labs(y = "Ratio of hybrids to natives", 
       x = "Ratio of neophytes to natives") -> fig2)

ggsave(dpi = 600, height = 3.68, width = 3.68, bg = "white",
       "Figures/hybrid_neo_ratio_geography.png")



# patchwork, combine all plots

p1 <- plot_grid(map_hybrids + theme(legend.justification = c(0,1)), 
          map_neophytes + theme(legend.justification = c(0,1)), 
          ncol = 1, labels = c("a", "b"), align = "v")

p2 <- plot_grid(fig1, fig2, ncol = 2, labels = c("c", "d"), align = "v")

plot_grid(p1, p2, nrow = 2)

ggsave(dpi = 600, height = 9.17, width = 7.53, bg = "white",
       "Figures/hybrid_neo_ratio_geography.png")
