source("R-code-clean/00_preamble.R")


# aim ---------------------------------------------------------------------

# compare distribution of neophyte native and hybrid species in the UK in
# relation to human footprint data-


# load dataset ----------------------------------------------------------

hyb_neo_nat_preston_complete <- read_excel("Data/hyb_neo_nat_preston_complete.xlsx")


# bring into long format
hyb_long <- bind_rows(
  hyb_neo_nat_preston_complete %>% 
    select(species_name_submitted, key_in_bsbi, taxonId.y, parents, name)%>% 
    rename(species_name = name) %>% 
    rename(taxonId = taxonId.y) %>% 
    rename(role = parents) %>% 
    distinct 
  ,
  hyb_neo_nat_preston_complete %>% 
    select(species_name_submitted, key_in_bsbi, taxonId.x, hybrid_species) %>% 
    mutate(role = "hybrid") %>% 
    distinct %>% 
    rename(species_name = hybrid_species) %>% 
    rename(taxonId = taxonId.x)
) %>% arrange(species_name_submitted)



# attach native status ----------------------------------------------------

# get native status of parents.
native_status <- read_delim("Data/results20220503144357.csv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

hyb_long <- left_join(hyb_long, native_status %>% 
                            select(taxonId, dataValue) %>% 
                            rename(native_status = dataValue)) 

# get different dataset to indicate native status
native_status_v2 <- read_csv("Data/BI_main.csv")

hyb_long <- left_join(hyb_long, native_status_v2 %>% 
                            select(taxon_name_binom, StaceIV_nativity)
                          , by = c("species_name" = "taxon_name_binom")) 


# fill rows amongst datasets
hyb_long <- hyb_long %>% 
  mutate(native_combined = ifelse(is.na(StaceIV_nativity) == F,
                                  StaceIV_nativity, native_status)) 
# recode native status
hyb_long$native_combined[hyb_long$native_combined == "?N"]  <- "native" 
hyb_long$native_combined[hyb_long$native_combined == "AN"]  <- "neophyte" 
hyb_long$native_combined[hyb_long$native_combined == "AR"]  <- "native"
hyb_long$native_combined[hyb_long$native_combined == "Arch-colonist"]  <- "native"
hyb_long$native_combined[hyb_long$native_combined == "Arch-cultd"]  <- "native"
hyb_long$native_combined[hyb_long$native_combined == "Arch-denizen"]  <- "native"
hyb_long$native_combined[hyb_long$native_combined == "N"]  <- "native"
hyb_long$native_combined[hyb_long$native_combined == "NE"]  <- "native"
hyb_long$native_combined[hyb_long$native_combined == "Neo-casual"]  <- "neophyte"
hyb_long$native_combined[hyb_long$native_combined == "Neo-natd"]  <- "neophyte"
hyb_long$native_combined[hyb_long$native_combined == "Neo-surv"]  <- "neophyte"
hyb_long$native_combined[hyb_long$native_combined == "Neonative"]  <- "neophyte"

hyb_long <- hyb_long %>% 
  select(-StaceIV_nativity, -native_status) %>% 
  distinct %>% 
  mutate(native_combined = ifelse(role == "hybrid", "hybrid", native_combined))
  
# fill in missing species
hyb_long <- hyb_long %>% 
  mutate(native_combined = ifelse(species_name == "Sorbus scalaria", "neophyte", native_combined)) %>% 
  mutate(native_combined = ifelse(species_name == "Conyza floribunda", "neophyte", native_combined)) %>% 
  mutate(native_combined = ifelse(species_name == "Euphrasia stricta", "neophyte", native_combined)) %>% 
  mutate(native_combined = ifelse(species_name == "Pilosella aurantiaca subsp. carpathicola", "neophyte", native_combined)) %>% 
  mutate(native_combined = ifelse(species_name == "Rumex frutescens", "neophyte", native_combined))
  
#write.xlsx(hyb_long, "hyb_long.xlsx", sheetName="Sheet1", showNA=FALSE)

# how many species in prestons list?
hyb_long %>% count(role)

# attach distribution -----------------------------------------------------

# load inventory data from UK
bsbi <- read_delim("Data/results20220503095052.csv", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

hyb_long <- left_join(hyb_long, bsbi %>% select(taxonId, dataValue)) 


# find human footprint ----------------------------------------------------

hyb_long <- hyb_long %>% 
  separate_rows(dataValue, sep = ",")

# load raster
hf_raster = rast("Data/wildareas-v3-2009-human-footprint.tif")
crs(hf_raster)

# load shape file of national grid
grid <- st_read("Data/OS-British-National-Grids-main/os_bng_grids.gpkg", 
                layer = "10km_grid")

grid <- st_transform(grid, crs(hf_raster))
crs(grid)

# crop raster to extent of polygons
hf_raster <- hf_raster %>% 
  terra::crop(grid)

# extract raster values
hf_per_tile <- terra::extract(hf_raster, vect(grid))  

# 128 is na value, replace with na
hf_mean <- 
  hf_per_tile %>% 
  mutate(hf  = ifelse(`wildareas-v3-2009-human-footprint` > 50, 
                      NA, `wildareas-v3-2009-human-footprint`)) %>% 
  group_by(ID) %>% 
  summarize(hf = mean(hf, na.rm = T))


# join tile averages for hf back to grid data frame, where tile names are indicated
hf_mean_sf <- 
  # back to sf
  st_as_sf(vect(grid)) %>% 
  # define id
  mutate(ID := seq_len(nrow(.))) %>% 
  # merge by id
  left_join(., hf_mean, by = "ID")


# now bring this data together with the species distribution data
hyb_long <- left_join(
  hyb_long,
  hf_mean_sf %>% st_drop_geometry(),
  by = c("dataValue" = "tile_name"))


# plot --------------------------------------------------------------------

# many hybrids probably have a small sample size (range), perhaps
# set a threshold
hyb_long <- full_join(hyb_long, 
                      hyb_long %>% 
                        group_by(species_name) %>% 
                        summarize(n = n_distinct(dataValue)))


hyb_long_threshold <- hyb_long %>% 
  group_by(species_name_submitted) %>% 
  filter(min(n, na.rm = TRUE) > 10)
  
hyb_long_threshold %>% select(species_name, native_combined, n) %>% 
  distinct %>% group_by(native_combined) %>% 
  summarise(median(n))

# this reduces the number considerably
length(unique(hyb_long_threshold$species_name_submitted))


# wrap species name to fit facet strips
hyb_long_threshold <- hyb_long_threshold %>% 
  mutate(clean_spec_name = str_replace(species_name_submitted, " \\s*\\([^\\)]+\\)", "")) %>% 
  mutate(clean_spec_name = str_wrap(clean_spec_name, 25))

# nas come from Ireland, for which I could not find a grid to match the hf to
# that means our analysis looks at the distribution of these taxa across
# Great Britain
hyb_long_threshold %>% filter(is.na(hf)) %>% View

(pp1 <- hyb_long_threshold %>% 
  ggplot(aes(y = fct_rev(native_combined), x = hf, color = native_combined, fill = native_combined)) +
  facet_wrap(.~ clean_spec_name) +
  ggridges::geom_density_ridges(
    quantile_lines = TRUE, quantiles = 2, 
    color = "black", alpha = .8
  ) +
  labs(x = "Human footprint index", y = "") +
  theme_ipsum(
    grid = "",
    axis_title_size = 12,
    plot_title_size = 12,
    axis_title_just = "mm",
    base_family = "Helvetica"
  ) +
  coord_cartesian(clip = "off") +
  theme(
    strip.text = element_text(face = "italic", size = 8),
    text = element_text(size = 14),
    plot.margin=unit(rep(0.2,4),"cm"),
    panel.spacing = unit(0.2, "lines"),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank() ,
    legend.key.height= unit(0.3, 'cm'),
    legend.key.width= unit(0.3, 'cm'),
    axis.line = element_line() ) +
  scale_fill_manual(values = c(hybrid = "#33e4ff" , native = "#007A99", neophyte = "#E6E6E6")))

shift_legend2 <- function(p) {
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  reposition_legend(p, 'center', panel=names)
}

pp1 <- shift_legend2(pp1)


# statistical model -------------------------------------------------------

b_mod_parent <- brm(hf ~ native_combined + (1|species_name_submitted), 
            data = hyb_long_threshold,
            family = gaussian,
            iter = 2000, warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
            file = "Models/triplet.RDS")

b_mod_parent <- readRDS("Models/triplet.RDS")

b_mod_parent

# plot means and posterior
b_mod_parent %>% 
  emmeans(~ native_combined,
          epred = TRUE) %>% 
  gather_emmeans_draws() %>% 
  ggplot(aes(x = .value, y = fct_rev(native_combined))) +
  stat_gradientinterval(aes(fill = native_combined), alpha =0.4, height =0.5) + 
  scale_fill_manual(values = c(hybrid = "#33e4ff" , native = "#007A99", neophyte = "#E6E6E6")) +
  theme_ipsum(
    grid = "",
    axis_title_size = 12,
    plot_title_size = 12,
    axis_title_just = "mm",
    base_family = "Helvetica"
  ) +
  theme(
    text = element_text(size = 14),
    plot.margin=unit(rep(0.2,4),"cm"),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line() ) +
  labs(x = "Human footprint index", 
       y = "") -> pp2


# plot contrasts
(b_mod_parent %>% 
  emmeans(~ native_combined,
          epred = TRUE) %>% 
    contrast(method = "pairwise") %>% 
    gather_emmeans_draws() %>% 
    ggplot(aes(x = .value,
               y = fct_rev(contrast)
               )) + 
    stat_gradientinterval(height =0.5, show.legend = F, fill = "#666666") +
    geom_vline(xintercept = 0, 
               linetype = 2) + 
    theme_ipsum(
      axis_title_size = 12,
      plot_title_size = 12,
      axis_title_just = "mm",
      grid = "",
      base_family = "Helvetica"
    ) +
    theme(
      text = element_text(size = 14),
      axis.text.y = element_text(vjust = 0, hjust = 0, face = 3),
      plot.margin=unit(rep(0.2,4),"cm"),
      legend.position = "none",
      legend.title = element_blank() ,
      axis.line.x = element_line() ) +
    labs(x = "Difference in human foot print index", 
         y = "") -> pp3)


# combine
plot_grid(pp1, pp2, pp3, nrow = 3, rel_heights = c(1, 0.3, 0.3) , labels = c("a", "b", "c"))
ggsave(dpi = 600, height = 6.85, width = 9.50, bg = "white",
       "Figures/hybrid_neo_parent_dis.png")



# how many hybrids have HFI values greater than those ever sampled for natives
# / neophytes
hyb_long_threshold %>% 
  group_by(clean_spec_name, species_name) %>% 
  summarize(role = first(native_combined),
            maxi = min(hf, na.rm = T),
            triplet = first(clean_spec_name)) %>% 
  group_by(clean_spec_name) %>% 
  arrange(clean_spec_name, desc(maxi)) %>% 
  mutate(rank=rank(maxi)) %>% View
  filter(role == "hybrid") %>% View

