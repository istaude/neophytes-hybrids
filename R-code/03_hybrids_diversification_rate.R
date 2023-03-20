source("R-code-clean/00_preamble.R")


# calculate diversification rates -----------------------------------------

nodes.info.1 %>% View
# calculate diversification rates
div_rates <- nodes.info.1 %>% 
  filter(level == "G") %>% 
  dplyr::select(family, genus, rn.bl)

# get data from kew: species richness per genera
kew_sp <- fread("Data/kew/checklist_names.txt")
kew_spn_per_gen <- kew_sp %>% 
  filter(taxon_status == "Accepted") %>% 
  filter(taxon_rank == "Species") %>% 
  filter(species_hybrid != "×") %>% 
  filter(species_hybrid != "+") %>%
  dplyr::select(taxon_name, genus) %>%
  distinct %>%
  count(genus) %>% 
  rename(nspp_in_genera = n)



# calculate speciation rates
div_rates <- left_join(div_rates, kew_spn_per_gen) %>% 
  mutate(div_rate_root_age_upper = log(nspp_in_genera*(1-0.5) + 0.5)/rn.bl,
         div_rate_root_age_lower = log(nspp_in_genera*(1-0.9) + 0.9)/rn.bl)


# identify which of these genera have hybrids
# flag genera that contain hybrid species
kew_hyb <- kew_sp %>% 
  filter(species_hybrid == "×" | species_hybrid == "+") %>% 
  filter(taxon_status == "Accepted") 
kew_hyb_gen <- kew_hyb %>% dplyr::select(genus) %>% distinct %>% mutate(hyb_present = "hybrid")

# same for neophytic species, identify which genera also produce neophytic species
# glonaf_higher_taxa_wcvp comes from the script 02_taxonomy..., needs to be loaded
glonaf_higher_taxa_wcvp <- fread("Data/glonaf_higher_taxa_wcvp.csv")

glonaf_neo_gen <- glonaf_higher_taxa_wcvp %>% dplyr::select(genus_acc) %>% distinct %>% 
  mutate(neo_present = "neophyte") %>% rename(genus = genus_acc)

# join
div_rates <- left_join(div_rates, kew_hyb_gen) %>% replace_na(list(hyb_present = "no hybrid"))
div_rates <- left_join(div_rates, glonaf_neo_gen) %>% replace_na(list(neo_present = "no neophyte"))


# filter genera that have no div rate information
div_rates <- div_rates %>% na.omit

# create grouping variable that indicates whether only hybrids and/or neophytes 
# are present
div_rates <-  div_rates %>% mutate(type = paste(hyb_present, neo_present, sep = " & "))


# full dataframe, including all info for analysis 4
write_csv(div_rates, "Data/hybrid_neophyte_analysis4.csv")

# filter extreme diversification rates that occur in super young lineages
# div_rates %>% arrange(desc(div_rate_root_age_lower )) %>% View
# div_rates <- div_rates %>% filter(div_rate_root_age_lower <= 20)

# avaiable for how many genera?
div_rates %>% nrow


div_rates <- div_rates %>% mutate(type = case_when(type == 'no hybrid & no neophyte' ~ 'none',
                                      type  == 'hybrid & no neophyte' ~ 'hybrids',
                                      type  == 'hybrid & neophyte' ~ 'hybrids and neophytes',
                                      type  == 'no hybrid & neophyte' ~ 'neophytes')) 
div_rates$type <-  fct_rev(factor(div_rates$type, levels = c("none",
                                                     "neophytes",
                                                     "hybrids", 
                                                     "hybrids and neophytes")))


# important to note that some genera have a speciation rate of 0, that happens
# when there is just one species in the genus

# otherwise the minimum diversificaiton rate is this here:
div_rates %>% filter(div_rate_root_age_lower != 0) %>% 
  summarize(min(div_rate_root_age_lower))

# first test: are genera carrying hybrids diversifying faster than those without?
# important to note, this plot excludes species with a diversificaiton rate of 0
# these species exclusively occur in the no hybrid categories
# also this plot excludes genera with extreme rates (above 20).
(ggplot(div_rates %>% filter(div_rate_root_age_lower <= 20), aes(x = type, y = div_rate_root_age_lower, fill = type)) + 
    geom_violinhalf(position = position_nudge(x = 0.2, y = 0)) +
    geom_point(aes(col = type), position = position_jitter(width = 0.1), alpha =0.8, size =0.1)   +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.2,
      alpha = 0.4
    ) +
    stat_n_text(size = 3, fontface = 3, y.pos = 1, family = "Helvetica") +
    scale_fill_manual(values = c(`none` = "#E6E6E6" , 
                                 `neophytes` = "#CCFFFF" , 
                                 `hybrids` = "#66E6FF" , 
                                 `hybrids and neophytes` = "#0099CC" )) +
    scale_color_manual(values = c(`none` = "#E6E6E6" , 
                                  `neophytes` = "#CCFFFF" , 
                                  `hybrids` = "#66E6FF" , 
                                  `hybrids and neophytes` = "#0099CC" )) +
    coord_flip()+
    theme_ipsum(
      grid = "",
      axis_title_size = 12,
      plot_title_size = 12,
      axis_title_just = "mm",
      base_family = "Helvetica"
    ) +
    scale_y_log10()+
    theme(
      text = element_text(size = 14),
      plot.margin=unit(rep(0.2,4),"cm"),
      legend.position = "none",
      axis.text.y = element_text(vjust = 0, hjust = 0, face = 3),
      legend.title = element_blank() ,
      axis.line.x = element_line() ) +
    labs(y = "Diversification rate (species per Myr)", 
         x = "") -> fig4a)


# statistical model
d_model_div <- div_rates %>% 
  select(div_rate_root_age_lower, type, nspp_in_genera) %>% 
  filter(div_rate_root_age_lower != 0) %>% 
  mutate(div_rate_root_age_lower_log = log10(div_rate_root_age_lower))


# for quick in text interpretation
modlm_div <- lm(data = d_model_div, 
    log10(div_rate_root_age_lower) ~ type - 1)
emmeans(modlm_div, "type", type = "response") %>% 
  contrast(method = "pairwise")


# for plotting
b_mod_div <- brm(data = d_model_div, 
             family = gaussian,
             div_rate_root_age_lower_log ~ type,
             iter = 2000, warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
             file = "Models/diversification_rate.rds")

b_mod_div <- readRDS("Models/diversification_rate.rds")

# plot means and posterior
(b_mod_div %>% 
  emmeans(~ type,
          epred = TRUE) %>% 
  gather_emmeans_draws() %>% 
ggplot(aes(x = 10^.value, y = type, fill = type)) +
  stat_gradientinterval(alpha =0.8, height =0.5) + 
    scale_fill_manual(values = c(`none` = "#E6E6E6" , 
                                 `neophytes` = "#CCFFFF" , 
                                 `hybrids` = "#66E6FF" , 
                                 `hybrids and neophytes` = "#0099CC" )) +
    theme_ipsum(
      grid = "",
      axis_title_size = 12,
      plot_title_size = 12,
      axis_title_just = "mm",
      base_family = "Helvetica"
    ) +
    scale_x_log10(breaks = scales::pretty_breaks(n = 5))+
    theme(
      text = element_text(size = 14),
      plot.margin=unit(rep(0.2,4),"cm"),
      legend.position = "none",
      axis.text.y = element_blank(),
      legend.title = element_blank() ,
      axis.line.x = element_line() ) +
  labs(x = "Diversification rate", 
       y = "") -> fig4b )


# lines are the 66% (thick), and 95% (slim) credible interval
# plot contrasts
(b_mod_div %>% 
  emmeans(~ type,
          epred = TRUE) %>% 
  contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
    ggplot(aes(x = 10^.value,
               y = fct_rev(factor(contrast, levels = c(
                 "hybrids and neophytes - none",
                 "hybrids and neophytes - neophytes",
                 "hybrids and neophytes - hybrids",
                 "hybrids - none",
                 "hybrids - neophytes",
                 "neophytes - none")))
    )) + 
    stat_gradientinterval(height =0.5, show.legend = F, fill = "#666666") +
    scale_x_log10()+
    geom_vline(xintercept = 1, 
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
    labs(x = "Contrasts (ratio) of diversification rates", 
         y = "") -> fig4c)

plot_grid(
plot_grid(fig4a, fig4b, nrow = 1, labels = c("a", "b"), rel_widths = c(1, 0.5)),
fig4c, nrow = 2, rel_heights = c(1, 0.7), labels = c("", "c"))


ggsave(dpi = 600, height = 5.38, width = 7.70, bg = "white",
       "Figures/diversification_rate_numbers_contrast.png")







# supplementary info ------------------------------------------------------

# including diversification rate of 0 and extreme trends
(ggplot(div_rates, aes(x = type, y = div_rate_root_age_lower + 0.01)) + 
    geom_violinhalf(position = position_nudge(x = 0.2, y = 0), fill = "#31cb69", alpha =0.4) +
    geom_point(position = position_jitter(width
                                          = 0.1), col = "#31cb69", alpha =0.1)   +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.2,
      alpha = 0.4,
      fill = "#31cb69"
    ) +
    stat_summary(fun=mean, geom="point", shape=21, size=3, 
                 color="black", fill="#cb3193") +
    coord_flip()+
    theme_ipsum(
      grid = "",
      axis_title_size = 12,
      plot_title_size = 12,
      axis_title_just = "mm",
      base_family = "Helvetica"
    ) +
    scale_y_log10()+
    theme(
      text = element_text(size = 14),
      plot.margin=unit(rep(0.2,4),"cm"),
      legend.position = "none",
      axis.text.y = element_text(vjust = 0, hjust = 0, face = 3),
      legend.title = element_blank() ,
      axis.line.x = element_line() ) +
    labs(y = "Diversification rate (species per Myr)", 
         x = "",
         title = "Including Diversification rates of\nzero and extreme ones (>20 spp. Myr)") -> fig_supp3a)


# including diversification rate of 0 and extreme trends
(ggplot(div_rates %>% filter(div_rate_root_age_upper <=20), aes(x = type, y = div_rate_root_age_upper )) + 
    geom_violinhalf(position = position_nudge(x = 0.2, y = 0), fill = "#31cb69", alpha =0.4) +
    geom_point(position = position_jitter(width
                                          = 0.1), col = "#31cb69", alpha =0.1)   +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.2,
      alpha = 0.4,
      fill = "#31cb69"
    ) +
    stat_summary(fun=mean, geom="point", shape=21, size=3, 
                 color="black", fill="#cb3193") +
    coord_flip()+
    theme_ipsum(
      grid = "",
      axis_title_size = 12,
      plot_title_size = 12,
      axis_title_just = "mm",
      base_family = "Helvetica"
    ) +
    scale_y_log10()+
    theme(
      text = element_text(size = 14),
      plot.margin=unit(rep(0.2,4),"cm"),
      legend.position = "none",
      axis.text.y = element_text(vjust = 0, hjust = 0, face = 3),
      legend.title = element_blank() ,
      axis.line.x = element_line() ) +
    labs(y = "Diversification rate (species per Myr)", 
         x = "",
         title = "Data as in main text,\nbut with e = 0.5") -> fig_supp3b)


mod_fig_supp3a <- lm(log10(div_rate_root_age_lower + 0.01) ~ type, 
   data = div_rates )


mod_fig_supp3b <- lm(log10(div_rate_root_age_upper) ~ type, 
   data = div_rates %>% 
     filter(div_rate_root_age_upper <=20) %>% 
     filter(div_rate_root_age_upper != 0))

emmeans(mod_fig_supp3a, "type", type = "response") %>% plot + 
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
  ylab("") + xlab("Mean diversification rate")-> fig_supp3c

emmeans(mod_fig_supp3b, "type", type = "response") %>% plot +
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
    axis.line.x = element_line() )+
  ylab("") + xlab("Mean diversification rate") -> fig_supp3d

ggarrange(
ggarrange(fig_supp3a, fig_supp3b, ncol = 2, labels = c("a", "b")),
ggarrange(fig_supp3c, fig_supp3d, ncol = 2, labels = c("", "")), nrow = 2)

ggsave(dpi = 600, height = 6.05, width = 9.05, bg = "white",
       "Figures/diversification_rate_supp.png")



# if we control for species richness in the model.
b_mod_div2 <- brm(data = d_model_div, 
                  family = gaussian,
                  div_rate_root_age_lower_log ~ type + nspp_in_genera,
                  iter = 2000, warmup = 500, chains = 4, cores = 4, backend = "cmdstanr")

b_mod_div2 %>% 
  emmeans(~ type,
          epred = TRUE) %>% 
  gather_emmeans_draws() %>% 
  ggplot(aes(x = 10^.value, y = type)) +
  stat_gradientinterval(fill = "#cb3193", alpha =0.4, height =0.5) + 
  theme_ipsum(
    grid = "X",
    axis_title_size = 12,
    plot_title_size = 12,
    axis_title_just = "mm"
  ) +
  scale_x_log10() +
  theme(
    text = element_text(size = 14),
    plot.margin=unit(rep(0.2,4),"cm"),
    legend.position = "none",
    legend.title = element_blank() ,
    axis.line = element_line() ) +
  labs(x = "Diversification rate (species per Myr)", 
       y = "")
# similar pattern, but questionable whether controlling for species richness 
# makes so much sense because it goes into the calculation of diversification rate
# anyway, this is at least robust: those genera with hybrids tend to have a higher
# diversification rate.
























# clade age ---------------------------------------------------------------
ggplot(div_rates, aes(x = type, y = rn.bl)) + 
  geom_violinhalf(position = position_nudge(x = 0.2, y = 0), fill = "#31cb69", alpha =0.4) +
  geom_point(position = position_jitter(width = 0.1), col = "#31cb69", alpha =0.1)   +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.2,
    alpha = 0.4,
    fill = "#31cb69"
  ) +
  stat_summary(fun=mean, geom="point", shape=21, size=3, 
               color="black", fill="#cb3193") +
  stat_n_text(size = 3, fontface = 3, hjust =0.1) +
  coord_flip()+
  theme_ipsum(
    grid = "X",
    axis_title_size = 12,
    plot_title_size = 12,
    axis_title_just = "mm"
  ) +
  scale_y_log10()+
  theme(
    text = element_text(size = 14),
    plot.margin=unit(rep(0.2,4),"cm"),
    legend.position = "none",
    legend.title = element_blank() ,
    axis.line = element_line() ) +
  labs(y = "Clade age (in Myr)", 
       x = "")


# -------------------------------------------------------------------------
# second test, are genera that produce many hybrids also diversifying faster?
div_rates_hyb <- left_join(kew_hyb_higher_taxa %>% count(genus) %>% rename(n_hyb = n),
          div_rates) %>% na.omit

# remove species with 0 as specciation rate
div_rates_hyb %>% filter(div_rate_root_age_lower == 0)
div_rates_hyb <- div_rates_hyb %>% filter(div_rate_root_age_lower != 0)

# unclear whether this is just an artifact, fit quantile regression: 
# this is as a continuous analogue to geom_boxplot().
div_rates_hyb %>% 
ggplot(aes(y = div_rate_root_age_lower, x = n_hyb)) +
  geom_point(alpha = 0.5, size = 2.4, pch = 21, fill = "magenta") +
  geom_quantile(quantiles = c(0.1, 0.5, 0.9), size = 2, aes(alpha = ..quantile..)) +
  scale_alpha(range = c(0.3, 0.7)) +
  #geom_smooth(col = "black", fill = "grey10", alpha = 0.2, method = "loess") +
  theme_ipsum(
    grid = "",
    axis_title_size = 12,
    plot_title_size = 12,
    axis_title_just = "mm"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  theme(
    text = element_text(size = 14),
    plot.margin=unit(rep(0.2,4),"cm"),
    legend.position = "none",
    aspect.ratio = 1,
    legend.title = element_blank() ,
    axis.line = element_line() ) +
  labs(y = "Diversification rate", x = "Number of hybrids in genus")


# statistical models
b_mod_div3_01 <- brm(bf(log10(div_rate_root_age_lower) ~ log10(n_hyb), quantile = 0.1), 
                     data = div_rates_hyb, 
                     family = asym_laplace(), 
                     iter = 2000, warmup = 500, chains = 4, cores = 4, 
                     backend = "cmdstanr")

b_mod_div3_05 <- update(b_mod_div3_01, 
                        formula = bf(log10(div_rate_root_age_lower) ~ log10(n_hyb), 
                                     quantile = 0.5))

b_mod_div3_09 <- update(b_mod_div3_01, 
                        formula = bf(log10(div_rate_root_age_lower) ~ log10(n_hyb), 
                                     quantile = 0.9))


# fitted regression lines for the different quantiles
fitted(b_mod_div3_01, dpar="mu") %>% 
  data.frame() %>% select(-Est.Error) %>% 
  rename(estimate_01 = Estimate, Q2.5_01 = Q2.5, Q97.5_01 = Q97.5) %>% 
  bind_cols(div_rates_hyb) %>% 
  bind_cols(
  fitted(b_mod_div3_05, dpar="mu") %>% 
  data.frame() %>% select(-Est.Error) %>%
  rename(estimate_05 = Estimate, Q2.5_05 = Q2.5, Q97.5_05 = Q97.5)) %>% 
  bind_cols(
  fitted(b_mod_div3_09, dpar="mu") %>% 
  data.frame() %>% select(-Est.Error) %>%
    rename(estimate_09 = Estimate, Q2.5_09 = Q2.5, Q97.5_09 = Q97.5)) %>% 

# plot commands  
  ggplot(aes(x = n_hyb)) +
  geom_point(aes(y = div_rate_root_age_lower), 
             alpha = 0.4, size = 2.4, pch = 21, fill = "#31cb69") +
  geom_line(aes(y = 10^estimate_01), lwd = 1, color = "black", alpha = .2) +
  geom_line(aes(y = 10^estimate_05), lwd = 1, color = "black", alpha = .5) +
  geom_line(aes(y = 10^estimate_09), lwd = 1, color = "black", alpha = .8) +
  
#  geom_smooth(aes(y = 10^estimate_01, ymin = 10^Q2.5_01, ymax = 10^Q97.5_01),
#              stat = "identity",
#              fill = "black", color = "black", alpha = 1/5) +
#  geom_smooth(aes(y = 10^estimate_05, ymin = 10^Q2.5_05, ymax = 10^Q97.5_05),
#             stat = "identity",
#              fill = "black", color = "black", alpha = 1/5) +
#  geom_smooth(aes(y = 10^estimate_09, ymin = 10^Q2.5_09, ymax = 10^Q97.5_09),
#              stat = "identity",
#              fill = "black", color = "black", alpha = 1/5) +
  theme_ipsum(
    axis_title_size = 12,
    plot_title_size = 12,
    axis_title_just = "mm"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  theme(
    text = element_text(size = 14),
    plot.margin=unit(rep(0.2,4),"cm"),
    legend.position = "none",
    aspect.ratio = 1,
    legend.title = element_blank() ,
    axis.line = element_line() ) +
  labs(y = "Diversification rate", x = "Number of hybrids in genus")

ggsave(dpi = 600, height = 3.68, width = 3.68, bg = "white",
       "Figures/diversification_rate_numbers.png")


# same for neophytic species?




