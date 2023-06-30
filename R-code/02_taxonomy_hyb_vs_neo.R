source("R-code/00_preamble.R")

# aim ---------------------------------------------------------------------
# check taxonomic overlap between neophytes and hybrids


# hybrids -----------------------------------------------------------------
kew_sp <- fread("Data/kew/checklist_names.txt", encoding = "UTF-8")

# kew species data
kew_hyb <- kew_sp %>% filter(species_hybrid == "×" | species_hybrid == "+")

# how many of these hybrids are accepted
kew_hyb %>% count(taxon_status)

# let's take only accepted hybrids
kew_hyb <- kew_hyb %>% filter(taxon_status == "Accepted") 

kew_hyb_higher_taxa <- kew_hyb %>% dplyr::select(plant_name_id, family, genus, taxon_name)


# neophytes ---------------------------------------------------------------
glonaf_sp <- read_delim("Data/GLONAF/Taxon_x_List_GloNAF_vanKleunenetal2018Ecology.csv", 
                        delim = "\t", escape_double = FALSE, 
                        locale = locale(encoding = "utf-16"), 
                        trim_ws = TRUE)

# extract genus
glonaf_higher_taxa <- glonaf_sp %>% 
  filter(hybrid == 0) %>%   
  filter(status == "naturalized") %>% 
  filter(name_status == "accepted") %>% 
  dplyr::select(standardized_name, family_tpl) %>% 
  distinct %>% 
  mutate(genus = sapply(strsplit(standardized_name, " "), "[", 1))

# genus names of glonaf may not match the accepted taxonomy of WCVP
# for example they list Spartina species but in fact they belong to
# the genus Sporobolus
glonaf_higher_taxa_wcvp <- left_join(
  glonaf_higher_taxa,
  kew_sp %>% 
    filter(taxon_rank == "Genus") %>% 
    dplyr::select(genus, taxon_status, accepted_plant_name_id) ) 

# now if a genus has more than 1 match...
glonaf_higher_taxa_wcvp %>% dplyr::select(genus, accepted_plant_name_id) %>% distinct %>% 
  count(genus) %>% arrange(desc(n))

# pick the accepted genus
glonaf_higher_taxa_wcvp <- bind_rows(
  left_join(glonaf_higher_taxa_wcvp,
            glonaf_higher_taxa_wcvp %>% dplyr::select(genus, accepted_plant_name_id) %>% distinct %>% 
              count(genus)) %>% filter(n == 1), 
  
  left_join(glonaf_higher_taxa_wcvp,
            glonaf_higher_taxa_wcvp %>% dplyr::select(genus, accepted_plant_name_id) %>% distinct %>% 
              count(genus)) %>% filter(n > 1) %>% filter(taxon_status == "Accepted") 
) %>% dplyr::select(-n)


# now retrieve the correct genus name
glonaf_higher_taxa_wcvp <- left_join(
  glonaf_higher_taxa_wcvp,
  kew_sp %>% 
    filter(taxon_rank == "Genus") %>%
    filter(taxon_status == "Accepted") %>% 
    dplyr::select(plant_name_id, genus) %>% 
    distinct %>% 
    rename(genus_acc = genus), by = c("accepted_plant_name_id" = "plant_name_id"))

# checks
glonaf_higher_taxa_wcvp %>% dplyr::select(accepted_plant_name_id, genus_acc) %>% distinct %>% 
  count(genus_acc) %>% arrange(desc(n))

# which species have nas
glonaf_higher_taxa_wcvp %>% filter(is.na(genus_acc))
glonaf_higher_taxa_wcvp <- glonaf_higher_taxa_wcvp %>% filter(!is.na(genus_acc))



write_csv(glonaf_higher_taxa_wcvp, "Data/glonaf_higher_taxa_wcvp.csv")
# join --------------------------------------------------------------------
tax_overlap_gen <- full_join(kew_hyb_higher_taxa %>% count(genus) %>% rename(n_hyb = n),
                             glonaf_higher_taxa_wcvp %>% count(genus_acc) %>% rename(n_neo = n), 
                             by = c("genus" = "genus_acc"))



# analyse -----------------------------------------------------------------

# first of all, how many genera that have neophytic species are also forming hybrids

# 1. number of genera in hybrids
tax_overlap_gen %>% filter(n_hyb >= 1) %>% dplyr::select(genus) %>% distinct %>% nrow

# 2. number of genera in neophytes
tax_overlap_gen %>% filter(n_neo >= 1) %>% dplyr::select(genus) %>% distinct %>% nrow

# 2.b how many plant genera are there in total
kew_sp %>% filter(taxon_rank == "Genus", 
                  taxon_status == "Accepted") %>% dplyr::select(genus) %>% distinct %>% nrow

# 2.c what is the percentage of genera that produces neophytes
tax_overlap_gen %>% filter(n_neo >= 1) %>% dplyr::select(genus) %>% distinct %>% nrow/
  kew_sp %>% filter(taxon_status == "Accepted") %>% dplyr::select(genus) %>% distinct %>% nrow

# 3. number of genera that have both
tax_overlap_gen %>% na.omit() %>% dplyr::select(genus) %>% distinct %>% nrow

# 4. what percentage of hybrid genera is also neophyte genera 
tax_overlap_gen %>% na.omit() %>% dplyr::select(genus) %>% distinct %>% nrow/
  tax_overlap_gen %>% filter(n_hyb >= 1) %>% dplyr::select(genus) %>% distinct %>% nrow

write_csv(tax_overlap_gen, "Data/tax_overlap_gen.csv")
# 74% of all genera that produce hybrids are also genera that produce neophytes

# as seen previously only 644 genera with neophytes also include hybrids,
# for what percentage of neophyte species do these genera account, and for how many
# naturalization events
tax_overlap_gen <- fread("Data/tax_overlap_gen.csv")

# how many neophyte species do these genera represent
tax_overlap_gen %>% na.omit() %>% 
  left_join(glonaf_higher_taxa_wcvp %>% 
              dplyr::select(-genus), by = c("genus" = "genus_acc")) %>% 
  nrow
glonaf_higher_taxa_wcvp %>% nrow

# how many naturalizations do these genera present
glonaf_kew_dis_comp <- fread("Data/glonaf_kew_dis_comp.csv")
relevant_species <- tax_overlap_gen %>% na.omit() %>% 
  left_join(glonaf_higher_taxa_wcvp %>% 
              dplyr::select(-genus), by = c("genus" = "genus_acc")) %>% 
  dplyr::select(standardized_name)

left_join(relevant_species, glonaf_kew_dis_comp) %>% nrow /
glonaf_kew_dis_comp %>% nrow

left_join(relevant_species, glonaf_kew_dis_comp) %>% na.omit %>% nrow /
  left_join(relevant_species, glonaf_kew_dis_comp) %>% nrow
# -------------------------------------------------------------------------

# second, in genera that produce many neophytes, are there also more hybrids
# we need to correct for the number of species per genera
kew_sp %>% dplyr::select(taxon_rank) %>%  distinct
# how many accepted species do we find per genus
kew_spn_per_gen <- kew_sp %>% 
  filter(taxon_status == "Accepted") %>% 
  filter(taxon_rank == "Species") %>% 
  filter(species_hybrid != "×") %>% 
  filter(species_hybrid != "+") %>%
  dplyr::select(taxon_name, genus) %>%
  distinct %>%
  count(genus) %>% 
  rename(nspp_in_genera = n)

# join to hybrid neophyte df
tax_overlap_gen <- left_join(tax_overlap_gen, kew_spn_per_gen)

# full data frame for analysis 2, including genus name, number of species, 
# number of hybrids and number of neophytes
write_csv(tax_overlap_gen, "Data/hybrid_neophytes_analysis2.csv")

# now we focus on genera that produce both hybrids and neophytic species.
# we are not interested in genera that do not produce hybrids, these can also
# have a high number of neophytic species. 
d_tax <- tax_overlap_gen %>%  
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>% 
  filter(n_hyb != 0) %>%
  filter(n_neo != 0) %>% 
  mutate(ratio_neo = n_neo / nspp_in_genera, ratio_hyb = n_hyb / nspp_in_genera)

#d_tax %>% filter(ratio_neo > 1)
#kew_sp %>% filter(genus == "Robinia") %>% filter(taxon_status == "Accepted") %>% View
#glonaf_higher_taxa_wcvp %>% filter(genus == "Robinia")
# level to 1
#d_tax <- d_tax %>% mutate(ratio_neo = ifelse(ratio_neo > 1, 1, ratio_neo))

# statistical model -------------------------------------------------------

# ratio
d_tax$ratio_hyb_log <- log10(d_tax$ratio_hyb)
d_tax$ratio_neo_log <- log10(d_tax$ratio_neo)
b_mod_tax <- brm(data = d_tax, 
             family = gaussian,
             ratio_hyb_log ~ ratio_neo_log,
             iter = 2000, warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
             file = "Models/taxonomy_overlap_ratio.rds")

b_mod_tax <- readRDS("Models/taxonomy_overlap_ratio.rds")
b_mod_tax 


# interpret the coefficient: 
(1.2^0.63 -1) * 100

fit_tax <- fitted(b_mod_tax) %>% 
  data.frame() %>% 
  bind_cols(d_tax)

draws_tax <- d_tax %>%
  modelr::data_grid(ratio_neo_log = modelr::seq_range(ratio_neo_log, n = 101)) %>%
  add_epred_draws(b_mod_tax, ndraws = 100)


ggplot(d_tax, aes(x = ratio_neo, y = ratio_hyb)) +
    geom_point(alpha = 0.3, size = 2, pch = 21, fill = "grey70") +
    scale_x_log10() +
    scale_y_log10() +
    geom_line(inherit.aes = F,
              data = draws_tax , 
              aes(y = 10^.epred, 
                  x = 10^ratio_neo_log, 
                  group = .draw), 
              col = "#33E4FFFF",
              alpha = .1) +
    geom_smooth(inherit.aes = F, 
                data = fit_tax,
                aes(x = 10^ratio_neo_log, 
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
      plot.margin = margin(1, 1,0,0, "cm"),
      legend.position = "none",
      legend.title = element_blank() ,
      axis.line = element_line() ) +
    labs(y = "Hybrids to sp. in genus", 
         x = "Neophytes to sp. in genus")  -> pa

# numbers
d_tax$n_hyb_log <- log10(d_tax$n_hyb)
d_tax$n_neo_log <- log10(d_tax$n_neo)

b_mod_tax_n <- brm(data = d_tax, 
                 family = gaussian,
                 n_hyb_log ~ n_neo_log,
                 iter = 2000, warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
                 file = "Models/taxonomy_overlap_numbers.rds")

b_mod_tax_n <- readRDS("Models/taxonomy_overlap_numbers.rds")
b_mod_tax_n


# interpret the coefficient: 
(1.2^0.63 -1) * 100

fit_tax_n <- fitted(b_mod_tax_n) %>% 
  data.frame() %>% 
  bind_cols(d_tax)

draws_tax_n <- d_tax %>%
  modelr::data_grid(n_neo_log = modelr::seq_range(n_neo_log, n = 101)) %>%
  add_epred_draws(b_mod_tax_n, ndraws = 100)


ggplot(d_tax, aes(x = n_neo, y = n_hyb)) +
  geom_point(alpha = 0.3, size = 2, pch = 21, fill = "grey70") +
  scale_x_log10() +
  scale_y_log10() +
  geom_line(inherit.aes = F,
            data = draws_tax_n , 
            aes(y = 10^.epred, 
                x = 10^n_neo_log, 
                group = .draw), 
            col = "#33E4FFFF",
            alpha = .1) +
  geom_smooth(inherit.aes = F, 
              data = fit_tax_n,
              aes(x = 10^n_neo_log, 
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
    plot.margin = margin(1, 1,0,0, "cm"),
    legend.position = "none",
    legend.title = element_blank() ,
    axis.line = element_line() ) +
  labs(y = "Number of neophytes", 
       x = "Number of hybrids") -> pb


# examples of genera with high proportions --------------------------------
d_tax %>% arrange(desc(ratio_hyb), ratio_neo) %>% View

View(d_tax[with(d_tax , order(-ratio_hyb, ratio_neo)), ])


# add phylogenetic tree ---------------------------------------------------

# load Molinia Venegas tree and have a look
load("Data/Trees/Tree.RData")

tip_states <- read.csv("Data/hybrid_neophytes_analysis2.csv", header = T)
tip_states$hyb_frac <- tip_states$n_hyb/tip_states$nspp_in_genera
tip_states$neo_frac <- tip_states$n_neo/tip_states$nspp_in_genera

tip_states_all <- read.csv("Data/hybrid_neophyte_analysis4.csv", header=T)
# add tip labels
tips <- as.data.frame(Tree$tip.label)
colnames(tips) <- "genus"
tip_states_all_Tree <- left_join(tips, tip_states_all, by="genus") %>% 
  left_join(tip_states[,-4], by = "genus")


#
tip_states_all_Tree <- tip_states_all_Tree %>% filter(!is.na(family))



# remove NAs and genera with no neophytes or hybrids
length(Tree$tip.label)
Tree <- drop.tip(Tree, 
        c(tip_states_all_Tree[tip_states_all_Tree$type == "no hybrid & no neophyte", 1]))
Tree <- drop.tip(Tree, 
        c(tip_states_all_Tree[which(is.na(tip_states_all_Tree$type) == T), 1]))
# numbers of hybrids and neophytes
nrow(tip_states_all_Tree)
tip_states_all_Tree <- tip_states_all_Tree %>% 
  filter(!is.na(type)) %>% 
  filter(type != "no hybrid & no neophyte")



tip_states_all_Tree_long_x <- reshape2::melt(tip_states_all_Tree[,c(1,13,14)], id.vars = "genus")  %>% 
  mutate(value = log10(value + 1)) %>% 
  mutate(value = ifelse(variable == "n_neo", - value, value))


#
tip_states_all_Tree_long_x[is.na(tip_states_all_Tree_long_x)] <- 0

# fractions
#tip_states_all_Tree_long_x <- reshape2::melt(tip_states_all_Tree[,c(1,15,16)], id.vars = "genus") %>% 
#  mutate(value = log10(value)) %>% 
#  mutate(value = ifelse(variable == "neo_frac", - value, value)) 

# visualize
(ggtree(Tree, layout="circular", size = 0.1) +
  geom_fruit(data=tip_states_all_Tree_long_x, 
             geom=geom_col, 
             mapping = aes(y=genus, x=value, fill = variable), 
             width= 10, orientation = "y", offset = 0.2) +
  scale_fill_manual(values = c("#33E4FFFF", "grey70"), 
                    labels = c("Number of hybrids", "Number of neophytes"),
                    name = "") +
  # geom_tiplab(size = .5) + 
  theme(legend.position = c(.0, .85),
        plot.margin = margin(0, 0,0,0, "cm"),
        text=element_text(family="Helvetica")) -> pc)

ggarrange(
  ggarrange(pb, pa, ncol = 2, labels = c("a", "b"), align = "v"),
  pc, nrow = 2, labels = c("","c"), heights = c(1.4,2))

ggsave(dpi = 600, height = 6.05, width = 6.02, bg = "white",
       "Figures/hybrid_neo_ratio_taxonomoy.png")



# phylo signal ------------------------------------------------------------

# data preparation
tips <- as.data.frame(Tree$tip.label)
colnames(tips) <- "genus"
tip_states_all_Tree <- dplyr::left_join(dplyr::left_join(tips, tip_states_all, by="genus"), tip_states[,-4], by = "genus")
str(tip_states_all_Tree)

# Phylogenetic signal for number of hybrids
traits <- tip_states_all_Tree$n_hyb
names(traits) <- tip_states_all_Tree$genus

phySig <- phylosig(Tree, tip_states_all_Tree$n_hyb, method = "K", test= T)
lambdaModel <- fitContinuous(Tree, traits, model = "lambda")
lambdaModel$opt$aicc

#compare to null model
nosigModel <- fitContinuous(Tree, traits, model="lambda", bounds=list(lambda=c(0, 0.0001)))
nosigModel$opt$aicc

p.lam <- pchisq( 2*(lambdaModel$opt$lnL - nosigModel$opt$lnL), 1, lower.tail=FALSE)
p.lam

# Phylogenetic signal for number of neophytes
traits <- tip_states_all_Tree$n_neo
names(traits) <- tip_states_all_Tree$genus

phySig <- phylosig(Tree, tip_states_all_Tree$n_neo, method = "K", test= T)
lambdaModel <- fitContinuous(Tree, traits, model = "lambda")
lambdaModel$opt$aicc

# compare to null model
nosigModel <- fitContinuous(Tree, traits, model="lambda", bounds=list(lambda=c(0, 0.0001)))
nosigModel$opt$aicc

p.lam <- pchisq( 2*(lambdaModel$opt$lnL - nosigModel$opt$lnL), 1, lower.tail=FALSE)
p.lam

# Phylogenetic signal for presence of hybrids or neophytes (binary character)
Tree$node.label <- c(1:length(Tree$node.label))
Tree3 <- di2multi(Tree)

tip_states_all$hyb_present <- gsub("no hybrid", 0, tip_states_all$hyb_present)
tip_states_all$hyb_present <- gsub("hybrid", 1, tip_states_all$hyb_present)

tip_states_all$neo_present <- gsub("no neophyte", 0 , tip_states_all$neo_present)
tip_states_all$neo_present <- gsub("neophyte", 1, tip_states_all$neo_present)

tips <- as.data.frame(Tree$tip.label)
colnames(tips) <- "genus"
tip_states_all_Tree <- dplyr::left_join(dplyr::left_join(tips, tip_states_all, by="genus"), tip_states[,-4], by = "genus")
str(tip_states_all_Tree)

#neo present
all_dat <- comparative.data(Tree3, tip_states_all_Tree[which(is.na(tip_states_all_Tree$neo_present) == F), ], genus, na.omit = F)
res_neo <- phylo.d(all_dat, binvar=neo_present, permut = 1000, rnd.bias=NULL)
plot(res_neo)

#hybrid present
all_dat_hyb <- comparative.data(Tree3, tip_states_all_Tree[which(is.na(tip_states_all_Tree$hyb_present) == F), ], genus, na.omit = F)
res_hyb <- phylo.d(all_dat_hyb, binvar=neo_present, permut = 1000, rnd.bias=NULL)
plot(res_hyb)

