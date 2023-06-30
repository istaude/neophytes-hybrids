source("R-code/00_preamble.R")

# aim ---------------------------------------------------------------------


# quantify how often a neophyte meets the same genus.
# a species being a neophyte in a region, may hybridize preferentially if there
# is also the same genus present in the region newly colonized. does that concur
# often or not so often?


# read and harmonize glonaf with powo -------------------------------------

# read glonaf data
glonaf_sp <- read_delim("Data/GLONAF/Taxon_x_List_GloNAF_vanKleunenetal2018Ecology.csv", 
                        delim = "\t", escape_double = FALSE, 
                        locale = locale(encoding = "utf-16"), 
                        trim_ws = TRUE)

# how many species?
glonaf_sp %>% select(standardized_name) %>% distinct %>% nrow


# join correct genus name, because otherwise this will be difficult to
# match with powo data
# extract genus
glonaf_higher_taxa <- glonaf_sp %>% 
  filter(hybrid == 0) %>% 
  filter(status == "naturalized") %>% 
  filter(name_status == "accepted") %>% 
  select(standardized_name, family_tpl) %>% 
  distinct %>% 
  mutate(genus = sapply(strsplit(standardized_name, " "), "[", 1))


# genus names of glonaf may not match the accepted taxonomy of WCVP
# for example they list Spartina species but in fact they belong to
# the genus Sporobolus
glonaf_higher_taxa_wcvp <- left_join(
  glonaf_higher_taxa,
  kew_sp %>% 
    filter(taxon_rank == "Genus") %>% 
    select(genus, taxon_status, accepted_plant_name_id) ) 

# now if a genus has more than 1 match...
glonaf_higher_taxa_wcvp %>% select(genus, accepted_plant_name_id) %>% distinct %>% 
  count(genus) %>% arrange(desc(n))

# pick the accepted genus
glonaf_higher_taxa_wcvp <- bind_rows(
  
  left_join(glonaf_higher_taxa_wcvp,
            glonaf_higher_taxa_wcvp %>% select(genus, accepted_plant_name_id) %>% distinct %>% 
              count(genus)) %>% filter(n == 1), 
  
  left_join(glonaf_higher_taxa_wcvp,
            glonaf_higher_taxa_wcvp %>% select(genus, accepted_plant_name_id) %>% distinct %>% 
              count(genus)) %>% filter(n > 1) %>% filter(taxon_status == "Accepted") 
) %>% select(-n)


# now retrieve the correct genus name
glonaf_higher_taxa_wcvp <- left_join(
  glonaf_higher_taxa_wcvp,
  kew_sp %>% 
    filter(taxon_rank == "Genus") %>%
    filter(taxon_status == "Accepted") %>% 
    select(plant_name_id, genus) %>% 
    distinct %>% 
    rename(genus_acc = genus), by = c("accepted_plant_name_id" = "plant_name_id"))


# checks
glonaf_higher_taxa_wcvp %>% select(accepted_plant_name_id, genus_acc) %>% distinct %>% 
  count(genus_acc) %>% arrange(desc(n))

glonaf_higher_taxa_wcvp <- glonaf_higher_taxa_wcvp %>% filter(!is.na(genus_acc))

# join with glonaf_sp
glonaf_sp <- left_join(glonaf_higher_taxa_wcvp %>% select(standardized_name, genus_acc) , glonaf_sp)


# how many species now?
glonaf_sp %>% select(standardized_name) %>% distinct %>% nrow



# join with distribution data ---------------------------------------------

# read region ids for l3 areas
glonaf_reg <- read_csv("Data/GLONAF/Region_GloNAF_vanKleunenetal2018Ecology.csv")

glonaf_full <- full_join(glonaf_sp, glonaf_reg)

glonaf_tdwg3 <- glonaf_full %>% select(standardized_name, genus_acc, tdwg3) %>% distinct

# there are 127,000 colonizations events that led to a neophyte tdwg3 combination

# read powo data ----------------------------------------------------------

# read species list
kew_sp <- fread("Data/kew/checklist_names.txt", encoding = "UTF-8")

# kew species data
kew_gen <- kew_sp %>% 
  filter(taxon_rank == "Genus") %>% 
  filter(taxon_status == "Accepted") %>% 
  select(genus, plant_name_id) %>% 
  distinct

# kew data on distribution
kew_dis <- fread("Data/kew/checklist_distribution.txt")
head(kew_dis)
nrow(kew_dis)

kew_gen_dis <- left_join(kew_gen, kew_dis)

kew_gen_dis <- kew_gen_dis %>% 
  filter(introduced == 0) %>% 
  select(genus, plant_name_id, area_code_l3)


# join powo and glonaf data -----------------------------------------------
glonaf_kew_dis_comp <- left_join(glonaf_tdwg3, kew_gen_dis, 
                                 by = c("genus_acc" = "genus", "tdwg3" = "area_code_l3"))

glonaf_kew_dis_comp %>% nrow

# get rid of the numeric tdwg3 regions, they are major ones
glonaf_kew_dis_comp <- glonaf_kew_dis_comp %>% 
  filter(!str_detect(tdwg3,"^\\s*[0-9]*\\s*$")) %>% 
  rename(tdwg3_neopyhte_occ = tdwg3, tdwg3_native_occ = plant_name_id)

# save this data frame, is useful for the analysis on how many of the genera
# with both neophytes and hybrids contribute to naturalizatione events.
write_csv(glonaf_kew_dis_comp, "Data/glonaf_kew_dis_comp.csv")

# how many species in this rather final df
glonaf_kew_dis_comp %>% select(standardized_name) %>% distinct %>% nrow

# what's the overall rate
# full number of neophyte events
glonaf_kew_dis_comp %>% nrow

# full number of neophytes event with potential to meet species of same genus
glonaf_kew_dis_comp %>% filter(!is.na(tdwg3_native_occ)) %>% nrow

# about 51% percent of the neophyte events meet same genus
glonaf_kew_dis_comp %>% filter(!is.na(tdwg3_native_occ)) %>% nrow/
glonaf_kew_dis_comp %>% nrow


# show per species info, how many species never met the same genus?
no_match <- glonaf_kew_dis_comp %>% 
  group_by(standardized_name) %>% 
  summarize(n_distinct(no_nats_events = tdwg3_neopyhte_occ), 
            no_nat_events_match = sum(!is.na(tdwg3_native_occ))) %>% 
  filter(no_nat_events_match == 0)

yes_match <- glonaf_kew_dis_comp %>% 
  group_by(standardized_name) %>% 
  summarize(n_distinct(no_nats_events = tdwg3_neopyhte_occ), 
            no_nat_events_match = sum(!is.na(tdwg3_native_occ))) %>% 
  filter(no_nat_events_match >= 1)

# 74% of the species found at least one time a match, 26% of neophytes never found
# the same genus.
