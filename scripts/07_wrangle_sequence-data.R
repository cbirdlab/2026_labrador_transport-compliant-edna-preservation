#### Wrangle Files for Analysis ####
## Kevin Labrador
## 2025-04-05

#### 01 INTRODUCTION ####
# Collate all input files for analysis in a single .rds object. 
# Create several iteration of the data files depending on features that will be retained for downstream analyses.

#### 02 INITIALIZE ####

#### 02a Housekeeping ####
# Clear global environment.
rm(list = ls())

# Set-up the working directory in the source file location: 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### 02b Load Libraries ####
pacman::p_load(
  janitor,
  ape,
  Biostrings,
  phyloseq,
  microViz,
  readxl,
  ggVennDiagram,
  ggnested,
  ggpubr,
  phyloseq,
  tidytree,
  tidyverse
)

#### 03 USER DEFINED VARIABLES ####

#### 03a Assign INPUT File Paths ####
path_field_sample_info <- 
  "../data/rme_sample-info.xlsx"

path_water_params <- 
  "../data/rme_water-params.xlsx"

path_site_info <- 
  "../data/rme_site-info.xlsx"

path_maynard_metadata <- 
  "../data/rme_maynard-metadata.xlsx"

path_seq_sample_info <- 
  "../data/seq_sample_info.csv"

path_zotu_table <- 
  "../data/rainbow-bridge/zotu_table.csv"

path_lca_table <- 
  "../data/rainbow-bridge/lowest_taxonomy_table.csv"

path_zotu_seqs <- 
  "../data/rainbow-bridge/zotu_sequences.fasta"


#### 03b Assign OUTPUT File Paths ####
path_compiled_files <- 
  "../data/rme_phyloseq.rds"

path_plot_16su_empirical <- 
  "../results/plot_16su-empricial.jpg"



#### 04 IMPORT FILES ####
site_info <- 
  path_site_info %>% 
  read_xlsx() %>% 
  clean_names() %>% 
  
  # Rename OKG to TLK
  mutate (
    site = 
      case_when (
        site == "Okgok" ~ "Talakhaya",
        T ~ site
      ),
    site_code =
      case_when (
        site_code == "OKG" ~ "TLK",
        T ~ site_code
      ),
    resiliency_category =
      case_when (
        site_code == "TLK" ~ "medium-low",
        T ~ resiliency_category
      )
  )


water_params <- 
  path_water_params %>% 
  read_xlsx() %>% 
  clean_names() %>% 
  
  # Rename OKG to TLK
  mutate (
    site = 
      case_when (
        site == "Okgok" ~ "Talakhaya",
        T ~ site
      ),
    site_code =
      case_when (
        site_code == "OKG" ~ "TLK",
        T ~ site_code
      )
  )
    

maynard_metadata <- 
  path_maynard_metadata %>% 
  read_xlsx() %>% 
  clean_names() %>% 
  # Remove OKG
  filter (site != "Okgok")

seq_info <- 
  path_seq_sample_info %>% 
  read_csv() %>% 
  clean_names() %>% 
  
  # Rneame okg to tlk
  mutate(across(where(is.character), ~ str_replace(.x, "okg", "tlk")))

zotu_table <- 
  path_zotu_table %>% 
  read_csv %>% 
  clean_names()

# Calculate number of reads
zotu_table %>% 
  select (
    contains ("nec"),
    contains ("ntc"),
    contains ("rme"),
    contains ("rmep")
  ) %>% 
  sum()


lca_table <- 
  path_lca_table %>% 
  read_csv %>% 
  clean_names()


seqs <- readDNAStringSet(path_zotu_seqs)

#### 05 PREPARE SAMPLE METADATA ####
water_params_mean <- 
  water_params %>% 
  select (-site_code) %>% 
  group_by (
    site
  ) %>% 
  summarize (
    across(
      .cols = where(is.numeric) & -replicate,
      .fns = list (
        mean = ~mean (.x, na.rm = T),
        sd = ~sd(.x, na.rm = T)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) 

sample_metadata <- 
  seq_info %>% 
  select (
    dna_extract_tube_id,
    sample_type,
    preservative,
    sample_id,
    pcr_quant_ng
  ) %>% 
  
  # Remove positive control
  filter (
    dna_extract_tube_id != "PC"
  ) %>% 
  
  # Assign site codes
  mutate (
    site_code = 
      case_when (
        sample_type == "sample" ~ substr(sample_id, 1, 3),
        grepl("nfc", sample_id) ~ substr(sample_id, 1, 3),
        grepl("nec", sample_id) ~ substr(sample_id, 1, 3),
        grepl("nec", dna_extract_tube_id) ~ "nec",
        grepl("PCR", dna_extract_tube_id) ~ "pcr",
        T ~ NA
      )
  ) %>%

  # Join sample info with site info
  left_join (
    site_info %>%
      select (
        site, 
        site_code, 
        resiliency_category,
        longitude,
        latitude,
        depth_ft
      ) %>% 
      mutate (
        site_code = tolower(site_code)
      )
  ) %>% 
  
  # Join with water params
  left_join(
    water_params_mean
  ) %>% 
  left_join(
    maynard_metadata %>% 
      mutate (
        site_code = tolower(site_code)
      )
  )


sample_id_decode <- 
  sample_metadata %>% 
  select (
    dna_extract_tube_id,
    sample_id
  ) %>% 
  mutate (
    sample_id =
      case_when (
        dna_extract_tube_id == "rme_15" ~ "nec_01",
        dna_extract_tube_id == "rme_24" ~ "nec_02",
        dna_extract_tube_id == "rme_35" ~ "nec_03",
        dna_extract_tube_id == "PCR-NTC1" ~ "ntc_01",
        dna_extract_tube_id == "PCR-NTC2" ~ "ntc_02",
        dna_extract_tube_id == "PCR-NTC3" ~ "ntc_03",
        dna_extract_tube_id == "PCR-NTC4" ~ "ntc_04",
        !is.na(sample_id) ~ str_remove(sample_id, "-\\d+"),
        T ~ tolower(dna_extract_tube_id)
      ) %>% 
      str_replace_all ("-", "_"),
    
    fraction = 
      case_when (
        grepl("rme_", dna_extract_tube_id) ~ "filter",
        grepl("rmep_|nec", dna_extract_tube_id) ~ "filtrate",
        grepl("ntc", sample_id) ~ "pcr",
        T ~ NA
      ),
    
    dna_extract_tube_id = 
      dna_extract_tube_id %>%
      tolower() %>% 
      str_replace_all ("-", "_"),
    
    new_id = 
      paste(sample_id, fraction, sep = "_")
    
  ) 

name_key <- 
  setNames(
    sample_id_decode$dna_extract_tube_id,
    sample_id_decode$new_id
  )

sample_metadata_for_phyloseq <-
  sample_metadata %>% 
  
  # Change sample id in metadata sheet
  mutate (sample_id = sample_id_decode$new_id) %>% 
  
  # Add fraction information
  mutate (fraction = sub(".*_", "", sample_id)) %>% 
  
  # # Create a new grouping variable (pooled sites, separate controls)
  # mutate (
  #   grp =
  #     factor (
  #       case_when (
  #         sample_type == "sample" ~ site_code,
  #         T ~ sample_id
  #       )
  #     )
  # ) %>%
  
  # Create a new grouping variable (pooled fractions, separate controls)
  mutate (
    grp =
      factor (
        case_when (
          sample_type == "sample" ~ str_remove(sample_id, "_(filter|filtrate)$"),
          T ~ sample_id
        )
      )
  ) %>%
  
  # Assign sample id as rownames
  column_to_rownames("sample_id") %>% 
  
  sample_data()

# Rename the zotu table
zotu_table_renamed <- 
  zotu_table %>% 
  rename(any_of(name_key))


# Quick stats
(data_stats <- 
    zotu_table_renamed %>% 
    select (
      zotu,
      taxid,
      taxid_rank,
      species,
      contains ("filter"),
      contains ("filtrate"),
      contains ("pcr")
    ) %>% 
    summarize (
      n_reads   = sum(across(matches("filter|filtrate|pcr")), na.rm = TRUE),
      n_zotus = n_distinct (zotu),
      n_taxid = n_distinct (taxid),
      n_species = n_distinct(species[taxid_rank == "species"], na.rm = TRUE)
    ) %>% 
    t() 
)
colnames(data_stats) <- "pre_classification"

#### 07 CLASIFY HABITATS (MARINE VS NON-MARINE) ####
# genus_filter
filter_marine_genera <- 
  c(
    "Abralia",
    "Abudefduf",
    "Acanthastrea",
    "Acanthurus",
    "Actinopyga",
    "Afrocucumis",
    "Amphipholis",
    "Anampses",
    "Anguilla",
    "Apagesoma",
    "Aphareus",
    "Atergatopsis",
    "Axiopsis",
    "Balistapus",
    "Beania",
    "Blenniella",
    "Bohadschia",
    "Calcinus",
    "Callogobius",
    "Calotomus",
    "Canarium",
    "Centropyge",
    "Cephalothrix",
    "Cerebratulus",
    "Cerithium",
    "Cetoscarus",
    "Chaetodon",
    "Cheilinus",
    "Cheilopogon",
    "Cherusius",
    "Chlorurus",
    "Chrysiptera",
    "Chthamalus",
    "Ciliopagurus",
    "Cirrhilabrus",
    "Cirripectes",
    "Comatella",
    "Conus",
    "Coralliophila",
    "Coris",
    "Ctenochaetus",
    "Cyphastrea",
    "Dendropoma",
    "Diadema",
    "Dichrometra",
    "Diodora",
    "Diploastrea",
    "Diplomma",
    "Dipolydora",
    "Echinothrix",
    "Elysia",
    "Encrasicholina",
    "Entomacrodus",
    "Epinephelus",
    "Euapta",
    "Euprymna",
    "Eviota",
    "Fragum",
    "Galathea",
    "Gymnothorax",
    "Halichoeres",
    "Holothuria",
    "Homophyllia",
    "Ircinia",
    "Istiblennius",
    "Kuhlia",
    "Kyphosus",
    "Labidodemas",
    "Labroides",
    "Lambis",
    "Liocarpilodes",
    "Liopropoma",
    "Lutjanus",
    "Macropharyngodon",
    "Macrophiothrix",
    "Madracis",
    "Melichthys",
    "Metalpheus",
    "Metapenaeopsis",
    "Microtis",
    "Monotaxis",
    "Moolgarda",
    "Mulloidichthys",
    "Naso",
    "Neoliomera",
    "Novaculichthys",
    "Octopus",
    "Ophiocoma",
    "Ophiocomella",
    "Ophiomastix",
    "Ophionereis",
    "Orbicella",
    "Oscarella",
    "Oxycheilinus",
    "Panulirus",
    "Parabetaeus",
    "Paracirrhites",
    "Paragoniastrea",
    "Paraxanthias",
    "Parupeneus",
    "Pempheris",
    "Percnon",
    "Periglypta",
    "Pomachromis",
    "Phasianella",
    "Phyrella",
    "Plectroglyphidodon",
    "Plesiastrea",
    "Pocillopora",
    "Pomacanthus",
    "Praealticus",
    "Pseudocheilinus",
    "Pseudochromis",
    "Ptereleotris",
    "Pterocaesio",
    "Rhabdoblennius",
    "Rhinecanthus",
    "Sabia",
    "Sargocentron",
    "Scarus",
    "Scomberoides",
    "Seriatopora",
    "Sistrum",
    "Stenella",
    "Stichopus",
    "Stylocoeniella",
    "Stylophora",
    "Sufflamen",
    "Thalamita",
    "Thalassoma",
    "Thysanopoda",
    "Timarete",
    "Trevathana",
    "Trimma",
    "Turbo",
    "Tylocarcinus",
    "Valenciennea",
    "Vasum",
    "Zanclus",
    "Zebrasoma"
  )


(families_with_dropped_genus <- 
    zotu_table_renamed %>%
    filter (genus == "LCA_dropped") %>% 
    select (family) %>% 
    distinct() %>% 
    arrange (family)
)

filter_marine_family <- 
  c(
    "Acanthuridae",
    "Astrocoeniidae",
    "Balistidae",
    "Cardiidae",
    "Delphinidae",
    "Diploastraeidae",
    "Exocoetidae",
    "Faviidae",
    "Lobophylliidae",
    "Merulinidae",
    "Ophidiidae",
    "Plesiastreidae",
    "Pocilloporidae"
  )


(order_with_dropped_families <- 
    zotu_table_renamed %>%
    filter (family == "LCA_dropped") %>% 
    select (order) %>% 
    distinct() %>% 
    arrange (order)
)

filter_marine_order <- "Anthozoa"


zotu_table_marine <- 
  zotu_table_renamed %>% 
  filter (
    (genus %in% filter_marine_genera) |
      (genus == "LCA_dropped" & family %in% filter_marine_family) |
      (genus == "LCA_dropped" & family == "LCA_dropped" & order %in% filter_marine_order) 
  ) %>% 
  mutate (
    habitat = "marine"
  ) 

marine_zotus <- zotu_table_marine$zotu

zotu_table_non_marine <- 
  zotu_table_renamed %>% 
  filter (
    !zotu %in% zotu_table_marine$zotu
  ) %>% 
  mutate (
    habitat = "non-marine"
  ) 

non_marine_zotus <- zotu_table_non_marine$zotu

zotu_table_with_habitat <- 
  rbind (
    zotu_table_marine,
    zotu_table_non_marine
  )

(data_stats <-
    cbind (
      data_stats,
      zotu_table_with_habitat %>% 
        select (
          zotu,
          taxid,
          taxid_rank,
          species,
          habitat,
          contains ("filter"),
          contains ("filtrate"),
          contains ("pcr")
        ) %>% 
        group_by(habitat) %>% 
        summarize (
          n_reads   = sum(across(matches("filter|filtrate|pcr")), na.rm = TRUE),
          n_zotus = n_distinct (zotu),
          n_taxid = n_distinct (taxid),
          n_species = n_distinct(species[taxid_rank == "species"], na.rm = TRUE)
        ) %>% 
        column_to_rownames("habitat") %>% 
        t()
    )
)


#### 08 PREPARE ZOTU METADATA ####
zotu_metadata <-
  zotu_table_with_habitat 

zotu_taxonomy <- 
  zotu_metadata %>% 
  select (
    taxid,
    zotu,
    kingdom,
    phylum,
    class,
    order,
    family,
    genus,
    species,
    habitat
  ) %>% 
  unique() %>% 
  mutate (
    taxid = trimws(taxid)
  )

zotu_taxonomy_for_phyloseq <-
  zotu_taxonomy %>% 
  column_to_rownames("zotu") %>% 
  as.matrix() %>% 
  tax_table()


#### 09 VISUALIZE TAXONOMIC COMPOSITION ####
# Assess taxonomic composition
(plot_sp_richness <- 
   zotu_metadata %>% 
   select (
     taxid,
     habitat,
     phylum,
     class,
     species
   ) %>% 
   filter (
     species != "LCA_dropped"
   ) %>% 
   distinct() %>% 
   ggnested(
     aes (
       x = fct_infreq(phylum),
       main_group = phylum,
       sub_group = class
     ),
     legend_title = "Phylum;\nClass"
   )+
   geom_bar() +
   theme_bw() +
   labs (
     x = "Phylum",
     y = "Species Richness",
     subtitle = "No. of unique species"
   ) + 
   facet_wrap(~habitat, scales = "free_x") + 
   theme (
     legend.position = "right",
     axis.text.x = element_text (angle = 45, hjust = 1, vjust = 1),
     legend.text = ggtext::element_markdown()
   )
)


(plot_rel_comp <- 
    zotu_table_with_habitat %>%
    select (
      taxid,
      habitat,
      phylum,
      class,
      species
    ) %>% 
    filter(
      species != "LCA_dropped"
    ) %>% 
    distinct() %>% 
    ggnested(
      aes (
        x = fct_infreq(phylum),
        main_group = phylum,
        sub_group = class
      ),
      legend_title = "Phylum;\nClass"
    )+
    geom_bar(position = "fill") +
    theme_bw() +
    labs (
      x = "Phylum",
      y = "Relative Species Richness (%)",
      subtitle = "Relative composition"
    ) + 
    facet_wrap(~habitat, scales = "free_x") + 
    theme (
      legend.position = "right",
      axis.text.x = element_text (angle = 45, hjust = 1, vjust = 1),
      legend.text = ggtext::element_markdown()
    ) 
)

(plot_16su_empirical <- 
    ggarrange(
      plot_sp_richness +
        theme (axis.text.x = element_blank(),
               axis.title.x = element_blank()),
      plot_rel_comp,
      labels = "AUTO",
      nrow = 2,
      common.legend = T,
      legend = "right",
      align = "v"
    )
)

#### 10 PREPARE ZOTU ABUNDANCE TABLE ####
community_matrix <- 
  zotu_table_with_habitat %>% 
  select (
    zotu,
    contains ("filter"),
    contains ("filtrate"),
    contains ("pcr")
  ) 


zotu_read_abundance <- 
  community_matrix %>% 
  mutate (
    n_reads = 
      rowSums (
        across (-c(zotu))
      )
  ) %>% 
  select (
    zotu,
    n_reads
  ) %>% 
  left_join(
    zotu_taxonomy
  ) %>% 
  arrange (n_reads)

(total_reads <- 
    zotu_read_abundance$n_reads %>% 
    sum()
)

zotu_matrix_for_phyloseq <- 
  community_matrix %>% 
  column_to_rownames("zotu") %>% 
  otu_table(
    taxa_are_rows = T
  )

#### 11 PREPARE FASTA FILE ####
old_names <- names(seqs)
new_names <- 
  str_extract(old_names, "^Zotu\\d+")


seqs_for_phyloseq <- seqs
names(seqs_for_phyloseq) <- new_names


#### 12 CREATE PHYLOSEQ FILE ####
dim(sample_metadata_for_phyloseq)
dim(zotu_taxonomy_for_phyloseq)
dim(zotu_matrix_for_phyloseq)
length(seqs_for_phyloseq)


# Create the phyloseq object
rme_phyloseq <- 
  phyloseq(
    zotu_matrix_for_phyloseq,
    zotu_taxonomy_for_phyloseq,
    sample_metadata_for_phyloseq,
    seqs_for_phyloseq
  ) %>% 
  # Remove zotus with NA in kingdom
  subset_taxa(., !is.na (kingdom)) %>% 
  tax_fix(
    min_length = 2,
    unknowns = c("LCA_dropped"),
    sep = "_", 
    anon_unique = TRUE,
    suffix_rank = "classified"
  )


# Validate phyloseq
phyloseq_validate (rme_phyloseq)


#### 13 EXPORT FILES ####
saveRDS(
  rme_phyloseq,
  path_compiled_files
)

ggsave (
  plot_16su_empirical,
  file = path_plot_16su_empirical,
  width = 12,
  height = 7,
  units = "in",
  dpi = 330
)


