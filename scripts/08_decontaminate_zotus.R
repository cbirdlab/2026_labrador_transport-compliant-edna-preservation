#### decontaminate zotus ####
## Kevin Labrador

#### 01 INTRODUCTION ####
# Assess negative controls before doing any further zotu filtering/clustering
# Decontaminate phyloseq object using the *decontam* package

#### 02 INITIALIZE ####

#### 02a Housekeeping ####
# Clear global environment.
rm(list = ls())

# Set-up the working directory in the source file location: 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### 02b Load Libraries ####
pacman::p_load(
  ggside,
  ggpubr,
  decontam,
  phyloseq,
  microViz,
  MiscMetabar,
  vegan,
  readxl,
  janitor,
  tidyverse
)


#### 03 USER DEFINED VARIABLES ####

#### 03a Assign INPUT File Paths ####
path_infiles <- 
  "../data/rme_phyloseq.rds"


#### 03b Assign OUTPUT File Paths ####
path_plot_ordination_preDecon <- 
  "../results/plot_ordination_pre-decon.png"

path_plot_ordination_postDecon <- 
  "../results/plot_ordination_post-decon.png"

path_plot_ordination_contam <- 
  "../results/plot_ordination_contam.png"

path_plot_ordination_without_marine <- 
  "../results/plot_ordination_post-decon-no-marine.png"

path_plot_ordination_samples_only <- 
  "../results/plot_neg_ctrl_ordination_post-decon-no-marine_samples-only.png"

path_permanova_samples <- 
  "../results/permanova_sample_ctrl.txt"

path_phyloseq_decontaminated <- 
  "../data/rme_phyloseq_decontaminated.rds"




#### 04 IMPORT FILES ####
phyloseq_raw <- 
  path_infiles %>% 
  readRDS() 

# Wrangle input phyloseq file
phyloseq <- 
  phyloseq_raw %>% 
  # Aggregate samples based on grouping
  merge_samples("grp") %>% 
  # Remove taxid on tax table to allow aggregation
  tax_mutate (taxid = NULL) %>% 
  # Add new info to sample data
  ps_mutate(
    sample_type = 
      case_when (
        grepl ("nec", rownames (sample_data(.))) ~ "extraction n.ctrl",
        grepl ("nfc", rownames (sample_data(.))) ~ "field n.ctrl",
        grepl ("ntc", rownames (sample_data(.))) ~ "pcr n.ctrl",
        T ~ "sample"
      ),
    sample_id = 
      rownames(sample_data(.)),
    fraction = 
      case_when (
        grepl ("filter", rownames(sample_data(.))) ~ "filter",
        grepl ("filtrate", rownames(sample_data(.))) ~ "filtrate",
        grepl ("pcr", rownames(sample_data(.))) ~ "pcr",
        T ~ "sample"
      ),
    preservative = 
      case_when (
        grepl ("_cs", rownames(sample_data(.))) ~ "No Pres",
        grepl ("_ds", rownames(sample_data(.))) ~ "DESS",
        grepl ("_tl", rownames(sample_data(.))) ~ "Buffer TL",
        T ~ "Lab Control"
      ) %>% 
      factor (
        levels = c(
          "No Pres",
          "DESS",
          "Buffer TL",
          "Lab Control"
        )
      ),
    
    plot_group = ifelse(
      sample_type == "sample",
      "sample",
      sample_type
    ) %>% 
      factor (
        levels = c(
          "sample",
          "field n.ctrl",
          "extraction n.ctrl",
          "pcr n.ctrl"
        )
      )
    
    # plot_group = ifelse(
    #   sample_type == "sample",
    #   as.character(preservative),
    #   sample_type
    # ) %>% 
    #   factor (
    #     levels = c(
    #       "No Pres",
    #       "DESS",
    #       "Buffer TL",
    #       "field n.ctrl",
    #       "extraction n.ctrl",
    #       "pcr n.ctrl"
    #     )
    #   )
  )


#### 05 VISUALIZE RAW DATASET ####

# Assign phyloseq object to analyze
ps <- phyloseq %>% 
  # Reorder factors
  ps_mutate (
    preservative = 
      factor (preservative, 
              levels = c(
                "No Pres",
                "DESS",
                "Buffer TL",
                "Lab Control"
              )
      ),
    sample_type = 
      factor (sample_type,
              levels = c(
                "sample",
                "extraction n.ctrl",
                "field n.ctrl",
                "pcr n.ctrl"
              )
      )
  )



# Plot library size vs n_zotus
df <- 
  ps %>% 
  sample_data %>% 
  as("matrix") %>% 
  as.data.frame() %>% 
  mutate (
    library_size = sample_sums (ps),
    n_zotu =  estimate_richness(ps, measures = "Observed")[[1]],
    preservative = 
      factor (preservative, 
              levels = c(
                "No Pres",
                "DESS",
                "Buffer TL",
                "Lab Control"
              )
      ),
    sample_type = 
      factor (sample_type,
              levels = c(
                "sample",
                "extraction n.ctrl",
                "field n.ctrl",
                "pcr n.ctrl"
              )
      )
  )

df_stats <- 
  df %>% 
  select (
    sample_type,
    library_size,
    n_zotu
  ) %>% 
  mutate (
    grp = 
      case_when (
        sample_type == "sample" ~ "sample",
        T ~ "ctrl"
      )
  )

wilcox.test (library_size ~ grp, data = df_stats, exact = F)  
wilcox.test (n_zotu ~ grp, data = df_stats)


(summary_stats <- 
    df_stats %>% 
    group_by (grp) %>% 
    get_summary_stats () 
)


(plot_reads_zotus <- 
    df %>% 
    ggplot(
      aes (
        x = library_size,
        y = n_zotu,
        fill = sample_type,
        color = sample_type,
        shape = preservative
      )
    ) +
    geom_point(
      size = 3, 
      alpha = 0.75
    )+
    # ggrepel::geom_text_repel(
    #   aes(label = sample_id,
    #       color = sample_type)
    # ) +
    scale_y_log10() + 
    scale_x_log10() +
    labs (
      x = expression("No. of reads ("~log[10]~"scale)"),
      y = expression("No. of zOTUs ("~log[10]~"scale)"),
      fill = "Sample Type",
      color = "Sample Type",
      shape = "Preservative"
    ) +
    theme_bw() + 
    scale_color_brewer (
      palette = "Dark2"
    ) + 
    scale_fill_brewer (
      palette = "Dark2"
    ) + 
    
    scale_shape_manual(
      values = c(22, 24, 25, 21)
    )
)


# Plot composition

## Assign Color palette
# taxa_df <- phyloseq@tax_table %>% as.data.frame()
# tax_levels <- taxa_df$class %>% unique()
# n_colors <- length(tax_levels)
# 
# set.seed (999)
# cols <- 
#   Polychrome::createPalette(n_colors, c("#FF0000", "#00FF00", "#0000FF")) %>% 
#   colorspace::darken(amount = 0.1)
# names(cols) <- tax_levels
# scales::show_col(cols)


(plot_composition <- 
    ps %>%
    comp_barplot (
      tax_level = "class",
      n_taxa = 4,
      other_name = "Other",
      palette = distinct_palette (n = 4, add = "grey50"),
      merge_other = T,
      facet_by = "plot_group")+ 
    coord_flip() + 
    labs (
      y = "Relative Composition",
      fill = "Class"
    ) 
)



# Plot PCA Ordination
(plot_pca <- 
    ps%>% 
    tax_transform(
      trans = "rclr", # robust centered log-ratio transformation
      #rank = "class"
    ) %>% 
    ord_calc(
      method = "PCA"
    ) %>% 
    ord_plot(
      fill = "sample_type",
      col = "sample_type",
      shape = "preservative",
      size = 2,
      #plot_taxa = 1:4,
      plot_samples = T,
    ) +
    #stat_ellipse(aes(colour = sample_type))+
    # ggrepel::geom_text_repel(
    #   aes(
    #     label = sample_id,
    #     col = sample_type
    #     ),
    #   size = 3
    # ) +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(
      values = c(22, 24, 25, 21)
    ) + 
    theme_bw() +
    labs (
      fill = "Sample Type",
      color = "Sample Type",
      shape = "Preservative"
    ) 
)


# Plot PCoA Ordination
(plot_pcoa <- 
    ps%>% 
    tax_transform(
      trans = "rclr", # robust centered log-ratio transformation
      #rank = "class"
    ) %>% 
    dist_calc("euclidean") %>% 
    ord_calc("PCoA") %>% 
    ord_plot(
      color = "sample_type",
      fill = "sample_type",
      shape = "preservative",
      size = 2
    ) +
    stat_ellipse(aes(colour = sample_type))+
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    scale_shape_manual(
      values = c(22, 24, 25, 21)
    ) + 
    theme_bw() +
    labs (
      fill = "Sample Type",
      color = "Sample Type",
      shape = "Preservative"
    ) 
)

# Merge plots
(plot_neg_ctrl_assessment_pre_decon <- 
    ggpubr::ggarrange (
      plot_composition + 
        theme (legend.position = "bottom") +
        guides(fill = guide_legend(nrow = 2)),
      ggarrange (
        plot_reads_zotus +
          theme (legend.position = "bottom") +
          guides(color = guide_legend(nrow = 2),
                 shape = guide_legend(nrow = 2)), 
        plot_pca + 
          theme (legend.position = "bottom") +
          guides(color = guide_legend(nrow = 2),
                 shape = guide_legend(nrow = 2)),
        nrow = 2,
        labels = c("B", "C"),
        common.legend = T,
        legend = "bottom"
      ),
      ncol = 2,
      labels = c("A", NA),
      widths = c(1, 0.75)
    )
)


#### 06 DECONTAMINATE ####
# Follow Callahan and Davis (2018) (https://benjjneb.github.io/decontam/vignettes/decontam_intro.html)

# Step 1. Assign phyloseq object to decontaminate
ps <- phyloseq_raw

# Step 2. Identify contaminants using both Frequency and Prevalence
sample_data(ps)$is.neg <- 
  sample_data(ps)$sample_type != "sample"

contam_df <- 
  isContaminant(
    ps,
    method = "either", 
    conc = "pcr_quant_ng",
    neg = "is.neg",
    threshold = c(0.1, 0.5) 
  )

table(contam_df$contaminant)
head(which(contam_df$contaminant))

# Visualize contamination based on frequency
plot_frequency(ps, taxa_names(ps)[c(2, 5, 1546, 1571)], conc="pcr_quant_ng") + 
  xlab("Total DNA (ng) in Pool") +
  theme_bw()


# Visualize contamination based on prevalence
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_type != "sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_type == "sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- 
  data.frame(
    pa.pos=taxa_sums(ps.pa.pos), 
    pa.neg=taxa_sums(ps.pa.neg),
    contaminant=contam_df$contaminant)

ggplot(
  data=df.pa, 
  aes(x=pa.neg, 
      y=pa.pos, 
      fill=contaminant)
) + 
  geom_jitter(
    pch = 21,
    col = "black",
    alpha = 0.25
  ) +
  labs (
    x = "Prevalence (Negative Controls)",
    y = "Prevalence (True Samples)",
    fill = "Contaminant") +
  theme_bw() + 
  scale_fill_manual (
    values = c("gray80", "red3")
  )


# Step 3 Prune contaminants out of the phyloseq object
ps_decontaminated <- 
  prune_taxa(!contam_df$contaminant, ps)

ps_contaminants <- 
  prune_taxa(contam_df$contaminant, ps)


#### 07 VISUALIZE DECONTAMINATED DATASET ####

# Assign phyloseq object to analyze
ps <- 
  ps_decontaminated %>% 
  # Aggregate samples based on grouping
  merge_samples("grp") %>% 
  # Remove taxid on tax table to allow aggregation
  tax_mutate (taxid = NULL) %>% 
  # Remove 0 observations
  prune_taxa(taxa_sums(.) > 0, .) %>% 
  prune_samples(sample_sums(.) > 0, .) %>% 
  # Add new info to sample data
  ps_mutate(
    sample_type = 
      case_when (
        grepl ("nec", rownames (sample_data(.))) ~ "extraction n.ctrl",
        grepl ("nfc", rownames (sample_data(.))) ~ "field n.ctrl",
        grepl ("ntc", rownames (sample_data(.))) ~ "pcr n.ctrl",
        T ~ "sample"
      ),
    sample_id = 
      rownames(sample_data(.)),
    fraction = 
      case_when (
        grepl ("filter", rownames(sample_data(.))) ~ "filter",
        grepl ("filtrate", rownames(sample_data(.))) ~ "filtrate",
        grepl ("pcr", rownames(sample_data(.))) ~ "pcr",
        T ~ "sample"
      ),
    preservative = 
      case_when (
        grepl ("_cs", rownames(sample_data(.))) ~ "No Pres",
        grepl ("_ds", rownames(sample_data(.))) ~ "DESS",
        grepl ("_tl", rownames(sample_data(.))) ~ "Buffer TL",
        T ~ "Lab Control"
      ) %>% 
      factor (
        levels = c(
          "No Pres",
          "DESS",
          "Buffer TL",
          "Lab Control"
        )
      ),
    plot_group = ifelse(
      sample_type == "sample",
      as.character(preservative),
      sample_type
    ) %>% 
      factor (
        levels = c(
          "No Pres",
          "DESS",
          "Buffer TL",
          "field n.ctrl",
          "extraction n.ctrl",
          "pcr n.ctrl"
        )
      )
  )


# Plot library size vs n_zotus
df <- 
  ps %>% 
  sample_data %>% 
  as("matrix") %>% 
  as.data.frame() %>% 
  mutate (
    library_size = sample_sums (ps),
    n_zotu =  estimate_richness(ps, measures = "Observed")[[1]],
    preservative = 
      factor (preservative, 
              levels = c(
                "No Pres",
                "DESS",
                "Buffer TL",
                "Lab Control"
              )
      ),
    sample_type = 
      factor (sample_type,
              levels = c(
                "sample",
                "extraction n.ctrl",
                "field n.ctrl",
                "pcr n.ctrl"
              )
      )
  )

df_stats <- 
  df %>% 
  select (
    sample_type,
    library_size,
    n_zotu
  ) %>% 
  mutate (
    grp = 
      case_when (
        sample_type == "sample" ~ "sample",
        T ~ "ctrl"
      )
  )

wilcox.test (library_size ~ grp, data = df_stats, exact = F)  
wilcox.test (n_zotu ~ grp, data = df_stats)


(summary_stats <- 
    df_stats %>% 
    group_by (grp) %>% 
    get_summary_stats () 
)


(plot_reads_zotus <- 
    df %>% 
    ggplot(
      aes (
        x = library_size,
        y = n_zotu,
        fill = sample_type,
        color = sample_type,
        shape = preservative
      )
    ) +
    geom_point(
      size = 3, 
      alpha = 0.75)
  +
    scale_y_log10() + 
    scale_x_log10() +
    labs (
      x = expression("No. of reads ("~log[10]~"scale)"),
      y = expression("No. of zOTUs ("~log[10]~"scale)"),
      fill = "Sample Type",
      color = "Sample Type",
      shape = "Preservative"
    ) +
    theme_bw() + 
    scale_color_brewer (
      palette = "Dark2"
    ) + 
    scale_fill_brewer (
      palette = "Dark2"
    ) + 
    
    scale_shape_manual(
      values = c(22, 24, 25, 21)
    )
)


# Plot composition

## Assign Color palette
# taxa_df <- phyloseq@tax_table %>% as.data.frame()
# tax_levels <- taxa_df$class %>% unique()
# n_colors <- length(tax_levels)
# 
# set.seed (999)
# cols <- 
#   Polychrome::createPalette(n_colors, c("#FF0000", "#00FF00", "#0000FF")) %>% 
#   colorspace::darken(amount = 0.1)
# names(cols) <- tax_levels
# scales::show_col(cols)

(plot_composition <- 
    ps %>%
    comp_barplot (
      tax_level = "class",
      n_taxa = 4,
      other_name = "Other",
      palette = distinct_palette (n = 4, add = "grey50"),
      merge_other = T,
      facet_by = "plot_group")+ 
    coord_flip() + 
    labs (
      y = "Relative Composition",
      fill = "Class"
    ) 
)



# Plot PCA Ordination
(plot_pca <- 
    ps%>% 
    tax_transform(
      trans = "rclr", # robust centered log-ratio transformation
      #rank = "class"
    ) %>% 
    ord_calc(
      method = "PCA"
    ) %>% 
    ord_plot(
      fill = "sample_type",
      col = "sample_type",
      shape = "preservative",
      size = 2,
      #plot_taxa = 1:4,
      plot_samples = T,
    ) +
    stat_ellipse(aes(colour = sample_type))+
    # ggrepel::geom_text_repel(
    #   aes(
    #     label = sample_id,
    #     col = sample_type
    #     ),
    #   size = 3
    # ) + 
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(
      values = c(22, 24, 25, 21)
    ) + 
    theme_bw() +
    labs (
      fill = "Sample Type",
      color = "Sample Type",
      shape = "Preservative"
    ) 
)


# Plot PCoA Ordination
(plot_pcoa <- 
    ps%>% 
    tax_transform(
      trans = "rclr", # robust centered log-ratio transformation
      #rank = "class"
    ) %>% 
    dist_calc("euclidean") %>% 
    ord_calc("PCoA") %>% 
    ord_plot(
      color = "sample_type",
      fill = "sample_type",
      shape = "preservative",
      size = 2
    ) +
    stat_ellipse(aes(colour = sample_type))+
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    scale_shape_manual(
      values = c(22, 24, 25, 21)
    ) + 
    theme_bw() +
    labs (
      fill = "Sample Type",
      color = "Sample Type",
      shape = "Preservative"
    ) 
)

# Merge plots
(plot_neg_ctrl_assessment_post_decon <- 
    ggpubr::ggarrange (
      plot_composition + 
        theme (legend.position = "bottom") +
        guides(fill = guide_legend(nrow = 2)),
      ggarrange (
        plot_reads_zotus +
          theme (legend.position = "bottom") +
          guides(color = guide_legend(nrow = 2),
                 shape = guide_legend(nrow = 2)), 
        plot_pca + 
          theme (legend.position = "bottom") +
          guides(color = guide_legend(nrow = 2),
                 shape = guide_legend(nrow = 2)),
        nrow = 2,
        labels = c("B", "C"),
        common.legend = T,
        legend = "bottom"
      ),
      ncol = 2,
      labels = c("A", NA),
      widths = c(1, 0.75)
    )
)




#### 08 VISUALIZE CONTAMINATION ####
# Assign phyloseq object to analyze
ps <- 
  ps_contaminants %>% 
  # Aggregate samples based on grouping
  merge_samples("grp") %>% 
  # Remove taxid on tax table to allow aggregation
  tax_mutate (taxid = NULL) %>% 
  # Remove 0 observations
  prune_taxa(taxa_sums(.) > 0, .) %>% 
  prune_samples(sample_sums(.) > 0, .) %>% 
  # Add new info to sample data
  ps_mutate(
    sample_type = 
      case_when (
        grepl ("nec", rownames (sample_data(.))) ~ "extraction n.ctrl",
        grepl ("nfc", rownames (sample_data(.))) ~ "field n.ctrl",
        grepl ("ntc", rownames (sample_data(.))) ~ "pcr n.ctrl",
        T ~ "sample"
      ),
    sample_id = 
      rownames(sample_data(.)),
    fraction = 
      case_when (
        grepl ("filter", rownames(sample_data(.))) ~ "filter",
        grepl ("filtrate", rownames(sample_data(.))) ~ "filtrate",
        grepl ("pcr", rownames(sample_data(.))) ~ "pcr",
        T ~ "sample"
      ),
    preservative = 
      case_when (
        grepl ("_cs", rownames(sample_data(.))) ~ "No Pres",
        grepl ("_ds", rownames(sample_data(.))) ~ "DESS",
        grepl ("_tl", rownames(sample_data(.))) ~ "Buffer TL",
        T ~ "Lab Control"
      ) %>% 
      factor (
        levels = c(
          "No Pres",
          "DESS",
          "Buffer TL",
          "Lab Control"
        )
      ),
    plot_group = ifelse(
      sample_type == "sample",
      as.character(preservative),
      sample_type
    ) %>% 
      factor (
        levels = c(
          "No Pres",
          "DESS",
          "Buffer TL",
          "field n.ctrl",
          "extraction n.ctrl",
          "pcr n.ctrl"
        )
      )
  )



# Plot library size vs n_zotus
df <- 
  ps %>% 
  sample_data %>% 
  as("matrix") %>% 
  as.data.frame() %>% 
  mutate (
    library_size = sample_sums (ps),
    n_zotu =  estimate_richness(ps, measures = "Observed")[[1]],
    preservative = 
      factor (preservative, 
              levels = c(
                "No Pres",
                "DESS",
                "Buffer TL",
                "Lab Control"
              )
      ),
    sample_type = 
      factor (sample_type,
              levels = c(
                "sample",
                "extraction n.ctrl",
                "field n.ctrl",
                "pcr n.ctrl"
              )
      )
  )

df_stats <- 
  df %>% 
  select (
    sample_type,
    library_size,
    n_zotu
  ) %>% 
  mutate (
    grp = 
      case_when (
        sample_type == "sample" ~ "sample",
        T ~ "ctrl"
      )
  )

wilcox.test (library_size ~ grp, data = df_stats, exact = F)  
wilcox.test (n_zotu ~ grp, data = df_stats)


(summary_stats <- 
    df_stats %>% 
    group_by (grp) %>% 
    get_summary_stats () 
)


(plot_reads_zotus <- 
    df %>% 
    ggplot(
      aes (
        x = library_size,
        y = n_zotu,
        fill = sample_type,
        color = sample_type,
        shape = preservative
      )
    ) +
    geom_point(
      size = 3, 
      alpha = 0.75)
  +
    scale_y_log10() + 
    scale_x_log10() +
    labs (
      x = expression("No. of reads ("~log[10]~"scale)"),
      y = expression("No. of zOTUs ("~log[10]~"scale)"),
      fill = "Sample Type",
      color = "Sample Type",
      shape = "Preservative"
    ) +
    theme_bw() + 
    scale_color_brewer (
      palette = "Dark2"
    ) + 
    scale_fill_brewer (
      palette = "Dark2"
    ) + 
    
    scale_shape_manual(
      values = c(22, 24, 25, 21)
    )
)


# Plot composition

## Assign Color palette
# taxa_df <- phyloseq@tax_table %>% as.data.frame()
# tax_levels <- taxa_df$class %>% unique()
# n_colors <- length(tax_levels)
# 
# set.seed (999)
# cols <- 
#   Polychrome::createPalette(n_colors, c("#FF0000", "#00FF00", "#0000FF")) %>% 
#   colorspace::darken(amount = 0.1)
# names(cols) <- tax_levels
# scales::show_col(cols)

(plot_composition <- 
    ps %>%
    comp_barplot (
      tax_level = "class",
      n_taxa = 4,
      other_name = "Other",
      palette = distinct_palette (n = 4, add = "grey50"),
      merge_other = T,
      facet_by = "plot_group")+ 
    coord_flip() + 
    labs (
      y = "Relative Composition",
      fill = "Class"
    ) 
)



# Plot PCA Ordination
(plot_pca <- 
    ps%>% 
    tax_transform(
      trans = "rclr", # robust centered log-ratio transformation
      #rank = "class"
    ) %>% 
    ord_calc(
      method = "PCA"
    ) %>% 
    ord_plot(
      fill = "sample_type",
      col = "sample_type",
      shape = "preservative",
      size = 2,
      #plot_taxa = 1:4,
      plot_samples = T,
    ) +
    stat_ellipse(aes(colour = sample_type))+
    # ggrepel::geom_text_repel(
    #   aes(
    #     label = sample_id,
    #     col = sample_type
    #     ),
    #   size = 3
    # ) + 
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(
      values = c(22, 24, 25, 21)
    ) + 
    theme_bw() +
    labs (
      fill = "Sample Type",
      color = "Sample Type",
      shape = "Preservative"
    ) 
)


# Plot PCoA Ordination
(plot_pcoa <- 
    ps%>% 
    tax_transform(
      trans = "rclr", # robust centered log-ratio transformation
      #rank = "class"
    ) %>% 
    dist_calc("euclidean") %>% 
    ord_calc("PCoA") %>% 
    ord_plot(
      color = "sample_type",
      fill = "sample_type",
      shape = "preservative",
      size = 2
    ) +
    stat_ellipse(aes(colour = sample_type))+
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    scale_shape_manual(
      values = c(22, 24, 25, 21)
    ) + 
    theme_bw() +
    labs (
      fill = "Sample Type",
      color = "Sample Type",
      shape = "Preservative"
    ) 
)

# Merge plots
(plot_neg_ctrl_assessment_contaminants <- 
    ggpubr::ggarrange (
      plot_composition + 
        theme (legend.position = "bottom") +
        guides(fill = guide_legend(nrow = 2)),
      ggarrange (
        plot_reads_zotus +
          theme (legend.position = "bottom") +
          guides(color = guide_legend(nrow = 2),
                 shape = guide_legend(nrow = 2)), 
        plot_pca + 
          theme (legend.position = "bottom") +
          guides(color = guide_legend(nrow = 2),
                 shape = guide_legend(nrow = 2)),
        nrow = 2,
        labels = c("B", "C"),
        common.legend = T,
        legend = "bottom"
      ),
      ncol = 2,
      labels = c("A", NA),
      widths = c(1, 0.75)
    )
)





#### 09 COMPARE COMPOSITION AFTER DECONTAMINATION AND REMOVING NON-MARINE ZOTUS ####



# Assign phyloseq object to analyze
ps <- 
  ps_decontaminated %>% 
  
  # Aggregate samples based on grouping
  merge_samples("grp") %>% 
  
  # Remove taxid on tax table to allow aggregation
  tax_mutate (taxid = NULL) %>% 
  
  # Remove non-marine samples
  subset_taxa(habitat == "marine") %>%
  
  # Remove 0 observations
  prune_taxa(taxa_sums(.) > 0, .) %>% 
  prune_samples(sample_sums(.) > 0, .) %>% 
  
  # Add new info to sample data
  ps_mutate(
    sample_type = 
      case_when (
        grepl ("nec", rownames (sample_data(.))) ~ "extraction n.ctrl",
        grepl ("nfc", rownames (sample_data(.))) ~ "field n.ctrl",
        grepl ("ntc", rownames (sample_data(.))) ~ "pcr n.ctrl",
        T ~ "sample"
      ) %>% factor (
        levels = c(
          "sample",
          "extraction n.ctrl",
          "field n.ctrl",
          "pcr n.ctrl"
        )
      ),
    sample_id = 
      rownames(sample_data(.)),
    fraction = 
      case_when (
        grepl ("filter", rownames(sample_data(.))) ~ "filter",
        grepl ("filtrate", rownames(sample_data(.))) ~ "filtrate",
        grepl ("pcr", rownames(sample_data(.))) ~ "pcr",
        T ~ "sample"
      ),
    preservative = 
      case_when (
        grepl ("_cs", rownames(sample_data(.))) ~ "No Pres",
        grepl ("_ds", rownames(sample_data(.))) ~ "DESS",
        grepl ("_tl", rownames(sample_data(.))) ~ "Buffer TL",
        T ~ "Lab Control"
      ) %>% 
      factor (
        levels = c(
          "No Pres",
          "DESS",
          "Buffer TL",
          "Lab Control"
        )
      ),
    plot_group = ifelse(
      sample_type == "sample",
      as.character(preservative),
      "Neg. Ctrls"
    ) %>% 
      factor (
        levels = c(
          "No Pres",
          "DESS",
          "Buffer TL",
          "Neg. Ctrls"
        )
      ),
    permanova_grp = 
      case_when (
        sample_type != "sample" ~ "ctrl",
        T ~ "sample"
      )
  ) 

# Rename samples
new_names <- gsub("_cs(?=$|_)", "_np", sample_names(ps), perl = TRUE)
ps <- rename_samples(ps, new_names)

# Plot library size vs n_zotus
df <- 
  ps %>% 
  sample_data %>% 
  as("matrix") %>% 
  as.data.frame() %>% 
  mutate (
    library_size = sample_sums (ps),
    n_zotu =  estimate_richness(ps, measures = "Observed")[[1]],
    preservative = 
      factor (preservative, 
              levels = c(
                "No Pres",
                "DESS",
                "Buffer TL",
                "Lab Control"
              )
      ),
    sample_type = 
      factor (sample_type,
              levels = c(
                "sample",
                "extraction n.ctrl",
                "field n.ctrl",
                "pcr n.ctrl"
              )
      )
  )

df_stats <- 
  df %>% 
  select (
    sample_type,
    library_size,
    n_zotu
  ) %>% 
  mutate (
    grp = 
      case_when (
        sample_type == "sample" ~ "sample",
        T ~ "ctrl"
      )
  )

wilcox.test (library_size ~ grp, data = df_stats, exact = F)  
wilcox.test (n_zotu ~ grp, data = df_stats)


(summary_stats <- 
    df_stats %>% 
    group_by (grp) %>% 
    get_summary_stats () 
)


(plot_reads_zotus <- 
    df %>% 
    ggplot(
      aes (
        x = library_size,
        y = n_zotu,
        fill = sample_type,
        color = sample_type,
        shape = preservative
      )
    ) +
    geom_point(
      size = 3, 
      alpha = 0.75)
  +
    scale_y_log10() + 
    scale_x_log10() +
    labs (
      x = expression("No. of reads ("~log[10]~"scale)"),
      y = expression("No. of zOTUs ("~log[10]~"scale)"),
      fill = "Sample Type",
      color = "Sample Type",
      shape = "Preservative"
    ) +
    theme_bw() + 
    scale_color_brewer (
      palette = "Dark2"
    ) + 
    scale_fill_brewer (
      palette = "Dark2"
    ) + 
    
    scale_shape_manual(
      values = c(22, 24, 25, 21)
    )
)


# Plot composition

## Assign Color palette
# taxa_df <- phyloseq@tax_table %>% as.data.frame()
# tax_levels <- taxa_df$class %>% unique()
# n_colors <- length(tax_levels)
# 
# set.seed (999)
# cols <- 
#   Polychrome::createPalette(n_colors, c("#FF0000", "#00FF00", "#0000FF")) %>% 
#   colorspace::darken(amount = 0.1)
# names(cols) <- tax_levels
# scales::show_col(cols)

(plot_composition <- 
    ps %>%
    comp_barplot (
      tax_level = "class",
      n_taxa = 4,
      other_name = "Other",
      palette = distinct_palette (n = 4, add = "grey50"),
      merge_other = T,
      facet_by = "plot_group")+ 
    coord_flip() + 
    labs (
      y = "Relative Composition",
      fill = "Class"
    ) 
)



# Plot PCA Ordination
(plot_pca <- 
    ps%>% 
    tax_transform(
      trans = "rclr", # robust centered log-ratio transformation
      #rank = "class"
    ) %>% 
    ord_calc(
      method = "PCA"
    ) %>% 
    ord_plot(
      fill = "sample_type",
      col = "sample_type",
      shape = "preservative",
      size = 2,
      #plot_taxa = 1:4,
      plot_samples = T,
    ) +
    stat_ellipse(aes(colour = sample_type))+
    # ggrepel::geom_text_repel(
    #   aes(
    #     label = sample_id,
    #     col = sample_type
    #     ),
    #   size = 3
    # ) + 
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(
      values = c(22, 24, 25, 21)
    ) + 
    theme_bw() +
    labs (
      fill = "Sample Type",
      color = "Sample Type",
      shape = "Preservative"
    ) 
)

# Merge plots
(plot_neg_ctrl_assessment_post_decon_without_marine <- 
    ggpubr::ggarrange (
      plot_composition + 
        theme (legend.position = "bottom") +
        guides(fill = guide_legend(nrow = 2)),
      ggarrange (
        plot_reads_zotus +
          theme (legend.position = "bottom") +
          guides(color = guide_legend(nrow = 2),
                 shape = guide_legend(nrow = 2)), 
        plot_pca + 
          theme (legend.position = "bottom") +
          guides(color = guide_legend(nrow = 2),
                 shape = guide_legend(nrow = 2)),
        nrow = 2,
        labels = c("B", "C"),
        common.legend = T,
        legend = "bottom"
      ),
      ncol = 2,
      labels = c("A", NA),
      widths = c(0.75, 1)
    )
)

# Do PERMANOVA
permanova_sample_vs_ctrl <- 
  ps %>% 
  tax_transform(
    trans = "rclr"
  ) %>% 
  dist_calc("euclidean") %>% 
  dist_permanova(
    variables = "permanova_grp",
    n_perms = 10000,
    by = "term")


# Do PERMDISP
grp <- 
  ps %>% 
  sample_data() %>% 
  pull (permanova_grp)

dist_sample_vs_ctrl <- 
  ps %>%   
  tax_transform(
    trans = "rclr"
  ) %>% 
  dist_calc("euclidean") %>% 
  dist_get() 

permdisp_sample_vs_ctrl <- 
  betadisper (
    d = dist_sample_vs_ctrl,
    group = grp
  )

# Permutation test
permutest_permdisp_sample_ctrl <- permutest(permdisp_sample_vs_ctrl, permutations = 10000)
plot(permdisp_sample_vs_ctrl)
boxplot(permdisp_sample_vs_ctrl)
  
  
  
  #### 10 COMPARE COMPOSITION AMONG PRESERVATION TREATMENTS ONLY ####

# Plot PCA Ordination 
(plot_pca_preservation_treatments <- 
   
   ps %>% 
   
   subset_samples(
     sample_type == "sample"
   ) %>% 
   
   tax_transform(
     trans = "rclr", # robust centered log-ratio transformation
     #rank = "class"
   ) %>% 
   ord_calc(
     method = "PCA"
   ) %>% 
   ord_plot(
     fill = "preservative",
     col = "preservative",
     #shape = "preservative",
     size = 2,
     #plot_taxa = 1:4,
     plot_samples = T,
   ) +
   stat_ellipse(aes(colour = preservative))+
   # ggrepel::geom_text_repel(
   #   aes(
   #     label = sample_id,
   #     col = sample_type
   #     ),
   #   size = 3
   # ) + 
   
   
   scale_color_manual(
     values = c(
       "Buffer TL" = "#E7298A",
       "DESS" = "#7570B3",
       "No Pres" = "#D95F02"
     )
   ) +
   
   scale_fill_manual(
     values = c(
       "Buffer TL" = "#E7298A",
       "DESS" = "#7570B3",
       "No Pres" = "#D95F02"
     )
   ) +
   theme_bw() +
   labs (
     fill = "Preservative",
     color = "Preservative"
   ) 
)


# Do PERMANOVA
permanova_preservatives <- 
  ps %>% 
  subset_samples(
    sample_type == "sample"
  ) %>% 
  tax_transform(
    trans = "rclr"
  ) %>% 
  dist_calc("euclidean") %>% 
  dist_permanova(
    variables = "preservative",
    n_perms = 10000,
    by = "term")

# Do PERMDISP
grp <- 
  ps %>% 
  subset_samples(
    sample_type == "sample"
  ) %>% 
  sample_data() %>% 
  pull (preservative)

dist_preservative <- 
  ps %>%   
  subset_samples(
    sample_type == "sample"
  ) %>% 
  tax_transform(
    trans = "rclr"
  ) %>% 
  dist_calc("euclidean") %>% 
  dist_get() 

permdisp_preservative <- 
  betadisper (
    d = dist_preservative,
    group = grp
  )

# Permutation test
permutest_permdisp_preservative <- permutest(permdisp_preservative, permutations = 10000)
plot(permdisp_preservative)
boxplot(permdisp_preservative)



#### 10 EXPORT FILES ####
saveRDS(
  ps_decontaminated,
  path_phyloseq_decontaminated
)

ggsave (
  plot_neg_ctrl_assessment_pre_decon,
  file = path_plot_ordination_preDecon,
  width = 10.5,
  height = 7,
  units = "in",
  dpi = 330
)

ggsave (
  plot_neg_ctrl_assessment_post_decon,
  file = path_plot_ordination_postDecon,
  width = 10.5,
  height = 7,
  units = "in",
  dpi = 330
)

ggsave (
  plot_neg_ctrl_assessment_contaminants,
  file = path_plot_ordination_contam,
  width = 10.5,
  height = 7,
  units = "in",
  dpi = 330
)


ggsave (
  plot_neg_ctrl_assessment_post_decon_without_marine,
  file = path_plot_ordination_without_marine,
  width = 12,
  height = 7,
  units = "in",
  dpi = 330
)

ggsave (
  plot_pca_preservation_treatments,
  file = path_plot_ordination_samples_only,
  width = 7,
  height = 5,
  units = "in",
  dpi = 330
)

sink(path_permanova_samples)
cat ("### Sample VS Control ###")
cat ("\n")
print (permanova_sample_vs_ctrl)
cat ("\n")
print (permutest_permdisp_sample_ctrl)
cat ("\n\n### VS Preservatives ###")
cat ("\n")
print (permanova_preservatives)
cat ("\n")
print (permutest_permdisp_preservative)
sink()