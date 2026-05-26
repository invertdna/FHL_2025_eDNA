# ======================================================================
# FUNCTIONAL ANNOTATION WORKFLOW
# Marker: COI
#
# Purpose:
#   1. Load eDNA index data and metadata
#   2. Validate and clean species names
#   3. Retrieve taxonomy from WoRMS
#   4. Retain marine taxa only
#   5. Assign trophic and vernacular groups
#
# Inputs:
#   - coi_cur_eDNAidx.rds
#   - coi_metadata.rds
#
# Outputs:
#   - coi_long_marine_with_taxonomy.rds
#   - taxonomy_lookup_details.csv
#   - coi_unknown_groups.csv
#
# Author: Beng
# Date: 2025-11-02
# ======================================================================


# ======================================================================
# LOAD REQUIRED PACKAGES
# ======================================================================

library(tidyverse)
library(worrms)
library(progress)


# ======================================================================
# CONFIGURATION
# ======================================================================

config <- list(
  abundance_file = "coi_cur_eDNAidx.rds",
  metadata_file  = "coi_metadata.rds",
  taxonomy_cache = "coi_taxonomy_cache.csv"
)


# ======================================================================
# SECTION 1 — LOAD AND PREPARE DATA
# ======================================================================

coi_data <- readRDS(config$abundance_file)
coi_meta <- readRDS(config$metadata_file)

# Ensure species names are stored in a dedicated column
if (!"Species" %in% colnames(coi_data)) {
  
  coi_data <- coi_data %>%
    as.data.frame() %>%
    rownames_to_column("Species")
}

# Convert abundance table to long format
coi_long <- coi_data %>%
  pivot_longer(
    cols      = -Species,
    names_to  = "Tube",
    values_to = "Proportion"
  )

# Standardize Tube IDs
coi_long <- coi_long %>%
  mutate(Tube = as.character(Tube))

coi_meta <- coi_meta %>%
  mutate(Tube = as.character(Tube))

# Merge metadata
coi_long <- coi_long %>%
  left_join(coi_meta, by = "Tube")

# Remove rows with missing metadata
coi_long <- coi_long %>%
  filter(
    !is.na(Day.Night),
    !is.na(Treatment)
  )


# ======================================================================
# SECTION 2 — CLEAN SPECIES NAMES
# ======================================================================

clean_species_name <- function(name) {
  
  cleaned <- name %>%
    str_to_lower() %>%
    str_remove("^uncultured ") %>%
    str_remove(" sp\\..*$") %>%
    str_remove(" cf\\..*$") %>%
    str_remove(" aff\\..*$") %>%
    str_remove(" [A-Z]{2,}[-_][0-9]+.*$") %>%
    str_remove(" [0-9]+$") %>%
    str_remove(" \\(.*\\)") %>%
    str_remove(" var\\..*$") %>%
    str_remove(" subsp\\..*$") %>%
    str_remove("\\?") %>%
    str_remove_all("[\\[\\]]") %>%
    str_squish()

  words <- str_split(cleaned, " ")[[1]]

  if (length(words) >= 2) {
    
    cleaned <- paste0(
      str_to_title(words[1]),
      " ",
      words[2]
    )
    
  } else {
    
    cleaned <- str_to_title(cleaned)
  }

  return(cleaned)
}


# ======================================================================
# SECTION 3 — TAXONOMY LOOKUP USING WoRMS
# ======================================================================

# Load cache if available
if (file.exists(config$taxonomy_cache)) {
  
  tax_cache <- read_csv(
    config$taxonomy_cache,
    show_col_types = FALSE
  )
  
} else {
  
  tax_cache <- tibble(
    OriginalName = character(),
    CleanedName  = character(),
    Family       = character(),
    Class        = character(),
    Phylum       = character(),
    Kingdom      = character(),
    LookupStatus = character()
  )
}

# Species requiring lookup
unique_species <- unique(coi_long$Species)

species_to_lookup <- setdiff(
  unique_species,
  tax_cache$OriginalName
)

# Progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent | :current/:total",
  total  = length(species_to_lookup),
  clear  = FALSE,
  width  = 60
)

# Function for WoRMS lookup
lookup_species <- function(sp) {

  pb$tick()

  cleaned <- clean_species_name(sp)

  result <- tryCatch({

    rec <- wm_records_name(
      name        = cleaned,
      fuzzy       = TRUE,
      marine_only = TRUE
    )

    if (is.data.frame(rec) && nrow(rec) > 0) {

      tibble(
        OriginalName = sp,
        CleanedName  = cleaned,
        Family       = rec$family[1],
        Class        = rec$class[1],
        Phylum       = rec$phylum[1],
        Kingdom      = rec$kingdom[1],
        LookupStatus = "match"
      )

    } else {

      tibble(
        OriginalName = sp,
        CleanedName  = cleaned,
        Family       = NA_character_,
        Class        = NA_character_,
        Phylum       = NA_character_,
        Kingdom      = NA_character_,
        LookupStatus = "not_found"
      )
    }

  }, error = function(e) {

    tibble(
      OriginalName = sp,
      CleanedName  = cleaned,
      Family       = NA_character_,
      Class        = NA_character_,
      Phylum       = NA_character_,
      Kingdom      = NA_character_,
      LookupStatus = "error"
    )
  })

  Sys.sleep(0.1)

  return(result)
}

# Query WoRMS
if (length(species_to_lookup) > 0) {

  new_taxonomy <- map_dfr(
    species_to_lookup,
    lookup_species
  )

  tax_cache <- bind_rows(
    tax_cache,
    new_taxonomy
  )

  write_csv(
    tax_cache,
    config$taxonomy_cache
  )
}


# ======================================================================
# SECTION 4 — FILTER MARINE TAXA
# ======================================================================

marine_taxa <- tax_cache %>%
  filter(LookupStatus == "match")

coi_long_marine <- coi_long %>%
  left_join(
    marine_taxa,
    by = c("Species" = "OriginalName")
  ) %>%
  filter(!is.na(Phylum))

message(
  paste(
    "Marine taxa retained:",
    n_distinct(coi_long_marine$Species)
  )
)


# ======================================================================
# SECTION 5 — TROPHIC GROUP ASSIGNMENT
# ======================================================================

assign_trophic_group <- function(sp, fam, cls, phy) {

  sp  <- tolower(ifelse(is.na(sp),  "", sp))
  fam <- tolower(ifelse(is.na(fam), "", fam))
  cls <- tolower(ifelse(is.na(cls), "", cls))
  phy <- tolower(ifelse(is.na(phy), "", phy))

  # Primary producers
  if (phy %in% c(
    "rhodophyta",
    "chlorophyta",
    "ochrophyta",
    "bacillariophyta",
    "haptophyta",
    "cryptophyta"
  )) {
    return("Primary producers")
  }

  # Apex predators
  if (
    fam %in% c(
      "salmonidae",
      "lamnidae",
      "scombridae"
    )
  ) {
    return("Apex predators")
  }

  # Filter feeders
  if (
    cls %in% c(
      "bivalvia",
      "ascidiacea"
    )
  ) {
    return("Filter feeders")
  }

  # Planktivores
  if (
    cls %in% c(
      "copepoda",
      "maxillopoda"
    )
  ) {
    return("Planktivores")
  }

  # Benthic invertebrates
  if (
    phy %in% c(
      "annelida",
      "echinodermata",
      "nemertea"
    )
  ) {
    return("Benthic invertebrates")
  }

  # Mesopredators
  if (
    cls %in% c(
      "actinopteri",
      "teleostei"
    )
  ) {
    return("Mesopredators")
  }

  return("Unknown")
}


# Apply trophic classification
coi_long_marine <- coi_long_marine %>%
  rowwise() %>%
  mutate(
    Trophic_Group = assign_trophic_group(
      Species,
      Family,
      Class,
      Phylum
    )
  ) %>%
  ungroup()


# ======================================================================
# SECTION 6 — VERNACULAR GROUP ASSIGNMENT
# ======================================================================

coi_long_marine <- coi_long_marine %>%
  mutate(
    Common_Group = case_when(

      # Fish
      Phylum == "Chordata" &
        Class %in% c(
          "Actinopteri",
          "Actinopterygii",
          "Teleostei"
        ) ~ "Ray-finned fishes",

      # Crustaceans
      Phylum == "Arthropoda" &
        Class == "Malacostraca" ~ "Crustaceans",

      # Polychaetes
      Phylum == "Annelida" ~ "Polychaetes",

      # Jellyfish
      Phylum == "Cnidaria" &
        Class %in% c(
          "Scyphozoa",
          "Hydrozoa"
        ) ~ "Jellyfish",

      # Phytoplankton
      Phylum %in% c(
        "Bacillariophyta",
        "Dinoflagellata",
        "Haptophyta",
        "Cryptophyta"
      ) ~ "Phytoplankton",

      # Macroalgae
      Phylum %in% c(
        "Rhodophyta",
        "Chlorophyta",
        "Ochrophyta"
      ) ~ "Macroalgae",

      TRUE ~ "Other"
    )
  )


# ======================================================================
# SECTION 7 — EXPORT RESULTS
# ======================================================================

# Save annotated dataset
saveRDS(
  coi_long_marine,
  "coi_long_marine_with_taxonomy.rds"
)

# Save taxonomy table
write_csv(
  tax_cache,
  "taxonomy_lookup_details.csv"
)

# Export unknown groups
unknown_groups <- coi_long_marine %>%
  filter(
    Trophic_Group == "Unknown" |
      Common_Group == "Other"
  ) %>%
  distinct(
    Species,
    Family,
    Class,
    Phylum,
    Trophic_Group,
    Common_Group
  )

write_csv(
  unknown_groups,
  "coi_unknown_groups.csv"
)


# ======================================================================
# FUNCTIONAL ANNOTATION WORKFLOW
# Marker: 12S (formerly MV1 / MarVer1)
#
# Purpose:
#   1. Load eDNA index data and metadata
#   2. Validate and clean species names
#   3. Retrieve taxonomy from WoRMS
#   4. Retain marine taxa only
#   5. Assign trophic and vernacular groups
#
# Inputs:
#   - mv1_cur_eDNAidx.rds
#   - mv1_metadata.rds
#
# Outputs:
#   - mv1_long_marine_with_taxonomy.rds
#   - taxonomy_lookup_details.csv
#   - mv1_unknown_groups.csv
#
# Author: Beng
# Date: 2025-11-02
# ======================================================================


# ======================================================================
# LOAD REQUIRED PACKAGES
# ======================================================================

library(tidyverse)
library(worrms)
library(progress)


# ======================================================================
# CONFIGURATION
# ======================================================================

config <- list(
  abundance_file = "mv1_cur_eDNAidx.rds",
  metadata_file  = "mv1_metadata.rds",
  taxonomy_cache = "mv1_taxonomy_cache.csv"
)


# ======================================================================
# SECTION 1 — LOAD AND PREPARE DATA
# ======================================================================

mv1_data <- readRDS(config$abundance_file)
mv1_meta <- readRDS(config$metadata_file)

# Ensure species names are stored in a dedicated column
if (!"Species" %in% colnames(mv1_data)) {
  
  mv1_data <- mv1_data %>%
    as.data.frame() %>%
    rownames_to_column("Species")
}

# Convert abundance table to long format
mv1_long <- mv1_data %>%
  pivot_longer(
    cols      = -Species,
    names_to  = "Tube",
    values_to = "Proportion"
  )

# Standardize Tube IDs
mv1_long <- mv1_long %>%
  mutate(Tube = as.character(Tube))

mv1_meta <- mv1_meta %>%
  mutate(Tube = as.character(Tube))

# Merge metadata
mv1_long <- mv1_long %>%
  left_join(mv1_meta, by = "Tube")

# Remove rows with missing metadata
mv1_long <- mv1_long %>%
  filter(
    !is.na(Day.Night),
    !is.na(Treatment)
  )


# ======================================================================
# SECTION 2 — CLEAN SPECIES NAMES
# ======================================================================

clean_species_name <- function(name) {
  
  cleaned <- name %>%
    str_to_lower() %>%
    str_remove("^uncultured ") %>%
    str_remove(" sp\\..*$") %>%
    str_remove(" cf\\..*$") %>%
    str_remove(" aff\\..*$") %>%
    str_remove(" [A-Z]{2,}[-_][0-9]+.*$") %>%
    str_remove(" [0-9]+$") %>%
    str_remove(" \\(.*\\)") %>%
    str_remove(" var\\..*$") %>%
    str_remove(" subsp\\..*$") %>%
    str_remove("\\?") %>%
    str_remove_all("[\\[\\]]") %>%
    str_squish()

  words <- str_split(cleaned, " ")[[1]]

  if (length(words) >= 2) {
    
    cleaned <- paste0(
      str_to_title(words[1]),
      " ",
      words[2]
    )
    
  } else {
    
    cleaned <- str_to_title(cleaned)
  }

  return(cleaned)
}


# ======================================================================
# SECTION 3 — TAXONOMY LOOKUP USING WoRMS
# ======================================================================

# Load cache if available
if (file.exists(config$taxonomy_cache)) {
  
  tax_cache <- read_csv(
    config$taxonomy_cache,
    show_col_types = FALSE
  )
  
} else {
  
  tax_cache <- tibble(
    OriginalName = character(),
    CleanedName  = character(),
    Family       = character(),
    Class        = character(),
    Phylum       = character(),
    Kingdom      = character(),
    LookupStatus = character()
  )
}

# Species requiring lookup
unique_species <- unique(mv1_long$Species)

species_to_lookup <- setdiff(
  unique_species,
  tax_cache$OriginalName
)

# Progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent | :current/:total",
  total  = length(species_to_lookup),
  clear  = FALSE,
  width  = 60
)

# Function for WoRMS lookup
lookup_species <- function(sp) {

  pb$tick()

  cleaned <- clean_species_name(sp)

  result <- tryCatch({

    rec <- wm_records_name(
      name        = cleaned,
      fuzzy       = TRUE,
      marine_only = TRUE
    )

    if (is.data.frame(rec) && nrow(rec) > 0) {

      tibble(
        OriginalName = sp,
        CleanedName  = cleaned,
        Family       = rec$family[1],
        Class        = rec$class[1],
        Phylum       = rec$phylum[1],
        Kingdom      = rec$kingdom[1],
        LookupStatus = "match"
      )

    } else {

      tibble(
        OriginalName = sp,
        CleanedName  = cleaned,
        Family       = NA_character_,
        Class        = NA_character_,
        Phylum       = NA_character_,
        Kingdom      = NA_character_,
        LookupStatus = "not_found"
      )
    }

  }, error = function(e) {

    tibble(
      OriginalName = sp,
      CleanedName  = cleaned,
      Family       = NA_character_,
      Class        = NA_character_,
      Phylum       = NA_character_,
      Kingdom      = NA_character_,
      LookupStatus = "error"
    )
  })

  Sys.sleep(0.1)

  return(result)
}

# Query WoRMS
if (length(species_to_lookup) > 0) {

  new_taxonomy <- map_dfr(
    species_to_lookup,
    lookup_species
  )

  tax_cache <- bind_rows(
    tax_cache,
    new_taxonomy
  )

  write_csv(
    tax_cache,
    config$taxonomy_cache
  )
}


# ======================================================================
# SECTION 4 — FILTER MARINE TAXA
# ======================================================================

marine_taxa <- tax_cache %>%
  filter(LookupStatus == "match")

mv1_long_marine <- mv1_long %>%
  left_join(
    marine_taxa,
    by = c("Species" = "OriginalName")
  ) %>%
  filter(!is.na(Phylum))

message(
  paste(
    "Marine taxa retained:",
    n_distinct(mv1_long_marine$Species)
  )
)


# ======================================================================
# SECTION 5 — TROPHIC GROUP ASSIGNMENT
# ======================================================================

assign_trophic_group <- function(sp, fam, cls, phy) {

  sp  <- tolower(ifelse(is.na(sp),  "", sp))
  fam <- tolower(ifelse(is.na(fam), "", fam))
  cls <- tolower(ifelse(is.na(cls), "", cls))
  phy <- tolower(ifelse(is.na(phy), "", phy))

  # Primary producers
  if (phy %in% c(
    "rhodophyta",
    "chlorophyta",
    "ochrophyta",
    "bacillariophyta",
    "haptophyta",
    "cryptophyta"
  )) {
    return("Primary producers")
  }

  # Apex predators
  if (
    fam %in% c(
      "salmonidae",
      "lamnidae",
      "scombridae"
    )
  ) {
    return("Apex predators")
  }

  # Filter feeders
  if (
    cls %in% c(
      "bivalvia",
      "ascidiacea"
    )
  ) {
    return("Filter feeders")
  }

  # Planktivores
  if (
    cls %in% c(
      "copepoda",
      "maxillopoda"
    )
  ) {
    return("Planktivores")
  }

  # Benthic invertebrates
  if (
    phy %in% c(
      "annelida",
      "echinodermata",
      "nemertea"
    )
  ) {
    return("Benthic invertebrates")
  }

  # Mesopredators
  if (
    cls %in% c(
      "actinopteri",
      "teleostei"
    )
  ) {
    return("Mesopredators")
  }

  return("Unknown")
}


# Apply trophic classification
mv1_long_marine <- mv1_long_marine %>%
  rowwise() %>%
  mutate(
    Trophic_Group = assign_trophic_group(
      Species,
      Family,
      Class,
      Phylum
    )
  ) %>%
  ungroup()


# ======================================================================
# SECTION 6 — VERNACULAR GROUP ASSIGNMENT
# ======================================================================

mv1_long_marine <- mv1_long_marine %>%
  mutate(
    Common_Group = case_when(

      # Fish
      Phylum == "Chordata" &
        Class %in% c(
          "Actinopteri",
          "Actinopterygii",
          "Teleostei"
        ) ~ "Ray-finned fishes",

      # Crustaceans
      Phylum == "Arthropoda" &
        Class == "Malacostraca" ~ "Crustaceans",

      # Polychaetes
      Phylum == "Annelida" ~ "Polychaetes",

      # Jellyfish
      Phylum == "Cnidaria" &
        Class %in% c(
          "Scyphozoa",
          "Hydrozoa"
        ) ~ "Jellyfish",

      # Phytoplankton
      Phylum %in% c(
        "Bacillariophyta",
        "Dinoflagellata",
        "Haptophyta",
        "Cryptophyta"
      ) ~ "Phytoplankton",

      # Macroalgae
      Phylum %in% c(
        "Rhodophyta",
        "Chlorophyta",
        "Ochrophyta"
      ) ~ "Macroalgae",

      TRUE ~ "Other"
    )
  )


# ======================================================================
# SECTION 7 — EXPORT RESULTS
# ======================================================================

# Save annotated dataset
saveRDS(
  mv1_long_marine,
  "mv1_long_marine_with_taxonomy.rds"
)

# Save taxonomy table
write_csv(
  tax_cache,
  "taxonomy_lookup_details.csv"
)

# Export unknown groups
unknown_groups <- mv1_long_marine %>%
  filter(
    Trophic_Group == "Unknown" |
      Common_Group == "Other"
  ) %>%
  distinct(
    Species,
    Family,
    Class,
    Phylum,
    Trophic_Group,
    Common_Group
  )

write_csv(
  unknown_groups,
  "mv1_unknown_groups.csv"
)
