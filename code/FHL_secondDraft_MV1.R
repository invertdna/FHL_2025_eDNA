library(tidyverse)
library(here)
library(vegan)
library(readxl)
library(microViz)
library(lubridate)
library(indicspecies)

a <- read.csv(here("data/FHL02_MV1_taxon_table_wide(2).csv")) #raw data
meta <- read.csv(here("data/Field_Data_Sheet_Digitized.csv")) #metadata

#Function for calculating the eDNA Index from Kelly et al. 2019
eDNA_index <- function(community_matrix){
  #community matrix with samples in cols, species in rows
  require(vegan)
  
  props <- decostand(community_matrix, method = "total", MARGIN = 2)
  
  idx <- decostand(props, method = "max", MARGIN = 1)
  idx <- as.data.frame(idx)
  colnames(idx) <- colnames(community_matrix)
  row.names(idx) <- row.names(community_matrix)
  
  return(idx)
}

#QC the raw data for min read depth and min occurrences in 10% of samples
b <- a %>% 
  dplyr::select(-2) %>% 
  drop_na() %>% 
  pivot_longer(cols = -BestTaxon) %>% 
  filter(value > 0) %>% 
  group_by(BestTaxon) %>% 
  add_tally() %>% 
  filter(n > 7) %>% ##FILTER for min occurrences
  dplyr::select(-n) %>% 
  pivot_wider(names_from = name, values_from = value, values_fill = 0) %>% 
  column_to_rownames("BestTaxon") %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(tot = rowSums(.)) %>% 
  filter(tot > 10000) %>% ##FILTER for total read-depth in sample
  dplyr::select(-tot) %>% 
  dplyr::select(which(colSums(.) > 0))

#get relevant metadata  
m <- meta %>% 
  dplyr::select(Tube.Number, Date, Time.Finish, Day.Night, Treatment..Light..NoLight..FieldNegative.) %>% 
  rename(Treatment = Treatment..Light..NoLight..FieldNegative.,
         Tube = Tube.Number,
         Time = Time.Finish) %>% 
  mutate(Date = mdy(Date),
         Time = hm(Time)) %>% 
  mutate(datetime = ymd_hms(paste(Date, Time)))
  
#make rownames consistent across data and metadata, etc
g <- b %>% 
  rownames_to_column("Tube") %>% 
  mutate(Tube = str_replace_all(Tube, pattern = "MV1\\.FHL", replacement = ""),
         Tube = str_replace_all(Tube, pattern = "_MV1", replacement = "")) %>% 
  mutate(Tube = as.numeric(Tube)) %>% 
  column_to_rownames("Tube")

#filter metadata to just samples we kept after QC
p <- m %>% 
  filter(Tube %in% rownames(g))
  p <- p[order(p$Tube),]
  p <- p %>% 
    unite(c(Date, Day.Night, Treatment), col = tmp, remove = F) %>% 
    mutate(time_idx = match(tmp, unique(tmp)))
  
#reorder to mirror metadata
g <- g[order(as.numeric(rownames(g))),] 

#check ordering
rownames(g) == p$Tube

#transform to eDNA index
g_idx <- g %>% 
  t() %>% 
  eDNA_index()



#calculate distance matrix
d <- g_idx %>% 
  t() %>% 
  vegdist(method = "bray")

##NMDS results; ordination of distance matrix just calculated
f <- d %>% 
  metaMDS(trymax = 1000) 

##plot NMDS
#check on biological triplicates
p %>% 
  # unite(c(Date, Day.Night, Treatment), col = tmp, remove = F) %>% 
  # mutate(time_idx = match(tmp, unique(tmp))) %>% 
  left_join(data.frame(Tube = as.numeric(colnames(g_idx)),
                       f$points), by = "Tube") %>%
  ggplot(aes(x = MDS1, y = MDS2, color = as.factor(time_idx))) +
  geom_point() +
  stat_chull(aes(fill = as.factor(time_idx)), alpha = 0.3) +
  theme_minimal()

#look for treatment effect
p %>% 
  left_join(data.frame(Tube = as.numeric(colnames(g_idx)),
                       f$points), by = "Tube") %>% 
  ggplot(aes(x = MDS1, y = MDS2, color = Day.Night, shape = Treatment)) +
    geom_point() +
    stat_chull(aes(fill = Day.Night), alpha = 0.3) +
    theme_minimal()


##PERMANOVA: what fraction of variance is explained by treatments?
res <- adonis2(as.matrix(d) ~ Day.Night + time_idx + Treatment %in% time_idx, by = "terms", data = p)
res

##Time Series
props <- g %>% 
  t() %>% 
  decostand(method = "total", MARGIN = 2) %>% 
  as.data.frame()

props %>% 
  rownames_to_column("BestTaxon") %>% 
  pivot_longer(-BestTaxon) %>% 
  mutate(name = as.numeric(name)) %>% 
  left_join(p, join_by(name == Tube)) %>%
  ggplot(aes(x = name, y = value, fill = BestTaxon)) +
    geom_col()  +
    xlab("Time point") + ylab("Proportion")
  # facet_wrap(~Day.Night)

##single-species analyses; time series plots
pdf("my_plots.pdf", width = 8, height = 6)
for (i in 1:nrow(g_idx)){
    p1 <- data.frame(p, 
               obs = unlist(g_idx[i,])) %>% 
      ggplot(aes(x = datetime, y = obs)) +
        geom_point() +
        ggtitle(rownames(g_idx)[i])
    print(p1)
}

for (i in 1:nrow(g_idx)){
  p2 <- data.frame(p, 
                   obs = unlist(g_idx[i,])) %>% 
    ggplot(aes(x = interaction(Day.Night, Treatment), y = obs)) +
    geom_boxplot() +
    ggtitle(rownames(g_idx)[i])
  print(p2)
}
dev.off()

##single-species analyses; IndVal; TODO need to look more into
#using eDNA index to start, then pres/abs
sc <- g %>% 
  # t() %>% 
  multipatt(cluster = interaction(p$Day.Night, p$Treatment), 
            func = "IndVal.g",
            permutations = 999
            )

summary(sc)
