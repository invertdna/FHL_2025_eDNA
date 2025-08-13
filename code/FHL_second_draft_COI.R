library(tidyverse)
library(here)
library(vegan)
library(readxl)
library(microViz)
library(lubridate)

a <- read.csv(here("data/FHL04_LoLe_taxon_table_wide.csv")) #raw data
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


#QC the raw data for min read depth
b <- a %>% 
  dplyr::select(-2) %>% 
  drop_na() %>% 
  pivot_longer(cols = -BestTaxon) %>% 
  filter(value > 0) %>% 
  group_by(BestTaxon) %>% 
  add_tally() %>% 
  filter(n > 5) %>% ##FILTER for min occurrences
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
  dplyr::select(Tube.Number, Date, Time.Finish, Day.Night, Treatment..Light..NoLight..FieldNegative., Temperature, Salinity) %>% 
  rename(Treatment = Treatment..Light..NoLight..FieldNegative.,
         Tube = Tube.Number,
         Time = Time.Finish) %>% 
  mutate(Date = mdy(Date),
         Time = hm(Time)) %>% 
  mutate(datetime = ymd_hms(paste(Date, Time)))
  
#make rownames consistent across data and metadata, etc
g <- b %>% 
  rownames_to_column("Tube") %>% 
  mutate(Tube = str_replace_all(Tube, pattern = "LoLe\\.FHL", replacement = ""),
         Tube = str_replace_all(Tube, pattern = "\\.COI.+", replacement = "")) %>% 
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

#make pres/absence dataset
bin = g
  bin[bin>0]<-1  
#calculate distance matrix, using robust aitchison distance to deal w compositionality
d <- bin %>% 
  vegdist(method = "jaccard", binary = T)

#visualize
  heatmap(as.matrix(d))

##NMDS results; ordination of distance matrix just calculated
f <- d %>% 
  metaMDS(trymax = 1000) 

##plot NMDS
#check on biological triplicates
p %>% 
  unite(c(Date, Day.Night, Treatment), col = tmp, remove = F) %>% 
  mutate(time_idx = match(tmp, unique(tmp))) %>% 
  bind_cols(f$points) %>% 
  ggplot(aes(x = MDS1, y = MDS2, color = as.factor(time_idx))) +
  geom_point() +
  stat_chull(aes(fill = as.factor(time_idx)), alpha = 0.3) +
  theme_minimal()

#look for treatment effect
p %>% 
  bind_cols(f$points) %>% 
  ggplot(aes(x = MDS1, y = MDS2, color = Day.Night, shape = Treatment)) +
    geom_point() +
    stat_chull(aes(fill = Day.Night), alpha = 0.3) +
    theme_minimal()


##PERMANOVA
res <- adonis2(as.matrix(d) ~ Day.Night + time_idx + Treatment %in% time_idx, by = "terms", data = p)
res



##Time Series
common <- b %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(tot = rowSums(.)) %>% 
  filter(tot > 500) %>% 
  dplyr::select(-tot) %>% 
  decostand(method = "total", MARGIN = 2)

common %>% 
  rownames_to_column("BestTaxon") %>% 
  pivot_longer(-BestTaxon) %>% 
  ggplot(aes(x = name, y = value, fill = BestTaxon)) +
    geom_col() + xlab("") + ylab("Proportion")

d_props <- common %>% 
  t() %>% 
  vegdist(method = "robust.aitchison")

#visualize
heatmap(as.matrix(d_props))

##NMDS results; ordination of distance matrix just calculated
f <- d_props %>% 
  metaMDS(trymax = 1000) 
p %>% 
  bind_cols(f$points) %>% 
  ggplot(aes(x = MDS1, y = MDS2, color = Day.Night, shape = Treatment)) +
  geom_point() +
  stat_chull(aes(fill = Day.Night), alpha = 0.3) +
  theme_minimal()



#visualize day/night differences in proportion for each species
# data.frame(p,
#            t(common)
# ) %>% 
#   pivot_longer(cols = -c(1:10)) %>% 
#   ggplot(aes(x = Day.Night, y = value)) +
#   geom_boxplot() +
#   facet_wrap(~name, scales = "free_y")

## if we drop Alitta, to see about other patterns
# common %>% 
#   rownames_to_column("BestTaxon") %>% 
#   filter(BestTaxon != "Alitta") %>% 
#   pivot_longer(-BestTaxon) %>% 
#   group_by(name) %>% 
#   mutate(value = value/sum(value)) %>%  #renormalize into proportions
#   ggplot(aes(x = name, y = value, fill = BestTaxon)) +
#   geom_col() 
#   # theme(legend.position = "none")


##Finally, we'll break down and do richness over time:
R <- rowSums(bin)
data.frame(R, p) %>% 
  ggplot(aes(x = datetime, y = R, color = Day.Night)) +
    geom_point() +
    geom_smooth()

data.frame(R, p) %>% 
  ggplot(aes(x = interaction(Day.Night, Treatment), y = R, color = Day.Night)) +
  geom_boxplot()

data.frame(R, p) %>% 
  ggplot(aes(x = Salinity, y = R, color = Day.Night)) +
    geom_point()

data.frame(R, p) %>% 
  ggplot(aes(x = Salinity, y = Temperature, color = Day.Night)) +
  geom_point()

data.frame(R, p) %>% 
  ggplot(aes(x = Day.Night, y = Salinity, color = Day.Night)) +
  geom_boxplot()

#temp and sal correlated, obv; TODO look at tide
lm(Temperature ~ Salinity, data = p) %>% summary()

##Group-specific trends:
w <- g %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("BestTaxon") %>% 
  left_join(a %>% dplyr::select(BestTaxon, Class)) %>% 
  relocate(c(1,ncol(.)))
  w[3:ncol(w)][w[3:ncol(w)]>0]<-1  #make binary

w %>% 
  drop_na() %>% 
  filter(Class != "") %>% 
  pivot_longer(-c(BestTaxon, Class), names_to = "Tube") %>% 
  dplyr::select(-1) %>% #drop species
  group_by(Class, Tube) %>% 
  summarise(R = sum(value)) %>% 
  mutate(Tube = as.numeric(Tube)) %>% 
  left_join(p) %>% 
  ggplot(aes(x = interaction(Day.Night, Treatment), y = R)) +
    geom_boxplot() +
    facet_wrap(~Class) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
plots1 <- list()
plots2 <- list()
richness <- list()
for (i in 1:length(unique(w$Class))){
  
  tmp <- w %>% 
    drop_na() %>% 
    filter(Class == unique(w$Class)[i]) %>% 
    dplyr::select(-2) %>% #drop Class
    column_to_rownames("BestTaxon") %>% 
    mutate(tot = rowSums(.)) %>% 
    filter(tot > 0) %>% 
    dplyr::select(-tot) %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(tot = rowSums(.)) %>% 
    filter(tot > 0) %>% 
    dplyr::select(-tot)
  
  if(ncol(tmp)>10){
  
  assign(paste0(unique(w$Class)[i],"_nmds"),
     tmp %>% 
      metaMDS(distance = "jaccard", trymax = 1000)
  )
  
  richness[[i]] <- m[match(rownames(tmp),m$Tube),] %>% 
    unite(c(Date, Day.Night, Treatment), col = tmp, remove = F) %>% 
    mutate(time_idx = match(tmp, unique(tmp))) 
  richness[[i]]$Richness = as.vector(rowSums(tmp))
  
  plots1[[i]] <-  richness[[i]] %>% 
    ggplot(aes(x = interaction(Day.Night, Treatment), y = Richness, fill = Day.Night)) +
      geom_boxplot() +
    theme_minimal() +
    ggtitle(unique(w$Class)[i])
    
  plots2[[i]] <-  m[match(rownames(tmp),m$Tube),] %>% 
    unite(c(Date, Day.Night, Treatment), col = tmp, remove = F) %>% 
    mutate(time_idx = match(tmp, unique(tmp))) %>% 
    bind_cols(eval(sym(paste0(unique(w$Class)[i],"_nmds")))$points) %>% 
    ggplot(aes(x = MDS1, y = MDS2, color = Day.Night, shape = Treatment)) +
    geom_point() +
    stat_chull(aes(fill = Day.Night), alpha = 0.3) +
    theme_minimal() +
    ggtitle(unique(w$Class)[i])
  }
  }
plots1
plots2

###some of the above analyses, but using eDNA Index 

g_idx <- g %>% 
  t() %>% 
  eDNA_index()

p %>% 
  bind_cols(t(g_idx)) %>% 
  ggplot(aes(x = datetime, y = Alitta)) +
    geom_point() +
    geom_line() +
    ggtitle("Polychaete Spawning")

d_idx <- g_idx %>% 
  t() %>% 
  vegdist(method = "bray")

#visualize
heatmap(as.matrix(d_idx))

##NMDS results; ordination of distance matrix just calculated
f_idx <- d_idx %>% 
  metaMDS(trymax = 1000) 
p %>% 
  bind_cols(f_idx$points) %>% 
  ggplot(aes(x = MDS1, y = MDS2, color = Day.Night, shape = Treatment)) +
  geom_point() +
  stat_chull(aes(fill = Day.Night), alpha = 0.3) +
  theme_minimal()
