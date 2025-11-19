# Libraries -----------------------------------------------------------------------------------
library(betareg)
library(vegan)
library(ggplot2)
library(dplyr)
library(here)
library(tibble)
library(tidyr)
library(ggpubr)

# Functions -----------------------------------------------------------------------------------
make_index <- function(data, index_variable, index_name) {
	data %>%
		group_by(across(all_of(index_variable))) %>%
		mutate(!!index_name := cur_group_id()) %>%
		ungroup()
}

# extract_param <- function (model, par) 
# {
# 	fit <- (methods::selectMethod("summary", signature = "stanfit"))(object = model, 
# 																																	 par = par)
# 	fit <- fit$summary
# 	return(fit %>% unlist() %>% as.data.frame %>% round(., 9))
# }


logit <- function(p) {
	if (any(p <= 0 | p >= 1)) stop("logit is only defined for 0 < p < 1")
	log(p / (1 - p))
}

inv_logit <- function(x) {
	1 / (1 + exp(-x))
}

estimate_length_of_pairs <- function(n){
	return((n*(n-1))/2)
}

design_gg_mat <- function(dist_matrix,metadata){
	
	dist_long <- dist_matrix %>% as.matrix() %>%
		as_tibble(rownames = "Tube_1") %>% 
		pivot_longer(-Tube_1, names_to = "Tube_2", values_to = "D") %>%
		mutate(Tube_1=as.numeric(Tube_1)) %>% 
		mutate(Tube_2=as.numeric(Tube_2)) %>% 
		filter(Tube_1 < Tube_2)
	
	pairs <- combn(metadata$Tube, 2)
	
	pairwise_df <- data.frame(
		Tube_1 = pairs[1, ],
		Tube_2 = pairs[2, ]
	)
	
	design_mat <- pairwise_df %>% 
		left_join(.,dist_long,by=c('Tube_1','Tube_2')) %>%
		left_join(.,
							metadata %>% select(Tube,Day.Night) %>% rename(time_1='Day.Night'),
							by=c('Tube_1'='Tube')) %>% 
		left_join(.,
							metadata %>% select(Tube,Day.Night) %>% rename(time_2='Day.Night'),
							by=c('Tube_2'='Tube')) %>% 
		mutate(time_cat=if_else(time_1==time_2,0,1)) %>% 
		left_join(.,
							metadata %>% select(Tube,Treatment) %>% rename(tre_1='Treatment'),
							by=c('Tube_1'='Tube')) %>%
		left_join(.,
							metadata %>% select(Tube,Treatment) %>% rename(tre_2='Treatment'),
							by=c('Tube_2'='Tube')) %>% 
		mutate(tre_cat=if_else(tre_1==tre_2,0,1)) %>% 
		left_join(.,
							metadata %>% select(Tube,Date) %>% make_index("Date", "date_idx") %>% rename(date_1='date_idx') %>% select(-Date),
							by=c('Tube_1'='Tube')) %>% 
		left_join(.,
							metadata %>% select(Tube,Date) %>% make_index("Date", "date_idx") %>% rename(date_2='date_idx') %>% select(-Date),
							by=c('Tube_2'='Tube')) %>% 
		mutate(delta_date=date_2-date_1) %>% 
		mutate(rep_cat=if_else(time_cat==0&tre_cat==0&delta_date==0,1,0)) %>% 
		filter(!(time_cat==1&tre_cat==1)) %>% #remove the combinations when Time of Day changes (Day v Night) and Treatment changes (Light v No light)
		mutate(beta_idx=ifelse(time_1=='night',2,1)) %>% 
		mutate(beta_idx=if_else(tre_cat==0,0,beta_idx)) %>% 
		mutate(light_day=if_else(beta_idx==1,1,0),
					 light_night=if_else(beta_idx==2,1,0)) %>% select(-beta_idx)
	# filter(!delta_date==7) 
	# select(Tube_1,Tube_2,D,time_cat,tre_cat,delta_date,rep_cat,beta_idx)
	
	return(design_mat)
}

pairwise_row_diff <- function(comm) {
	n <- nrow(comm)
	out <- list()
	
	for (i in 1:(n - 1)) {
		for (j in (i + 1):n) {
			diff_vec <- (comm[i, ] - comm[j, ])
			name <- paste(rownames(comm)[i], rownames(comm)[j], sep = "_vs_")
			out[[name]] <- diff_vec
		}
	}
	
	return(out)
}

# Load data -----------------------------------------------------------------------------------
coi_eDNAidx <- readRDS(here('Output','coi_cur_eDNAidx.rds'))
coi_bin <- readRDS(here('Output','coi_cur_pa.rds'))
coi_prop <- readRDS(here('Output','coi_cur_prop.rds'))
coi_prop_comm <- readRDS(here('Output','coi_cur_prop_common.rds'))
coi_meta <- readRDS(here('Output','coi_metadata.rds'))

mfu_eDNAidx <- readRDS(here('Output','mfu_cur_eDNAidx.rds'))
mfu_bin <- readRDS(here('Output','mfu_cur_pa.rds'))
mfu_prop <- readRDS(here('Output','mfu_cur_prop.rds'))
mfu_prop_comm <- readRDS(here('Output','mfu_cur_prop_common.rds'))
mfu_meta <- readRDS(here('Output','mfu_metadata.rds'))

mv1_eDNAidx <- readRDS(here('Output','mv1_cur_eDNAidx.rds'))
mv1_bin <- readRDS(here('Output','mv1_cur_pa.rds'))
mv1_prop <- readRDS(here('Output','mv1_cur_prop.rds'))
mv1_prop_comm <- readRDS(here('Output','mv1_cur_prop_common.rds'))
mv1_meta <- readRDS(here('Output','mv1_metadata.rds'))


dist_coi_eDNAidx <- coi_eDNAidx %>% t() %>% vegdist(method = "bray", binary = T)
dist_coi_bin <- coi_bin %>% t() %>% vegdist(method = "jaccard", binary = T)
dist_coi_prop <- coi_prop %>% t() %>% vegdist(method = "bray", binary = T)
dist_coi_prop_comm <- coi_prop_comm %>% t() %>% vegdist(method = "bray", binary = T)

dist_mfu_eDNAidx <- mfu_eDNAidx %>% t() %>% vegdist(method = "bray", binary = T)
dist_mfu_bin <- mfu_bin %>% t() %>% vegdist(method = "jaccard", binary = T)
dist_mfu_prop <- mfu_prop %>% t() %>% vegdist(method = "bray", binary = T)
dist_mfu_prop_comm <- mfu_prop_comm %>% t() %>% vegdist(method = "bray", binary = T)

dist_mv1_eDNAidx <- mv1_eDNAidx %>% t() %>% vegdist(method = "bray", binary = T)
dist_mv1_bin <- mv1_bin %>% t() %>% vegdist(method = "jaccard", binary = T)
dist_mv1_prop <- mv1_prop %>% t() %>% vegdist(method = "bray", binary = T)
dist_mv1_prop_comm <- mv1_prop_comm %>% t() %>% vegdist(method = "bray", binary = T)


# Create design matrix ------------------------------------------------------------------------
design_mat_coi_eDNAidx <- design_gg_mat(dist_matrix = dist_coi_eDNAidx,metadata = coi_meta)
design_mat_coi_bin <- design_gg_mat(dist_matrix = dist_coi_bin,metadata = coi_meta)
design_mat_coi_prop <- design_gg_mat(dist_matrix = dist_coi_prop,metadata = coi_meta)
design_mat_coi_prop_comm <- design_gg_mat(dist_matrix = dist_coi_prop_comm,metadata = coi_meta)

design_mat_mfu_eDNAidx <- design_gg_mat(dist_matrix = dist_mfu_eDNAidx,metadata = mfu_meta)
design_mat_mfu_bin <- design_gg_mat(dist_matrix = dist_mfu_bin,metadata = mfu_meta)
design_mat_mfu_prop <- design_gg_mat(dist_matrix = dist_mfu_prop,metadata = mfu_meta)
design_mat_mfu_prop_comm <- design_gg_mat(dist_matrix = dist_mfu_prop_comm,metadata = mfu_meta)

design_mat_mv1_eDNAidx <- design_gg_mat(dist_matrix = dist_mv1_eDNAidx,metadata = mv1_meta)
design_mat_mv1_bin <- design_gg_mat(dist_matrix = dist_mv1_bin,metadata = mv1_meta)
design_mat_mv1_prop <- design_gg_mat(dist_matrix = dist_mv1_prop,metadata = mv1_meta)
design_mat_mv1_prop_comm <- design_gg_mat(dist_matrix = dist_mv1_prop_comm,metadata = mv1_meta)

estimate_length_of_pairs(ncol(coi_eDNAidx)) #The total amount of pairwise comparisons between samples
nrow(design_mat_coi_eDNAidx) #The amount of pairwise comparisons retained and factorized (given estimable factors)

#The amount of pairwise comparisons removed
estimate_length_of_pairs(ncol(coi_eDNAidx))-nrow(design_mat_coi_eDNAidx)

#Create a list for all the design matricies so that it is easier to do for-loops downstream
design_mat_list <- list(design_mat_coi_eDNAidx,
												design_mat_coi_bin,
												design_mat_coi_prop,
												# design_mat_coi_prop_comm,
												design_mat_mfu_eDNAidx,
												design_mat_mfu_bin,
												design_mat_mfu_prop,
												# design_mat_mfu_prop_comm,
												design_mat_mv1_eDNAidx,
												design_mat_mv1_bin,
												design_mat_mv1_prop #design_mat_mv1_prop_comm
)


# Analysis ------------------------------------------------------------------------------------

TrophicLevel_COI <- read.csv(here('Output','COI_TAX.csv'))
valid_species_COI <- TrophicLevel_COI %>% filter(!is.na(TrophicLevel)) %>% pull(OriginalName) 
comm <- coi_eDNAidx %>% filter(rownames(.)%in%valid_species_COI) %>% #select marine species
	mutate(total_occurrences=rowSums(. > 0)) %>% #Filter only species that have more than 10 occurrences
	filter(total_occurrences>10) %>% select(-total_occurrences)

# Compute the species-specific pairwise difference
pair_sp_diff <- pairwise_row_diff(t(comm))

# Transform the table in a long format and join the trophic level
pair_sp_diff_long <- bind_cols(pair_sp_diff) %>% 
	cbind(TrophicLevel_COI %>% filter(OriginalName%in%rownames(comm)) %>% select(CleanedName,TrophicLevel)) %>% 
	rename(Species='CleanedName') %>% 
	pivot_longer(cols = -c('Species','TrophicLevel')) %>% 
	separate(name, into = c("Tube_1", "Tube_2"), sep = "_vs_") %>% 
	mutate(Tube_1=as.numeric(Tube_1)) %>% 
	mutate(Tube_2=as.numeric(Tube_2)) 

# Old plot
# pair_sp_diff_long_meta <- pair_sp_diff_long %>% 
# 	left_join(.,design_mat_coi_eDNAidx,by=c('Tube_1','Tube_2')) %>% 
# 	filter(!is.na(D)) %>% 
# 	mutate(ToD=paste0(time_1,'_',time_2)) %>% 
# 	rename(Day_v_night='time_cat') %>% 
# 	rename(Date_difference='delta_date') %>% 
# 	rename(Light_v_Nolight='tre_cat') %>% 
# 	rename(Biological_replciate='rep_cat') 
# 
# pair_sp_diff_long_meta %>% 
# 	group_by(ToD,Light_v_Nolight,light_day,light_night,Biological_replciate,Species,TrophicLevel) %>% 
# 	summarize(eDNA_idx_change=mean(value)) %>% 
# 	ungroup() %>% 
# 	filter(Biological_replciate==0) %>% select(-Biological_replciate) %>% 
# 	filter(eDNA_idx_change< -0.05|eDNA_idx_change> 0.05) %>%
# 	arrange(TrophicLevel) %>% as.data.frame()
# 
# sp_sig <- pair_sp_diff_long_meta %>% 
# 	# filter(TrophicLevel=='Planktivores') %>%
# 	mutate(Treatment=if_else(Day_v_night==1,paste0(time_1,'_',time_2),NA)) %>% 
# 	mutate(Treatment=if_else(Light_v_Nolight==1,paste0(time_1,'-',tre_1,'_',tre_2),Treatment)) %>% 
# 	mutate(value=if_else(Treatment=='night-NoLight_Light',value*-1,value)) %>% 
# 	mutate(Treatment=if_else(Treatment=='night-Light_NoLight','night-NoLight_Light',Treatment)) %>% 
# 	mutate(value=if_else(Treatment=='day-NoLight_Light',value*-1,value)) %>% 
# 	mutate(Treatment=if_else(Treatment=='day-Light_NoLight','day-NoLight_Light',Treatment)) %>% 
# 	group_by(Treatment,Species) %>% 
# 	summarize(eDNA_idx_change=mean(value)) %>% 
# 	filter(eDNA_idx_change< -0.1|eDNA_idx_change> 0.1) %>% as.data.frame() %>% 
# 	pull(Species) %>% unique()
# 
# ch_data <- pair_sp_diff_long_meta %>%
# 	# filter(TrophicLevel=='Planktivores') %>%
# 	filter(Species%in%sp_sig) %>%
# 	mutate(Treatment=if_else(Day_v_night==1,paste0(time_1,'-',time_2),NA)) %>% 
# 	mutate(Treatment=if_else(Light_v_Nolight==1,paste0(tre_1,'-',tre_2,'(',time_1,')'),Treatment)) %>% 
# 	mutate(value=if_else(Treatment=='NoLight-Light(night)',value*-1,value)) %>% 
# 	mutate(Treatment=if_else(Treatment=='NoLight-Light(night)','Light-NoLight(night)',Treatment)) %>% 
# 	mutate(value=if_else(Treatment=='NoLight-Light(day)',value*-1,value)) %>% 
# 	mutate(Treatment=if_else(Treatment=='NoLight-Light(day)','Light-NoLight(day)',Treatment)) %>% 
# 	mutate(Treatment=if_else(Treatment=='night-day','Night-Day',Treatment)) %>% 
# 	mutate(Treatment=if_else(Treatment=='day-night','Day-Night',Treatment)) %>% 
# 	filter(!is.na(Treatment)) %>%
# 	mutate(Treatment = factor(Treatment,
# 														levels = c(
# 															'Light-NoLight(night)',
# 															'Night-Day',
# 															'Day-Night',
# 															'Light-NoLight(day)'
# 														))) %>% 
# 	filter(Date_difference<2)


# New plot
ch_data <-
	pair_sp_diff_long %>% 
	left_join(.,design_mat_coi_eDNAidx,by=c('Tube_1','Tube_2')) %>% 
	filter(!is.na(D)) %>% 
	rename(tod_1='time_1',tod_2='time_2') %>% 
	rename(l_1='tre_1',l_2='tre_2') %>% 
	mutate(Treatment=paste0(tod_1,'-',tod_2)) %>% 
	mutate(treatment=if_else(tod_1!=tod_2,paste0(tod_1,'-',tod_2),NA)) %>% 
	mutate(treatment=if_else(tod_1==tod_2&l_1!=l_2,paste0(l_1,'-',l_2,'(',tod_1,')'),treatment)) %>% 
	filter(delta_date<2) %>% 
	mutate(value=if_else(treatment=='night-day',value*-1,value)) %>% 
	mutate(treatment=if_else(treatment=='night-day','day-night',treatment)) %>% 
	mutate(value=if_else(treatment=='NoLight-Light(day)',value*-1,value)) %>% 
	mutate(treatment=if_else(treatment=='NoLight-Light(day)','Light-NoLight(day)',treatment)) %>% 
	mutate(value=if_else(treatment=='NoLight-Light(night)',value*-1,value)) %>% 
	mutate(treatment=if_else(treatment=='NoLight-Light(night)','Light-NoLight(night)',treatment)) %>% 
	filter(!is.na(treatment)) %>% 
	filter(!(tod_1!=tod_2&l_1!=l_2)) %>% 
	filter(!(delta_date==1&tod_1==tod_2)) %>% 
	mutate(treatment=if_else(treatment=='Light-NoLight(night)','Light-No_Light(night)',treatment),
				 treatment=if_else(treatment=='Light-NoLight(day)','Light-No_Light(day)',treatment),
				 treatment=if_else(treatment=='day-night','Day-Night',treatment)) %>% 
	rename(Day_v_night='time_cat') %>% 
	rename(Date_difference='delta_date') %>% 
	rename(Light_v_Nolight='tre_cat') %>% 
	rename(Biological_replciate='rep_cat') %>% 
	select(-Treatment) %>% rename(Treatment='treatment')

sp_sig <- ch_data %>% 
	group_by(Species,Treatment) %>% 
	summarise(eDNA_idx_change=mean(value)) %>% 
	filter(eDNA_idx_change< -0.1|eDNA_idx_change> 0.1) %>% 
	ungroup() %>% distinct(Species) %>% 
	pull(Species)


p0 <-
	ch_data %>%
	filter(Species%in%sp_sig) %>% 
	ggline(x='Treatment',y='value',add = c('mean_se'),color='TrophicLevel')+
	labs(y = 'eDNA index change')+
	geom_abline(intercept=0,slope=0,lty=2)+
	annotate("text",x = Inf,y = 0,label = "no change",hjust = 1.1,vjust = -0.5)
ggsave(here('plots','sp_changes','Trophic Level response.jpg'),p0)

p1 <- ch_data %>%
	filter(Species%in%sp_sig) %>% 
	filter(TrophicLevel=='Planktivores') %>% 
	ggline(x='Treatment',y='value',add = c('mean_se'),color='Species')+
	labs(y = 'eDNA index change')+
	ggtitle('Planktivores')+
	geom_abline(intercept=0,slope=0,lty=2)+
	annotate("text",x = Inf,y = 0,label = "no change",hjust = 1.1,vjust = -0.5)
ggsave(here('plots','sp_changes','Planktivores.jpg'),p1)

p2 <- ch_data %>%
	filter(Species%in%sp_sig) %>% 
	filter(TrophicLevel=='Gelatinous zooplankton') %>% 
	ggline(x='Treatment',y='value',add = c('mean_se'),color='Species')+
	labs(y = 'eDNA index change')+
	ggtitle('Gelatinous zooplankton')+
	geom_abline(intercept=0,slope=0,lty=2)+
	annotate("text",x = Inf,y = 0,label = "no change",hjust = 1.1,vjust = -0.5)
ggsave(here('plots','sp_changes','Gelatinous zooplankto.jpg'),p2)


p3 <- ch_data %>%
	filter(Species%in%sp_sig) %>% 
	filter(TrophicLevel=='Filter feeders') %>% 
	ggline(x='Treatment',y='value',add = c('mean_se'),color='Species')+
	labs(y = 'eDNA index change')+
	ggtitle('Filter feeders')+
	geom_abline(intercept=0,slope=0,lty=2)+
	annotate("text",x = Inf,y = 0,label = "no change",hjust = 1.1,vjust = -0.5)
ggsave(here('plots','sp_changes','Filter feeders.jpg'),p3)

p4 <- ch_data %>%
	filter(Species%in%sp_sig) %>% 
	filter(TrophicLevel=='Benthic invertebrates') %>% 
	ggline(x='Treatment',y='value',add = c('mean_se'),color='Species')+
	labs(y = 'eDNA index change')+
	ggtitle('Benthic invertebrates')+
	geom_abline(intercept=0,slope=0,lty=2)+
	annotate("text",x = Inf,y = 0,label = "no change",hjust = 1.1,vjust = -0.5)
ggsave(here('plots','sp_changes','Benthic invertebrates.jpg'),p4)

p5 <- ch_data %>%
	filter(Species%in%sp_sig) %>% 
	filter(TrophicLevel=='Primary producers') %>% 
	ggline(x='Treatment',y='value',add = c('mean_se'),color='Species')+
	labs(y = 'eDNA index change')+
	ggtitle('Primary producers')+
	geom_abline(intercept=0,slope=0,lty=2)+
	annotate("text",x = Inf,y = 0,label = "no change",hjust = 1.1,vjust = -0.5)
ggsave(here('plots','sp_changes','Primary producers.jpg'),p5)

p6 <- ch_data %>%
	filter(Species%in%sp_sig) %>% 
	filter(TrophicLevel=='Forage species') %>% 
	ggline(x='Treatment',y='value',add = c('mean_se'),color='Species')+
	labs(y = 'eDNA index change')+
	ggtitle('Forage species')+
	geom_abline(intercept=0,slope=0,lty=2)+
	annotate("text",x = Inf,y = 0,label = "no change",hjust = 1.1,vjust = -0.5)
ggsave(here('plots','sp_changes','Forage species.jpg'),p6)

p7 <- ch_data %>%
	filter(Species%in%sp_sig) %>% 
	filter(TrophicLevel=='Parasites') %>% 
	ggline(x='Treatment',y='value',add = c('mean_se'),color='Species')+
	labs(y = 'eDNA index change')+
	geom_abline(intercept=0,slope=0,lty=2)+
	ggtitle('Parasites')+
	annotate("text",x = Inf,y = 0,label = "no change",hjust = 1.1,vjust = -0.5)
ggsave(here('plots','sp_changes','Parasites.jpg'),p7)

p8 <- ch_data %>%
	filter(Species%in%sp_sig) %>% 
	filter(TrophicLevel=='Pathogens') %>% 
	ggline(x='Treatment',y='value',add = c('mean_se'),color='Species')+
	labs(y = 'eDNA index change')+
	ggtitle('Pathogens')+
	geom_abline(intercept=0,slope=0,lty=2)+
	annotate("text",x = Inf,y = 0,label = "no change",hjust = 1.1,vjust = -0.5)
ggsave(here('plots','sp_changes','Pathogens.jpg'),p8)


