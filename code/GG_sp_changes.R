# Libraries -----------------------------------------------------------------------------------
library(betareg)
library(vegan)
library(ggplot2)
library(dplyr)
library(here)
library(tibble)
library(tidyr)
library(ggpubr)
# library(MoMAColors)

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
coi_raw <- readRDS(here('Output','coi_cur_raw.rds'))
coi_eDNAidx <- readRDS(here('Output','coi_cur_eDNAidx.rds'))
coi_bin <- readRDS(here('Output','coi_cur_pa.rds'))
coi_prop <- readRDS(here('Output','coi_cur_prop.rds'))
coi_prop_comm <- readRDS(here('Output','coi_cur_prop_common.rds'))
coi_meta <- readRDS(here('Output','coi_metadata.rds'))

mfu_raw <- readRDS(here('Output','mfu_cur_raw.rds'))
mfu_eDNAidx <- readRDS(here('Output','mfu_cur_eDNAidx.rds'))
mfu_bin <- readRDS(here('Output','mfu_cur_pa.rds'))
mfu_prop <- readRDS(here('Output','mfu_cur_prop.rds'))
mfu_prop_comm <- readRDS(here('Output','mfu_cur_prop_common.rds'))
mfu_meta <- readRDS(here('Output','mfu_metadata.rds'))

mv1_raw <- readRDS(here('Output','mv1_cur_raw.rds'))
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
TrophicLevel_MV1 <- read.csv(here('Output','MV1_TAX.csv'))
valid_species_COI <- TrophicLevel_COI %>% filter(!is.na(TrophicLevel)) %>% pull(OriginalName) 
valid_species_MV1 <- TrophicLevel_MV1 %>% filter(!is.na(TrophicLevel)) %>% pull(OriginalName) 
comm_coi <- coi_eDNAidx %>% filter(rownames(.)%in%valid_species_COI) %>% #select marine species
	mutate(total_occurrences=rowSums(. > 0)) %>% #Filter only species that have more than 10 occurrences
	filter(total_occurrences>10) %>% select(-total_occurrences)
comm_mv1 <- mv1_eDNAidx %>% filter(rownames(.)%in%valid_species_MV1) %>% #select marine species
	mutate(total_occurrences=rowSums(. > 0)) %>% #Filter only species that have more than 10 occurrences
	filter(total_occurrences>10) %>% select(-total_occurrences)

# Compute the species-specific pairwise difference
pair_sp_diff_coi <- pairwise_row_diff(t(comm_coi))
pair_sp_diff_mv1 <- pairwise_row_diff(t(comm_mv1))

# Transform the table in a long format and join the trophic level
pair_sp_diff_long_coi <- bind_cols(pair_sp_diff_coi) %>% 
	cbind(TrophicLevel_COI %>% filter(OriginalName%in%rownames(comm_coi)) %>% select(CleanedName,TrophicLevel)) %>% 
	rename(Species='CleanedName') %>% 
	pivot_longer(cols = -c('Species','TrophicLevel')) %>% 
	separate(name, into = c("Tube_1", "Tube_2"), sep = "_vs_") %>% 
	mutate(Tube_1=as.numeric(Tube_1)) %>% 
	mutate(Tube_2=as.numeric(Tube_2)) 

pair_sp_diff_long_mv1 <- bind_cols(pair_sp_diff_mv1) %>% 
	cbind(TrophicLevel_MV1 %>% filter(OriginalName%in%rownames(comm_mv1)) %>% select(CleanedName,TrophicLevel)) %>% 
	rename(Species='CleanedName') %>% 
	pivot_longer(cols = -c('Species','TrophicLevel')) %>% 
	separate(name, into = c("Tube_1", "Tube_2"), sep = "_vs_") %>% 
	mutate(Tube_1=as.numeric(Tube_1)) %>% 
	mutate(Tube_2=as.numeric(Tube_2)) 


# New plot
ch_data_coi <-
	pair_sp_diff_long_coi %>% 
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

ch_data_mv1 <-
	pair_sp_diff_long_mv1 %>% 
	left_join(.,design_mat_mv1_eDNAidx,by=c('Tube_1','Tube_2')) %>% 
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

sp_sig_coi <- ch_data_coi %>% 
	group_by(Species,Treatment) %>% 
	summarise(eDNA_idx_change=mean(value)) %>% 
	filter(eDNA_idx_change< -0.1|eDNA_idx_change> 0.12) %>% 
	ungroup() %>% distinct(Species) %>% 
	pull(Species)

sp_sig_mv1 <- ch_data_mv1 %>% 
	group_by(Species,Treatment) %>% 
	summarise(eDNA_idx_change=mean(value)) %>% 
	filter(eDNA_idx_change< -0.1|eDNA_idx_change> 0.12) %>% 
	ungroup() %>% distinct(Species) %>% 
	pull(Species)

# 
# response_group <- ch_data %>% 
# 	group_by(Species,Treatment) %>% 
# 	summarise(eDNA_idx_change=mean(value)) %>% 
# 	filter(eDNA_idx_change< -0.1|eDNA_idx_change> 0.12)
# 
# # sp_sig <- response_group %>% ungroup() %>% distinct(Species) %>% 
# # 	pull(Species)
# 
# 
# resp_sp <- response_group %>%
# 	mutate(response_g=if_else(Treatment=='Day-Night'&eDNA_idx_change>0.1,'1',NA)) %>% 
# 	mutate(response_g=if_else(Treatment=='Day-Night'&eDNA_idx_change< -0.1,'2',response_g)) %>% 
# 	mutate(response_g=if_else(Treatment=='Light-No_Light(night)'&eDNA_idx_change> 0.1,'3',response_g)) %>% 
# 	mutate(response_g=if_else(Treatment=='Light-No_Light(night)'&eDNA_idx_change< -0.1,'4',response_g)) %>% 
# 	mutate(response_g=if_else(Treatment=='Light-No_Light(day)'&eDNA_idx_change> 0.1,'5',response_g)) %>% 
# 	mutate(response_g=if_else(Treatment=='Light-No_Light(day)'&eDNA_idx_change< -0.1,'6',response_g)) %>% 
# 	group_by(Species) %>% 
# 	summarise(all_responses = paste(response_g, collapse = ", ")) %>% 
# 	ungroup() %>% 
# 	mutate(Response=if_else(all_responses=='1','Day-thriving',NA)) %>% 
# 	mutate(Response=if_else(all_responses=='1, 3','Day-thriving & light-sensitive',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='1, 6, 4','Day-thriving & light-sensitive',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='1, 6, 3','Day-thriving & light-sensitive',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='1, 5, 4','Day-thriving & light-sensitive',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='1, 5, 3','Day-thriving & light-sensitive',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='1, 5','Day-thriving & light-sensitive',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='1, 4','Day-thriving & light-sensitive',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='1, 6','Day-thriving & light-sensitive',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='4','Light-shy',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='6','Light-shy',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='2','Night-thriving',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='2, 3','Night-thriving',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='2, 4','Night-thriving',Response)) %>% 
# 	mutate(Response=if_else(all_responses=='2, 6, 4','Night-thriving',Response)) 
# 
# 
# plot_data <- ch_data %>% left_join(.,resp_sp,by='Species') %>% 
# 	ungroup() %>% 
# 	filter(!is.na(TrophicLevel)) %>% 
# 	filter(Species%in%sp_sig) %>%
# 	filter(!is.na(Response)) %>%
# 	group_by(Species,TrophicLevel,Response,Treatment) %>% 
# 	summarise(mean_eDNA_ch=mean(value)) %>% 
# 	mutate(Treatment=factor(Treatment,levels = c('Light-No_Light(night)','Day-Night','Light-No_Light(day)'))) %>% 
# 	ungroup() %>% 
# 	mutate(TrophicLevel=if_else(!(TrophicLevel%in%trophic_order),'Other',TrophicLevel))
# 
# 
# trophic_order <- c("Benthic invertebrates", 
# 									 "Gelatinous zooplankton",
# 									 # "Forage species",
# 									 "Planktivores",
# 									 "Filter feeders",
# 									 "Primary producers",
# 									 'Other')
# 
# tl_sp <- plot_data %>% 
# 	mutate(TrophicLevel=if_else(!(TrophicLevel%in%trophic_order),'Other',TrophicLevel)) %>% 
# 	distinct(Species,TrophicLevel) %>%
# 	mutate(TrophicLevel = factor(TrophicLevel, levels = trophic_order))
# 
# # Create the species order grouped by trophic level
# species_order <- tl_sp %>%
# 	arrange(TrophicLevel, Species) %>%
# 	pull(Species)
# 
# tl_sp %>%
# 	arrange(TrophicLevel, Species) %>% 
# 	count(TrophicLevel)
# 
# tl_summ <- tl_sp %>%
# 	arrange(TrophicLevel, Species) %>% 
# 	mutate(row_num = row_number()) %>%
# 	group_by(TrophicLevel) %>%
# 	summarise(
# 		st_row = min(row_num),
# 		end_row = max(row_num),
# 		n_species = n()
# 	) 
# 
# 
# #Benthic
# benthic_color <- sample(moma.colors('Lupi',n=300)[61:110],tl_summ$n_species[tl_summ$TrophicLevel=='Benthic invertebrates'])
# 	# c("#AF4B66","#BA6177", "#BD687C", "#C06F82", "#C47688","#CE8C98")
# 
# #Gelatinous zooplankton 3
# gel_color <- sample(moma.colors('Lupi',n=300)[271:300],tl_summ$n_species[tl_summ$TrophicLevel=='Gelatinous zooplankton'])
# 	# c("#E37447", "#DD6745", "#D85A44")
# 
# #Forage species 2
# oth_color <- sample(moma.colors('Lupi',n=300)[250:260],tl_sp %>% count(TrophicLevel) %>% filter(TrophicLevel=='Other') %>% pull(n))
# 	# c("#E3B15B", "#FEB550")
# 
# #Planktivores 2
# plank_color <- sample(moma.colors('Lupi',n=300)[205:222],tl_sp %>% count(TrophicLevel) %>% filter(TrophicLevel=='Planktivores') %>% pull(n))
# 	# c("#3292A0", "#4F9794")
# 
# #Filter feeders 1
# filt_color <- sample(moma.colors('Lupi',n=300)[1:45],tl_sp %>% count(TrophicLevel) %>% filter(TrophicLevel=='Filter feeders') %>% pull(n))
# 	# c("#61BEA4")
# 
# #Primary producers 34
# prim_color <- sample(moma.colors('Lupi',n=300)[c(135:201,226:245)],tl_sp %>% count(TrophicLevel) %>% filter(TrophicLevel=='Primary producers') %>% pull(n))
# 	# c("#C0A588" ,"#B0A572" ,"#AEA56E" ,"#ACA56B" ,"#A9A568" ,"#A7A564" ,"#A5A561" ,"#A2A55E" ,"#A0A55B" ,"#9EA557" ,"#9BA554" ,"#99A551" ,"#96A44F" ,"#93A452" ,"#8FA355" ,"#8BA258" ,"#88A25B" ,"#84A15E" ,"#80A061" ,"#7CA064" ,"#799F67" ,"#759E6A" ,"#719E6D" ,"#6D9D6F" ,"#6A9C72" ,"#669C75" ,"#769E85" ,"#7D9F83" ,"#85A180" ,"#8CA27D" ,"#93A37A" ,"#9BA477" ,"#A2A674" ,"#A9A771" )
# 
# sp_color <- c('#AF4B66',benthic_color,
# 							'#E37447',gel_color,
# 							"#3292A0",plank_color,
# 							"#61BEA4",filt_color,
# 							"#9BA554",prim_color,
# 							'#E3B15B',oth_color)
# 
# color_breaks <- c(trophic_order[1],species_order[tl_summ$st_row[tl_summ$TrophicLevel==trophic_order[1]]:tl_summ$end_row[tl_summ$TrophicLevel==trophic_order[1]]],
# 									trophic_order[2],species_order[tl_summ$st_row[tl_summ$TrophicLevel==trophic_order[2]]:tl_summ$end_row[tl_summ$TrophicLevel==trophic_order[2]]],
# 									trophic_order[3],species_order[tl_summ$st_row[tl_summ$TrophicLevel==trophic_order[3]]:tl_summ$end_row[tl_summ$TrophicLevel==trophic_order[3]]],
# 									trophic_order[4],species_order[tl_summ$st_row[tl_summ$TrophicLevel==trophic_order[4]]:tl_summ$end_row[tl_summ$TrophicLevel==trophic_order[4]]],
# 									trophic_order[5],species_order[tl_summ$st_row[tl_summ$TrophicLevel==trophic_order[5]]:tl_summ$end_row[tl_summ$TrophicLevel==trophic_order[5]]],
# 									trophic_order[6],species_order[tl_summ$st_row[tl_summ$TrophicLevel==trophic_order[6]]:tl_summ$end_row[tl_summ$TrophicLevel==trophic_order[6]]])
# 
# plot_data %>% ungroup() %>% 
# ggplot()+
# 	geom_point(aes(x=Treatment,y=mean_eDNA_ch,colour = Species),alpha=0.5)+
# 	geom_line(aes(x=Treatment,y=mean_eDNA_ch,colour = Species,group = Species), linewidth = 0.4, alpha = 0.6,lty=2)+
# 	geom_point(data=plot_data %>% group_by(TrophicLevel,Response,Treatment) %>%
# 						 	summarise(mean_eDNA_ch=mean(mean_eDNA_ch)),
# 						 aes(x=Treatment,y=mean_eDNA_ch,colour = TrophicLevel),size=4)+
# 	geom_line(data=plot_data %>% group_by(TrophicLevel,Response,Treatment) %>%
# 						 	summarise(mean_eDNA_ch=mean(mean_eDNA_ch)),
# 						 aes(x=Treatment,y=mean_eDNA_ch,colour = TrophicLevel,group = TrophicLevel),size=1)+
# 	scale_color_manual(values=sp_color, breaks = color_breaks)+
# 	facet_grid(Response~TrophicLevel)+
# 	ylab('eDNA index change')+
# 	guides(color = guide_legend(ncol = 3))+
# 	geom_abline(intercept=0,slope=0,lty=2)+
# 	theme_bw()+
# 	theme(
# 		axis.text.x = element_text(angle = 50, hjust = 1)
# 	)



trophic_color_coi <- c('#AF4B66',#benthic_color,
									 '#E3B15B',#for_color,
									 '#E37447',#gel_color,
									 "#3292A0",#plank_color,
									 "#61BEA4",#filt_color,
									 "#9BA554"#prim_color,
									 )

names(trophic_color_coi) <- c('Benthic invertebrates',
													'Forage species',
													'Gelatinous zooplankton',
													'Planktivores',
													'Filter feeders',
													'Primary producers')

plot_data_coi <- ch_data_coi %>% 
	filter(!is.na(TrophicLevel)) %>% 
	group_by(Species,TrophicLevel,Treatment) %>% 
	summarise(mean_eDNA_ch=mean(value)) %>% 
	mutate(Treatment=factor(Treatment,levels = c('Light-No_Light(night)','Day-Night','Light-No_Light(day)'))) %>% 
	ungroup() %>% 
	filter(TrophicLevel%in%names(trophic_color_coi))

trophic_color_mv1 <- c('#AF4B66',#Apex predators,
									 '#E37447',#Forage species,
									 "#3292A0"#Mesopredators
									 )

names(trophic_color_mv1) <- c('Apex predators',
													'Forage species',
													'Mesopredators')

plot_data_2_mv1 <- ch_data_mv1 %>% 
	filter(!is.na(TrophicLevel)) %>% 
	group_by(Species,TrophicLevel,Treatment) %>% 
	summarise(mean_eDNA_ch=mean(value)) %>% 
	mutate(Treatment=factor(Treatment,levels = c('Light-No_Light(night)','Day-Night','Light-No_Light(day)'))) %>% 
	ungroup() %>% 
	filter(TrophicLevel%in%names(trophic_color_mv1))

p1 <- ggplot()+
	geom_boxplot(data=plot_data_coi,
						 aes(x=Treatment,y=mean_eDNA_ch,colour = TrophicLevel),outlier.shape = NA,size=0.8,coef = 0)+
	geom_point(data=plot_data_coi %>% filter(!(Species%in%sp_sig_coi)),aes(x=Treatment,y=mean_eDNA_ch),color='grey',alpha=0.5)+
	geom_line(data=plot_data_coi %>% filter(!(Species%in%sp_sig_coi)),aes(x=Treatment,y=mean_eDNA_ch,group = Species),color='grey',alpha=0.6)+
	geom_point(data=plot_data_coi %>% filter(Species%in%sp_sig_coi),
						 aes(x=Treatment,y=mean_eDNA_ch,colour = TrophicLevel),alpha=0.5)+
	geom_line(data=plot_data_coi %>% filter(Species%in%sp_sig_coi),
						aes(x=Treatment,y=mean_eDNA_ch,colour = TrophicLevel,group = Species),alpha=0.6)+
	facet_wrap(~TrophicLevel,ncol=2)+
	scale_color_manual(values=trophic_color_coi)+
	labs(y='eDNA index change (mean of all pairwise comparisons)',x='Pairwise comparison (treatment)')+
	geom_abline(intercept=0,slope=0,lty=2)+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
				legend.position = 'none',
				axis.title = element_text(size = 15),
				axis.text = element_text(size = 14),
				strip.text = element_text(size=15))
# ggsave(here('plots','Figure 2_coi.jpg'),p1,width=8,heigh=12)

p2 <-
ggplot()+
	geom_boxplot(data=plot_data_2_mv1,
						 aes(x=Treatment,y=mean_eDNA_ch,colour = TrophicLevel),outlier.shape = NA,size=0.8,coef = 0)+
	geom_point(data=plot_data_2_mv1 %>% filter(!(Species%in%sp_sig_mv1)),aes(x=Treatment,y=mean_eDNA_ch),color='grey',alpha=0.5)+
	geom_line(data=plot_data_2_mv1 %>% filter(!(Species%in%sp_sig_mv1)),aes(x=Treatment,y=mean_eDNA_ch,group = Species),color='grey',alpha=0.6)+
	geom_point(data=plot_data_2_mv1 %>% filter(Species%in%sp_sig_mv1),
						 aes(x=Treatment,y=mean_eDNA_ch,colour = TrophicLevel),alpha=0.5)+
	geom_line(data=plot_data_2_mv1 %>% filter(Species%in%sp_sig_mv1),
						aes(x=Treatment,y=mean_eDNA_ch,colour = TrophicLevel,group = Species),alpha=0.6)+
	facet_wrap(~TrophicLevel,ncol=1)+
	scale_color_manual(values=trophic_color_mv1)+
	labs(y='eDNA index change (mean of all pairwise comparisons)',x='Pairwise comparison (treatment)')+
	geom_abline(intercept=0,slope=0,lty=2)+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
				legend.position = 'none',
				axis.title = element_text(size = 15),
				axis.text = element_text(size = 14),
				strip.text = element_text(size=15))
# ggsave(here('plots','Figure 2_mv1.jpg'),p2,width=6,heigh=12)




trophic_color_comb <- c('#AF4B66',#benthic,
											 '#E37447',#gel,
											 "#3292A0",#plank,
											 '#E3B15B',#for_co1,
											 "#61BEA4",#filt,
											 "#9BA554",#prim,
											 '#E3B15B',#for_mv1,
											 '#053060',#apex
											 '#7bb7d5'#meso
)


names(trophic_color_comb) <- c('Benthic invertebrates (COI)',
															'Gelatinous zooplankton (COI)',
															'Planktivores (COI)',
															'Forage species (COI)',
															'Filter feeders (COI)',
															'Primary producers (COI)',
															'Forage species (12S)',
															'Apex predators (12S)',
															'Mesopredators (12S)')

plot_data_3 <- ch_data_coi %>% 
	filter(!is.na(TrophicLevel)) %>% 
	group_by(Species,TrophicLevel,Treatment) %>% 
	summarise(mean_eDNA_ch=mean(value)) %>% 
	mutate(Treatment=if_else(Treatment=='Light-No_Light(night)','Light - NoLight (night)',Treatment)) %>% 
	mutate(Treatment=if_else(Treatment=='Day-Night','Day - Night',Treatment)) %>% 
	mutate(Treatment=if_else(Treatment=='Light-No_Light(day)','Light - NoLight (day)',Treatment)) %>% 
	mutate(Treatment=factor(Treatment,levels = c('Light - NoLight (night)','Day - Night','Light - NoLight (day)'))) %>% 
	# mutate(Treatment=factor(Treatment,levels = c('Light-No_Light(night)','Day-Night','Light-No_Light(day)'))) %>% 
	ungroup() %>% 
	filter(TrophicLevel%in%names(trophic_color_coi)) %>% mutate(Marker='COI') %>% 
	rbind(.,
				ch_data_mv1 %>% 
					filter(!is.na(TrophicLevel)) %>% 
					group_by(Species,TrophicLevel,Treatment) %>% 
					summarise(mean_eDNA_ch=mean(value)) %>% 
					mutate(Treatment=if_else(Treatment=='Light-No_Light(night)','Light - NoLight (night)',Treatment)) %>% 
					mutate(Treatment=if_else(Treatment=='Day-Night','Day - Night',Treatment)) %>% 
					mutate(Treatment=if_else(Treatment=='Light-No_Light(day)','Light - NoLight (day)',Treatment)) %>% 
					mutate(Treatment=factor(Treatment,levels = c('Light - NoLight (night)','Day - Night','Light - NoLight (day)'))) %>% 
					# mutate(Treatment=factor(Treatment,levels = c('Light-No_Light(night)','Day-Night','Light-No_Light(day)'))) %>% 
					ungroup() %>% 
					filter(TrophicLevel%in%names(trophic_color_mv1)) %>% mutate(Marker='12S')
	) %>% 
	mutate(TrophicLevel=paste0(TrophicLevel,' (',Marker,')')) %>% 
	mutate(TrophicLevel=factor(TrophicLevel,levels = names(trophic_color_comb)))

p3 <-
ggplot()+
	geom_boxplot(data=plot_data_3,
							 aes(x=Treatment,y=mean_eDNA_ch,colour = TrophicLevel),outlier.shape = NA,size=0.8,coef = 0)+
	geom_point(data=plot_data_3 %>% filter(!(Species%in%c(sp_sig_mv1,sp_sig_coi))),aes(x=Treatment,y=mean_eDNA_ch),color='grey',alpha=0.5)+
	geom_line(data=plot_data_3 %>% filter(!(Species%in%c(sp_sig_mv1,sp_sig_coi))),aes(x=Treatment,y=mean_eDNA_ch,group = Species),color='grey',alpha=0.6)+
	geom_point(data=plot_data_3 %>% filter(Species%in%c(sp_sig_mv1,sp_sig_coi)),
						 aes(x=Treatment,y=mean_eDNA_ch,colour = TrophicLevel),alpha=0.5)+
	geom_line(data=plot_data_3 %>% filter(Species%in%c(sp_sig_mv1,sp_sig_coi)),
						aes(x=Treatment,y=mean_eDNA_ch,colour = TrophicLevel,group = Species),alpha=0.6)+
	facet_wrap(~TrophicLevel,ncol=3)+
	scale_color_manual(values=trophic_color_comb)+
	labs(y='eDNA index change (mean of all pairwise comparisons)',x='Pairwise comparison (treatment)')+
	geom_abline(intercept=0,slope=0,lty=2)+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
				legend.position = 'none',
				axis.title = element_text(size = 19),
				axis.text = element_text(size = 15),
				strip.text = element_text(size=16))

p3_l <-
	data.frame(Treatment = NA, mean_eDNA_ch = NA, 
						 label = c("+ Day / − Night", "+ Light / − No Light")) %>% 
	ggplot()+
	geom_hline(aes(yintercept=0.1, linetype='within ± 0.1 eDNA\nindex change\n(insignificant change)'), colour='grey50') +
	geom_hline(aes(yintercept=-0.1, linetype='within ± 0.1 eDNA\nindex change\n(insignificant change)'), colour='grey50') +
	geom_point(aes(x=Treatment, y=mean_eDNA_ch, shape=label), alpha=0)+
	scale_shape_manual(
		name="eDNA index change",
		values=c(NA, NA),
		labels=c("+ Day  − Night", "+ Light  − No Light")) +
	scale_linetype_manual(
		name=NULL,
		values=c('within ± 0.1 eDNA\nindex change\n(insignificant change)'='solid')) +
	guides(
		shape    = guide_legend(order=1, override.aes=list(alpha=0.7)),
		linetype = guide_legend(order=2)) +
	theme(legend.position="right",
				legend.key=element_blank(),
				legend.text  = element_text(size=14),
				legend.title = element_text(size=15))

p3_legend <- cowplot::get_legend(p3_l)

p3_f <- cowplot::plot_grid(p3,p3_legend,rel_widths = c(8,2))
ggsave(here('plots','Figure 2_cmbxxxx.jpg'),p3_f,width=12,heigh=12)



eDNA_idx_change< -0.1|eDNA_idx_change> 0.1


# ggsave(here('plots','Figure 2_cmb.jpg'),p3,width=11,heigh=12)




x <-
plot_data_coi %>% 
	# plot_data_2_mv1 %>% 
	pivot_wider(names_from = Treatment,
							values_from = mean_eDNA_ch) %>% 
	# filter(TrophicLevel=='Benthic invertebrates') %>% 
	# filter(!(Species%in%sp_sig_coi)) %>% 
	arrange(TrophicLevel,Species) %>% 
	as.data.frame() %>% 
	mutate(`Day-Night`=round(`Day-Night`,3)) %>% 
	mutate(`Light-No_Light(day)`=round(`Light-No_Light(day)`,3)) %>% 
	mutate(`Light-No_Light(night)`=round(`Light-No_Light(night)`,3))

# write.csv(x,here('coi.csv'),row.names = F)
# write.csv(x,here('mv1.csv'),row.names = F)

plot_data_coi %>% 
	pivot_wider(names_from = Treatment,
							values_from = mean_eDNA_ch) %>% 
	filter(TrophicLevel=='Benthic invertebrates') %>%
	arrange(TrophicLevel,Species) %>% 
	print(n=100)



# Sp accumulation curves and rarefaction curves -----------------------------------------------

# COI
sac <- specaccum(coi_bin %>% t(), method = "random")

sp_acc <- data.frame(Sites = sac$sites, Richness = sac$richness, sd = sac$sd)

pp1 <- sp_acc %>%
	mutate(lo=Richness-(2*sd)) %>% 
	mutate(up=Richness+(2*sd)) %>% 
	ggplot()+
	geom_errorbar(aes(x=Sites,ymin = lo,ymax = up))+
	geom_point(aes(x=Sites,y=Richness))+
	labs(x='Sampling effort (no of samples)',y='ASV richness')+
	theme_bw()+
	theme(axis.text = element_text(size=15),
				axis.title = element_text(size=19),
				strip.text = element_text(size=16))



rare_list <- rarecurve(
	coi_raw,
	step = 100,
	# sample = min(rowSums(coi_bin %>% t())),
	tidy = TRUE
)

pp2 <- rare_list %>% 
	ggplot(aes(x = Sample, y = Species, group = Site)) +
	geom_line(alpha = 0.6) +
	labs(
		x = "Sequencing depth (reads)",
		y = "ASV richness") +
	theme_bw()+
	theme(axis.text = element_text(size=15),
				axis.title = element_text(size=19),
				strip.text = element_text(size=16))

p1 <- cowplot::plot_grid(pp1,pp2)
ggsave(here('plots','coi_acc_rare_curves.jpg'),p1,width = 12,height = 6)


# MV1
sac <- specaccum(mv1_bin %>% t(), method = "random")

sp_acc <- data.frame(Sites = sac$sites, Richness = sac$richness, sd = sac$sd)

pp1 <- sp_acc %>%
	mutate(lo=Richness-(2*sd)) %>% 
	mutate(up=Richness+(2*sd)) %>% 
	ggplot()+
	geom_errorbar(aes(x=Sites,ymin = lo,ymax = up))+
	geom_point(aes(x=Sites,y=Richness))+
	labs(x='Sampling effort (no of samples)',y='ASV richness')+
	theme_bw()+
	theme(axis.text = element_text(size=15),
				axis.title = element_text(size=19),
				strip.text = element_text(size=16))



rare_list <- rarecurve(
	mv1_raw,
	step = 100,
	# sample = min(rowSums(coi_bin %>% t())),
	tidy = TRUE
)

pp2 <- rare_list %>% 
	ggplot(aes(x = Sample, y = Species, group = Site)) +
	geom_line(alpha = 0.6) +
	labs(
		x = "Sequencing depth (reads)",
		y = "ASV richness") +
	theme_bw()+
	theme(axis.text = element_text(size=15),
				axis.title = element_text(size=19),
				strip.text = element_text(size=16))

p2 <- cowplot::plot_grid(pp1,pp2)
ggsave(here('plots','mv1_acc_rare_curves.jpg'),p2,width = 12,height = 6)
