# Libraries -----------------------------------------------------------------------------------
library(betareg)
library(vegan)
library(ggplot2)
library(dplyr)
library(here)
library(tibble)
library(tidyr)


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

#Create a list for all the names of the design matricies so that it is easier to know which matrix belongs to which in the for-loops downstream
names_of_iterations <- c('coi_eDNAidx',
												 'coi_presence_absence',
												 'coi_prop',
												 # 'coi_relative_abundance_common_sp',
												 'mfu_eDNAidx',
												 'mfu_presence_absence',
												 'mfu_relative_abundance',
												 # 'mfu_relative_abundance_common_sp',
												 'mv1_eDNAidx',
												 'mv1_presence_absence',
												 'mv1_relative_abundance' #'mv1_relative_abundance_common_sp'
)

# This is an example (here is for COI and eDNA index) of a design matrix
design_mat_coi_eDNAidx %>% 
	count(time_cat,tre_cat,light_day,light_night,delta_date,rep_cat) %>% 
	rename(Day_v_night='time_cat') %>% 
	rename(Date_difference='delta_date') %>% 
	rename(Light_v_Nolight='tre_cat') %>% 
	rename(Biological_replciate='rep_cat')

# Create a list for capturing all the model outputs
s <- vector('list')

for (k in 1:length(names_of_iterations)) {
	design_mat <- design_mat_list[[k]]
	
	# model <- betareg(D ~ delta_date + time_cat + tre_cat + beta_idx + rep_cat, data = design_mat, link = "logit")
	model <- betareg(D ~ delta_date + time_cat + light_day + light_night + rep_cat, data = design_mat, link = "logit")
	
	s[[k]] <- summary(model)$coefficients$mean %>% as.data.frame() %>% rownames_to_column('param') %>% 
		mutate(param=if_else(param=='delta_date','Time',param)) %>% 
		mutate(param=if_else(param=='time_cat','Samp_period',param)) %>% 
		mutate(param=if_else(param=='light_day','Light_day',param)) %>% 
		mutate(param=if_else(param=='light_night','Light_night',param)) %>% 
		mutate(param=if_else(param=='rep_cat','Bio_rep',param)) %>% 
		mutate(param=if_else(param=='(Intercept)','Intercept',param)) %>% 
		mutate(param = factor(param,
													levels = c("Intercept", "Time", "Samp_period", "Light_day", "Light_night", "Bio_rep")))
	
}

# Plot fig 1 (make sure you have selected the correct model output s[[1]]=COI on eDNAidx);
# Check names_of_iterations[k] to know which data has been inputed in the model

fig_1_a <-
		s[[1]] %>%
	mutate(marker='COI') %>% 
	bind_rows(	s[[4]] %>% mutate(marker='MV1')) %>% 
	mutate(param=if_else(param=='Time','Time (days)',param)) %>% 
	mutate(param=if_else(param=='Samp_period','Sampling period',param)) %>% 
	mutate(param=if_else(param=='Light_day','Light effect during day',param)) %>% 
	mutate(param=if_else(param=='Light_night','Light effect during night',param)) %>% 
	mutate(param=if_else(param=='Bio_rep','Biological replicates',param)) %>% 
	mutate(param=factor(param,levels=c('Intercept',
																		 'Time (days)',
																		 'Sampling period',
																		 'Light effect during day',
																		 'Light effect during night',
																		 'Biological replicates'
																		 ))) %>% 
		mutate(lo=Estimate-`Std. Error`) %>% 
		mutate(up=Estimate+`Std. Error`) %>% 
		ggplot()+
		geom_point(aes(x=param,y=Estimate))+
		geom_errorbar(aes(x=param,ymin=lo,ymax=up),width=0.2) +
		theme_bw()+
	facet_grid(~marker)+
		labs(y='Dissimilarity effect (logit-scale)',
				 x= 'Parameters')+
		theme(axis.text = element_text(size=15),
					axis.title = element_text(size=19),
					strip.text = element_text(size=16),
					axis.text.x = element_text(angle = 45, 
																		 vjust = 1.0, hjust = 1))



	# ggsave(here('plots',paste0('Fig_1_',names_of_iterations[k],'.jpg')),fig_1,width = 10,height = 8)
	
# Creating the data that will be used for figure 3
# First creating the COI eDNA idx
output <- s[[1]]
plot_dat1 <- design_mat_list[[1]] %>%
		group_by(time_cat,tre_cat,delta_date,rep_cat,light_day,light_night) %>%
		summarise(D_mean=mean(D)) %>%
		mutate(I=output$Estimate[output$param=='Intercept']) %>%
		mutate(gamma=output$Estimate[output$param=='Time']) %>%
		mutate(alpha=output$Estimate[output$param=='Samp_period']) %>%
		mutate(beta_1=output$Estimate[output$param=='Light_day']) %>%
		mutate(beta_2=output$Estimate[output$param=='Light_night']) %>%
		mutate(theta=output$Estimate[output$param=='Bio_rep']) %>%
		mutate(beta=if_else(light_day==1,beta_1,beta_2)) %>%
		mutate(beta=if_else(tre_cat==0,0,beta)) %>%
		mutate(D_est_logit=
					 	I+
					 	gamma*delta_date+
					 	alpha * time_cat+
					 	beta * tre_cat) %>%
		# theta * rep_cat) %>%
		mutate(D_est=inv_logit(D_est_logit)) %>%
		mutate(group=if_else(time_cat==0,'Same_time_of_day',NA)) %>%
		mutate(group=if_else(time_cat==1,'Diff_time_of_day',group)) %>%
		mutate(group=if_else(time_cat==0&tre_cat==1&light_day==1,'Day_light_effect',group)) %>%
		mutate(group=if_else(time_cat==0&tre_cat==1&light_night==1,'Night_light_effect',group)) %>% 
	ungroup() %>% 
	select(delta_date,D_est,D_mean,group)

# Then creating the MV1 eDNA idx
output <- s[[4]]
plot_dat2 <- design_mat_list[[4]] %>%
		group_by(time_cat,tre_cat,delta_date,rep_cat,light_day,light_night) %>%
		summarise(D_mean=mean(D)) %>%
		mutate(I=output$Estimate[output$param=='Intercept']) %>%
		mutate(gamma=output$Estimate[output$param=='Time']) %>%
		mutate(alpha=output$Estimate[output$param=='Samp_period']) %>%
		mutate(beta_1=output$Estimate[output$param=='Light_day']) %>%
		mutate(beta_2=output$Estimate[output$param=='Light_night']) %>%
		mutate(theta=output$Estimate[output$param=='Bio_rep']) %>%
		mutate(beta=if_else(light_day==1,beta_1,beta_2)) %>%
		mutate(beta=if_else(tre_cat==0,0,beta)) %>%
		mutate(D_est_logit=
					 	I+
					 	gamma*delta_date+
					 	alpha * time_cat+
					 	beta * tre_cat) %>%
		# theta * rep_cat) %>%
		mutate(D_est=inv_logit(D_est_logit)) %>%
		mutate(group=if_else(time_cat==0,'Same_time_of_day',NA)) %>%
		mutate(group=if_else(time_cat==1,'Diff_time_of_day',group)) %>%
		mutate(group=if_else(time_cat==0&tre_cat==1&light_day==1,'Day_light_effect',group)) %>%
		mutate(group=if_else(time_cat==0&tre_cat==1&light_night==1,'Night_light_effect',group)) %>% 
	ungroup() %>% 
	select(delta_date,D_est,D_mean,group)

# Then joiing the two data together to be able to plot them together
fig_1_b <- 
	plot_dat1 %>% mutate(marker='COI') %>% rbind(.,plot_dat2 %>% mutate(marker='MV1')) %>% 
		ggplot()+
		geom_line(aes(x=delta_date,y=D_est,colour = factor(group)),size=1,lty=2)+
		geom_point(aes(x=delta_date,y=D_mean,colour = factor(group)),size=3)+
		# ggtitle(names_of_iterations[[k]])+
		scale_color_manual(
			values = c("Same_time_of_day" = "#1f78b4",
								 "Diff_time_of_day" = "#e31a1c",
								 "Day_light_effect" = "#33a02c",
								 "Night_light_effect" = "#6a3d9a"),
			name = "Group",
			labels = c("Light effect during day", "Different time of day", "Light effect during night", "Same time of day"))+
		labs(x='Time (days)',
				 y='Community dissimilarity')+
	facet_wrap(~marker)+
		theme_bw()+
	theme(axis.text = element_text(size=15),
				axis.title = element_text(size=19),
				strip.text = element_text(size=16),
				legend.text = element_text(size=14),
				legend.title = element_blank(),
				legend.position = 'bottom')

fig_1 <- cowplot::plot_grid(fig_1_a,fig_1_b,nrow = 2,align = 'v')
ggsave(here('plots',paste0('Fig_2','.jpg')),fig_1,width = 10,height = 12)

# c("#1f78b4", "#e31a1c", "#33a02c", "#6a3d9a")
