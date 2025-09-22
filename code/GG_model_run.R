# Libraries -----------------------------------------------------------------------------------
library(vegan)
library(rstan);options(mc.cores = parallel::detectCores())
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

extract_param <- function (model, par) 
{
	fit <- (methods::selectMethod("summary", signature = "stanfit"))(object = model, 
																																	 par = par)
	fit <- fit$summary
	return(fit %>% unlist() %>% as.data.frame %>% round(., 9))
}


logit <- function(p) {
	if (any(p <= 0 | p >= 1)) stop("logit is only defined for 0 < p < 1")
	log(p / (1 - p))
}

inv_logit <- function(x) {
	1 / (1 + exp(-x))
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
		filter(!(time_cat==1&tre_cat==1)) %>% 
		mutate(beta_idx=ifelse(time_1=='night',2,1))
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

for (k in 1:length(names_of_iterations)) {
design_mat <- design_mat_list[[k]]

stan_data <- list(
	N = nrow(design_mat),         
	D = design_mat$D,
	delta_date = design_mat$delta_date,
	time_cat = design_mat$time_cat,
	tre_cat = design_mat$tre_cat,
	rep_cat = design_mat$rep_cat,
	beta_idx = design_mat$beta_idx,
	J = max(design_mat$beta_idx)
)


fit <- stan(
	file = "/Users/gledguri/Local/GitHub/FHL_2025_eDNA/code/model.stan",
	data = stan_data,
	chains = 4,
	iter = 2000,
	warmup = 1000)


output <- extract_param(model = fit,par = c('gamma','alpha','beta','theta','intercept','phi')) %>% 
	rownames_to_column('param') %>% 
	mutate(param=if_else(param=='beta[1]','Light_day',param)) %>% 
	mutate(param=if_else(param=='beta[2]','Light_night',param)) %>% 
	mutate(param=if_else(param=='gamma','Date',param)) %>% 
	mutate(param=if_else(param=='alpha','Time (Day_v_Night)',param)) %>% 
	mutate(param=if_else(param=='theta','Bio_rep',param)) %>% 
	mutate(param=if_else(param=='intercept','Intercept',param)) %>% 
	filter(!param=='phi') %>% 
	mutate(param = factor(param, 
												levels = c("Intercept", "Date", "Time (Day_v_Night)", "Light_day", "Light_night", "Bio_rep")))

fig_1 <- output %>% 
	ggplot()+
	geom_point(aes(x=param,y=mean))+
	geom_errorbar(aes(x=param,ymin=`2.5%`,ymax=`97.5%`),width=0.2) +
	theme_bw()+
	labs(y='Dissimilarity effect (logit-scale)',
			 x= 'Parameters')+
	ggtitle(names_of_iterations[k])+
	theme(axis.text = element_text(size=14),
				axis.title = element_text(size=16),
				axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1))

ggsave(here('plots',paste0('Fig_1_',names_of_iterations[k],'.jpg')),fig_1,width = 10,height = 8)

fig_2 <- design_mat %>% 
	group_by(time_cat,tre_cat,delta_date,rep_cat) %>% 
	mutate(D_mean=mean(D)) %>% 
	mutate(I=output$mean[output$param=='Intercept']) %>%
	mutate(gamma=output$mean[output$param=='Date']) %>% 
	mutate(alpha=output$mean[output$param=='Time (Day_v_Night)']) %>% 
	mutate(beta_1=output$mean[output$param=='Light_day']) %>% 
	mutate(beta_2=output$mean[output$param=='Light_night']) %>% 
	mutate(theta=output$mean[output$param=='Bio_rep']) %>% 
	mutate(beta=if_else(beta_idx==2,beta_2,beta_1)) %>% 
	mutate(D_est_logit=
				 	I+
				 	gamma*delta_date+
				 	alpha * time_cat+
				 	beta * tre_cat+
				 	theta * rep_cat) %>% 
	mutate(D_est=inv_logit(D_est_logit)) %>% 
	ggplot()+
	geom_point(aes(x=D,y=D_est),alpha=0.3,pch=4)+
	geom_point(aes(x=D_mean,y=D_est),size=4,color='red')+
	theme_bw()+
	ggtitle(names_of_iterations[[k]])+
	geom_abline(intercept = 0,slope=1,lty=2)+
	labs(x='Estimated Dissimilarity',
			 y='Predicted Dissimilarity')
ggsave(here('plots',paste0('Fig_2_',names_of_iterations[k],'.jpg')),fig_2,width = 10,height = 8)


fig_3 <- design_mat %>% 
	mutate(beta_idx=if_else(tre_cat==0,0,beta_idx)) %>%
	group_by(time_cat,tre_cat,delta_date,rep_cat,beta_idx) %>%
	summarise(D_mean=mean(D)) %>% 
	mutate(I=output$mean[output$param=='Intercept']) %>%
	mutate(gamma=output$mean[output$param=='Date']) %>% 
	mutate(alpha=output$mean[output$param=='Time (Day_v_Night)']) %>% 
	mutate(beta_1=output$mean[output$param=='Light_day']) %>% 
	mutate(beta_2=output$mean[output$param=='Light_night']) %>% 
	mutate(theta=output$mean[output$param=='Bio_rep']) %>% 
	mutate(beta=if_else(beta_idx==2,beta_2,beta_1)) %>% 
	mutate(D_est_logit=
				 	I+
				 	gamma*delta_date+
				 	alpha * time_cat+
				 	beta * tre_cat) %>% 
				 	# theta * rep_cat) %>%
	mutate(D_est=inv_logit(D_est_logit)) %>% 
	mutate(group=if_else(time_cat==0,'Same_time_of_day',NA)) %>% 
	mutate(group=if_else(time_cat==1,'Diff_time_of_day',group)) %>% 
	mutate(group=if_else(time_cat==0&tre_cat==1&beta_idx==1,'Day_light_effect',group)) %>% 
	mutate(group=if_else(time_cat==0&tre_cat==1&beta_idx==2,'Night_light_effect',group)) %>% 
	ggplot()+
	geom_line(aes(x=delta_date,y=D_est,colour = factor(group)),size=1,lty=2)+
	geom_point(aes(x=delta_date,y=D_mean,colour = factor(group)),size=3)+
	ggtitle(names_of_iterations[[k]])+
	scale_color_manual(
		values = c("Same_time_of_day" = "#1f78b4",
							 "Diff_time_of_day" = "#e31a1c",
							 "Day_light_effect" = "#33a02c",
							 "Night_light_effect" = "#6a3d9a"),
		name = "Group",
		labels = c("Light effect during day", "Different ToD", "Light effect during night", "Same ToD"))+
	labs(x='Time',
				 y='Community dissimilarity')+
	theme_bw()
ggsave(here('plots',paste0('Fig_3_',names_of_iterations[k],'.jpg')),fig_3,width = 12,height = 8)
}
# c("#1f78b4", "#e31a1c", "#33a02c", "#6a3d9a")
