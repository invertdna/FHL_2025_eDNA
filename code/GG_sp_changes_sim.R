library(dplyr)
library(ggpubr)
library(tibble)
library(tidyr)

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


# sp1 loves to be in the day but not at night
# sp2 loves to be in the night but not at day
# sp3 loves to be in the night but more so when there's light
# sp4 is hates the light (the day and the bright light)

meta <- data.frame(tube=as.character(c(1:8)),
									 tr=rep(rep(c('day','night'),each=2),2),
									 date=rep(c(1:2),each=4),
									 l=rep(c('light','no_light'),4))
comm <- data.frame(sp_1=c(3,4,1,0,6,4,0,1),
									 sp_2=c(1,2,8,6,1,0,6,8),
									 # sp_3=c(1,2,9,6,1,0,10,3),
									 sp_4=c(0,1,1,6,1,0,1,7)) %>% 
	t() %>% as.data.frame() %>% setNames(as.character(1:8))


pair_sp_diff <- pairwise_row_diff(t(comm)) %>% bind_cols(.) %>% 
	rownames_to_column('species') %>% 
	pivot_longer(cols = -species) %>% 
	separate(name, into = c("Tube_1", "Tube_2"), sep = "_vs_") 

ch_data <-pair_sp_diff %>%
	left_join(.,meta %>% rename(tod_1='tr',
															date_1='date',
															l_1='l'),by=c('Tube_1'='tube')) %>% 
	left_join(.,meta %>% rename(tod_2='tr',
															date_2='date',
															l_2='l'),by=c('Tube_2'='tube')) %>% 
	mutate(treatment=if_else(tod_1!=tod_2,paste0(tod_1,'-',tod_2),NA)) %>% 
	mutate(treatment=if_else(tod_1==tod_2&l_1!=l_2,paste0(l_1,'-',l_2,'(',tod_1,')'),treatment)) %>% 
	mutate(delta_date=abs(date_1-date_2)) %>% 
	filter(delta_date<2)  %>% 
	mutate(value=if_else(treatment=='night-day',value*-1,value)) %>% 
	mutate(treatment=if_else(treatment=='night-day','day-night',treatment)) %>% 
	mutate(value=if_else(treatment=='no_light-light(day)',value*-1,value)) %>% 
	mutate(treatment=if_else(treatment=='no_light-light(day)','light-no_light(day)',treatment)) %>% 
	mutate(value=if_else(treatment=='no_light-light(night)',value*-1,value)) %>% 
	mutate(treatment=if_else(treatment=='no_light-light(night)','light-no_light(night)',treatment)) %>% 
	filter(!is.na(treatment)) %>% 
	filter(!(tod_1!=tod_2&l_1!=l_2)) %>% 
	filter(!(delta_date==1&tod_1==tod_2)) 

meta %>% cbind(comm %>% t() %>% as.data.frame() %>% select('sp_1'))
ch_data %>% group_by(species,treatment) %>% summarise(ch=mean(value))
# ch_data %>% group_by(species) %>% summarise(ch=mean(value))

ch_data %>% 
	ggline(x='treatment',y='value',color='species',add=c('mean_se'))+
	geom_abline(intercept=0,slope=0,lty=2)+
	annotate("text",x = Inf,y = 0,label = "no change",hjust = 1.1,vjust = -0.5)
