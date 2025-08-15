```{r StackedBars}
#COI

COI_night_minTubevec <- c(1,15,37,51,67,82)
COI_night_maxTubevec <- c(5,20,42,56,70,84)

plotCOI <- COI_common_prop %>% 
	rownames_to_column("Tube") %>% 
	pivot_longer(-Tube) %>% 
	left_join(COI_p %>% dplyr::select(Tube, datetime) %>% mutate(Tube = as.character(Tube))) %>% 
	ggplot(aes(x = as.numeric(Tube), y = value, fill = name)) +
	geom_col() +
	#theme(legend.position = "none") +
	
	
	# Add alternating day/night shading rectangles first
	geom_rect(data = data.frame(
		xmin = 0.5,
		xmax = 90.5,
		ymin = -Inf, ymax = -0.05
	), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
	fill = "lightblue", alpha = 0.3, inherit.aes = FALSE) +
	
	geom_rect(data = data.frame(
		xmin = COI_night_minTubevec+0.5,
		xmax = COI_night_maxTubevec+0.5,
		ymin = -Inf, ymax = -0.05
	), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
	fill = "darkblue", alpha = 0.3, inherit.aes = FALSE) +
	
	geom_col() + xlab("") + ylab("Proportion") +
	theme_minimal() +
	theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
		legend.position = "top",
		legend.text = element_text(size = 7,
															 margin = margin(t = 1, r = 10, b = 1, l = 10, unit = "pt")),
		legend.title = element_blank()
	) +
	
	scale_x_continuous(labels = unique(COI_p$Date), 
										 breaks = COI_p$Tube[match(unique(COI_p$Date), COI_p$Date)]) +
	# Expand the plot area to show the shading below
	coord_cartesian(ylim = c(-0.1, max(COI_common_prop, na.rm = TRUE) * 1.05), 
									expand = FALSE)

ggsave(filename = "figures/timeseries_COI.pdf",
			 plot = plotCOI, 
			 width = 12,
			 height = 7)


##MV1
MV1_night_minTubevec <- c(1,15,37,51,65,80,92)
MV1_night_maxTubevec <- c(6,20,42,56,70,84,97)

plotMV1 <- MV1_common_prop %>% 
	rownames_to_column("Tube") %>% 
	pivot_longer(-Tube) %>% 
	left_join(MV1_common_p %>% dplyr::select(Tube, datetime) %>% mutate(Tube = as.character(Tube))) %>% 
	ggplot(aes(x = as.numeric(Tube), y = value, fill = name)) +
	geom_col() +
	#theme(legend.position = "none") +
	
	
	# Add alternating day/night shading rectangles first
	geom_rect(data = data.frame(
		xmin = 0.5,
		xmax = 74,
		ymin = -Inf, ymax = -0.05
	), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
	fill = "lightblue", alpha = 0.3, inherit.aes = FALSE) +
	
	geom_rect(data = data.frame(
		xmin = MV1_night_minTubevec+0.5,
		xmax = MV1_night_maxTubevec+0.5,
		ymin = -Inf, ymax = -0.05
	), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
	fill = "darkblue", alpha = 0.3, inherit.aes = FALSE) +
	
	geom_col() + xlab("") + ylab("Proportion") +
	theme_minimal() +
	theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
		legend.position = "top",
		legend.text = element_text(size = 7,
															 margin = margin(t = 1, r = 10, b = 1, l = 10, unit = "pt")),
		legend.title = element_blank()
	) +
	
	scale_x_continuous(labels = unique(MV1_p$Date), 
										 breaks = MV1_p$Tube[match(unique(MV1_p$Date), MV1_p$Date)]) +
	# Expand the plot area to show the shading below
	coord_cartesian(ylim = c(-0.1, max(MV1_common_prop, na.rm = TRUE) * 1.05), 
									expand = FALSE)

ggsave(filename = "figures/timeseries_MV1.pdf",
			 plot = plotMV1, 
			 width = 12,
			 height = 7)



##MFU