#temporary analysis script -- COI day vs night comparison

tmp <- data.frame(COI_p, COI_idx)
plots <- list()

pdf("figures/myplots_species.pdf")  
for (i in 11:ncol(tmp)){
plots[[i-10]] <-  tmp %>% 
    ggplot(aes(x = Day.Night, y = !!sym(names(tmp)[i]))) +
    geom_boxplot() +
    ggtitle(names(tmp)[i]) +
    ylab(paste0(names(tmp)[i], " Index"))
}
print(plots)
dev.off()

data.frame(
tot = rowSums(COI_04),
COI_p) %>% 
  group_by(Day.Night) %>% 
  summarize(sum(tot))


tmp %>% 
  pivot_longer(cols = 11:ncol(.)) %>% 
  group_by(name, Day.Night) %>% 
  summarise(m = mean(value)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Day.Night, values_from = m) %>% 
  mutate(diff = day - night) %>% 
  arrange(desc(diff)) %>% 
  print(n = 100)

