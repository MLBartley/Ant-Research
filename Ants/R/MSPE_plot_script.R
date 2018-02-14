
load("./Ants/output/MSPE_pen2_table.Rda")
load("./Ants/output/MSPE_pen3_table.Rda")


#if models ran at penalty tuning parameter multiple times, average together


MSPE_results2_summ <- MSPE_results %>% 
  group_by(penalty) %>% 
  summarise(MSPE = mean(MSPE))


MSPE_results3_summ <- MSPE_results_THREE %>% 
  group_by(penalty) %>% 
  summarise(MSPE = mean(MSPE))

#plot separately 

compare_pen2 <- ggplot(data = MSPE_results2_summ, aes(x = penalty, y = MSPE)) +
  geom_point(aes(x = penalty)) 

plot_pen2 <- compare_pen2 + 
  geom_vline( xintercept=-23, col = 'red') + 
  xlab("Penalty (e^value)") + 
  xlim(-30, 30) +
  # theme_bw(base_size = 13, base_family = "Helvetica") + 
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        axis.line = element_line(colour = "black")) 

print(plot_pen2)

#want to combine data into one - need column for number of states

MSPE_results <- MSPE_results %>% 
  
  
  #remove any low outliers from unsuccessful runs
  
  
  #plot with color by model + add vert lines 
  
  
  
  #Woo done, now save, send to Ephraim, and add to manuscript
  