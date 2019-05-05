###############################################################################
## This script aims visualize the varying MSPE values
##
## Created: April 15, 2019
## Updated 1:
###############################################################################

#outline:
## load mspe estimates
## plot graph

load("./NIMBLE/data-prepped/MSPE_simple.Rdata")


MSPE_results_summary <- MSPE_results %>%
  group_by(penalty) %>%
  summarise(MSPE = mean(MSPE))


compare_pen <- ggplot(data = MSPE_results_summary,
                      aes(x = penalty, y = MSPE)) +
  geom_point(aes(x = penalty))

plot_pen <- compare_pen +
  geom_line() +
  # geom_vline( xintercept=-6, col = 'red') +
  # labs(x = expression("Penalty (" ~tau~ "= e^value)")) +
  # ylim(0.017875, 0.017915) +
  # xlim(-15, 0) +
  # theme_bw(base_size = 13, base_family = "Helvetica") +
  theme(text = element_text(size=23),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

plot_pen
