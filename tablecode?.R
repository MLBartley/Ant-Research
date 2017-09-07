plot(sim$start, 1:length(sim$start_time), main = "Penalized HMM", xlab = "Seconds", ylab = "Cumulative     
      Interaction Count", 
  xlim = c(0, 1*60*60), cex.main = 4.5, cex.lab = 1.2, cex.axis = 1.2)
states <- results[[4]]$X.est
rr <- rle(states[, 1])
rr$values <- round(rr$values, digits = 0)
embedded.chain <- rr$values
cs <- c(0, cumsum(rr$lengths)) * 1 - 1
cols <- c("#bc535644", "#538bbc44")

for (j in 1:length(embedded.chain)) {
  rect(cs[j], 0, cs[j + 1], length(sim$start_time), col = cols[embedded.chain[j]], 
    density = NA)
}

points(sim$start, 1:length(sim$start))

legend(x = 100, y = 150, legend = c("Low Rate of Trophallaxis", "High Rate of Trophallaxis"), pch  = c(15, 15), col = c("#bc535699", "#538bbc99"))
