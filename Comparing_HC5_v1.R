
### Loading necessary packages ###
package.list = c("openxlsx","ssdtools","tidyverse", "EnvStats", "ggplot2", "gridExtra", "reshape2") 
tmp.install = which(lapply(package.list, require, character.only = TRUE) == FALSE)
if(length(tmp.install)>0) install.packages(package.list[tmp.install])
lapply(package.list, require, character.only = TRUE)
# If all responses show "TRUE", all packages are successfully loaded.


### Data to pick up which chemicals to be analyzed ###
dat01 <- read.xlsx("TableS1_v1.xlsx", sheet = "Table_S1")
StudyChemicals <- dat01$original.CAS
length(StudyChemicals) # number of chemicals analyzed


### Read the data used for the analysis ###
# Please get the data from the EnviroTox database (https://envirotoxdatabase.org/)
# Click the "Search" button without entering any search terms, and then click "Download as Excel File"
# Please change the file name if it is different from the one provided below.
EnviroTox_test <- read.xlsx("envirotox_20211110145304.xlsx", sheet="test")
EnviroTox_chem <- read.xlsx("envirotox_20211110145304.xlsx", sheet="substance")
EnviroTox_taxo <- read.xlsx("envirotox_20211110145304.xlsx", sheet="taxonomy")
head(EnviroTox_test)


### Select and process the data ###
EnviroTox_test_selected1 <- EnviroTox_test %>%
  filter(Test.statistic == "EC50" & Test.type == "A" | Test.statistic == "LC50" & Test.type == "A") %>% # Select only acute data
  filter(Effect.is.5X.above.water.solubility == "0") %>% # Exclude irrelevant data in terms of water solubility
  filter(original.CAS %in% StudyChemicals) %>% # Select study chemicals
  mutate(original.CAS = EnviroTox_chem[match(.$CAS, EnviroTox_chem$CAS), "original.CAS"]) %>%
  mutate_at(vars(Effect.value), as.numeric) %>% # Make "Effect.value" numeric values
  mutate(Effect.value = replace(.$Effect.value, !is.na(.$Effect.value), .$Effect.value * 10^3)) %>%  # Transform unit (mg/L to µg/L)
  mutate(Unit = replace(Unit, Unit == "mg/L", "µg/L")) %>%
  mutate(Substance = EnviroTox_chem[match(.$original.CAS, EnviroTox_chem$original.CAS), "Chemical.name"]) %>%
  separate(Substance, into = c("Short_name"), sep = ";", extra = "drop")
head(EnviroTox_test_selected1)


### Remove unnecessary columns ###
EnviroTox_test_selected2 <- EnviroTox_test_selected1 %>%
  select(-Source, -version, -Reported.chemical.name, -Effect.is.5X.above.water.solubility)

### Calculate geometric mean and select chemicals analyzed based on the number of species ###
EnviroTox_test_selected3 <- aggregate(EnviroTox_test_selected2$Effect.value,
                                      by = list(original.CAS = EnviroTox_test_selected2$original.CAS,
                                                Test.type = EnviroTox_test_selected2$Test.type,
                                                Latin.name = EnviroTox_test_selected2$Latin.name),
                                      function(x) geoMean(x)) %>%
  rename(Effect.value = x) %>%
  mutate(Trophic.Level = EnviroTox_taxo[match(.$Latin.name, EnviroTox_taxo$Latin.name), "Trophic.Level"]) %>%
  mutate(Substance = EnviroTox_chem[match(.$original.CAS, EnviroTox_chem$original.CAS), "Chemical.name"]) %>%
  separate(Substance, into = c("Short_name"), sep = ";", extra = "drop") %>%
  group_by(original.CAS, Test.type) %>%
  filter(n() >= 5)


### Calculate reference HC5 values and bimodality coefficients ###
result <- EnviroTox_test_selected3 %>%
  group_by(original.CAS) %>%
  summarise(reference_HC5 = quantile(Effect.value, probs = 0.05),
            BC = mousetrap::bimodality_coefficient(log10(Effect.value)))


### Gathering information ###
summary_info <- result
# Change the order based on the number of species
order_index <- match(dat01$original.CAS, summary_info$original.CAS)
summary_info <- summary_info[order_index, ]
summary_info$use_group_updated <- dat01$use_group_updated # Use group
summary_info$Short_name <- dat01$Short_name # Chemical name
summary_info$No_species <- dat01$No_species_acute # Number of species


### Select the best model for each chemical based on the complete dataset ###
n.original.CAS <- 35 
colnames(EnviroTox_test_selected3)[4] <- "Conc"
best.dist <- rep(NA, n.original.CAS)

# Model selection
for (i in 1:n.original.CAS) {
  t.data1 <- EnviroTox_test_selected3 %>%
    filter(original.CAS == summary_info$original.CAS[i])
  
  t.fits <- ssd_fit_dists(t.data1, dists = c("lnorm", "llogis", "gamma", "weibull", "burrIII3"))
  best.dist[i] <- ssd_gof(t.fits)$dist[ssd_gof(t.fits)$aic == min(ssd_gof(t.fits)$aic)]
}

# Replace distribution names 
full.names <- best.dist
full.names <- gsub("lnorm", "Log-normal", full.names)
full.names <- gsub("weibull", "Weibull", full.names)
full.names <- gsub("llogis", "Log-logistic", full.names)
full.names <- gsub("gamma", "Gamma", full.names)
full.names <- gsub("burrIII3", "Burr III", full.names)
summary_info$best.dist <- full.names


### Derive HC5 estimates based on the single-distribution approach ###
# Note this calculation will take a long time (e.g., 6 hours) 
# Suggest using results available in "result_list_save_sp15.obj" to proceed (see below)
# Settings
n.sim <- 1000 # Number of simulated datasets
n.toxdata <- 15 # Number of species (= toxicity data) selected
result_list <- vector("list", length = n.original.CAS) # Create an object to store results
colnames(EnviroTox_test_selected3)[4] <- "Conc"
### Derive HC5 estimates based on the single-distribution approach ###
# Note this calculation will take a long time (e.g., 6 hours)
# Settings
n.sim <- 1000 # Number of simulated datasets
n.toxdata <- 15 # Number of species (= toxicity data) selected
result_list <- vector("list", length = n.original.CAS) # Create an object to store results
colnames(EnviroTox_test_selected3)[4] <- "Conc"

# Start the simulation
set.seed(101)
for (i in 1:n.original.CAS) {
  t.data1 <- EnviroTox_test_selected3 %>%
    filter(original.CAS == summary_info$original.CAS[i])
  
  # Empty matrix to store results for i
  res_i <- matrix(NA, nrow = n.sim, ncol = 6)

  for (j in 1:n.sim) {
    # Initial sampling
    t.data2 <- sample_n(t.data1, size = n.toxdata, replace = FALSE)
    
    # Ensure t.data2 has rows with three different values for Trophic.Level
    while (length(unique(t.data2$Trophic.Level)) != 3) {
      t.data2 <- sample_n(t.data1, size = n.toxdata, replace = FALSE)
    }
    
    # Estimate SSDs
    fits <- ssd_fit_dists(t.data2, dists = c("lnorm", "llogis", "gamma", "weibull", "burrIII3"))
    hc.ind <- ssd_hc(fits, average = FALSE, delta = 1000) # Added "delta = 1000" not to drop results from models with larger AICc values
    hc.ave <- ssd_hc(fits, delta = 1000)
    
    # Store results for i
    res_i[j, 1] <- hc.ave$est
    try(res_i[j, 2] <- hc.ind$est[hc.ind$dist == "lnorm"], silent = TRUE)
    try(res_i[j, 3] <- hc.ind$est[hc.ind$dist == "llogis"], silent = TRUE)
    try(res_i[j, 4] <- hc.ind$est[hc.ind$dist == "burrIII3"], silent = TRUE)
    try(res_i[j, 5] <- hc.ind$est[hc.ind$dist == "weibull"], silent = TRUE)
    try(res_i[j, 6] <- hc.ind$est[hc.ind$dist == "gamma"], silent = TRUE)
  }
  
  # Store results in list
  result_list[[i]] <- res_i
  print(i)
}

result_list_save <- result_list
# saveRDS(result_list_save, "result_list_save_sp15.obj")

result_list <- readRDS("result_list_save_sp15.obj")
# result_list <- result_list[order_index]


### Make "Figure S1"-like figures ###
plots3 <- list()

for (i in 1:35) {  
  t.comp1 <- data.frame(result_list[[i]])
  colnames(t.comp1) <- c("Averaging", "Log-normal", "Log-logistic", "Burr III", "Weibull", "Gamma")
  
  medians <- apply(t.comp1, 2, median, na.rm = TRUE)
  
  t.comp1_long <- melt(t.comp1)
  
  column_order <- c("Averaging", "Log-normal", "Log-logistic", "Burr III", "Weibull", "Gamma")

  # Set the upper and lower limits
  y_max <- log10(max(t.comp1_long$value, na.rm = TRUE)) + 2
  y_min <- log10(min(t.comp1_long$value, na.rm = TRUE)) - 0.5

  p <- ggplot(t.comp1_long, aes(x = factor(variable, levels = column_order), y = log10(value), fill = variable)) +
    geom_violin(trim = FALSE, width = 1) +
    geom_point(data = data.frame(variable = names(medians), value = log10(medians)), aes(x = factor(variable, levels = column_order), y = value), color = "black", size = 3, shape = 17) +
    scale_fill_brewer(palette = 'Set2') + 
    labs(title = summary_info$Short_name[i], x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x = element_line(color = "black"),
          panel.border = element_rect(fill = NA, color = "black", size = 0.5),
          plot.title = element_text(size = 12),
          legend.position = "none",
          panel.spacing.y = unit(2, "mm"))

  p <- p + geom_hline(yintercept = log10(summary_info$reference_HC5[i]), color = grey(0.2, 0.5), size = 2) +
      annotate("text", x = 0.5, y = y_max * 0.9, label = summary_info$best.dist[i], vjust = 1, hjust = 0, size = 3, color = 1) +
      annotate("text", x = 2.5, y = y_max * 0.9, label = paste("n =", summary_info$No_species[i]), vjust = 1, hjust = 0, size = 3, color = 1) +
      annotate("text", x = 4.5, y = y_max * 0.9, label = paste("BC =", round(summary_info$BC[i], 2)), vjust = 1, hjust = 0, size = 3, color = 1) +
      ylim(y_min, y_max)

  plots3[[i]] <- p
}

# See plots
combined_plot3_1 <- grid.arrange(grobs = plots3[1:12], ncol = 3, nrow = 4) 
combined_plot3_2 <- grid.arrange(grobs = plots3[13:24], ncol = 3, nrow = 4)
combined_plot3_3 <- grid.arrange(grobs = plots3[25:35], ncol = 3, nrow = 4) 
# Save
# ggsave("Combined_plot_diff_sp15_page1.pdf", plot = combined_plot3_1, width = 10, height = 10)
# ggsave("Combined_plot_diff_sp15_page2.pdf", plot = combined_plot3_2, width = 10, height = 10)
# ggsave("Combined_plot_diff_sp15_page3.pdf", plot = combined_plot3_3, width = 10, height = 10)


### Make "Figure 2" ###
# Calculate the deviations and plot them
# Prepare matrices
t.res.mean <- matrix(-9999, ncol = 6, nrow = 35)
t.res.upper <- matrix(-9999, ncol = 6, nrow = 35)
t.res.lower <- matrix(-9999, ncol = 6, nrow = 35)

# Calculate the deviations
for (i in 1:35) {
  result_list_diff <- log10(result_list[[i]]) - log10(summary_info$reference_HC5[i]) # Reference HC5 value - HC5 of each distribution
  t.res.mean[i, ] <- apply(result_list_diff, 2, median, na.rm = TRUE) # Calculate the median deviation for each distribution
  t.res.upper[i, ] <- apply(result_list_diff, 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE)) # Calculate the 97.5 percentile deviation for each distribution
  t.res.lower[i, ] <- apply(result_list_diff, 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE)) # Calculate the 2.5 percentile deviation for each distribution
}


# Creating a flag for bimodality
summary_info$Cate_BC <- summary_info$BC
summary_info <- summary_info %>%
  dplyr::mutate(Cate_BC = ifelse(BC > 0.555, 1, 2))

# Make Figure 2
png("Comparison_sp15.png",  width = 1000, height = 750)
par(mfrow = c(1, 1), cex = 3, mar = c(2, 2.2, 0.5, 0.5), mgp = c(2, 0.5, 0))

plot(1:6, rep(-9999, 6), cex = 0, xlab = "", ylab = "", ylim = c(-8.1, 4), xlim = c(1, 6.5), xaxt = 'n')
lines(c(0, 7), c(0, 0), lty = 2)

# Median
mean.median <- apply(t.res.mean, 2, mean)
sd.median <- apply(t.res.mean, 2, sd)
arrows(1:6, mean.median - sd.median, 1:6, mean.median + sd.median, length = 0, lwd = 2)
points(1:6, mean.median, pch = 21, bg = grey(0.2, 0.8), col = 0, cex = 2)

for (i in 1:6) {
  points(jitter(rep(i + 0.15, 35), factor = 0.2), t.res.mean[, i], pch = c(4, 21)[summary_info$Cate_BC], lwd = c(2, 1)[summary_info$Cate_BC],
         bg = c(grey(1, 0.5), grey(0.2, 0.5))[summary_info$Cate_BC], col = c(1, 0)[summary_info$Cate_BC], cex = 0.8)
}

# Upper
mean.upper <- apply(t.res.upper, 2, mean)
sd.upper <- apply(t.res.upper, 2, sd)
arrows(1:6 + 0.3, mean.upper - sd.upper, 1:6 + 0.3, mean.upper + sd.upper, length = 0, lwd = 2, col = rgb(1, 0.08, 0.58, 0.9))
points(1:6 + 0.3, mean.upper, pch = 22, bg = rgb(1, 0.08, 0.58, 0.9), col = 0, cex = 2)

for (i in 1:6) {
  points(jitter(rep(i + 0.3 + 0.15, 35), factor = 0.2), t.res.upper[, i], pch = c(4, 22)[summary_info$Cate_BC], lwd = c(2, 1)[summary_info$Cate_BC],
         bg = c(grey(1, 0.5), rgb(1, 0.08, 0.58, 0.3))[summary_info$Cate_BC], col = c(rgb(1, 0.08, 0.58, 1), rgb(1, 0.08, 0.58, 0.3))[summary_info$Cate_BC], cex = 0.8)
}

# Lower
mean.lower <- apply(t.res.lower, 2, mean)
sd.lower <- apply(t.res.lower, 2, sd)
points(1:6 + 0.3, mean.lower, pch = 23, col = 0, cex = 2, bg = rgb(30/255, 144/255, 255/255, alpha = 0.9))
arrows(1:6 + 0.3, mean.lower - sd.lower, 1:6 + 0.3, mean.lower + sd.lower, length = 0, lwd = 2, col = rgb(30/255, 144/255, 255/255, alpha = 0.9))

for (i in 1:6) {
  points(jitter(rep(i + 0.3 + 0.15, 35), factor = 0.2), t.res.lower[, i], pch = c(4, 23)[summary_info$Cate_BC], lwd = c(2, 1)[summary_info$Cate_BC],
         bg = c(grey(1, 0.5), rgb(30/255, 144/255, 255/255, alpha = 0.3))[summary_info$Cate_BC], col = c(rgb(30/255, 144/255, 255/255, alpha = 1), 0)[summary_info$Cate_BC], cex = 0.8)
}

axis(1, 1:6 + 0.15, label = FALSE)
dev.off()














