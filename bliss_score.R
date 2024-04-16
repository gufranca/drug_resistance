

setwd("path/")

load("modules_kuramochi.rda")

kura_ola_cb_data = read.table("kura_all_ola_cb_norm_dmso.tsv", header = T,
                              row.names = 1)

palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]  # Extract color info
palette3_all <- unlist(mapply(brewer.pal,                     # Create vector with all colors
                              palette3_info$maxcolors,
                              rownames(palette3_info)))
set.seed(333223)                                             # Set random seed
palette3 <- sample(palette3_all, 7)                    # Sample colors
palette3                                     


bliss_excess <- function(combo, A, B){
  # combo: fraction of inhibition of the combination AB
  # A: fraction of inhibition of drug A at dose x
  # B: fraction of inhibition of drug B at dose y
  
  bliss = combo - (A + B - (A*B))
  return(bliss)
}


bliss_average_all <- function(kura_ola_cb_data) {
  # Calculates bliss score for all the doses and replicates, then take the
  # average between replicates to generate a matrix for heatmap plot
  
  # to interate the rows in pairs
  indices = seq(1,nrow(kura_ola_cb_data), 2)
  bliss_scores_all = c()
  samples = c()
  for(i in indices) {
    ind = i + 1 # to subset pairs of rows that contain drug A and Combo
    subset_data = kura_ola_cb_data[i:ind, ]
    # transform back to fraction and calculate the inhibition (1 - viability)
    # for values that are negative (proliferation, set inhibition = 0)
    subset_data = 1 - (subset_data / 100)
    subset_data[subset_data < 0] <- 0
    sample_name = rownames(subset_data[2,])
    samples = c(samples, sample_name)
    # drug B only inhibition
    B = subset_data[2,1]
    B_vector = rep(B, 9)
    # combo inhibition
    combo = as.numeric(subset_data[2, 2:ncol(subset_data)])
    #print(combo)
    # drug A inhibition
    A_vector = as.numeric(subset_data[1, 2:ncol(subset_data)])
    # calculate scores for all doses
    bliss_scores =  mapply(bliss_excess, combo, A_vector, B_vector)
    bliss_scores_all = rbind(bliss_scores_all, bliss_scores)
    rownames(bliss_scores_all) <- samples  
  }
  
  
  # calculate means across the three replicates
  indices_rep = seq(1,nrow(bliss_scores_all), 3)
  bliss_scores_mean = c()
  samples = c()
  for(i in indices_rep) {
    ind_rep = i + 2
    subset_data_rep = bliss_scores_all[i:ind_rep, ]
    sample_means = colMeans(subset_data_rep)
    bliss_scores_mean = rbind(bliss_scores_mean, sample_means)
    sample_name = rownames(subset_data_rep)[1] # just take the last one
    sample_name = paste(strsplit(sample_name, split = "_")[[1]][1:2],
                        collapse = "_")
    samples = c(samples, sample_name)
    rownames(bliss_scores_mean) <- samples
    colnames(bliss_scores_mean) <- c("1", "2.5", "5", "10","20","40", "80",
                                     "160", "320")
  }
  
  # reorder the rows for heatmap
  bliss_scores_mean = bliss_scores_mean[c(1,5,2,4,3), ]
  return(bliss_scores_mean)
}  
  
bliss_scores_all = bliss_average_all(kura_ola_cb_data)
# transform with *100
bliss_scores_all = bliss_scores_all*100

library(ComplexHeatmap)
library(RColorBrewer)

Heatmap(bliss_scores_all, cluster_rows = F, cluster_columns = F,
        col = colorRampPalette(rev(brewer.pal(10, "RdBu")))(20), border = T,
        rect_gp = gpar(col = "black", lwd = 1))


