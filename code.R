# Install and load required R packages from GitHub and CRAN
install.packages("ReMFPCA")
library(ReMFPCA)      # For multivariate functional PCA
library(expm)         # For matrix exponential functions
library(fda)          # For functional data analysis
library(ggplot2)      # For data visualization
library(scales)       # For formatting scales (e.g., percentages)
library(factoextra)   # For cluster visualization
library(cluster)      # For clustering algorithms
library(aricode)      # For clustering evaluation metrics
library(mclust)       # For model-based clustering
library(gridExtra)    # For arranging multiple plots

# Create a B-spline basis with 120 basis functions over [0,1]
bs2 <- create.bspline.basis(c(0, 1), 120)
mdbs2 <- Basismfd(bs2)  # Convert to a multivariate functional data basis

# Initialize multivariate functional data objects from motion sense data
mfd1 = Mfd(X = motion_sense_data$user_acceleration, mdbs = mdbs2, method = "data")
mfd2 = Mfd(X = motion_sense_data$pitch_attitude, mdbs = mdbs2, method = "data")
mvmfd_obj_ori = Mvmfd(list(mfd1, mfd2))  # Combine into a multivariate object

# Evaluate the functional data at 1000 points over [0,1]
D1 = mvmfd_obj_ori$eval(list(seq(0, 1, length.out = 1000), seq(0, 1, length.out = 1000)))[[1]]
D2 = mvmfd_obj_ori$eval(list(seq(0, 1, length.out = 1000), seq(0, 1, length.out = 1000)))[[2]]

# Compute weights based on covariance to normalize the data
w_1 = sqrt(1 / sum(diag(cov(t(D1))) / 1000))
w_2 = sqrt(1 / sum(diag(cov(t(D2))) / 1000))
D_1_w = D1 * w_1  # Apply weight to user acceleration data
D_2_w = D2 * w_2  # Apply weight to pitch attitude data

# Create new multivariate functional data objects with weighted data
mfd_v1 <- Mfd(X = D_1_w, mdbs = mdbs2)
mfd_v2 <- Mfd(X = D_2_w, mdbs = mdbs2)
mvmfd_obj <- Mvmfd(mfd_v1, mfd_v2)

# Define tuning parameters for regularization
alpha1 = list(a1 = 2^seq(-35, 5, length.out = 10), a2 = 2^seq(-35, 5, length.out = 10))
sparse_tuning = as.integer(seq(1, 40, length.out = 20))

# Perform Regularized Multivariate Functional PCA (ReMFPCA) with different settings
result_nonsparse <- Remfpca(mvmfd_obj = mvmfd_obj, ncomp = 3, smooth_tuning = list(0, 0), 
                            sparse_tuning = 0, sparse_type = "hard", smoothing_type = "coefpen")
result_sparse <- Remfpca(mvmfd_obj = mvmfd_obj, ncomp = 3, smooth_tuning = list(0, 0), 
                         sparse_tuning = c(1:90), sparse_type = "hard", smoothing_type = "coefpen")
result_sparse_smooth <- Remfpca(mvmfd_obj = mvmfd_obj, ncomp = 3, smooth_tuning = alpha1, 
                                sparse_tuning = sparse_tuning, sparse_type = "hard", smoothing_type = "coefpen")
result_smooth <- Remfpca(mvmfd_obj = mvmfd_obj, ncomp = 3, smooth_tuning = alpha1, 
                         sparse_tuning = 0, sparse_type = "hard", smoothing_type = "coefpen")
result_joint = Remfpca(mvmfd_obj = mvmfd_obj, ncomp = 3, smooth_tuning = alpha1, 
                       sparse_tuning = 0, sparse_type = "hard", alpha_orth = TRUE, smoothing_type = "coefpen")

###################### Figure 7.1 ######################
par(mfrow = c(2, 3))
m <- matrix(c(1, 2, 3, 4, 5, 6, 7, 7, 7), nrow = 3, ncol = 3, byrow = TRUE)
layout(mat = m, heights = c(0.4, 0.4, 0.2))
yname = c("user acceleration", "pitch angle change rate ")
variance_percentage = ((result_smooth$values) / sum((result_smooth$values)))
for (p in 1:2) {
  for (j in 1:3) {
    if (j == 1) {
      # Compute inner products for comparison across methods
      cond0 = ReMFPCA::inprod_mvmfd(mvmfd_obj[1, ], result_smooth$pc_mfd[j, ])
      cond1 = ReMFPCA::inprod_mvmfd(mvmfd_obj[1, ], result_sparse_smooth$pc_mfd[j, ])
      cond2 = ReMFPCA::inprod_mvmfd(mvmfd_obj[1, ], result_nonsparse$pc_mfd[j, ])
      cond3 = ReMFPCA::inprod_mvmfd(mvmfd_obj[1, ], result_joint$pc_mfd[j, ])
      if (p == 1) {
        par(mai = c(0.4, 0.6, 0.4, 0))
        ReMFPCA::plot(as.numeric(cond2 / abs(cond2)) * result_nonsparse$pc_mfd[j, p], col = "grey", lwd = 2, 
                      main = paste0("PC ", j, " ( ", percent(variance_percentage[j], accuracy = 0.01), " )"), 
                      ylab = yname[p], cex.lab = 2, cex.main = 2, cex.axis = 1.8, xaxt = "n")
        legend("bottomright",
               legend = c(bquote(paste(bold(alpha)[seq], " = ", .(round(result_sparse_smooth$smooth_tuning[[j]][[p]], 3)))),
                          bquote(paste(bold(alpha)[joint], " = ", .(round(result_joint$smooth_tuning[[p]], 3))))),
               lty = 0, cex = 2, bty = "n")
      } else {
        par(mai = c(0.7, 0.6, 0.1, 0))
        ReMFPCA::plot(as.numeric(cond0 / abs(cond0)) * result_nonsparse$pc_mfd[j, p], col = "grey", lwd = 2, 
                      ylab = yname[p], cex.lab = 2, cex.main = 2, xlab = "Time", cex.axis = 1.8)
        legend("bottomright",
               legend = c(bquote(paste(bold(alpha)[seq], " = ", .(round(result_sparse_smooth$smooth_tuning[[j]][[p]], 3)))),
                          bquote(paste(bold(alpha)[joint], " = ", .(round(result_joint$smooth_tuning[[p]], 3))))),
               lty = 0, cex = 2, bty = "n")
      }
      ReMFPCA::plot(as.numeric(cond1 / abs(cond1)) * result_sparse_smooth$pc_mfd[j, p], lwd = 3, add = T, col = "red")
      ReMFPCA::plot(as.numeric(cond0 / abs(cond0)) * result_smooth$pc_mfd[j, p], lwd = 3, add = T, col = "black", lty = 2)
      ReMFPCA::plot(as.numeric(cond3 / abs(cond3)) * result_joint$pc_mfd[j, p], lwd = 3, add = T, col = "deepskyblue", lty = 6)
    } else {
      if (j == 2) {
        if (p == 1) {
          par(mai = c(0.4, 0.4, 0.4, 0.2))
        } else {
          par(mai = c(0.7, 0.4, 0.1, 0.2))
        }
        cond3 = cond3
      }
      if (j == 3) {
        if (p == 1) {
          par(mai = c(0.4, 0.2, 0.4, 0.4))
        } else {
          par(mai = c(0.7, 0.2, 0.1, 0.4))
        }
        cond3 = -1 * cond3
      }
      cond0 = 1 * cond0
      cond1 = 1 * cond1
      cond2 = 1 * cond2
      if (p == 1) {
        ReMFPCA::plot(as.numeric(cond2 / abs(cond2)) * result_nonsparse$pc_mfd[j, p], col = "grey", lwd = 2, 
                      main = paste0("PC ", j, " ( ", percent(variance_percentage[j], accuracy = 0.01), " )"), 
                      cex.lab = 2, cex.main = 2, cex.axis = 1.8, xaxt = "n")
        legend("bottomright",
               legend = c(bquote(paste(bold(alpha)[seq], " = ", .(round(result_sparse_smooth$smooth_tuning[[j]][[p]], 3)))),
                          bquote(paste(bold(alpha)[joint], " = ", .(round(result_joint$smooth_tuning[[p]], 3))))),
               lty = 0, cex = 2, bty = "n")
      } else {
        ReMFPCA::plot(as.numeric(cond2 / abs(cond2)) * result_nonsparse$pc_mfd[j, p], col = "grey", lwd = 2, 
                      cex.lab = 2, cex.main = 2, xlab = "Time", cex.axis = 1.8)
        legend("bottomright",
               legend = c(bquote(paste(bold(alpha)[seq], " = ", .(round(result_sparse_smooth$smooth_tuning[[j]][[p]], 3)))),
                          bquote(paste(bold(alpha)[joint], " = ", .(round(result_joint$smooth_tuning[[p]], 3))))),
               lty = 0, cex = 2, bty = "n")
      }
      ReMFPCA::plot(as.numeric(cond1 / abs(cond1)) * result_sparse_smooth$pc_mfd[j, p], lwd = 3, add = T, col = "red")
      ReMFPCA::plot(as.numeric(cond0 / abs(cond0)) * result_smooth$pc_mfd[j, p], lwd = 3, add = T, col = "black", lty = 2)
      ReMFPCA::plot(as.numeric(cond3 / abs(cond3)) * result_joint$pc_mfd[j, p], lwd = 3, add = T, col = "deepskyblue", lty = 6)
    }
  }
}

# Define colors for the legend
colorset = c("grey", "black", "red", "deepskyblue")
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
plot_colors <- colorset
legend(x = "center", inset = 0,
       legend = c("Non Regularized MFPCA", "Sequential Smooth MFPCA", "Smooth and Sparse MFPCA", "Joint Smooth MFPCA"), 
       col = plot_colors, lwd = c(2, 3, 3, 3), cex = 3, horiz = FALSE, lty = c(1, 2, 1, 6), ncol = 2)

###################### Figure 7.2 ######################
par(mfrow = c(1, 2))
plot(result_smooth$lsv[, 1], col = c(rep(1, 24), rep(2, 24), rep(3, 24), rep(4, 24)), pch = 19, 
     main = "Smooth MFPCA", ylab = "PC 1 scores", xlab = "", cex.lab = 1.4)
legend("topleft", legend = c("Jogging", "Standing", "Walking", "Sitting"), 
       col = c(1, 2, 3, 4), pch = 19, bty = "n")
plot(result_sparse_smooth$lsv[, 1], col = c(rep(1, 24), rep(2, 24), rep(3, 24), rep(4, 24)), pch = 19, 
     main = "Smooth and Sparse MFPCA", ylab = "PC 1 scores", xlab = "", cex.lab = 1.4)

# Scale PC scores by eigenvalues
coef_non = sweep(result_nonsparse$lsv, 2, sqrt((96) * result_nonsparse$values), `*`)
coef_smooth = sweep(result_smooth$lsv, 2, sqrt((96) * result_nonsparse$values), `*`)
coef_sparse = sweep(result_sparse_smooth$lsv, 2, sqrt((96) * result_nonsparse$values), `*`)

par(mfrow = c(1, 2))
plot(coef_smooth[, 1], col = c(rep(1, 24), rep(2, 24), rep(3, 24), rep(4, 24)), pch = 19, 
     main = "Smooth MFPCA", ylab = "PC 1 scores", xlab = "", cex.lab = 1.4)
legend("topleft", legend = c("Jogging", "Standing", "Walking", "Sitting"), 
       col = c(1, 2, 3, 4), pch = 19, bty = "n")
plot(coef_sparse[, 1], col = c(rep(1, 24), rep(2, 24), rep(3, 24), rep(4, 24)), pch = 19, 
     main = "Smooth and Sparse MFPCA", ylab = "PC 1 scores", xlab = "", cex.lab = 1.4)

# Silhouette plots to determine optimal number of clusters
S1 = fviz_nbclust(scale(coef_sparse[, 1]), FUNcluster = kmeans, method = "silhouette") + 
  theme_classic() + ggtitle(paste0("PC score ", 1)) + theme_minimal() + labs(x = "k") + 
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

S2 = fviz_nbclust(scale(coef_sparse[, 2]), FUNcluster = kmeans, method = "silhouette") + 
  theme_classic() + ggtitle(paste0("PC score ", 2)) + theme_minimal() + labs(x = "k") + 
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

S12 = fviz_nbclust(scale(coef_sparse[, 1:2]), FUNcluster = kmeans, method = "silhouette") + 
  theme_classic() + ggtitle(paste0("PC score ", 1, " & ", 2)) + theme_minimal() + labs(x = "k") + 
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

grid.arrange(S1, S2, S12, widths = c(5, 5, 5, 0.5), nrow = 1)

# Perform k-means clustering on scaled PC scores
set.seed(120)
cluster_sparse_pc1 = kmeans(x = scale(coef_sparse[, 1]), 3, nstart = 100)
set.seed(120)
cluster_sparse_pc2 = kmeans(x = scale(coef_sparse[, 2]), 4, nstart = 100)
set.seed(120)
cluster_sparse_pc12 = kmeans(x = scale(coef_sparse[, 1:2]), 4, nstart = 100)

cluster_sparse_pc1
cluster_sparse_pc2
cluster_sparse_pc12
# Prepare data frames for cluster visualization
labels <- c("jogging", "stand", "walk", "sit")
fpc1_cluster_sparse <- data.frame(index = rep(1:24, 4), label = cluster_sparse_pc1$cluster, y = rep(labels, each = 24))
fpc2_cluster_spars <- data.frame(index = rep(1:24, 4), label = cluster_sparse_pc2$cluster, y = rep(labels, each = 24))
fpc12_cluster_spars <- data.frame(index = rep(1:24, 4), label = cluster_sparse_pc12$cluster, y = rep(labels, each = 24))

###################### Figure 7.3 ######################
plot_pc1_sparse = ggplot(fpc1_cluster_sparse, aes(x = index, y = y, fill = factor(label))) + 
  geom_tile() +
  scale_fill_manual(values = c('red', 'forestgreen', 'deepskyblue'), labels = c('Cluster 1', 'Cluster 2', 'Cluster 3'), name = NULL) +
  theme(
    legend.text = element_text(size = 25),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 25),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1, 'cm')
  ) + xlab("") + ylab("") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, 24.5)) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5), linetype = "dashed", color = "black", lwd = 0.5)

plot_pc2_sparse = ggplot(fpc2_cluster_spars, aes(x = index, y = y, fill = factor(label))) + 
  geom_tile() +
  scale_fill_manual(values = c('red', "deepskyblue", 'forestgreen', "pink"), labels = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'), name = NULL) +
  theme(
    legend.text = element_text(size = 25),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 25),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1, 'cm')
  ) + xlab("") + ylab("") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, 24.5)) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5), linetype = "dashed", color = "black", lwd = 0.5)

plot_pc12_sparse = ggplot(fpc12_cluster_spars, aes(x = index, y = y, fill = factor(label))) + 
  geom_tile() +
  scale_fill_manual(values = c('red', 'forestgreen', 'deepskyblue', "pink"), labels = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'), name = NULL) +
  theme(
    legend.text = element_text(size = 25),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 25),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1, 'cm')
  ) + xlab("") + ylab("") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, 24.5)) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5), linetype = "dashed", color = "black", lwd = 0.5)

# Display Figure 7.3 plots
plot_pc1_sparse
plot_pc2_sparse
plot_pc12_sparse

###################### Figure C.4 ######################
bs_smooth <- create.bspline.basis(c(0, 1), 20)
mdbs_smooth <- Basismfd(bs_smooth)

ori_data = mvmfd_obj$eval(list(seq(0, 1, length.out = 1000), seq(0, 1, length.out = 1000)))
mfd1_smooth = Mfd(X = D1, mdbs = mdbs_smooth)
mfd2_smooth = Mfd(X = D2, mdbs = mdbs_smooth)
mvmfd_obj_smooth = Mvmfd(list(mfd1_smooth, mfd2_smooth))
ori_data_smooth = mvmfd_obj_smooth$eval(list(seq(0, 1, length.out = 1000), seq(0, 1, length.out = 1000)))

mfd1_smooth_f = Mfd(X = ori_data_smooth[[1]], mdbs = mdbs2)
mfd2_smooth_f = Mfd(X = ori_data_smooth[[2]], mdbs = mdbs2)
mvmfd_obj_smooth_f = Mvmfd(list(mfd1_smooth_f, mfd2_smooth_f))

yname = c("user acceleration", "pitch angle change rate ")
par(mfrow = c(2, 3))
variance_percentage = ((result_sparse_smooth$values) / sum((result_sparse_smooth$values)))
maxv = c(0, 0)
minv = c(100000, 1000000)
min5 = seq(0, 1, length.out = 300)
weights_svd = c(w_1, w_2)
for (p in 1:2) {
  for (j in 1:3) {
    sigma = sqrt(result_sparse_smooth$values[j])
    f1 = ReMFPCA::mean(mvmfd_obj_smooth_f[, p])
    f2 = f1 + (2.5 * sigma) * result_sparse_smooth$pc_mfd[j, p] * (1 / weights_svd[p])
    f3 = f1 - (2.5 * sigma) * result_sparse_smooth$pc_mfd[j, p] * (1 / weights_svd[p])
    maxv[p] = max(maxv[p], max(f2$eval(min5), f3$eval(min5)))
    minv[p] = min(minv[p], min(f2$eval(min5), f3$eval(min5)))
  }
}
for (p in 1:2) {
  x_grids <- seq(range(min5)[1], range(min5)[2], length.out = 100)
  for (j in 1:3) {
    sigma = sqrt(result_sparse_smooth$values[j])
    f1 = ReMFPCA::mean(mvmfd_obj_smooth_f[, p])
    if (j == 1) {
      f2 = f1 - (2 * sigma) * result_sparse_smooth$pc_mfd[j, p] * (1 / weights_svd[p])
      f3 = f1 + (2 * sigma) * result_sparse_smooth$pc_mfd[j, p] * (1 / weights_svd[p])
    } else if (j == 2) {
      f2 = f1 + (2 * sigma) * result_sparse_smooth$pc_mfd[j, p] * (1 / weights_svd[p])
      f3 = f1 - (2 * sigma) * result_sparse_smooth$pc_mfd[j, p] * (1 / weights_svd[p])
    } else {
      f2 = f1 + (2 * sigma) * result_sparse_smooth$pc_mfd[j, p] * (1 / weights_svd[p])
      f3 = f1 - (2 * sigma) * result_sparse_smooth$pc_mfd[j, p] * (1 / weights_svd[p])
    }
    if (j == 1) {
      if (p == 1) {
        par(mai = c(0.4, 0.6, 0.4, 0))
        ReMFPCA::plot(f1, ylim = c(minv[p], maxv[p]), 
                      main = paste0("PC ", j, " ( ", percent(variance_percentage[j], accuracy = 0.01), " )"), 
                      ylab = yname[p], cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.8)
      } else {
        par(mai = c(0.7, 0.6, 0.1, 0))
        ReMFPCA::plot(f1, ylim = c(minv[p], maxv[p]), ylab = yname[p], 
                      cex.lab = 1.5, cex.main = 1.5, xlab = "Time", cex.axis = 1.8)
      }
    } else {
      if (p == 1) {
        if (j == 2) {
          par(mai = c(0.4, 0.4, 0.4, 0.2))
        } else {
          par(mai = c(0.4, 0.2, 0.4, 0.4))
        }
        ReMFPCA::plot(f1, ylim = c(minv[p], maxv[p]), 
                      main = paste0("PC ", j, " ( ", percent(variance_percentage[j], accuracy = 0.01), " )"), 
                      ylab = yname[p], cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.8)
      } else {
        if (j == 2) {
          par(mai = c(0.7, 0.4, 0.1, 0.2))
        } else {
          par(mai = c(0.7, 0.2, 0.1, 0.4))
        }
        ReMFPCA::plot(f1, ylim = c(minv[p], maxv[p]), ylab = yname[p], 
                      cex.lab = 1.5, cex.main = 1.5, xlab = "Time", cex.axis = 1.8)
      }
    }
    points(x_grids, f2$eval(x_grids), pch = "+", col = 3)
    points(x_grids, f3$eval(x_grids), pch = "-", col = 2)
  }
}

###################### Figure C.6 ######################
bs_smooth <- create.bspline.basis(c(0, 1), 150)
mdbs_smooth <- Basismfd(bs_smooth)
ori_data = mvmfd_obj_ori$eval(list(seq(0, 1, length.out = 1000), seq(0, 1, length.out = 1000)))
mfd1_smooth = Mfd(X = ori_data[[1]], mdbs = mdbs_smooth)
mfd2_smooth = Mfd(X = ori_data[[2]], mdbs = mdbs_smooth)
mvmfd_obj_smooth = Mvmfd(list(mfd1_smooth, mfd2_smooth))
ori_data_smooth = mvmfd_obj_smooth$eval(list(seq(0, 1, length.out = 1000), seq(0, 1, length.out = 1000)))

mfd1_smooth_f = Mfd(X = ori_data_smooth[[1]], mdbs = mdbs2)
mfd2_smooth_f = Mfd(X = ori_data_smooth[[2]], mdbs = mdbs2)
mvmfd_obj_smooth_f = Mvmfd(list(mfd1_smooth_f, mfd2_smooth_f))

mvmfd_obj_sample_s = mvmfd_obj_smooth_f
mvmfd_obj_sample = mvmfd_obj
cluster_index = list()
cluster_index[[1]] = which(cluster_sparse_pc12$cluster == 3)
cluster_index[[2]] = which(cluster_sparse_pc12$cluster == 4)
cluster_index[[3]] = which(cluster_sparse_pc12$cluster == 2)
cluster_index[[4]] = which(cluster_sparse_pc12$cluster == 1)

colorset <- c("red", "darkblue", "darkgreen", "orange")
m <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9), nrow = 3, ncol = 4, byrow = TRUE)
layout(mat = m, heights = c(0.4, 0.4, 0.2))

# Define y-axis ticks for the log scale
log_y_ticks <- c(0.001, 0.01, 0.1, 1, 3)

for (p in 1:2) {
  for (j in 1:4) {
    if (p == 1) {
      if (j == 1) {
        par(mar = c(0.5, 5, 2, 0.2))
      } else if (j == 2) {
        par(mar = c(0.5, 3, 2, 2.2))
      } else if (j == 3) {
        par(mar = c(0.5, 1, 2, 4.2))
      } else {
        par(mar = c(0.5, 1, 2, 4.2))
      }
      ReMFPCA::plot(mvmfd_obj_smooth[cluster_index[[j]], p],
                    col = "grey", lwd = 1, main = paste0("cluster ", j, " ( n = ", length(cluster_index[[j]]), " )"), 
                    cex.lab = 2, cex.main = 2, log = "y", ylim = c(0.001, 3), yaxt = "n", cex.axis = 1.8)
      axis(2, at = log_y_ticks, labels = log_y_ticks, cex.axis = 1.8)
    } else {
      if (j == 1) {
        par(mar = c(0.5, 5, 2, 0.2))
      } else if (j == 2) {
        par(mar = c(0.5, 3, 2, 2.2))
      } else if (j == 3) {
        par(mar = c(0.5, 1, 2, 4.2))
      } else {
        par(mar = c(0.5, 1, 2, 4.2))
      }
      ReMFPCA::plot(mvmfd_obj_smooth[cluster_index[[j]], p],
                    col = "grey", lwd = 1, cex.lab = 2, cex.main = 2, ylim = c(-2, 0.5), cex.axis = 1.8)
    }
    ReMFPCA::plot(ReMFPCA::mean(mvmfd_obj_smooth[cluster_index[[j]], p]), add = TRUE, col = colorset[j], lwd = 2)
    ReMFPCA::plot(ReMFPCA::mean(mvmfd_obj_smooth[((j - 1) * 24 + 1):(j * 24), p]), add = TRUE, col = colorset[j], lty = 2, lwd = 2)
  }
}

# Add legend for Figure C.6
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
plot_colors <- colorset
legend(x = "center", inset = 0,
       legend = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "Jogging", "Standing", "Walking", "Sitting"), 
       col = plot_colors, lwd = 3, cex = 2, horiz = TRUE, lty = c(rep(1, 4), rep(2, 4)))


