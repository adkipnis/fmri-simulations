#================================================#
#          fmri-sim results evaluation           #
#                 AK, 06-10-2020                 #
#================================================#

rm(list = ls()) #clear environment
library(tidyverse)
library(afex)
library(cowplot)
library(ggExtra)
library(ggridges)
library(e1071)

scale = 1 #scaling factor for figures
show_plot = F
# mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# Get filelist
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
files = list.files(getwd())
result_files = files[grepl("results", files)]
result_files = result_files[grepl("csv", result_files)]

# Load each csv and adjust permutation index
perm_per_csv = 2
p = 0
df = read.csv(result_files[1])
for (path in result_files[2:length(result_files)]){
  tmp = read.csv(path)
  unique = max(tmp['perm_num'])
  tmp['perm_num'] = tmp['perm_num']+p
  df = rbind(df, tmp)
  p = p+unique
  
}
rm(tmp) 

# Tidy up df
df$recovered = df$recovered*100
df$oe_average = (df$sub.01_data_rdm_oe_rel + df$sub.02_data_rdm_oe_rel +
                   df$sub.03_data_rdm_oe_rel + df$sub.04_data_rdm_oe_rel +
                   df$sub.05_data_rdm_oe_rel)/5
df$GT = ordered(as.factor(df$GT), levels = c("V1_left", "V1_right",
                                             "V2_left", "V2_right",
                                             "V3_left", "V3_right",
                                             "V4_left", "V4_right",
                                             "VVC_left", "VVC_right",
                                             "VMV_left", "VMV_right",
                                             "PHA_left", "PHA_right",
                                             "FFC_left", "FFC_right",
                                             "IT_left", "IT_right",
                                             "MT+MST_left", "MT+MST_right"))
#df$GT[which(is.na(df$GT))] = "MT+MST_right" 
df$GT = recode(df$GT, "V1_left" = "V1-left",
               "V1_right" = "V1-right",
               "V2_left" = "V2-left",
               "V2_right" = "V2-right",
               "V3_left" = "V3-left",
               "V3_right" = "V3-right",
               "V4_left" = "V4-left",
               "V4_right" = "V4-right",
               "VVC_left" = "VVC-left",
               "VVC_right" = "VVC-right",
               "VMV_left" = "VMV-left",
               "VMV_right" = "VMV-right",
               "PHA_left" = "PHA-left",
               "PHA_right" = "PHA-right",
               "FFC_left" = "FFC-left",
               "FFC_right" = "FFC-right",
               "IT_left" = "IT-left",
               "IT_right" = "IT-right",
               "MT+MST_left" = "MT/MST-left",
               "MT+MST_right" = "MT/MST-right")
df$pattern_subset = as.factor(df$pattern_subset)
df$n_runs = as.factor(df$n_runs)
df$prec_type = ordered(as.factor(df$prec_type),
                       levels = c("none", "res-univariate", "res-total", "instance-based"))
df$prec_type = recode(df$prec_type, "none" = "None",
                      "res-univariate" = "Univariate",
                      "res-total" = "Time-based",
                      "instance-based" = "Instance-based")
df$snr_rel = as.factor(df$snr_rel)
df$percent_better = df$n_sig_better/19 * 100
df_best <- df %>% filter(n_runs == 32, pattern_subset == 50)
df$nc_mean = (df$nc_low + df$nc_high)/2

# Data frame for exporting values to latex 
df_summary = setNames(data.frame(matrix(ncol =2 , nrow = 0)), c("var", "num"))
df_summary[nrow(df_summary) + 1,] = c("n_full", nrow(df))
df_summary[nrow(df_summary) + 1,] = c("n_best", nrow(df_best))
df_summary[nrow(df_summary) + 1,] = c("n_datasets", length(result_files))
df_summary[nrow(df_summary) + 1,] = c("n_datasets_3", 3*length(result_files))
###############################################################################
# ROI size
###############################################################################
df_rs = df[,c("GT", "sub.01_GT_roi_size", "sub.02_GT_roi_size", "sub.03_GT_roi_size",
              "sub.04_GT_roi_size", "sub.05_GT_roi_size")]

df_rs <- gather(df_rs, sub, roi_size, c("sub.01_GT_roi_size", "sub.02_GT_roi_size", "sub.03_GT_roi_size",
                                          "sub.04_GT_roi_size", "sub.05_GT_roi_size"), factor_key=TRUE)
df_rs$sub <- recode(df_rs$sub, "sub.01_GT_roi_size" = "1",
                  "sub.02_GT_roi_size" = "2", "sub.03_GT_roi_size" = "3",
                  "sub.04_GT_roi_size" = "4", "sub.05_GT_roi_size" = "5")
df_roi_size <- df_rs %>% group_by(GT) %>% summarise(roi_size = mean(roi_size))

# roi_sizes = rbind(df$sub.01_GT_roi_size, df$sub.02_GT_roi_size,
#                   df$sub.03_GT_roi_size, df$sub.04_GT_roi_size,
#                   df$sub.05_GT_roi_size)
# oe_rels = rbind(df$sub.01_data_rdm_oe_rel, df$sub.02_data_rdm_oe_rel,
#                 df$sub.03_data_rdm_oe_rel, df$sub.04_data_rdm_oe_rel,
#                 df$sub.05_data_rdm_oe_rel)
# cor(roi_sizes, oe_rels)



###############################################################################
# Reliability
###############################################################################

###############################################################################
# OE Reliabilities per ROI
df_oe <- df %>% group_by(GT) %>%
  filter(n_runs == 32, pattern_subset == 50) %>%
  summarise(Mean_OE = mean(oe_average, na.rm =  TRUE),
            MSE_OE = sd(oe_average, na.rm =  TRUE)/sqrt(n()))

oe_range = round(range(df_oe$Mean_OE), digits = 2)
df_summary[nrow(df_summary) + 1,] = c("oe_range_1", oe_range[1])
df_summary[nrow(df_summary) + 1,] = c("oe_range_2", oe_range[2])

# Early vs late areas
oe_EVC = round(mean(df_oe$Mean_OE[1:8]), digits = 2)
oe_rest = round(mean(df_oe$Mean_OE[9:20]), digits = 2)
df_summary[nrow(df_summary) + 1,] = c("oe_EVC", oe_EVC)
df_summary[nrow(df_summary) + 1,] = c("oe_rest", oe_rest)

# Left vs right areas
oe_lr_corr = cor(df_oe$Mean_OE[seq(1, nrow(df_oe), 2)],
                 df_oe$Mean_OE[seq(2, nrow(df_oe), 2)])
df_summary[nrow(df_summary) + 1,] = c("oe_lr_corr", round(oe_lr_corr, digits = 2))


plot_oe_rel <-  ggplot(df_best, aes(x=GT, y=oe_average)) + 
  geom_boxplot(color="black", fill = "lightsteelblue", outlier.alpha = 0.1)+
  coord_fixed(ratio = 7, ylim=c(0, 1), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  # labs(caption = "*for RDMs from 32 runs with 50 conditions")+
  # xlab("ROI")+
  ylab(expression(atop(bold("Odd-Even Reliability*"), "[Pearson's r]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 10, angle = 70, hjust=1.1),
        axis.text.y = element_text(size = 12))
if(show_plot){plot_oe_rel}
# ggsave("results_bp_ROI_OE.png", device = "png",  scale = scale, dpi = 320)

###############################################################################
# Noise Ceilings (best case)
df_nc <- df %>% group_by(GT) %>%
  filter(n_runs == 32, pattern_subset == 50) %>%
  summarise(Lower_bound = mean(nc_low), Upper_bound = mean(nc_high)) 
df_nc$roi_size = df_roi_size$roi_size 
  
# Range
nc_range = range(df_nc$Lower_bound)
df_summary[nrow(df_summary) + 1,] = c("nc_range_1", round(nc_range[1], digits = 2))
df_summary[nrow(df_summary) + 1,] = c("nc_range_2", round(nc_range[2], digits = 2))

# ROI size vs noise ceiling
nc_rs_corr = cor(df_nc$Lower_bound, df_nc$roi_size) 
df_summary[nrow(df_summary) + 1,] = c("nc_rs_corr", format(round(nc_rs_corr, digits = 2), nsmall = 2))


# Left vs right areas
nc_lr_corr = cor(df_nc$Lower_bound[seq(1, nrow(df_nc), 2)],
                 df_nc$Lower_bound[seq(2, nrow(df_nc), 2)]) 
df_summary[nrow(df_summary) + 1,] = c("nc_lr_corr", format(round(nc_lr_corr, digits = 2), nsmall = 2))

# Early vs late areas
nc_EVC = mean(df_nc$Lower_bound[1:8])
nc_rest = mean(df_nc$Lower_bound[9:20])
df_summary[nrow(df_summary) + 1,] = c("nc_EVC", round(nc_EVC, digits = 2))
df_summary[nrow(df_summary) + 1,] = c("nc_rest", round(nc_rest, digits = 2))

# Correlation with OE
nc_oe_dif = mean(df_nc$Lower_bound - df_oe$Mean_OE)
nc_oe_cor = cor(df_oe$Mean_OE, df_nc$Lower_bound)
df_summary[nrow(df_summary) + 1,] = c("nc_oe_dif", round(nc_oe_dif, digits = 2))
df_summary[nrow(df_summary) + 1,] = c("nc_oe_cor", round(nc_oe_cor, digits = 2))

plot_nc <- ggplot(df_nc) +
  aes(ymin = Lower_bound,
      ymax = Upper_bound,
      xmin = as.numeric(GT) - .3,
      xmax = as.numeric(GT) + .3,
      x = GT) +
  geom_rect(position=position_dodge(), col = 'black', fill = 'lightgray', stat="identity") + 
  geom_text(aes(x=GT, y=Upper_bound+0.03, label=roi_size), size=2.5) +
  coord_fixed(ratio = 7, ylim=c(0, 1), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  labs(caption = "*for RDMs from 32 runs with 50 conditions")+
  xlab("ROI")+
  ylab(expression(atop(bold("Noise Ceiling*"), "[Pearson's r]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 10, angle = 70, hjust=1.1),
        axis.text.y = element_text(size = 12))
if(show_plot){plot_nc}
# ggsave("results_ROI_NC.png", device = "png",  scale = scale, dpi = 320)

###############################################################################
# OE-NC Grid Plot 
plot_grid(plot_oe_rel, plot_nc, labels=c("A", "B"), ncol = 1, nrow = 2)
ggsave("results_grid_OE_NC.png", device = "png",  scale = scale, dpi = 320,
       width = 9, height = 9)


###############################################################################
# OE Reliabilities per ROI and SNR

AOV_oe_snr <- aov_ez("GT", "oe_average", df_best,
                     within = "snr_rel", na.rm = TRUE)
AOV_oe_snr_table <- nice(AOV_oe_snr)
df_summary[nrow(df_summary) + 1,] = c("eta_OE", round(as.numeric(AOV_oe_snr_table[1,5]), digits = 2))


plot_oe_rel_snr <- ggplot(df_best) +
  aes(y = oe_average, x = GT,  fill = snr_rel) +
  geom_boxplot(color="black", outlier.alpha = 0.1)+
  coord_fixed(ratio = 7, ylim=c(0, 1), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  labs(fill = "Relative SNR")+
  scale_fill_brewer() +
  xlab("ROI")+
  ylab(expression(atop(bold("Odd-Even Reliability*"), "[Pearson's r]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 10, angle = 70, hjust=1.1),
        axis.text.y = element_text(size = 12))+
  theme(legend.position="top")
if(show_plot){plot_oe_rel_snr}
# ggsave("results_bp_ROI_OE_SNR.png", device = "png",  scale = scale, dpi = 320)


###############################################################################
# Noise Ceilings per SNR (best case)
df_nc_snr <- df %>% group_by(GT, snr_rel) %>%
  filter(n_runs == 32,pattern_subset == 50) %>%
  summarise(Lower_bound = mean(nc_low), Upper_bound = mean(nc_high)) 

AOV_nc_snr <- aov_ez("GT", "Lower_bound", df_nc_snr,
                     within = "snr_rel")
AOV_nc_snr_table <- nice(AOV_nc_snr)
df_summary[nrow(df_summary) + 1,] = c("eta_NC", round(as.numeric(AOV_nc_snr_table[1,5]), digits = 2))


plot_nc_snr <- ggplot(df_nc_snr) +
  aes(ymin = Lower_bound,
      ymax = Upper_bound,
      xmin = as.numeric(GT) - .3,
      xmax = as.numeric(GT) + .3,
      x = GT, fill = snr_rel) +
  geom_rect(position=position_dodge(), col = 'black', stat="identity") + 
  labs(fill = "relative SNR")+
  scale_fill_brewer() +
  coord_fixed(ratio = 7, ylim=c(0, 1), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  labs(caption = "*for RDMs from 32 runs with 50 conditions",
       fill = "Relative SNR")+
  xlab("ROI")+
  ylab(expression(atop(bold("Noise Ceiling*"), "[Lower bound; Pearson's r]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 10, angle = 70, hjust=1.1),
        axis.text.y = element_text(size = 12))+
  theme(legend.position="none")
if(show_plot){plot_nc_snr}
# ggsave("results_ROI_NC_SNR.png", device = "png",  scale = scale, dpi = 320)

###############################################################################
# OE-NC SNR Grid Plot 
plot_grid(plot_oe_rel_snr, plot_nc_snr, labels=c("A", "B"), ncol = 1, nrow = 2)
ggsave("results_grid_OE_NC_SNR.png", device = "png",  scale = scale, dpi = 320,
       width = 9, height = 9)


###############################################################################
# Point estimates
###############################################################################

###############################################################################
# Point estimates: Per GT and SNR
df_pe <- df %>% group_by(GT) %>%
  summarise(Mean_Perf = mean(point_est), MSE_Perf = sd(point_est)/sqrt(n()),
            Mean_StD = mean(std), MSE_StD = sd(std)/sqrt(n())) 

# Left vs right areas
perf_lr_corr = cor(df_pe$Mean_Perf[seq(1, nrow(df_nc), 2)],
                   df_pe$Mean_Perf[seq(2, nrow(df_nc), 2)]) 
df_summary[nrow(df_summary) + 1,] = c("perf_lr_corr", round(perf_lr_corr, digits = 2))


df_pe_snr <- df %>% group_by(GT, snr_rel) %>%
  summarise(Mean_Perf = mean(point_est), MSE_Perf = sd(point_est)/sqrt(n()),
            Mean_StD = mean(std), MSE_StD = sd(std)/sqrt(n())) 

# Range
perf_range = range(df_pe_snr$Mean_Perf)
df_summary[nrow(df_summary) + 1,] = c("perf_range_1", round(perf_range[1], digits = 2))
df_summary[nrow(df_summary) + 1,] = c("perf_range_2", round(perf_range[2], digits = 2))


est_snr0.5 = df %>% filter(snr_rel == 0.5)
est_snr1 = df %>% filter(snr_rel == 1)
est_snr2 = df %>% filter(snr_rel == 2)
kurt_0.5 = kurtosis(est_snr0.5$point_est)
kurt_1 = kurtosis(est_snr1$point_est)
kurt_2 = kurtosis(est_snr2$point_est)

df_summary[nrow(df_summary) + 1,] = c("kurt_0.5", round(kurt_0.5, digits = 2))
df_summary[nrow(df_summary) + 1,] = c("kurt_1", round(kurt_1, digits = 2))
df_summary[nrow(df_summary) + 1,] = c("kurt_2", format(round(kurt_2, digits = 2), nsmall = 2))



AOV_perf_snr <- aov_ez("GT", "point_est", df,
                     within = "snr_rel")
AOV_perf_snr_table <- nice(AOV_nc_snr)
df_summary[nrow(df_summary) + 1,] = c("eta_perf", round(as.numeric(AOV_perf_snr_table[1,5]), digits = 2))

plot_sim <- ggplot(df) +
  aes(y = point_est, x = GT, fill = snr_rel) +
  geom_boxplot(color="black", outlier.alpha = 0.1)+
  coord_fixed(ratio = 7, ylim=c(0, 1), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  labs(fill = "Relative SNR")+
  ylab(expression(atop(bold("Model Performance"), "[Similarity with Data RDM; Pearson's r]")))+
  scale_fill_brewer() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 10, angle = 70, hjust=1.1),
        axis.text.y = element_text(size = 12))+
  theme(legend.position="top")
if(show_plot){plot_sim}
# ggsave("results_bp_ROI_cor.png", device = "png",  scale = scale, dpi = 320)


###############################################################################
# Point estimate variability: Per GT and SNR
std_skew = skewness(df$std)
std_median = median(df$std)
df_summary[nrow(df_summary) + 1,] = c("std_skew", format(round(std_skew, digits = 2), nsmall = 2))
df_summary[nrow(df_summary) + 1,] = c("std_median", round(std_median, digits = 2))

AOV_perf_snr_var <- aov_ez("GT", "std", df,
                       within = "snr_rel")
AOV_perf_snr_var_table <- nice(AOV_perf_snr_var)
df_summary[nrow(df_summary) + 1,] = c("eta_var", round(as.numeric(AOV_perf_snr_var_table[1,5]), digits = 2))



plot_sim_var <- ggplot(df) +
  aes(y=std, x = GT, fill = snr_rel) +
  geom_boxplot(color="black", outlier.alpha = 0.1)+
  coord_fixed(ratio = 14, ylim=c(0, 0.5), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,0.5,0.05)) +
  scale_fill_brewer() +
  labs(fill = "Relative SNR")+
  xlab("ROI")+
  ylab(expression(atop(bold("Model Uncertainty"), "[Standard Deviation of Performance]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 10, angle = 70, hjust=1.1),
        axis.text.y = element_text(size = 12))+
  theme(legend.position="none")
if(show_plot){plot_sim_var}

# ggsave("results_bp_ROI_cor_SD.png", device = "png",  scale = scale, dpi = 320)

###############################################################################
# Point-Estimates Grid Plot 
plot_grid(plot_sim, plot_sim_var, labels=c("A", "B"), ncol = 1, nrow = 2)
ggsave("results_grid_sim_sim_var.png", device = "png",  scale = scale, dpi = 320,
       width = 9, height = 9)

###############################################################################
# Scatter plot for point estimates per NC (lower bound)

nc_perf_cor = cor(df$nc_low, df$point_est)
df_summary[nrow(df_summary) + 1,] = c("nc_perf_cor", round(nc_perf_cor, digits = 2))

scatterplot <- ggplot(df) + 
  aes(x = nc_low, y = point_est, color=snr_rel) + 
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_vline(xintercept=0, linetype="dashed", color = "gray") +
  geom_point(alpha=0.5) +
  scale_color_brewer(type='seq')+
  coord_cartesian(xlim=c(-1, 1), ylim=c(-0.1, 1), expand = FALSE)+
  scale_x_continuous(breaks=seq(-1,1,0.1)) +
  scale_y_continuous(breaks=seq(-0.1,1,0.1)) +
  labs(color = "Relative SNR")+
  xlab(expression(atop(bold("Noise Ceiling of Data"), "[Lower bound; Pearson's r]")))+
  ylab(expression(atop(bold("Model Performance"), "[Similarity with Ground Truth; Pearson's r]")))+
  theme_bw()+
  theme(axis.title.x = element_text(vjust = -2, size = 14))+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 11))+
  theme(legend.position=c(0.15, 0.8),
        legend.box.background = element_rect(colour = "black"))

  
scatter_marginal <- ggMarginal(
  scatterplot,
  type = "density",
  margins = "both",
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)
if(show_plot){scatter_marginal}

ggsave("results_scatter_NC_perf.png", plot = scatter_marginal, device = "png",
       scale = scale, dpi = 320, width = 9, height = 6)

# ###############################################################################
# # Scatter plot for point estimates per OE-rel
# 
# cor(df_best$oe_average, df_best$point_est)
# 
# 
# scatterplot2 <- ggplot(df_best) +
#   aes(x = oe_average, y = point_est, color=snr_rel) +
#   geom_hline(yintercept=0, linetype="dashed", color = "gray") +
#   geom_vline(xintercept=0, linetype="dashed", color = "gray") +
#   geom_point() +
#   scale_color_brewer(type='seq')+
#   coord_cartesian(xlim=c(-1, 1), ylim=c(0, 1), expand = FALSE)+
#   scale_x_continuous(breaks=seq(-1,1,0.1)) +
#   scale_y_continuous(breaks=seq(-1,1,0.1)) +
#   labs(color = "Relative SNR")+
#   xlab(expression(atop(bold("Noise Ceiling of Data"), "[Lower bound; Pearson's r]")))+
#   ylab(expression(atop(bold("Similarity with Ground Truth"), "[Pearson's r]")))+
#   theme_bw()+
#   theme(axis.title.x = element_text(vjust = -2, size = 14))+
#   theme(axis.title.y = element_text(
#     size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
#   theme(axis.text.x = element_text(size = 11),
#         axis.text.y = element_text(size = 11))+
#   theme(legend.position=c(0.15, 0.8),
#         legend.box.background = element_rect(colour = "black"))
# 
# 
# scatter_marginal2 <- ggMarginal(
#   scatterplot2,
#   type = "density",
#   margins = "both",
#   size = 5,
#   groupColour = TRUE,
#   groupFill = TRUE
# )
# 
# ggsave("results_scatter_OE_perf.png", plot = scatter_marginal, device = "png",
#        scale = scale, dpi = 320)
###############################################################################
# Model recovery
###############################################################################

###############################################################################
# Model Recovery per ROI

df_rec <- df %>% group_by(GT) %>%
  summarise(Mean_Acc = mean(recovered), MSE_Acc = sd(recovered)/sqrt(n()),
            better = round(mean(n_sig_better)/19, digits=2))
rec_range = range(df_rec$Mean_Acc)
df_summary[nrow(df_summary) + 1,] = c("rec_range_1", round(rec_range[1], digits = 2))
df_summary[nrow(df_summary) + 1,] = c("rec_range_2", round(rec_range[2], digits = 2))

# plot_rec <- ggplot(df_rec, aes(y=Mean_Acc, x = GT)) +
#   geom_bar(stat="identity", width = 0.9, color="black", fill = "palegreen3", position=position_dodge())+
#   geom_errorbar(aes(ymin = Mean_Acc-MSE_Acc, ymax = Mean_Acc+MSE_Acc), width=.2,
#                 position=position_dodge(0.9))+
#   geom_text(aes(label=better), position=position_dodge(width=0.9), vjust=-1)+
#   scale_fill_brewer() +
#   coord_fixed(ratio = 1/9, ylim=c(45, 100), expand = FALSE)+
#   # coord_cartesian(ylim=c(45, 100), expand = FALSE)+
#   scale_y_continuous(breaks=seq(50,100,5)) +
#   labs(fill = "Relative SNR")+
#   # xlab("ROI")+
#   ylab(expression(atop(bold("Model Recovery"), "[Percentage correct]")))+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5, size = 16))+
#   theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
#   theme(axis.title.x=element_blank())+
#   theme(axis.title.y = element_text(
#     size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
#   theme(axis.text.x = element_text(size = 10, angle = 70, hjust=1.1),
#         axis.text.y = element_text(size = 12))+
#   theme(legend.position="top")
# plot_rec

###############################################################################
# Model Recovery per ROI and SNR
AOV_rec_snr <- aov_ez("GT", "recovered", df,
                           within = "snr_rel")
AOV_rec_snr_table <- nice(AOV_rec_snr)
df_summary[nrow(df_summary) + 1,] = c("eta_rec_snr", round(as.numeric(AOV_rec_snr_table[1,5]), digits = 2))

df_rec_snr <- df %>% group_by(GT, snr_rel) %>%
  summarise(Mean_Acc = mean(recovered), MSE_Acc = sd(recovered)/sqrt(n()),
            better = mean(n_sig_better)) 

plot_rec_snr <- ggplot(df_rec_snr, aes(y=Mean_Acc, x = GT, fill = snr_rel)) +
  geom_bar(stat="identity", width = 0.9, color="black", position=position_dodge())+
  geom_errorbar(aes(ymin = Mean_Acc-MSE_Acc, ymax = Mean_Acc+MSE_Acc), width=.2,
                position=position_dodge(0.9))+
  scale_fill_brewer() +
  # coord_fixed(ratio = 1/7, ylim=c(45, 100), expand = FALSE)+
  coord_cartesian(ylim=c(0, 100), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,100,10)) +
  labs(fill = "Relative SNR")+
  xlab("ROI")+
  ylab(expression(atop(bold("Model Recovery"), "[Percentage correct]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 10, angle = 70, hjust=1.1),
        axis.text.y = element_text(size = 12))+
  theme(legend.position="top")

if(show_plot){plot_rec_snr}
ggsave("results_ROI_acc_SNR.png", device = "png",  scale = scale, dpi = 320)


# 
# ###############################################################################
# # Model Recovery Grid Plot 
# plot_grid(plot_rec, plot_rec_snr, labels=c("A", "B"), ncol = 1, nrow = 2)
# ggsave("results_grid_rec_rec_var.png", device = "png",  scale = scale, dpi = 320,
#        width = 9, height = 9)

###############################################################################
# Noise Normalization
###############################################################################

###############################################################################
# Noise normalization: Recovery

AOV_rec_nn <- aov_ez("GT", "recovered", df,
                     within = "prec_type")
AOV_rec_nn_table <- nice(AOV_rec_nn)
df_summary[nrow(df_summary) + 1,] = c("eta_rec_nn", round(as.numeric(AOV_rec_nn_table[1,5]), digits = 2))


df_nn <- df %>% group_by(prec_type, snr_rel) %>%
  summarise(Mean_Acc = mean(recovered), MSE_Acc = sd(recovered)/sqrt(n()),
            Mean_StD = mean(std), MSE_StD = sd(std)/sqrt(n())) 

df_summary[nrow(df_summary) + 1,] = c("perf_unn", round(mean(df_nn$Mean_Acc[1:6]), digits = 2))
df_summary[nrow(df_summary) + 1,] = c("perf_mnn", round(mean(df_nn$Mean_Acc[7:12]), digits = 2))

plot_nn <- ggplot(df_nn) +
  aes(y = Mean_Acc, x = snr_rel, fill = prec_type) +
  geom_bar(stat="identity", width = 0.9, color="black", position=position_dodge())+
  geom_errorbar(aes(ymin = Mean_Acc-MSE_Acc, ymax = Mean_Acc+MSE_Acc), width=.2,
                position=position_dodge(0.9))+
  coord_cartesian(ylim=c(0, 100), expand = FALSE)+
  # coord_fixed(ratio = 1/50, ylim=c(45, 100), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,100,10)) +
  scale_fill_brewer(palette = "OrRd") +
  labs(fill = "Noise Normalization")+
  ylab(expression(atop(bold("Model Recovery"), "[Percentage correct]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  theme(legend.position="top")

if(show_plot){plot_nn}
# ggsave("results_NN_Acc.png", device = "png",  scale = scale, dpi = 320)

###############################################################################
# Noise normalization: Estimate Variability 

AOV_rec_nn_var <- aov_ez("GT", "std", df,
                         within = "prec_type")
AOV_rec_nn_var_table <- nice(AOV_rec_nn_var)
df_summary[nrow(df_summary) + 1,] = c("eta_rec_nn_var", round(as.numeric(AOV_rec_nn_var_table[1,5]), digits = 2))


plot_nn_var <- ggplot(df) +
  aes(y = std, x = snr_rel, fill = prec_type) +
  geom_boxplot(color="black", outlier.alpha = 0.1)+
  # coord_fixed(ratio = 2, ylim=c(0, 0.5), expand = FALSE)+
  coord_cartesian(ylim=c(0, .5), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,0.5,0.05)) +
  scale_fill_brewer(palette = "OrRd") +
  labs(fill = "Noise Normalization")+
  xlab("Relative SNR")+
  ylab(expression(atop(bold("Model Uncertainty"), "[Standard Deviation of Performance]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  theme(legend.position="blank")

if(show_plot){plot_nn_var}
# ggsave("results_bp_NN_SD.png", device = "png",  scale = scale, dpi = 320)

###############################################################################
# Noise Normalization Grid Plot 
plot_grid(plot_nn, plot_nn_var, labels=c("A", "B"), ncol = 1, nrow = 2)
ggsave("results_grid_nn_nn_var.png", device = "png",  scale = scale, dpi = 320,
       width = 9, height = 9)



###############################################################################
# Noise Normalization and rc
###############################################################################

###############################################################################
# Noise normalization: Recovery


df_nn_rc <- df %>% group_by(prec_type, n_runs, pattern_subset) %>%
  summarise(Mean_Acc = mean(recovered), MSE_Acc = sd(recovered)/sqrt(n()),
            Mean_StD = mean(std), MSE_StD = sd(std)/sqrt(n())) 

plot_nn_rc <- ggplot(df_nn_rc) +
  aes(y = Mean_Acc, x = n_runs, fill = prec_type) +
  geom_bar(stat="identity", width = 0.9, color="black", position=position_dodge())+
  geom_errorbar(aes(ymin = Mean_Acc-MSE_Acc, ymax = Mean_Acc+MSE_Acc), width=.2,
                position=position_dodge(0.9))+
  coord_cartesian(ylim=c(0, 100), expand = FALSE)+
  facet_wrap(facets = "pattern_subset", drop = T)+
  # coord_fixed(ratio = 1/50, ylim=c(45, 100), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,100,10)) +
  scale_fill_brewer(palette = "OrRd") +
  labs(fill = "Noise Normalization")+
  ylab(expression(atop(bold("Model Recovery"), "[Percentage correct]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  theme(legend.position="top")

if(show_plot){plot_nn_rc}
# ggsave("results_NN_Acc.png", device = "png",  scale = scale, dpi = 320)

###############################################################################
# Noise normalization: Estimate Variability 


plot_nn_rc_var <- ggplot(df) +
  aes(y = std, x = n_runs, fill = prec_type) +
  geom_boxplot(color="black", outlier.alpha = 0.1)+
  # coord_fixed(ratio = 2, ylim=c(0, 0.5), expand = FALSE)+
  coord_cartesian(ylim=c(0, .5), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,0.5,0.05)) +
  scale_fill_brewer(palette = "OrRd") +
  facet_wrap(facets = "pattern_subset", drop = T)+
  labs(fill = "Noise Normalization")+
  xlab("Number of Runs")+
  ylab(expression(atop(bold("Model Uncertainty"), "[Standard Deviation of Performance]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  theme(legend.position="blank")

if(show_plot){plot_nn_rc_var}
# ggsave("results_bp_NN_SD.png", device = "png",  scale = scale, dpi = 320)

###############################################################################
# Noise Normalization Grid Plot 
plot_grid(plot_nn_rc, plot_nn_rc_var, labels=c("A", "B"), ncol = 1, nrow = 2)
ggsave("results_grid_nn_rc_nn_rc_var.png", device = "png",  scale = scale, dpi = 320,
       width = 9, height = 12)



###############################################################################
# Runs vs. Conditions
###############################################################################
###############################################################################
# Runs vs. Conditions: reliability
# plot_rc_rel <- ggplot(df) +
#   aes(y = oe_average, x = n_runs, fill = pattern_subset) +
#   geom_boxplot(color="black", outlier.alpha = 0.1)+
#   coord_cartesian(ylim=c(0, 1), expand = FALSE)+
#   scale_y_continuous(breaks=seq(0,1,0.1)) +
#   scale_fill_brewer(palette = "Purples") +
#   labs(fill = "Number of Conditions")+
#   xlab("Number of Runs")+
#   ylab(expression(atop(bold("Odd-Even Reliability"), "[Pearson's r]")))+
#   facet_wrap(facets = "pattern_subset_type", drop = T) +
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5, size = 16))+
#   theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
#   theme(axis.title.x = element_text(size = 14))+
#   theme(axis.title.y = element_text(
#     size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
#   theme(axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12))+
#   theme(legend.position="top")
# 
# if(show_plot){plot_rc_rel}


###############################################################################
# Runs vs. Conditions: Recovery
df_rc <- df %>% group_by(n_runs, pattern_subset) %>%
  summarise(Mean_Acc = mean(recovered), MSE_Acc = sd(recovered)/sqrt(n()),
            Mean_StD = mean(std), MSE_StD = sd(std)/sqrt(n())) 

mean_rec_5p = mean(df_rc[df_rc$pattern_subset == 5,]$Mean_Acc)
mean_rec_10p = mean(df_rc[df_rc$pattern_subset == 10,]$Mean_Acc)
df_summary[nrow(df_summary) + 1,] = c("mean_rec_5p", round(mean_rec_5p, digits = 2))
df_summary[nrow(df_summary) + 1,] = c("mean_rec_10p", round(mean_rec_10p, digits = 2))

AOV_rec_rc <- aov_ez("GT", "recovered", df,
                         within = "pattern_subset")
AOV_rec_rc_table <- nice(AOV_rec_rc)
df_summary[nrow(df_summary) + 1,] = c("eta_rec_rc", round(as.numeric(AOV_rec_rc_table[1,5]), digits = 2))


plot_rc <- ggplot(df_rc) +
  aes(y = Mean_Acc, x = n_runs, fill = pattern_subset) +
  geom_bar(stat="identity", width = 0.9, color="black", position=position_dodge())+
  geom_errorbar(aes(ymin = Mean_Acc-MSE_Acc, ymax = Mean_Acc+MSE_Acc), width=.2,
                position=position_dodge(0.9))+
  coord_cartesian(ylim=c(0, 100), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,100,10)) +
  scale_fill_brewer(palette = "Purples") +
  labs(fill = "Number of Conditions")+
  xlab("Number of Runs")+
  ylab(expression(atop(bold("Model Recovery"), "[Percentage correct]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  theme(legend.position="top")

if(show_plot){plot_rc}
# ggsave("results_Runs_Conditions_Acc.png", device = "png",  scale = scale, dpi = 320)

###############################################################################
# Runs vs. Conditions: Estimate Variability
mean_var_5p = mean(df_rc[df_rc$pattern_subset == 5,]$Mean_StD)
mean_var_10p = mean(df_rc[df_rc$pattern_subset == 10,]$Mean_StD)
mean_var_25p = mean(df_rc[df_rc$pattern_subset == 25,]$Mean_StD)
difference_10_25 =  mean_var_10p - mean_var_25p
df_summary[nrow(df_summary) + 1,] = c("mean_var_5p", round(mean_var_5p, digits = 2))
df_summary[nrow(df_summary) + 1,] = c("mean_var_10p", round(mean_var_10p, digits = 2))
df_summary[nrow(df_summary) + 1,] = c("difference_10_25", round(difference_10_25, digits = 2))

AOV_rec_rc_var <- aov_ez("GT", "std", df,
                     within = "pattern_subset")
AOV_rec_rc_var_table <- nice(AOV_rec_rc_var)
df_summary[nrow(df_summary) + 1,] = c("eta_rec_rc_var", round(as.numeric(AOV_rec_rc_var_table[1,5]), digits = 2))

plot_rc_var <- ggplot(df) +
  aes(y = std, x = n_runs, fill = pattern_subset) +
  geom_boxplot(color="black", outlier.alpha = 0.1)+
  coord_cartesian(ylim=c(0, .5), expand = FALSE)+
  scale_y_continuous(breaks=seq(0,0.5,0.05)) +
  scale_fill_brewer(palette = "Purples") +
  labs(fill = "Number of Conditions")+
  xlab("Number of Runs")+
  ylab(expression(atop(bold("Model Uncertainty"), "[Standard Deviation of Performance]")))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  theme(legend.position="none")

if(show_plot){plot_rc_var}
# ggsave("results_bp_Runs_Conditions_SD.png", device = "png",  scale = scale, dpi = 320)

###############################################################################
# Run-Condition Grid Plot 
plot_grid(plot_rc, plot_rc_var, labels=c("A", "B"), ncol = 1, nrow = 2)
ggsave("results_grid_rc_rc_var.png", device = "png",  scale = scale, dpi = 320,
       width = 9, height = 9)

###############################################################################
# Scatter plot for point estimates per NC (lower bound) per runs

scatterplot <- ggplot(df) + 
  aes(x = nc_low, y = point_est, color=pattern_subset) + 
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_vline(xintercept=0, linetype="dashed", color = "gray") +
  geom_point(alpha=0.5) +
  scale_color_brewer(type='seq', palette = "Purples")+
  coord_cartesian(xlim=c(-1, 1), ylim=c(-0.1, 1), expand = FALSE)+
  scale_x_continuous(breaks=seq(-1,1,0.1)) +
  scale_y_continuous(breaks=seq(-0.1,1,0.1)) +
  labs(color = "Number of Conditions")+
  xlab(expression(atop(bold("Noise Ceiling of Data"), "[Lower bound; Pearson's r]")))+
  ylab(expression(atop(bold("Model Performance"), "[Similarity with Ground Truth; Pearson's r]")))+
  theme_bw()+
  theme(axis.title.x = element_text(vjust = -2, size = 14))+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 11))+
  theme(legend.position=c(0.15, 0.8),
        legend.box.background = element_rect(colour = "black"))


scatter_marginal_rc <- ggMarginal(
  scatterplot,
  type = "density",
  margins = "both",
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)
if(show_plot){scatter_marginal_rc}

ggsave("results_scatter_NC_perf_rc.png", plot = scatter_marginal_rc, device = "png",
       scale = scale, dpi = 320, width = 9, height = 6)


###############################################################################
# Model recovery significant model comparisons

# df_serial <- df %>% filter(pattern_subset_type == 'orderly* sampled conditions')

plot_sig_comp <- ggplot(df) +
  aes(x = percent_better) + 
  geom_histogram(color="black",  alpha=0.5, position="identity", bins = 20)+
  facet_wrap(facets = "pattern_subset", drop = T)+
  coord_cartesian(ylim=c(0, 1000), expand = FALSE)+
  scale_x_continuous(breaks=seq(0,100,20))+
  xlab(expression(atop(bold("Model Persuasiveness"), "[Percentage models sig. outperformed]")))+
  ylab(expression(atop(bold("Count of occurences"), "[out of N observations]")))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 14, vjust = -1))+
  theme(axis.title.y = element_text(
    size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
if(show_plot){plot_sig_comp}


# ggplot(df_serial, aes(x = percent_better, y = pattern_subset)) +
#   geom_density_ridges(aes(height=..density..), scale= 0.95, stat="density")+
#   facet_wrap(facets = "n_runs", drop = T)


ggsave("results_conditions_hist.png", device = "png",  scale = scale, dpi = 320)

################################################################################
write.table(df_summary, "summary.txt", col.names = F, row.names=F, sep=" = ", quote = FALSE)
