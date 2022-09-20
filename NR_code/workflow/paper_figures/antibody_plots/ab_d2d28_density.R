library(tidyverse)

ab_titer_joined <- readRDS("ab_data/reshaped/ab_titer_joined.rds")

plot_dat <- ab_titer_joined %>%
  filter(timepoint %in% c("D2_4weeks")) %>%
  select(`Code ID`, S2P_num, RBD_num, Readout, timepoint, date_samp) %>%
  gather(key = antigen, value = titer, -c(`Code ID`, Readout, timepoint, date_samp)) %>%
  mutate(antigen = gsub("_num", "", antigen)) %>%
  mutate(antigen = gsub("S2P", "S-2P", antigen)) %>%
  mutate(antigen_isotype = paste(antigen, Readout, sep = " ")) %>%
  group_by(`Code ID`, antigen_isotype) %>%
  summarise(titer = mean(titer, na.rm = TRUE)) %>%
  ungroup()

plot_dat <- plot_dat %>%
  filter(!grepl("IgM", antigen_isotype))

mean_cv_dat <- plot_dat %>%
  group_by(antigen_isotype) %>%
  summarize(mean = mean(titer),
            sd = sd(titer)) %>%
  mutate(cv = sd /mean)

mean_cv_dat <- mean_cv_dat %>%
  mutate(annotation = paste0("mean = ", formatC(mean, format = "e", digits = 2), "\n cv = ", round(cv, 2)))

plot_dat <- plot_dat %>%
  left_join(mean_cv_dat)

pdf("plots/paper_figures/ab_titer_density_d2D28_nonparametric_bw_selection.pdf", height = 4, width = 4)
#https://aakinshin.net/posts/kde-bw/
p <- ggplot(plot_dat, aes(x = titer)) +
  stat_bin(aes(y=..density..), bins = 20, fill = "grey63") +
  stat_density(color = "red", bw = "SJ", fill = NA) +
  geom_text(data = mean_cv_dat, x = -Inf, y = Inf, hjust = -.1, vjust = 1.2, aes(label = annotation), size = 3) +
  facet_wrap(~antigen_isotype) +
  scale_x_log10() +theme_bw() +
  geom_rug() +
  xlab("Au/ml") +
  ylab("Density") +
  ggtitle("Sheather-Jones Bandwidth selection")
print(p)

p <- ggplot(plot_dat, aes(x = titer)) +
  stat_bin(aes(y=..density..), bins = 20, fill = "grey63") +
  stat_density(color = "red", bw = "ucv", fill = NA) +
  geom_text(data = mean_cv_dat, x = -Inf, y = Inf, hjust = -.1, vjust = 1.2, aes(label = annotation), size = 3) +
  facet_wrap(~antigen_isotype) +
  scale_x_log10() +theme_bw() +
  geom_rug()+
  xlab("Au/ml") +
  ylab("Density") +
  ggtitle("unbiased cross validation bandwidth selection")
print(p)

p <- ggplot(plot_dat, aes(x = titer)) +
  stat_bin(aes(y=..density..), bins = 20, fill = "grey63") +
  stat_density(color = "red", bw = "bcv", fill = NA) +
  geom_text(data = mean_cv_dat, x = -Inf, y = Inf, hjust = -.1, vjust = 1.2, aes(label = annotation), size = 3) +
  #geom_label(data = mean_cv_dat, x = -Inf, y = Inf, hjust = -.1, vjust = 1.2, aes(label = annotation), size = 3) +
  facet_wrap(~antigen_isotype) +
  scale_x_log10() +theme_bw() +
  #annotation_logticks()+
  geom_rug()+
  xlab("Au/ml") +
  ylab("Density") +
  ggtitle("Biased cross validation bandwidth selection")
print(p)
dev.off()

file.remove("plots/paper_figures/ab_titer_density_d2D28.pdf")

pdf("plots/paper_figures/ab_titer_density_d2D28.pdf", height = 4, width = 4)


p <- ggplot(plot_dat, aes(x = titer)) +
  geom_density(adjust = .8)+
  facet_wrap(~antigen_isotype) +
  scale_x_log10() +theme_bw() +
  geom_rug()+
  ggtitle("d2D28 Titer")
print(p)

p <- ggplot(plot_dat, aes(x = titer)) +
  geom_histogram(bins = 20)+
  geom_rug()+
  facet_wrap(~antigen_isotype) +
  scale_x_log10() +theme_bw() +
  ggtitle("d2D28 Titer")
print(p)


dev.off()




pdf("plots/paper_figures/ab_titer_density_d2D28_vary_bandwidth.pdf", height = 4, width = 4)
for(bw in seq(.6, 2, by = .2)){
  nbins <- round((20/bw))
  p <- ggplot(plot_dat, aes(x = titer)) +
    stat_bin(aes(y=..density..), bins = nbins, fill = "blue2") +
    geom_density(adjust = bw, color = "red")+
    geom_text(x = -Inf, y = Inf, hjust = -.1, vjust = 1.2, aes(label = annotation), size = 3) +
    facet_wrap(~antigen_isotype) +
    scale_x_log10() +theme_bw() +
    geom_rug()+
    ggtitle(paste("d2D28 Titer", "bandwdith =", bw, "\nn hist bins = ", nbins))
  print(p)
  
}


dev.off()


pdf("plots/paper_figures/ab_titer_density_d2D28_vary_bandwidth_nolog.pdf", height = 4, width = 4)
for(bw in seq(.6, 2, by = .2)){
  nbins <- round((20/bw))
  p <- ggplot(plot_dat, aes(x = titer)) +
    stat_bin(aes(y=..density..), bins = nbins, fill = "blue2") +
    geom_density(adjust = bw, color = "red")+
    facet_wrap(~antigen_isotype) +
    theme_bw() +
    geom_rug()+
    ggtitle(paste("d2D28 Titer", "bandwdith =", bw, "\nn hist bins = ", nbins))
  print(p)
  
}


dev.off()


