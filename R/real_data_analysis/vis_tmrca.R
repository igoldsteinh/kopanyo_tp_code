# visualize median tmrca
library(tidyverse)

snp5_tmrca <- read_csv(here::here("data", "snp5_median_tmrca.csv"))

snp10_tmrca <- read_csv(here::here("data", "snp10_median_tmrca.csv"))


snp5_hist <- snp5_tmrca %>% 
             ggplot() +
             aes(x = median_tmrca) +
             geom_histogram(binwidth = 0.5) +
             theme_bw() +
             ggtitle("Median TMRCA (SNP 5 Analysis)") +
             xlab("Median TMRCA (Years)") +
             ylab("Count") +
             scale_x_continuous(breaks = c(0, 5, 10, 15, 20),
                                labels = c("0", "5", "10", "15", "20"))

snp10_hist <- snp10_tmrca %>% 
  ggplot() +
  aes(x = median_tmrca) +
  geom_histogram(binwidth = 0.5) +
  theme_bw() +
  ggtitle("Median TMRCA (SNP 10 Analysis)") +
  xlab("Median TMRCA (Years)") +
  ylab("Count") +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20))
