
# script to reproduce the analysis and figures from Bulgarelli et al., 2015
#
# if you use any of the following code, please cite:
#
# Davide Bulgarelli, Ruben Garrido-Oter, Philipp C. Münch, Aaron Weiman,
# Johannes Dröge, Yao Pan, Alice C. McHardy, Paul Schulze-Lefert
# Structure and Function of the Bacterial Root Microbiota in Wild
# and Domesticated Barley. Cell Host and Microbe, 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# plot figures

system("mkdir -p ./figures")
system("rm -f ./figures/*pdf")

source("figure_2_and_s2.R")
source("figure_4.R")
source("figure_5.R")
source("figure_6A1.R")
source("figure_6A2.R")
source("figure_6B.R")

figure6_A1 <- create_figure_6A1()
figure6_A2 <- create_figure_6A2()
figure6_B <- create_figure_6B()

pdf("./figures/figure_6.pdf")
  print(figure6_B)
  print(figure6_A1)
  print(figure6_A2)
dev.off()

