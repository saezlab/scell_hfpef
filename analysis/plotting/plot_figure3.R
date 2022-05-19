## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-05-18
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## composition plots for fig 3
library(cowplot)


marker.hmap= readRDS("output/figures/main/Fig3/marker.comp.rds")


marker.hmap= unify_axis(marker.hmap)+
  theme(axis.text.x = element_text(angle= 0),
        axis.title= element_blank())+
  coord_equal()

p.marker
p.marker= unify_axis(p.marker)+
  theme(axis.text.x = element_text(angle= 0))

p.leftpart= plot_grid(p.marker, plot_grid(marker.hmap, NULL, NULL, NULL, ncol = 1), rel_widths = c(1,1))
p.leftpart

pdf("output/figures/main/Fig3/left.part.pdf",
    width = 8,
    height= 7)
p.leftpart
dev.off()
