

umap1 = readRDS( "output/figures/main/umap_fibs.rds")


ecm= readRDS("output/figures/main/Fig2/ecm.ora.rds")
ecm =unify_axis(ecm)+
  coord_equal()+
  theme(axis.title = element_blank())

progeny= readRDS("output/figures/main/Fig2/progeny.rds")

library(cowplot)
x= grid.grabExpr(draw(progeny[[2]][[2]]))

unify_axis(ecm)
plot_grid(unify_axis(umap1)+
            theme(#axis.title = element_blank(),
                  axis.text= element_blank()),
          ecm,
          ncol = 2, align= "h", axis= "t")
