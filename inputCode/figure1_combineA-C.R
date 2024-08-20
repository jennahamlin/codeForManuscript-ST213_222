library(cowplot)

pdf(here::here('02_Manuscript/outputFigures', "figure1A-C_RatesandCases.pdf"),
     width = 12, height = 10)

cowplot::plot_grid(pBig222, pBig1, countMap, nrow = 3, 
                   labels = c("A", "B", "C"), label_size = 20)
dev.off()
