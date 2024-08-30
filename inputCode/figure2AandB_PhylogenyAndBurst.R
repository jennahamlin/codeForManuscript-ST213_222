library(here)
library(dplyr)
library(ape)
library(ggplot2)
library(ggforce)
library(tidytree)
library(ggtree)

df <- data.frame(x=c(1, 2, 2, 3, 3, 4, 8),
                 y=c(2, 4, 5, 4, 7, 9, 10))

#create scatter plot with circle
ccPlot <- 
  ggplot(data = df, aes(x, y)) +
  theme_void()+
  geom_circle(aes(x0=4, y0=8, r=0.54), fill='#e41a1c', linetype = "blank",
              inherit.aes=FALSE) +
  geom_circle(aes(x0=2, y0=8, r=0.54), fill='#3c6ff8', linetype = "blank",
              inherit.aes=FALSE) +
  geom_circle(aes(x0=1.5, y0=9, r=0.27), fill='#4daf4a', linetype = "blank",
              inherit.aes=FALSE) +
  geom_circle(aes(x0=2.5, y0=9, r=0.27), linetype='dashed', color='black',
              fill='white', lwd=1, inherit.aes=FALSE) +
  geom_circle(aes(x0=3.5, y0=9, r=0.27), linetype='dashed', color='black',
              fill='white', lwd=1, inherit.aes=FALSE) +
  geom_circle(aes(x0=4.5, y0=9, r=0.27), linetype='dashed', color='black',
              fill='white', lwd=1, inherit.aes=FALSE) +
  geom_circle(aes(x0=3.5, y0=6.5, r=0.27), fill = 'white', linetype = "dashed",
              lwd = 1, inherit.aes=FALSE) +
  geom_circle(aes(x0=4.5, y0=6.5, r=0.27), fill = '#ffd966', linetype = 'blank',
              inherit.aes=FALSE) +
  geom_circle(aes(x0=5, y0=7, r=0.27), fill = '#984ea3', linetype='blank',
              inherit.aes=FALSE) +
  geom_circle(aes(x0=3.5, y0=5.5, r=0.27), linetype='dashed', color='black',
              fill='white', lwd=1, inherit.aes=FALSE) +
  geom_circle(aes(x0=5.75, y0=6.5, r=0.27), fill = '#9fc5e8', linetype='blank',
              inherit.aes=FALSE) +
  annotate("text", x = 2, y = 8, label = "ST222", size = 8, fontface = "bold") +
  annotate("text", x = 4, y = 8, label = "ST213", size = 8, fontface = "bold") +
  annotate("text", x = 1.5, y = 9, label = "ST1742", size = 4.5, fontface = "bold") +
  annotate("text", x = 2.5, y = 9, label = "ST2497", size = 4.5, fontface = "bold") +
  annotate("text", x = 3.5, y = 9, label = "ST2728", size = 4.5, fontface = "bold") +
  annotate("text", x = 4.5, y = 9, label = "ST1601", size = 4.5, fontface = "bold") +
  annotate("text", x = 3.5, y = 6.5, label = "ST2517", size = 4.5, fontface = "bold") +
  annotate("text", x = 4.5, y = 6.5, label = "ST227", size = 4.5, fontface = "bold") +
  annotate("text", x = 5, y = 7, label = "ST289", size = 4.5, fontface = "bold") +
  annotate("text", x = 3.5, y = 5.5, label = "ST2519", size = 4.5, fontface = "bold") +
  annotate("text", x = 5.75, y = 6.5, label = "ST276", size = 4.5, fontface = "bold") +
  annotate("segment", x = 1.95 , xend = 1.65, y = 8.54, yend = 8.8, 
           colour='black', lwd = 1.5) +
  annotate("segment", x = 2 , xend = 2.35, y = 8.54, yend = 8.8, 
           colour='black', lwd = 1.5) +
  annotate("segment", x = 3.95 , xend = 3.65, y = 8.54, yend = 8.8, 
           colour='black', lwd = 1.5) +
  annotate("segment", x = 4 , xend = 4.35, y = 8.54, yend = 8.8, 
           colour='black', lwd = 1.5) +
  annotate("segment", x = 3.5 , xend = 3.8, y = 6.76, yend = 7.5, 
           colour='black', lwd = 1.5) +
  annotate("segment", x = 4.2 , xend = 4.45, y = 7.5, yend = 6.79, 
           colour='black', lwd = 1.5) +
  annotate("segment", x = 4.4 , xend = 4.85, y = 7.65, yend = 7.2, 
           colour='black', lwd = 1.5) +
  annotate("segment", x = 5.22 , xend = 5.55, y = 6.89, yend = 6.65, 
           colour='black', lwd = 1.5) +
  annotate("segment", x = 3.5 , xend = 3.5, y = 5.75, yend = 6.22, 
           colour='black', lwd = 1.5) +
  annotate("segment", x = 2.52 , xend = 3.46, y = 8, yend = 8, 
           colour='black', lwd = 1.5) +
  coord_fixed()
ccPlot

################################################################################
## Figure 2b phylogeny of st222/213 respective of other sequence types

## ran gubbins and then ran raxml to build phylogeny from gubbins output
## https://github.com/nickjcroucher/gubbins/issues/272 See here for 
## tree scale questions
core <- ape::read.tree(here::here("02_Manuscript/inputFiles", 
                                 "figure2a-RAxML_bipartitions.tre"))
trePos <- core
## Rename sampleID to a consistent naming structure
sampleID <- c("AZ00033059",	"AZ00086253",	"D8604", "05LE000011",	"M2018015271",
              "LP000099",	"C214-S",	"C184-S",	"M2019018415",	"CL19-203613",
              "C284-S",	"C342-S",	"CL19-203028",	"C166-O", "C195",	
              "D7111",	"IDR1600008620",	"HUM-2014017122",	"AZ00073479",	
              "D9224",	"D8504", "C2015013836",	"U303-12",	"ISO-2-MA",	"D9064",
              "CL18-202329",	"C2008002767",	"LP000119",	"LP000271", 
              "LP000324",	"CDPH-031",	"C2006001573",	"80548A", "Reference")

phyloID <- c("C400", "C401", "C402",  "C403", "C404", "C405",	"C214",	"C184",
             "C286","C269", "C284", "C342", "C319", "C166",	"C195", "C406",	
             "C407", "C408", "C409", "C410", "C411",	"C412",	"C413",	"C414",
             "C415", "C416", "C417", "C418", "C419", "C420", "C421", "C422",	
             "C423", "Toronto-2005")

## Toronto-2005 as reference (GCA_001592705.1)
newNames <- as.data.frame(cbind(oldLabel=sampleID, newLabel=phyloID))

## Change the names 
trePos$tip.label<-newNames[[2]][match(trePos$tip.label, newNames[[1]])]
trePos$tip.label<-sapply(trePos$tip.label, function(x) parse(text=x))

included <- list('ST1' = c('C400', 'C401', 'C402'),
                'ST62' = c('C403', 'C404', 'C405'), 
                'ST289' = c("C184", 'C214'),
                'ST213' = c('C286', 'C269', 'C284'),
                'ST222' = c('C342', 'C319', 'C195', 'Toronto - 2005'),
                'ST1742' = c('C166'),
                'ST187' = c('C406', 'C407', 'C408'),
                'ST242' = c('C409', 'C410', 'C411'),
                'ST36' = c('C412', 'C413', 'C414'),
                'ST94' = c('C415', 'C416', 'C417'),
                'ST731' = c('C418', 'C419', 'C420'),
                'ST42' = c('C421', 'C422', 'C423'))

stTree <- ggtree::groupOTU(trePos, included, group_name = 'included')

treeP <- ggtree::ggtree(stTree, size = 1, color = '#797a7d', ladderize = TRUE) + 
  ggtree::geom_tiplab(ggplot2::aes(color = included), fontface = 'bold', 
                      size = 6) +
  ggplot2::theme(legend.position='none') +
  ggplot2::coord_cartesian(clip='off') +
  ggplot2::scale_color_manual(values=c('ST1' = "#00188f", 'ST62' = "#ff8c00",
                              'ST222' = "#3f83f1", 'ST213' = "#e81123",
                              'ST187' = "#ec008c", 'ST242' = "#68217a", 
                              'ST36' = "#009e49", 'ST94' = "#00b294", 
                              'ST731' = "#00bcf2", 'ST42' = "#bad80a", 
                              'ST289' = '#984ea3', 'ST1742' = '#4daf4a')) + 
  ggplot2::theme(plot.margin = ggplot2::margin(t = 0.5,
                                               r = 4,
                                               b = 0.5,
                                               l = 0,
                                               unit= 'cm')) +
  ggtree::geom_nodepoint(ggplot2::aes(
                                  subset = 
                                    !is.na(as.numeric(stTree$node.label)) &
                                  as.numeric(stTree$node.label) > 85), 
                         size = 2.5) +
  ggtree::geom_treescale(label = 'nucleotide substitutions per site', 
                         linesize = 2,
                         color = '#797a7d', fontsize = 7) 


combined <- cowplot::plot_grid(ccPlot, treeP, labels = c("A", "B"),
                               label_size = 20,  nrow = 1)
combined <- combined +
  cowplot::draw_text('ST62', y = 0.10, x = 0.67, size = 20,
                     color = '#ff8c00', fontface = 'bold') + 
  cowplot::draw_text('ST213', y = 0.26, x = 0.84, size = 20,
                    color = '#e81123', fontface = 'bold') +
  cowplot::draw_text('ST289', y = 0.31, x = 0.84, size = 20,
                    color = '#984ea3', fontface = 'bold') +
  cowplot::draw_text('ST222', y = 0.34, x = 0.84, size = 20,
                     color = '#3f83f1', fontface = 'bold') +
  cowplot::draw_text('ST1742', y = 0.38, x = 0.84, size = 20,
                     color = '#4daf4a', fontface = 'bold') +
  cowplot::draw_text('ST1', y = 0.45, x = 0.67, size = 20,
                     color = '#00188f', fontface = 'bold') +
  cowplot::draw_text('ST731', y = 0.54, x = 0.84, size = 20,
                     color = '#00bcf2', fontface = 'bold') +
  cowplot::draw_text('ST94', y = 0.615, x = 0.95, size = 20,
                     color = '#00b294', fontface = 'bold') +
  cowplot::draw_text('ST42', y = 0.69, x = 0.95, size = 20,
                     color = '#bad80a', fontface = 'bold') +
  cowplot::draw_text('ST36', y = 0.78, x = 0.95, size = 20,
                     color = '#009e49', fontface = 'bold') +
  cowplot::draw_text('ST187', y = 0.86, x = 0.95, size = 20,
                     color = '#ec008c', fontface = 'bold') +
  cowplot::draw_text('ST242', y = 0.93, x = 0.95, size = 20,
                     color = '#68217a', fontface = 'bold') +
  cowplot::draw_text('Clonal\ncomplex', y = 0.33, x = 0.96, size = 20,
                     fontface = 'bold') +
  cowplot::draw_line(x = 0.90, y = c(0.24, 0.39), color = "black", size = 2)

combined
 
ggplot2::ggsave(here::here("02_Manuscript/outputFigures", 'figure2AandB-CCandPv3.pdf'),
                combined,  device = 'pdf', width = 45, height = 30, units= 'cm')

 

