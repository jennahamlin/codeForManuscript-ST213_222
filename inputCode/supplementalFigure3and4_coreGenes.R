## This script takes as input two files: 1) meta data file and 2) genes 
## presence/absence csv. It will filter the data and then determine the count
## of genes (core and accessory) along with the percent of genes (core (100%),
## core-1 (99-99.9%), softcore (95-98.9%), shell (15-94.9%), cloud (14.9% or two
## genomes), singleton (single genome)). The output is parsed by binned
## years(2001-2005, 2006-2010, 2011-2015, 2016-2022) and by region (C, ENC, NE). 

library(magrittr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(grid)
library(here)

getGeneCounts <- function(dfMeta, dfGene, var, core=T){
  
  ## File prep to make gene counts------------------------ 
  ## Inputs: 
  ##        - dfMeta: name of meta data file 'supplementalTable2Metadata-ST222'
  ##        - dfGene: name of gene presence absence data file 'coreGenePresenceAbsence.csv'
  ##        - var: region bins (e.g., 'ENC') or year bins (e.g., '2001-2005') 
  ##        - core: get sum of genes (core = T) or percent of genes (core = F) 
  
  ## Read in meta data'
  dfMeta=c("supplementalTable2Metadata-ST222.csv")
  meta <- read.csv(here::here("02_Manuscript", "inputFiles", dfMeta), 
                   header = T, stringsAsFactors = F)
  
  ## Sort meta data on isolateAccession column to match up when binding to data
  meta <- meta[order(meta$isolateAccession),]

  ## Exclude Canada isolate (C180) and single/double locus variants of ST222
  ## or ST213
  meta <- meta %>%
    dplyr::filter(State != 'Canada') %>%
    dplyr::filter(elGatoSTCalls == '222' | elGatoSTCalls == '213') %>%
    dplyr::filter(isolateAccession != 'U12') %>%
    dplyr::filter(isolateAccession != 'U10')
    
   
  ## Read in presence absence csv for gene counts
  dfGene=c("coreGenePresenceAbsence.csv")
  data <- read.csv(here::here("02_Manuscript", "inputFiles/pangenome", dfGene),
                   check.names = F,  na.strings = "")
  
  ## must remove the non ST222 or ST213 isolates from core/accessory counts
  data <- data[,!names(data) %in%
                 c('C166-O_OH-2013-OC-ST1742', 'C184-S_IN-2006-CS-ST289',
                   'C185-S_NY-2006-CS-ST289', 'C198-S_CT-2005-CS-ST276',
                   'C214-S_NY-2015-CS-ST289')]

  ## Get row names which are gene names
  geneNames <- data[,1]

  ## Subset dataframe to include only pres/abs matrix
  data = data[,15:ncol(data)]

  ## Add gene names back to subset data matrix
  #row.names(data) <- geneNames

  ## Transpose data and then convert to data frame to bind meta data
  data <- as.data.frame(t(data))

  ## Get isolate names for checking data transformation
  isoNames <- data[,1]

  ## Bind data and isolate names
  data <- cbind(isoNames, data)
  
  ## Bind meta data and data to allow filtering based on year or region
  data <- cbind(meta, data)

  ## Parse data and make counts------------------------- 
  ## Inputs: 
  ##        - var: region bins (e.g., 'ENC') or year bins (e.g., '2001-2005')  
  ##        - incore: get sum genes (core = T) or percent genes (core = F)
  
  dataCount <- function (var, incore=core) {

  if(var %in% c('2001-2005', '2006-2010', '2011-2015', '2016-2020')){
    dataOut <- data %>%
    dplyr::filter(panGenomeYear == var ) %>%
    dplyr::filter(elGatoSTCalls != 'NF')
    subGroup <- dataOut$panGenomeYear[1]
    } else if (var %in% c('C', 'NE', 'NC')){
      dataOut <- data %>%
      dplyr::filter(panGenomeRegion == 'Include' ) %>%
      dplyr::filter(Region == var) %>%
      dplyr::filter(elGatoSTCalls != 'NF')
      subGroup <- dataOut$Region[1]
      } else{
        print("You are requesting data that does not exist in the dataframe")
      }
    ## Remove the binded metadata and re-transpose to have isolates on columns
    ## and gene names as rows. If error of differing rows, check this by 
    ## adjusting the number of columns
    
    # dataOut <- data %>%
    # dplyr::filter(panGenomeYear == '2001-2005' )
    # 
    
    dataOut <- dataOut[,41:ncol(dataOut)]
    dataOut <- as.data.frame(t(dataOut))

    ## Convert text to 0 (absence) or 1 (presence)
    dataOut[!is.na(dataOut)] = 1
    dataOut[is.na(dataOut)] = 0
    
    ## Make rownames (genes) as a column
    dataOut <- tibble::rownames_to_column(dataOut, 'genes')
    
    ## Determine the number of columns, which is # of isolates
    isolateCount <- ncol(dataOut[,-1])

    ## Count the occurrences of 1 (presence) for a gene row by column count
    dataOut$geneCounts <- rowSums(dataOut == "1")
    
    ## Output presence/absence for each subgroup
    ## If geneCounts == total isolate count then gene is present in all 
    ## samples and assigned 1. If not, then assigned 0
    ## remove the first row, which is the header
    dataOutTest <- dataOut[,-1]
    dataOutTest <- cbind(geneNames, dataOutTest)
    dataOutTest <- dataOutTest[,c("geneNames", "geneCounts")]
    dataOutTest <- dataOutTest %>% dplyr::mutate (Total = 
                                             dplyr::case_when(
                                               geneCounts == isolateCount ~ 1,
                                               geneCounts != isolateCount ~ 0) )
    names(dataOutTest)[names(dataOutTest) == 'Total'] <- subGroup

    write.csv(dataOutTest, here::here("02_Manuscript/inputFiles/pangenome", paste0("dataOut", subGroup,".csv")))
    
    ## Core genes are those found in all isolates
    coreOut <- dataOut %>%
      dplyr::filter(geneCounts == isolateCount) %>%
      dplyr::count()
    
    ## Accessory genes are those that in any number of isolates less then 
    ## the total number of isolates
    accOut <- dataOut %>%
      dplyr::filter(geneCounts > 0 & geneCounts < isolateCount) %>%
      dplyr::count()
    
    ## Total genes is the sum of core genes and accessory genes
    totOut <- coreOut + accOut

  if(incore == TRUE) {
    
    print("Determing the number of core and accessory genes in the sample")
    
    ## Return the variables as a vector 
    vecOut <- c(subGroup, isolateCount, coreOut, accOut, totOut)

    } else {
    
        print("Determining the breakdown of accessory genes(e.g., singleton, cloud, shell, soft core)")
      
      ## Determine percent of genes for breakdown
      dataOut <- dataOut %>%
        dplyr::mutate(percent = geneCounts/isolateCount)

      core1 <-sum(dataOut$percent >= 0.99  & dataOut$percent < 1)
      softcore <- sum(dataOut$percent >= 0.95  & dataOut$percent < 0.99)
      shell <- sum(dataOut$percent >= 0.15  & dataOut$percent < 0.95 & dataOut$geneCounts != 1)
      cloud <- sum( dataOut$geneCounts != 1 & dataOut$geneCounts != 0 & dataOut$percent < .15 )
      singleton <- sum(dataOut$geneCounts == 1 )

      coreOut <- format(round(coreOut*100/totOut))
      core1 <- format(round(core1*100/totOut))
      softcore <- format(round(softcore*100/totOut))
      shell <- format(round(shell*100/totOut))
      cloud <- format(round(cloud*100/totOut))
      singleton <- format(round(singleton*100/totOut))

      ## Return the variables as a vector 
      vecOut <- c(subGroup, singleton, cloud, shell, softcore, core1, coreOut)
               
    }
    
    return(vecOut)
    
  }

## Run dataCount function   
output <- dataCount(var)
return(output)
}

## Create dataframe for count with associated column headers
countDF <- data.frame('Group' = character(),
                     'Num.of.isolates' = character(),
                     'Num.of.core.genes' = character(),
                     'Num.of.accessory.genes' = character(),
                     'Num.of.total.genes' = character(),
                     stringsAsFactors = F)

## Create dataframe for percent with associated column headers
percentDF <- data.frame('Group' = character(),
                     'Singleton (%)' = character(),
                     'Cloud (%)' = character(),
                     'Shell (%)' = character(),
                     'Soft Core (%)' = character(),
                     'Core-1 (%)' = character(),
                     'Core (%)' = character(),
                     stringsAsFactors = F, 
                     check.names = F)

## Parse data and make counts------------------------- 
## Inputs: 
##        - varlist: region bins (e.g., 'ENC') or year bins (e.g., '2001-2005')  
##        - df: csv file of gene presence/absence
##        - core: get sum genes (core = T) or percent genes (core = F)

makeCount <- function(varlist, df, core){
  getwd()
  datalist <- list()
  if(core == T){
  for (i in varlist){
    i_out <- getGeneCounts('supplementalTable2Metadata-ST222.csv', df, i, core)
    i_out <- as.data.frame(i_out)
    colnames(i_out) <- colnames(countDF)
    datalist[[i]] <- i_out
    } 
    } else {
      for (i in varlist){
    i_out <- getGeneCounts('supplementalTable2Metadata-ST222.csv', df, i, core)
    i_out <- as.data.frame(i_out)
    colnames(i_out) <- colnames(percentDF)
    datalist[[i]] <- i_out
      }
  }
  return(datalist)
}

## Make 3 figures per year/region and combine into one figure for manuscript
## 1) Count of gene categories (year- or region percent) as a bar plot
## 2) Plot of presence/absence of core genes by region/year
## 3) Rarefaction curve for region/year

### By Year
genesYear <- c('2001-2005', '2006-2010', '2011-2015', '2016-2020')

yearCount <- makeCount(genesYear, 'coreGenePresenceAbsence.csv', core=T )
yearCount <- do.call(rbind,yearCount)
rownames(yearCount) <- NULL
yearCount

yearPercent<- makeCount(genesYear, 'coreGenePresenceAbsence.csv', core=F )
yearPercent <- do.call(rbind,yearPercent)
rownames(yearPercent) <- NULL
yearPercent

## 1) Count of gene categories
yearLong <- tidyr::gather(yearPercent, Type, Percent,
                          'Singleton (%)':'Core (%)', factor_key = TRUE)

yearLong$Percent<-as.numeric(yearLong$Percent)

## N values come from yearCount column Num.of.isolates
myXlabs <- c('2001-2005\n(n=14)', '2006-2010\n(n=20)', '2011-2015\n(n=48)',
             '2016-2022\n(n=105)')

yearSubA <- 
  ggplot(yearLong, aes(x = Group, y = Percent, group = Type, fill = Type )) +
  ggplot2::theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y.left = element_line(linewidth = 1, color = "black"),
        axis.text=element_text(color="black", face = 'bold'),
        legend.title = element_text(face = "bold"),
        plot.margin = margin(10,2,5,15)) +
  geom_bar(stat='identity') +
  ylab("Percent of Pan-Genome") +
  theme(axis.title.y = element_text( face="bold", colour = "black"),
        axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold") )+
        #,
        #text = element_text(linewidth = 14)) +
  xlab("")+
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     expand=c(0,0)) +
  scale_x_discrete(labels=myXlabs) +
  scale_fill_manual(values=c("#8756d5",
                             "#42a5f5",
                             "#388e5a",
                             "#ffee58",
                             "#ffb74d",
                             "#e53935"),
                    name="Gene Categories",
                    breaks=c('Singleton (%)', 'Cloud (%)', 'Shell (%)',
                             'Soft Core (%)', 'Core-1 (%)', 'Core (%)'),
                    labels=c('Singleton - single genome',
                             'Cloud - 2 genomes (14.9%)',
                             'Shell (15 - 94.9%)',
                             'Soft Core (95 - 98.9%)',
                             'Core-1 (99 - 99.9%)',
                             'Core (100%)'))

## 2) Presence/Absence of Core Genes
## Parsing data for the bar plot. I could not figure out how to export this 
## data with out writing to a csv. But this was also good for checking

year0105 <- read.csv(here::here("02_Manuscript/inputFiles/pangenome",
                                  'dataOut2001-2005.csv'), check.names = F,
                       na.strings = "")

year0610 <- read.csv(here::here("02_Manuscript/inputFiles/pangenome",
                                'dataOut2006-2010.csv'), check.names = F,
                     na.strings = "")

year1115 <- read.csv(here::here("02_Manuscript/inputFiles/pangenome", 
                                  'dataOut2011-2015.csv'), check.names = F,
                       na.strings = "")

year1620 <- read.csv(here::here("02_Manuscript/inputFiles/pangenome", 
                                  'dataOut2016-2020.csv'), check.names = F,
                       na.strings = "")

## Grab year counts, which is 1 if total isolates have gene presence and
## zero if all isolates do not. 
yr0105 <- year0105[,4]
yr0610 <- year0610[,4]
yr1115 <- year1115[,4]
yr1620 <- year1620[,4]

combinedyears <- cbind(yr0105, yr0610, yr1115, yr1620)
combinedyears <- as.data.frame(combinedyears)
names(combinedyears)[1] <- 'y0105' 
names(combinedyears)[2] <- 'y0610'
names(combinedyears)[3] <- 'y1115'
names(combinedyears)[4] <- 'y1620'

combinedyears <- combinedyears %>%
   dplyr::mutate(Total = rowSums(combinedyears)/ncol(combinedyears))

genes <- year0105[,2]
combinedyears <- cbind(genes, combinedyears)
combinedyears <- subset (combinedyears, Total!= 0)
combinedyears <- combinedyears[-c(6)]

combinedyears$sum <- rowSums(combinedyears[,2:5])

yearLong <- tidyr::gather(combinedyears, year, presence, y0105:y1620, factor_key = T) # line 335

yearLong %>% 
  dplyr::group_by(year) %>%
  dplyr::summarise(count = n(), count1 = sum(presence ==0))  # line 337-339

yearLong$genes <- factor(yearLong$genes, 
                         levels = unique((yearLong$genes)[order(yearLong$presence,
                                                                yearLong$sum,
                                                                decreasing = TRUE)]))
yearSubB <- ggplot2::ggplot(yearLong) +
  ggplot2::geom_tile(width = 1, aes(x = genes, y = year,
                                    fill = factor(presence)), show.legend = F) +
  scale_fill_manual(values = c('#95989D', '#09204E')) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(50,50,50,90)) +
  annotation_custom(grob = grid::textGrob(label = '2001 - 2005 \n(n = 14)',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 12)),
                    ymin = 1,
                    ymax = 1,
                    xmin = -225,
                    xmax = -225) +
  annotation_custom(grob = grid::textGrob(label = 'Presence',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 16)),
                    ymin = 4.7,
                    ymax = 4.7,
                    xmin = 950,
                    xmax = 950) +
  annotation_custom(grob = grid::rectGrob(gp = gpar(fill = '#09204E',
                                                    col = '#09204E')),
                    ymin = 4.6, 
                    ymax = 4.9, 
                    xmin = 575,
                    xmax = 650) +
  annotation_custom(grob = grid::textGrob(label = 'Absence',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 16)),
                    ymin = 4.7,
                    ymax = 4.7,
                    xmin = 2150,
                    xmax = 2150) +
  annotation_custom(grob = grid::rectGrob(gp = gpar(fill = '#95989D',
                                                    col = '#95989D')),
                    ymin = 4.6,
                    ymax = 4.9,
                    xmin = 1775,
                    xmax = 1850) +
  annotation_custom(grob = grid::textGrob(label = '2006 - 2010\n(n = 20)',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 12)),
                    ymin = 2,
                    ymax = 2,
                    xmin = -225, 
                    xmax = -225) +
  annotation_custom(grob = grid::textGrob(label = '2011 - 2015\n(n = 48)',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 12)),
                    ymin = 3,
                    ymax = 3,
                    xmin = -225,
                    xmax = -225) +
  annotation_custom(grob = grid::textGrob(label = '2016 - 2020\n(n = 104)',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 12)),
                    ymin = 4,
                    ymax = 4,
                    xmin = -225,
                    xmax = -225) +
    coord_cartesian(clip = 'off') 

# ## 3) Rarefaction curves
subgroupList = c("2001-2005", "2006-2010", "2011-2015", "2016-2020", "allIsolates")
genes= NULL

for (i in 1:length(subgroupList)){
  conserved = colMeans(read.table(here::here(
    "02_Manuscript/inputFiles/pangenome",
    paste0(subgroupList[i], "_number_of_consevered_genes.Rtab"))))
  total = colMeans(read.table(here::here(
    "02_Manuscript/inputFiles/pangenome",
    paste0(subgroupList[i],"_number_of_genes_in_pan_genome.Rtab"))))
  if (is.null(genes)){
    genes = data.frame(genes_to_genomes = c(conserved,total),
                       genomes = c(c(1:length(conserved)),
                                   c(1:length(conserved))),
                       Key = c(rep("Conserved genes",length(conserved)),
                               rep("Total genes",length(total))),
                       group = paste0(subgroupList[i]))

  } else {
    new_genes =  data.frame(genes_to_genomes = c(conserved,total),
                            genomes = c(c(1:length(conserved)),
                                        c(1:length(conserved))),
                            Key = c(rep("Conserved genes",length(conserved)),
                                    rep("Total genes",length(total))),
                            group = paste0(subgroupList[i]))
    genes = rbind(genes, new_genes)
  }
}

yearSubC <-
  ggplot(data = genes[genes$Key == 'Conserved genes',],
       aes(x = genomes, y = genes_to_genomes)) +
  geom_line(linewidth = 1, aes(color = group), show.legend = FALSE,
            position = position_jitter(w=2, h=4)) +
  geom_line(data = genes[genes$Key == 'Total genes',],
            size = 0.5, linetype = 'dashed', aes(color = group),
            show.legend = FALSE,
            position = position_jitter(w=0.5, h=0.5))+
  ylim(c(0,6000)) +
  xlim(c(1,max(genes$genomes))) +
  xlab("Number of genomes") +
  ylab("Number of genes")+
  theme_bw() +
  theme(aspect.ratio=1/2,
        axis.title.x = element_text(size=14, face="bold", colour = "black"),
        axis.title.y = element_text(size=14, face="bold", colour = "black"),
        axis.text.x = element_text(face="bold", size = 12),
        axis.text.y = element_text(face="bold", size = 12),
        plot.margin = margin(2,300,2,10)) +
  scale_color_manual(values=c("#ffa500",
                             "#9689df",
                             "#588cfc",
                             "#6aa91c",
                             "#1d6501"))+
  annotation_custom(grob = grid::textGrob(label = 'Total\nGenes',
                                          gp=gpar(fontface='bold')),
                    ymin = 4400,
                    ymax = 4400,
                    xmin = 210,
                    xmax = 210) +
  annotation_custom(grob = grid::textGrob(label = 'Core\nGenes',
                                          gp=gpar(fontface='bold')),
                    ymin = 4400,
                    ymax = 4400,
                    xmin = 230,
                    xmax = 230
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#ffa500",
                                                   lty = 'dashed', lwd = 2)),
                    ymin = 3500,
                    ymax = 3500,
                    xmin = 200,
                    xmax = 220
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#ffa500", lwd = 2)),
                    ymin = 3500,
                    ymax = 3500,
                    xmin = 225,
                    xmax = 245
                    ) +
  annotation_custom(grob = grid::textGrob(label = '2001 - 2005',
                                          gp=gpar(fontface='bold')),
                    ymin = 3550,
                    ymax = 3550,
                    xmin = 270,
                    xmax = 270
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#9689df",
                                                   lty = 'dashed', lwd = 2)),
                    ymin = 3100,
                    ymax = 3100,
                    xmin = 200,
                    xmax = 220
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#9689df", lwd = 2)),
                    ymin = 3100,
                    ymax = 3100,
                    xmin = 225,
                    xmax = 245
                    ) +
  annotation_custom(grob = grid::textGrob(label = '2006 - 2010',
                                          gp=gpar(fontface='bold')),
                    ymin = 3150,
                    ymax = 3150,
                    xmin = 270,
                    xmax = 270
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#588cfc",
                                                   lty = 'dashed', lwd = 2)),
                    ymin = 2700,
                    ymax = 2700,
                    xmin = 200,
                    xmax = 220
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#588cfc", lwd = 2)),
                    ymin = 2700,
                    ymax = 2700,
                    xmin = 225,
                    xmax = 245
                    ) +
  annotation_custom(grob = grid::textGrob(label = '2011 -2015',
                                          gp=gpar(fontface='bold')),
                    ymin = 2800,
                    ymax = 2800,
                    xmin = 270,
                    xmax = 270
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#6aa91c",
                                                   lty = 'dashed', lwd = 2)),
                    ymin = 2300,
                    ymax = 2300,
                    xmin = 200,
                    xmax = 220
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#6aa91c", lwd = 2)),
                    ymin = 2300,
                    ymax = 2300,
                    xmin = 225,
                    xmax = 245
                    ) +
  annotation_custom(grob = grid::textGrob(label = '2016 - 2020',
                                          gp=gpar(fontface='bold')),
                    ymin = 2350,
                    ymax = 2350,
                    xmin = 270,
                    xmax = 270
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#1d6501",
                                                   lty = 'dashed', lwd = 2)),
                    ymin = 1950,
                    ymax = 1950,
                    xmin = 200,
                    xmax = 220
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#1d6501", lwd = 2)),
                    ymin = 1950,
                    ymax = 1950,
                    xmin = 225,
                    xmax = 245
                    ) +
  annotation_custom(grob = grid::textGrob(label = 'All Isolates',
                                          gp=gpar(fontface='bold')),
                    ymin = 1950,
                    ymax = 1950,
                    xmin = 270,
                    xmax = 270
                    ) +
   coord_cartesian(clip = 'off')

                             
ggpubr::ggarrange(yearSubA, yearSubC,
                  nrow = 2,
                  labels = c('A', 'B'))
                  
### By Region
genesRegion <- c('C', 'NC', 'NE')

regionCount <- makeCount(genesRegion,'coreGenePresenceAbsence.csv', core=T )
regionCount <- do.call(rbind,regionCount)
rownames(regionCount) <- NULL
regionCount

regionPercent<- makeCount(genesRegion, 'coreGenePresenceAbsence.csv', core=F )
regionPercent <- do.call(rbind,regionPercent)
rownames(regionPercent) <- NULL
regionPercent

## 1) Gene category counts
## N values parsed by regions is output from makeCount function
regionLong <- tidyr::gather(regionPercent, Type, Percent, 
                          'Singleton (%)':'Core (%)', factor_key = TRUE)

regionLong$Percent<-as.numeric(regionLong$Percent)

## N values come from yearCount column Num.of.isolates
myXlabs <- c('Northeast\n(n=61)', 
             'North Central\n(n=90)',
             'Central\n(n=36)')

## make order be the same between plot a and b
regionLong$Group <- factor(regionLong$Group, levels = c("NE", "NC", C))

regionSubA <- 
  ggplot(regionLong, aes(x = Group, y = Percent, group = Type, fill = Type )) +
  ggplot2::theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y.left = element_line(linewidth = 1, color = "black"),
        axis.text=element_text(color="black", face = 'bold'),
        legend.title = element_text(face = "bold"),
        plot.margin = margin(20,0,5,15)) +
  geom_bar(stat='identity') +
  ylab("Percent of Pan-Genome") +
  theme(axis.title.y = element_text(size=16, face="bold", colour = "black"),
        axis.text.x = element_text(face="bold", size = 12),
        axis.text.y = element_text(face="bold", size = 12),
        text = element_text(size = 14)) +
  xlab("")+
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     expand=c(0,0)) +
  scale_x_discrete(labels=myXlabs) +
  scale_fill_manual(values=c("#8756d5",
                             "#42a5f5", 
                             "#388e5a",
                             "#ffee58",
                             "#ffb74d",
                             "#e53935"),
                    name="Gene Categories",
                    breaks=c('Singleton (%)', 'Cloud (%)', 'Shell (%)', 
                             'Soft Core (%)', 'Core-1 (%)', 'Core (%)'),
                    labels=c('Singleton - single genome', 
                             'Cloud - 2 genomes (14.9%)',
                             'Shell (15 - 94.9%)',
                             'Soft Core (95 - 98.9%)',
                             'Core-1 (99 - 99.9%)',
                             'Core (100%)')) 

## 2)
## Parsing data for the bar plot. I could not figure out how to export this 
## data with out writing to a csv. But this was also good for checking

regionC <- read.csv(here::here("02_Manuscript/inputFiles/pangenome", 'dataOutC.csv'),
                   check.names = F,  na.strings = "")

regionNE <- read.csv(here::here("02_Manuscript/inputFiles/pangenome", 'dataOutNE.csv'),
                   check.names = F,  na.strings = "")

regionNC <- read.csv(here::here("02_Manuscript/inputFiles/pangenome", 'dataOutNC.csv'),
                   check.names = F,  na.strings = "")

## Grab region counts, which is 1 if total isolates have gene presence and
## zero if all isolates do not. 
regC <- regionC[,4]
regNE <- regionNE[,4]
regNC <- regionNC[,4]
 
combinedRegions <- cbind(regC,regNE,regNC)
combinedRegions <- as.data.frame(combinedRegions)
names(combinedRegions)[1] <- 'C' 
names(combinedRegions)[2] <- 'NE'
names(combinedRegions)[3] <- 'NC'

combinedRegions <- combinedRegions %>%
  dplyr::mutate(Total = rowSums(combinedRegions)/ncol(combinedRegions))

genes <- regionC[,2]
combinedRegions <- cbind(genes, combinedRegions)

## Remove any genes that are not present in all regions
combinedRegions <- subset (combinedRegions, Total!= 0)

combinedRegions <- combinedRegions[-c(5)]

combinedRegions$sum <- rowSums(combinedRegions[,2:4])

regionLong <- tidyr::gather(combinedRegions, region, presence, C:NC, factor_key =T )

regionLong %>% 
  dplyr::group_by(region) %>%
  dplyr::summarise(count = n(), count1 = sum(presence == 0))  # line 337-339

regionLong$genes <- factor(regionLong$genes, 
                         levels = unique((regionLong$genes)[order(regionLong$presence,
                                                                regionLong$sum,
                                                                decreasing = TRUE)]))
regionSubB <- ggplot2::ggplot(regionLong) +
  ggplot2::geom_tile(width = 1, aes(x = genes, y = region,
                                    fill = factor(presence)), show.legend = F) +
  scale_fill_manual(values = c('#95989D', '#09204E')) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(50,50,50,90)) +
  annotation_custom(grob = grid::textGrob(label = 'Central\n(C, n = 36)',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 12)),
                    ymin = 1,
                    ymax = 1,
                    xmin = -205,
                    xmax = -205) +
  annotation_custom(grob = grid::textGrob(label = 'Presence',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 16)),
                    ymin = 3.7,
                    ymax = 3.7,
                    xmin = 1000,
                    xmax = 1000) +
  annotation_custom(grob = grid::rectGrob(gp = gpar(fill = '#09204E',
                                                    col = '#09204E')),
                    ymin = 3.65, 
                    ymax = 3.75, 
                    xmin = 650,
                    xmax = 750) +
  annotation_custom(grob = grid::textGrob(label = 'Absence',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 16)),
                    ymin = 3.65,
                    ymax = 3.75,
                    xmin = 2000,
                    xmax = 2000) +
  annotation_custom(grob = grid::rectGrob(gp = gpar(fill = '#95989D',
                                                    col = '#95989D')),
                    ymin = 3.65,
                    ymax = 3.75,
                    xmin = 1650,
                    xmax = 1750) +
  annotation_custom(grob = grid::textGrob(label = 'North Central\n(NC, n = 90)',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 12)),
                    ymin = 2,
                    ymax = 2,
                    xmin = -205, 
                    xmax = -205) +
  annotation_custom(grob = grid::textGrob(label = 'Northeast\n(NE, n = 61)',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 12)),
                    ymin = 3,
                    ymax = 3,
                    xmin = -205,
                    xmax = -205) +
    coord_cartesian(clip = 'off') 

################ Decided to not include rarefaction curve for year bins
## 3) Rarefaction curve by region. The for loop was written by ABiL and the
## input files were generated via Roary output
subgroupList = c("allIsolates", "C", "NC","NE")
genes= NULL

for (i in 1:length(subgroupList)){
  conserved = colMeans(read.table(here::here(
    "02_Manuscript/inputFiles/pangenome",
    paste0(subgroupList[i], "_number_of_consevered_genes.Rtab"))))
  total = colMeans(read.table(here::here(
    "02_Manuscript/inputFiles/pangenome",
    paste0(subgroupList[i],"_number_of_genes_in_pan_genome.Rtab"))))
  if (is.null(genes)){
    genes = data.frame(genes_to_genomes = c(conserved,total),
                       genomes = c(c(1:length(conserved)),
                                   c(1:length(conserved))),
                       Key = c(rep("Conserved genes",length(conserved)),
                               rep("Total genes",length(total))),
                       group = paste0(subgroupList[i]))

  } else {
    new_genes =  data.frame(genes_to_genomes = c(conserved,total),
                            genomes = c(c(1:length(conserved)),
                                        c(1:length(conserved))),
                            Key = c(rep("Conserved genes",length(conserved)),
                                    rep("Total genes",length(total))),
                            group = paste0(subgroupList[i]))
    genes = rbind(genes, new_genes)
  }
}

regionSubC <-
  ggplot(data = genes[genes$Key == 'Conserved genes',],
       aes(x = genomes, y = genes_to_genomes)) +
  geom_line(linewidth = 1, aes(color = group), show.legend = FALSE) +
  geom_line(data = genes[genes$Key == 'Total genes',],
            size = 1, linetype = 'dashed', aes(color = group),
            show.legend = FALSE)+
  ylim(c(0,6000)) +
  xlim(c(1,max(genes$genomes))) +
  xlab("Number of genomes") +
  ylab("Number of genes")+
  theme_bw() +
  theme(aspect.ratio=1/2,
        axis.title.x = element_text(size=16, face="bold", colour = "black"),
        axis.title.y = element_text(size=16, face="bold", colour = "black"),
        axis.text.x = element_text(face="bold", size = 12),
        axis.text.y = element_text(face="bold", size = 12),
        plot.margin = margin(2,300,2,10)) +
  scale_color_manual(values=c("#42a5f5",
                             "#388e5a",
                             "#ffb74d",
                             "#e53935")) +
  annotation_custom(grob = grid::textGrob(label = 'Total\nGenes',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 16)),
                    ymin = 4400,
                    ymax = 4400,
                    xmin = 210,
                    xmax = 210) +
  annotation_custom(grob = grid::textGrob(label = 'Core\nGenes',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 16)),
                    ymin = 4400,
                    ymax = 4400,
                    xmin = 235,
                    xmax = 235
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#42a5f5",
                                                   lty = 'dashed', lwd = 2)),
                    ymin = 3500,
                    ymax = 3500,
                    xmin = 200,
                    xmax = 220
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#42a5f5", lwd = 2)),
                    ymin = 3500,
                    ymax = 3500,
                    xmin = 225,
                    xmax = 245
                    ) +
  annotation_custom(grob = grid::textGrob(label = 'All Isolates',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 14)),
                    ymin = 3550,
                    ymax = 3550,
                    xmin = 275,
                    xmax = 275
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#e53935",
                                                   lty = 'dashed', lwd = 2)),
                    ymin = 3100,
                    ymax = 3100,
                    xmin = 200,
                    xmax = 220
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#e53935", lwd = 2)),
                    ymin = 3100,
                    ymax = 3100,
                    xmin = 225,
                    xmax = 245
                    ) +
  annotation_custom(grob = grid::textGrob(label = 'Northeast',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 14)),
                    ymin = 3150,
                    ymax = 3150,
                    xmin = 275,
                    xmax = 275
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#ffb74d",
                                                   lty = 'dashed', lwd = 2)),
                    ymin = 2700,
                    ymax = 2700,
                    xmin = 200,
                    xmax = 220
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#ffb74d", lwd = 2)),
                    ymin = 2700,
                    ymax = 2700,
                    xmin = 225,
                    xmax = 245
                    ) +
  annotation_custom(grob = grid::textGrob(label = 'North central',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 14)),
                    ymin = 2800,
                    ymax = 2800,
                    xmin = 278,
                    xmax = 278
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#388e5a",
                                                   lty = 'dashed', lwd = 2)),
                    ymin = 2300,
                    ymax = 2300,
                    xmin = 200,
                    xmax = 220
                    ) +
  annotation_custom(grob = grid::linesGrob(gp=gpar(col = "#388e5a", lwd = 2)),
                    ymin = 2300,
                    ymax = 2300,
                    xmin = 225,
                    xmax = 245
                    ) +
  annotation_custom(grob = grid::textGrob(label = 'Central',
                                          gp=gpar(fontface='bold',
                                                  fontsize = 14)),
                    ymin = 2350,
                    ymax = 2350,
                    xmin = 275,
                    xmax = 275
                    ) +
   coord_cartesian(clip = 'off')

## Combine the three subplots and export to pdf as 16 x 12 inches
ggpubr::ggarrange(regionSubA, regionSubB, regionSubC,
                  nrow = 3,
                  labels = c('A', 'B', 'C'))
