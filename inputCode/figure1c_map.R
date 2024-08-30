library(rnaturalearth)
library(sf)
library(dplyr)
library(ggplot2)
library(magrittr)

usa <- rnaturalearth::ne_download(type='states', returnclass = 'sf')

# Download the lake data
lakes <- rnaturalearth::ne_download(type = 'lakes', category = 'physical',
                                   returnclass = 'sf') %>% dplyr::select(name)

# Subset for only great lakes
greatLakes <- lakes %>%
  filter(name == 'Lake Huron'
         | name == 'Lake Ontario'
         | name == 'Lake Michigan'
         | name == 'Lake Erie'
         | name == 'Lake Superior')

## Add column 'count' with 0 in all states
usa$Count <- 0

## Read in the meta data with counts
isoLoc <- read.csv(here::here('02_Manuscript/inputFiles',
                              'supplementalTable2Metadata-ST213andST222.csv'),
                   header = T, stringsAsFactors = F)

isoState <- isoLoc %>%
  dplyr::group_by(State) %>%
  dplyr::summarize(stateCount = n())

## Add total counts per state for all years
## There is most definitely a better way to do this. 

usa$Count[9] <- 5 #CO
usa$Count[24] <- 4 #CT
usa$Count[34] <- 13 #IL
usa$Count[35] <- 9 #IN
usa$Count[25] <- 6 #MA
usa$Count[45] <- 6 #MD
usa$Count[49] <- 1 #ME
usa$Count[50] <- 56 #MI
usa$Count[1] <- 32 #MN
usa$Count[18] <- 1  #MO
usa$Count[46] <- 1 #NJ
usa$Count[47] <- 37 #NY
usa$Count[38] <- 26 #OH
usa$Count[48] <- 7 #PA
usa$Count[27] <- 6 #RI
usa$Count[28] <- 3 #VT
usa$Count[41] <- 8 #WI

## Center of states, to be used to place state abbreviation
## Center of states, to be used to place state abbreviation
#centroids <- data.frame(name = datasets::state.name, 
#                        Longitude = datasets::state.center$x,
#                        Latitude = datasets::state.center$y)
## remove DC for binding
usa <- usa[-44,]

## order for binding
usa <- usa[order(usa$name),]

#usa <- cbind(usa, centroids)

## remove AK and HI
usa <- usa[-2,]
usa <- usa[-10,]

## Count of: 0, 1, 3, 4, 5, 6, 7, 8, 9, 13, 27, 32, 37, 56
myColors <-  c("#eeeeee", '#EAEAFF', '#E7D4FF', '#F8BFFF', '#FFABF3', '#FF98C2',
               '#FF8A87', '#FFB177', '#E8DA85', '#C0D291', '#A6BD9D',"#86d780",
               "#00898a", "#2a4858")
               
               
               '#A8A8A8')
  
  
myColors <-  c("#eeeeee", '#fff68f', "#fafa6e", "#d1f072", "#aae479", "#86d780", "#64c987",
               "#44b98d", "#23aa8f", "#00998f", "#00898a", "#007882", "#176877",
               "#245769", "#2a4858")


myColors <-  c("#eeeeee", '#f3f299', "#f2d84f", '#b8b82a',
               "#d1f072", "#86d780", "#64c987", "#4caf50",
               "#599e30", "#23aa8f",'#5c959f', "#007882",
               "#396778", "#2a4858")


## Adjust abbreviations for: MD, DE, NJ, CT, RI, NH, VT, MS
#usa$Longitude[7]  <- -73.9#DE
#usa$Longitude[28]  <- -73.2 #NJ
#usa$Latitude[37]  <- 41 #RI
#usa$Longitude[37]  <- -69.7 #RI

countMap <- ggplot(data = usa) +
  geom_sf(data = usa, aes(fill = factor(Count)), col = 'black') +
  geom_sf(data = greatLakes, fill = 'lightblue1', alpha = 0.3) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_fill_manual(values=myColors, breaks=c('0', '1', '3', '4', '5', '6', '7',
                                              '8', '9', '13','26', '32',
                                              '37', '56')) +
  geom_segment(aes(x = -74.4, y = 38.7, xend = -75.4, yend = 38.7), col ='#76747e') + #DE
  geom_segment(aes(x = -73.7, y = 39.9, xend = -74.5, yend = 39.9), col ='#76747e') + #NJ
  geom_segment(aes(x = -71.5, y = 41.6, xend = -70.1, yend = 41), col ='#76747e') + #RI
  coord_sf(xlim = c(-125, -69), ylim = c(25, 50), expand = TRUE) +
  theme_classic() +
  theme(legend.position = "right") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank()) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))+
  theme(legend.text=element_text(size=10)) +
  theme(legend.title=element_text(size=10, face="bold")) +
  guides(fill=guide_legend(nrow=3, byrow=T, title = "ST222/213 Total Case Count from 1992-2020"))

countMap
