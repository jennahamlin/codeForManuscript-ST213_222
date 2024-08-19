##This code randomly selects 3 isolates associated with sequence types of 
#interest from the Bionumerics meta data file. 
library(dplyr)
library(stringr)

isoBioFull <- read.csv(
  here::here('01_Development/04_InputFiles', 
             'P01_st1-metadata_Bionumerics_v01_2021-03-17_JH.csv'),
  header = T)

isoBioFull$ST <-  stringr::str_replace(isoBioFull$ST, '\\*', '')

#need to remove 'ST1 Study - new accessions' as these isolates are duplicates
isoBioFullCount <- isoBioFull %>%
  dplyr::filter(Source != 'Environmental' & Source != 'Environmental ') %>%
  dplyr::filter(Source != 'EN' & Source != 'OE') %>%
  dplyr::filter(Comments != 'ST1 Study - new accessions' & Comments != 'EU') %>%
  dplyr::distinct(., Sample_ID, .keep_all = T) %>%
  dplyr::group_by(ST) %>%          
  dplyr::summarize(caseCount = n())

#order the grouped by data, to evalute and select sequence types
isoBioFullCount <- isoBioFullCount[order(isoBioFullCount$caseCount,
                                        decreasing=TRUE), ]

print(isoBioFullCount, n = 25)
# A tibble: 15 x 2
# ST      caseCount
# <chr>       <int>
# ""           1123
# "1"           179
# "Novel"       141
# "1*"          129
# "222*"        112
# "213*"         69
# "187*"         58
# "731*"         58
# "242*"         54
# "36*"          48
# "62*"          40
# "42*"          38
# "94*"          36
# "224*"         30
# "?"            25


#Based on the above, I will randomly pick 3 isolates from these sequences types:
#1/1*, 222*/213*, 187/187*, 731/731*, 242/242*, 36/36*, 62/62*, 42/42*, 94/94*, 
#and 242/242* 731/731*

#parse down full bionumeric metadata file for selection of random isolates
isoBioFullRand <- isoBioFull %>%
  dplyr::filter(Source != 'Environmental' & Source != 'Environmental ') %>%
  dplyr::filter(Source != 'EN' & Source != 'OE') %>%
  dplyr::filter(Comments != 'ST1 Study - new accessions' & Comments != 'EU')%>%
  dplyr::filter(Source != "" & Source != 'Unknown')%>%
  dplyr::filter(ST != "" & Year != "")%>%
  dplyr::distinct(., Sample_ID, .keep_all = T)

#Function for selecting 3 random isolates using sample, in which the selection 
#of sequence types can either be 1 or 1*. This reflects the mechanism in which 
#the sequence types have been verified (e.g. sanger or pipeline). The input
#parameters are dataframe (df) and the sequence type (st) of interest. 
randomSelect <- function(df, st){
  stWild<-  paste(st, "*", sep="")
  df[ sample( which( df$ST == st | df$ST == stWild) , 3 ) , ]
  
}

#test to confirm function works and pulls both 1 and 1*
randomSelect(isoBioFullRand, "1")

# Key Alt.Identifier Year Outbreak Sample_ID   Source    SRR_ID ST State Comments
# 333   WGMLST_LEG0002926   MA-APHL-2017 2010             ISO-59 Clinical           1*    MA         
# 104 WGLEGIONELLA0001873      ST1-Study 1994              D3933       CS   Central  1    OH       US
# 111 WGLEGIONELLA0001877      ST1-Study 1992              D3375       CS Southeast  1    VA       US

#create vector for the sequence types of interest - identified above
x <- c("1", "36", "42", "62", "94",  "187", "222", "213", "242", '731')



#list apply over the randomSelect function 
randSelectOut <- lapply(x, function(x) randomSelect(isoBioFullRand, x))

#convert list of dataframes to one dataframe
dataOut <- dplyr::bind_rows(randSelectOut, .id = 'Key')

#write to csv file
write.csv(dataOut, here::here('01_Development/04_InputFiles', 
                              'random3ForFLASH.csv'), row.names = F)
