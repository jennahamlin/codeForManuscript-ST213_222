library(dplyr)
library(ggplot2)
library(here)
library(grid)

#function to add the linear regression equation and r2 value
linear = function(k) {
  z <- list(xx = format(as.numeric(coef(k)[1]), digits = 2),
            yy = format(as.numeric(abs(coef(k)[2])), digits = 2),
            r2 = format(summary(k)$adj.r.squared, digits = 3));
  if (coef(k)[2] >= 0)  {
    eq <- substitute(italic(hat(y)) == xx + yy %.% 
                       italic(x)*","~~italic(r)^2~"="~r2,z)
  } else {
    eq <- substitute(italic(hat(y)) == xx - yy %.% 
                       italic(x)*","~~italic(r)^2~"="~r2,z)  
  }
  as.character(as.expression(eq));              
}

isoBioFull <- read.csv(
  here::here('02_Manuscript/inputFiles', 
             'supplementalTable3MetaData-ST1.csv'),
  header = T, stringsAsFactors = F)

isoBioFullCount <- isoBioFull %>%
  dplyr::filter(State != "") %>%
  dplyr::rename(year = Year) %>%
  dplyr::group_by(year) %>%          
  dplyr::summarize(caseCount = n()) 

#data frame which is the original count data
isoBioFullCountOG <- isoBioFullCount

#convert year to numeric to bind future prediction rows
isoBioFullCount$year <- as.numeric(isoBioFullCount$year)

#generate years for which you want to predict the rate of case counts
new.year <- data.frame(year=c(2021:2030))

#bind data frames, and NA will be placed for future years in the case counts
isoBioFullCount <- dplyr::bind_rows(isoBioFullCount, new.year)

#look at historical data to see if there is a relationship
#plot only for years that have data (1982 - 2020)
plot(caseCount ~ year, data = isoBioFullCount[1:29,])

#dealing with the years with zero cases; we will combine years to give us rate
#of cases per year not counts. So for 1996 and 1998 with 1 and 4 cases
#respectively, the rate will be 1.667 (n cases = 5, n years = 3, 5/3)
#and be the values assigned to each year

isoBioFullCount <- isoBioFullCount %>% # remove 1996 and 1998
  dplyr::slice(-c(5:6))

isoBioFullCount <- rbind(isoBioFullCount, c(1996, 5/3))
isoBioFullCount <- rbind(isoBioFullCount, c(1997, 5/3))
isoBioFullCount <- rbind(isoBioFullCount, c(1998, 5/3))

isoBioFullCount <- isoBioFullCount[order(isoBioFullCount$year),]

#data frame which is rates of cases per year rather than counts; for plotting
isoBioFullCountRate <- isoBioFullCount

#take the log of casecount column
isoBioFullCount$caseCount <- log(isoBioFullCount$caseCount)

#rename data frame to reflect that we took the log on
isoBioFullCountLog <- isoBioFullCount

#check isoBioFullCountRate by plot; this is now rates of cases per year on the Y 
p1 <- ggplot2::ggplot(isoBioFullCountRate[1:29,], aes(x = year, y = caseCount,
                                                      group = 1))+ #add group = 1 to aes for line
  geom_line()+
  geom_point()+
  xlab("Year - in Past") + ylab("Frequency of ST1")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme(axis.text.x = element_text(angle=45, hjust=1))

lModel <- lm(caseCount ~ year, data = isoBioFullCountLog[1:29,])

st1 <- summary(lModel)
#Call:
#  lm(formula = caseCount ~ year, data = isoBioFullCountLog[1:29, 
#  ])
#
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-1.5083 -0.5140  0.1598  0.5019  1.3315 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -137.19323   32.31803  -4.245 0.000231 ***
#  year           0.06921    0.01611   4.296 0.000201 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.7259 on 27 degrees of freedom
#Multiple R-squared:  0.406,	Adjusted R-squared:  0.384 
#F-statistic: 18.46 on 1 and 27 DF,  p-value: 0.0002014

###########using tigerstats for additional check on equation of the line
#library(tigerstats)
# lmGC(caseCount ~ year, data = isoBioFullCountLog)
# Linear Regression
# 
# Correlation coefficient r =   0.6372  
# 
# Equation of Regression Line:
#   
#   caseCount = -137.1932 + 0.0692 * year 
# 
# Residual Standard Error:	s = 0.7259  
# R^2 (unadjusted):		R^2 = 0.406 

#exponentiate the predict value along with prediction intervals for 1992 - 2020
#along with the future (2021 - 2030) this is to generate the fit and prediction
#intervals to plot
isoExpHis <- exp(predict(lModel, list(year = isoBioFullCountRate$year), 
                         interval = 'prediction'))

#convert the prediction values from the regression line equation to dataframe
isoExpHis <- as.data.frame(isoExpHis)

#combine with timeValues (aka years) for plotting
isoExpHisTime <- cbind(isoBioFullCountRate$year, isoExpHis)

#rename columns to year and casecount
isoExpHisTime <- isoExpHisTime %>%
  rename(year = `isoBioFullCountRate$year`) %>%
  rename(caseCount = fit)

#estimated case count for next 10 years with intervals
#    futTime      fit      lwr      upr
#30 2021 14.671698 2.9802667  72.228008
#31 2022 15.723135 3.1601983  78.228310
#32 2023 16.849922 3.3489596  84.778528
#33 2024 18.057459 3.5468868  91.931843
#34 2025 19.351534 3.7543273  99.746731
#35 2026 20.738348 3.9716412 108.287495
#36 2027 22.224547 4.1992014 117.624859
#37 2028 23.817253 4.4373944 127.836624
#38 2029 25.524099 4.6866211 139.008382
#39 2030 27.353265 4.9472973 151.234311

pBig1 <- ggplot(isoBioFullCountRate[1:29,], aes(x = year, y = caseCount,
                                                group = 1))+ #add group = 1 to aes for line
  geom_line(lwd = 1)+
  geom_point(size = 3)+
  theme(text=element_text(family="sans"))+
  xlab("") + ylab("Frequency of ST1")+
  theme_bw()+
  ylim(0, 300)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme(axis.text.x = element_text(size = 12, angle=45, hjust=1))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(colour = "black", size = 16, face = "bold"))+
  geom_line(data = isoExpHisTime[1:29,], aes(year, caseCount),
            col = 'blue', linetype = 'dashed', lwd = 2)+
  geom_ribbon(data = isoExpHisTime[1:29,], aes(ymin = lwr, ymax = upr), alpha = 0.2)+
  geom_line(data = isoExpHisTime[30:39,], aes(year, caseCount),
            col = 'red', linetype = 'dashed', lwd = 2)+
  geom_ribbon(data = isoExpHisTime[30:39,], aes(ymin = lwr, ymax = upr), alpha = 0.2)+
  scale_x_continuous(breaks = seq(min(isoBioFullCountRate$year),
                                  max(isoBioFullCountRate$year), by = 2))

#add equation of the line and r^2 value
pBig1 <- pBig1 + annotate("text", x = 1996, y = 292,
                          label = linear(lm(caseCount ~ year,
                                            data = isoBioFullCountLog[1:29,])),
                          colour="black", size = 3.5, parse=TRUE)+
  theme(text=element_text(family="sans"))

#add time line information and dotted line
gline = linesGrob(gp = gpar(col = "black", lty = "dotted" ))

pInset1 <- ggplot(isoBioFullCountRate[1:29,], aes(x = year, y = caseCount,
                                                  group = 1))+ #add group = 1 to aes for line
  geom_line()+
  geom_point()+
  theme(text = element_text(family="sans"),
        axis.title.y = element_text(size = 8))+
  xlab("") + ylab("")+
  theme_bw()+
  ylim(0, 30)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme(axis.title.y = element_text(colour = "black", size = 15))+
  geom_line(data = isoExpHisTime[1:29,], aes(year, caseCount),
            col = 'blue', linetype = 'dashed', lwd = 1)+
  theme(axis.text.x = element_text(size = 8, angle=45, hjust=1, vjust = 1))+
  theme(axis.text.y = element_text(size = 8))

#combine big plot with inset
pBig1 <- pBig1 + ggplot2::annotation_custom(ggplotGrob(pInset1),
                                            xmin = 1990,  xmax = 2000,
                                            ymin = 78, ymax = 278)
pBig1 <- pBig1 +
  theme(plot.margin = unit(c(0.5, 1, 0.25, 1), "cm"))

pBig1 <- pBig1 + coord_cartesian(clip="off")

#don't forget to dev off
#dev.off()
 
## Check residuals using the function checkresiduals from forecast package
## Add years with no isolates collected for time series analysis
isoBioFullCountOGTest <- rbind(isoBioFullCountOG, c(1997, 0))

## Order the data by ascending year
isoBioFullCountCR <- isoBioFullCountOGTest[order(isoBioFullCountOGTest$year),]

## Add 1 to all values so we can take the log
isoBioFullCountCR <- isoBioFullCountCR%>%
   dplyr::mutate(isoBioFullCountCR[,2] + 1)

isoBioFullCountCR$caseCount <- log(isoBioFullCountCR$caseCount)

## Convert to time series data first
timeSeriesCount <- ts(isoBioFullCountCR$caseCount, start = 1992, end = 2020)
tseries::adf.test(timeSeriesCount)

#Augmented Dickey-Fuller Test
# data:  timeSeriesCount
# Dickey-Fuller = -4.1691, Lag order = 3, p-value = 0.01654
# alternative hypothesis: stationary

#This indicates that the time series is stationary due to the p-value being 
#less than 0.05 and thus rejecting the null H0 and this is what we want so we 
#can do regression analysis

###############################Dickey-Fuller Test background and stationary data
#For regression analysis to be performed, data has to be stationary. Or the 
#equation has to be rewritten in such a form that indicates a relationship 
#among stationary variables.

#What is stationary? Stationarity refers to a situation where the underlying
#stochastic process that generates the data is invariant with respect to time.
#On the other hand, if the characteristics over the time changes we call it 
#a non-stationary process.

#A time series is said to be “stationary” if it has no trend, exhibits constant 
#variance over time, and has a constant autocorrelation structure over time.

#One way to test whether a time series is stationary is to perform an augmented
#Dickey-Fuller test, which uses the following null and alternative hypotheses:

#H0: The time series is non-stationary. In other words, it has some 
#time-dependent structure and does not have constant variance over time.

#HA: The time series is stationary.

#If the p-value from the test is less than some significance level 
#(e.g. α = .05), then we can reject the null hypothesis and conclude that the
#time series is stationary.

isoTime <- forecast::tslm(timeSeriesCount ~ time(timeSeriesCount))
forecast::checkresiduals(isoTime) # supplemental figure 2

ggplot2::ggsave(here::here("02_Manuscript/outputFigures", 'supplementalFigure2CheckresidualsST1.pdf'),
                pBig1,  device = 'pdf', width = 45, height = 30, units= 'cm')

# Breusch-Godfrey test for serial correlation of order up to 6
# 
# data:  Residuals from Linear regression model
# LM test = 6.4948, df = 6, p-value = 0.3701

#See figure 1a code for more explanation. But we conclude there is no
#autocorrelation because the p-value for the Breusch-Godfrey test is 0.3701, 
#which is what we want