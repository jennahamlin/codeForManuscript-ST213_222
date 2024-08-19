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
  xlab("Year - in Past") + ylab("Rates of Cases per Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme(axis.text.x = element_text(angle=45, hjust=1))

lModel <- lm(caseCount ~ year, data = isoBioFullCountLog[1:29,])

st1 <- summary(lModel)
# 
# Call:
#   lm(formula = caseCount ~ year, data = isoBioFullCountLog[1:29, 
#   ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.5327 -0.5119  0.1447  0.5422  1.3398 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -140.30787   32.74696  -4.285 0.000208 ***
#   year           0.07077    0.01632   4.335 0.000181 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7355 on 27 degrees of freedom
# Multiple R-squared:  0.4104,	Adjusted R-squared:  0.3886 
# F-statistic:  18.8 on 1 and 27 DF,  p-value: 0.0001814

###########using tigerstats for additional check on equation of the line
#library(tigerstats)
# lmGC(caseCount ~ year, data = isoBioFullCountLog)
# Linear Regression
# 
# Correlation coefficient r =  0.6406 
# 
# Equation of Regression Line:
#   
#   caseCount = -140.3079 + 0.0708 * year 
# 
# Residual Standard Error:	s   = 0.7355 
# R^2 (unadjusted):		R^2 = 0.4104 

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
# 2021 14.878836 2.9252698  75.678409
# 2022 16.005782 3.1138454  82.272894
# 2023 17.218085 3.3126693  89.493528
# 2024 18.522210 3.5222005  97.402821
# 2025 19.925111 3.7429165 106.069705
# 2026 21.434270 3.9753147 115.570203
# 2027 23.057735 4.2199130 125.988177
# 2028 24.804164 4.4772508 137.416145
# 2029 26.682870 4.7478903 149.956191
# 2030 28.703872 5.0324174 163.720973

pBig1 <- ggplot(isoBioFullCountRate[1:29,], aes(x = year, y = caseCount,
                                                group = 1))+ #add group = 1 to aes for line
  geom_line(lwd = 1)+
  geom_point(size = 3)+
  theme(text=element_text(family="sans"))+
  xlab("") + ylab("Frequency of ST1 per year")+
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

#firstUSA <- grid::textGrob("ST1 - 1982 \nCA & IN, USA",
#                      gp=gpar(fontsize=8, col = 'gray'))

#xTitle1 <- textGrob("Past and Future Rate: 1992-2030",
#                    gp=gpar(fontsize=16, col = 'black', fontfamily = "sans",
#                            fontface = "bold"))

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
#+
#annotation_custom(firstUSA, xmin=1982.5, xmax=1982.5, ymin = 70, ymax = 70)+
#annotation_custom(gline, xmin=1982.5, xmax=1982.5, ymin = 10, ymax= 40)+
# annotation_custom(xTitle1,xmin=2005, xmax=2005, ymin = 335, ymax = 335)

pBig1 <- pBig1 + coord_cartesian(clip="off")

#don't forget to dev off
#dev.off()
# 
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
# Dickey-Fuller = -4.2888, Lag order = 3, p-value = 0.01227
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
forecast::checkresiduals(isoTime)

# Breusch-Godfrey test for serial correlation of order up to 6
# 
# data:  Residuals from Linear regression model
# LM test = 6.4948, df = 6, p-value = 0.3701

#See figure 1a code for more explanation. But we conclude there is no
#autocorrelation because the p-value for the Breusch-Godfrey test is 0.3701, 
#which is what we want