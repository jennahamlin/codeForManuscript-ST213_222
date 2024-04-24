#nice summary of prediction and confidence intervals in R
#https://rpubs.com/aaronsc32/regression-confidence-prediction-intervals

#nice summary ofInterpreting Log Transformations in a Linear Model
#https://library.virginia.edu/data/articles/interpreting-log-transformations-in-a-linear-model

#must load library as summarize(caseCount = n()); does not run
library(dplyr)
library(ggplot2)
library(magrittr)
library(here)
library(grid)

#function to add the linear regression equation and r2 value
linear = function(k) {
  z <- list(xx = format(as.numeric(coef(k)[1]), digits = 3),
            yy = format(as.numeric(abs(coef(k)[2])), digits = 3),
            r2 = format(summary(k)$adj.r.squared, digits = 4));
  if (coef(k)[2] >= 0)  {
    eq <- substitute(italic(hat(y)) == xx + yy %.%
                       italic(x)*","~~italic(r)^2~"="~r2,z)
  } else {
    eq <- substitute(italic(hat(y)) == xx - yy %.%
                       italic(x)*","~~italic(r)^2~"="~r2,z)
  }
  as.character(as.expression(eq));
}

#read in the data
isoLoc <- read.csv(here::here('02_Manuscript/inputFiles',
                              'supplementalTable2Metadata-ST222.csv'),
                   header = T, stringsAsFactors = F)

#summarize case counts by years
isoCount <- isoLoc %>%
  dplyr::filter(Source != 'Environmental' & Source != 'Unknown') %>%
  dplyr::select(Year, Type, State) %>%        #select these two columns
  na.omit %>%                                 #remove NA
  dplyr::filter(State != "Canada") %>%
  dplyr::rename(year = Year) %>%
  dplyr::group_by(year) %>%
  dplyr::summarize(caseCount = n()) #count summary and rename column header

#data frame which is the original count data
isoCountOG <- isoCount

#remove years which are on either side of the years with zero cases (1992-1996)
isoCount <- isoCount %>%
  dplyr::slice(-c(1:2))

#remove years which are on either side of the years with zero cases (1998-2000)
isoCount <- isoCount %>%
  dplyr::slice(-c(2:3))

#convert year to numeric to bind future prediction rows
isoCount$year <- as.numeric(isoCount$year)

#generate years for which you want to predict the rate of case counts
new.year <- data.frame(year=c(2021:2030))

#bind data frames, and NA will be placed for future years in the case counts
isoCount <- dplyr::bind_rows(isoCount, new.year)

# #look at historical data to see if there is a relationship
# plot(caseCount ~ year, data = isoCount[1:21,])

#dealing with the years with zero cases; we will assign the calculated rate for
#the years which were removed due to zero cases. So 1992-1996 will result in a
#rate of 0.04 (i.e., n cases = 2, n years = 5; 2/5) for each year. 1998-2000
#will each have a rate of 0.667 (i.e., n cases = 2, n years = 3, 2/3).

isoCount <- rbind(isoCount, c(1992, 2/5))
isoCount <- rbind(isoCount, c(1993, 2/5))
isoCount <- rbind(isoCount, c(1994, 2/5))
isoCount <- rbind(isoCount, c(1995, 2/5))
isoCount <- rbind(isoCount, c(1996, 2/5))

isoCount <- rbind(isoCount, c(1998, 2/3))
isoCount <- rbind(isoCount, c(1999, 2/3))
isoCount <- rbind(isoCount, c(2000, 2/3))

isoCount <- isoCount[order(isoCount$year),]

#data frame which is rates of cases per year rather than counts; for plotting
isoCountRate222 <- isoCount

#take the log of casecount column
isoCount$caseCount <- log(isoCount$caseCount)

# #rename data frame to reflect that we took the log on
isoCountLog <- isoCount

# #check isoCountRate222 by plot; this is now rates of cases per year on the Y
# p1 <- ggplot(isoCountRate222[1:29,], aes(x = year, y = caseCount, group = 1)) + #add group = 1 to aes for line
#   geom_line() +
#   geom_point() +
#   xlab("Year - in Past") + ylab("Rates of Cases per Year") +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   theme(axis.text.x = element_text(angle=45, hjust=1))
#
# p1
#
# #generate a linear model to estimate slope (or regression coefficient, as year)
# #and intercept for log values
lModel <- lm(caseCount ~ year, data = isoCountLog[1:29,])

st222 <- summary(lModel)

#this tells us that the year has a significant relationship with rate of case
#count per year. By that I mean, the value of the slope is 0.14726 with a
#p-value of  1.17e-09. The intercept of the line is -294.27833

# Call:
#   lm(formula = caseCount ~ year, data = isoCountLog[1:29, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.48981 -0.34875  0.02398  0.41904  1.37421 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -294.27833   32.65852  -9.011 1.26e-09 ***
#   year           0.14726    0.01628   9.045 1.17e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7335 on 27 degrees of freedom
# Multiple R-squared:  0.7519,	Adjusted R-squared:  0.7427 
# F-statistic: 81.82 on 1 and 27 DF,  p-value: 1.168e-09
 
# ###########using tigerstats for additional check on equation of the line
# library(tigerstats)
# lmGC(caseCount ~ year, data = isoCountLog)
#
# Linear Regression
# 
# Correlation coefficient r =  0.8671 
# 
# Equation of Regression Line:
#   
#   caseCount = -294.2783 + 0.1473 * year 
# 
# Residual Standard Error:	s   = 0.7335 
# R^2 (unadjusted):		R^2 = 0.7519 

# #exponentiate the predicted value for the past and present years (1992 - 2020)
# #this is to generate the fit and prediction intervals to plot
isoExpHis <- exp(predict(lModel, list(year = isoCountRate222$year),
                          interval = 'prediction'))
# #convert the prediction values from the regression line equation to dataframe
isoExpHis <- as.data.frame(isoExpHis)

# #combine with timeValues (aka years) for plotting
isoExpHisTime222 <- cbind(isoCountRate222$year, isoExpHis)

 # #rename columns to year and casecount
isoExpHisTime222 <- isoExpHisTime222 %>%
  rename(year = `isoCountRate222$year`) %>%
  rename(caseCount = fit)

# #estimated case counts with intervals for next ten years
# #   futTime      fit       lwr       upr
# 30 2021  27.9444025  5.5818294 139.898511
# 31 2022  32.3778637  6.3985551 163.837937
# 32 2023  37.5147066  7.3302852 191.991605
# 33 2024  43.4665247  8.3926465 225.118353
# 34 2025  50.3626161  9.6033273 264.116074
# 35 2026  58.3527927 10.9823428 310.047545
# 36 2027  67.6106341 12.5523356 364.171100
# 37 2028  78.3372592 14.3389134 427.977073
# 38 2029  90.7656948 16.3710301 503.231091
 
# #specify different height and width to include the increased margin 
#  pdf(here::here('02_Figures/',
#                 "currentAndPredictedWTimelineST222.pdf"), 
#      width = 9, height = 6)
# 
#  pdf(here::here('02_Figures',
#                  "currentAndPredictedWTimelineST222TESt.pdf"),
#       width = 8, height = 5)
 
# #big plot for st222
pBig222 <- ggplot(isoCountRate222[1:29,], aes(x = year, y = caseCount,
                                              group = 1))+ #add group = 1 to aes for line
  geom_line(lwd = 1)+
  geom_point(size = 3)+
  theme(text=element_text(family="sans"))+
  xlab("") + ylab("Frequency of ST222/213 per year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme(axis.text.x = element_text(size = 12, angle=45, hjust=1, vjust = 1))+
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(colour = "black", size = 16, face="bold"))+
  scale_x_continuous(breaks = seq(min(isoExpHisTime222$year),
                                  max(isoExpHisTime222$year), by = 2))+
  geom_line(data = isoExpHisTime222[1:29,], aes(year, caseCount),
            col = 'blue', linetype = 'dashed', lwd = 2)+
  geom_ribbon(data = isoExpHisTime222[1:29,], aes(ymin = lwr, ymax = upr),
              alpha = 0.2)+
  geom_line(data = isoExpHisTime222[30:39,], aes(year, caseCount),
            col = 'red', linetype = 'dashed', lwd = 2)+
  geom_ribbon(data = isoExpHisTime222[30:39,], aes(ymin = lwr, ymax = upr),
              alpha = 0.2)

#add linear equation
pBig222 <- pBig222 + annotate("text", x = 1995.5, y = 570,
                              label = linear(lm(caseCount ~ year,
                                                data = isoCountLog[1:29,])),
                              colour="black", size = 3.5, parse=TRUE)+
  theme(text=element_text(family="sans"))


#add time line information and dotted line
gline = linesGrob(gp = gpar(col = "black", lty = "dotted" ))

first213 <- textGrob("ST213 - 1992\n OH, USA",
                     gp=gpar(fontsize = 8, col = 'gray'))

first222 <- textGrob("ST222 - 1998\n PA, USA",
                     gp=gpar(fontsize = 8, col = 'gray'))

firstEU <- textGrob("ST222 - 2010\n Germany, EU",
                    gp=gpar(fontsize = 8, col = 'gray'))

xTitle222 <- textGrob("Past and Future Frequency from 1992 to 2030",
                      gp=gpar(fontsize = 16, col = 'black', fontfamily = "sans",
                              fontface = "bold"))

#add inset with zoom in of st2222
pInset222 <- ggplot(isoCountRate222[1:29,],
                    aes(x = year, y = caseCount,
                        group = 1))+ #add group = 1 to aes for line
  geom_line()+
  geom_point()+
  theme(text=element_text(family="sans"))+
  xlab("") + ylab("")+
  theme_bw()+
  ylim(0,40)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  geom_line(data = isoExpHisTime222[1:29,], aes(year, caseCount),
            col = 'blue', linetype = 'dashed', lwd = 1) +
  theme(axis.text.x = element_text(size = 8, angle=45, hjust=1, vjust = 1))+
  theme(axis.text.y = element_text(size = 8))

#combine big plot with inset
pBig222 <- pBig222 + ggplot2::annotation_custom(ggplotGrob(pInset222),
                                                xmin = 1990,  xmax = 2000,
                                                ymin = 125, ymax = 535)

#increase inner plot margin and add text annotation of the time line
pBig222 <- pBig222 +
  theme(plot.margin = unit(c(1, 1, 0, 1), "cm")) +
  annotation_custom(first213, xmin=1992, xmax=1992, ymin = 90, ymax = 90)+
  annotation_custom(gline, xmin=1992, xmax=1992, ymin= 20, ymax= 50)+
  annotation_custom(first222, xmin=1998, xmax=1998, ymin = 90, ymax = 90)+
  annotation_custom(gline, xmin=1998, xmax=1998, ymin= 20, ymax= 50) +
  annotation_custom(firstEU, xmin=2010, xmax=2010, ymin = 90, ymax = 90)+
  annotation_custom(gline, xmin=2010, xmax=2010, ymin= 20, ymax= 50)+
  annotation_custom(xTitle222,xmin=2010, xmax=2010, ymin = 675, ymax = 675)

pBig222 <- pBig222 + coord_cartesian(clip="off")

# # # #don't forget to dev off
# # # dev.off()
# # # 
# 
# ##check residuals using the function checkresiduals from forecast package
# #add years with no isolates collected for time series analysis
isoCountCR <- isoCount %>%
  dplyr::mutate(isoCount[,2] + 1)

#convert to time series data first
timeSeriesCount <- ts(isoCountCR$caseCount, start = 1992, end = 2020)

tseries::adf.test(timeSeriesCount)

# Augmented Dickey-Fuller Test
# 
# data:  tsCountLogDiff
# Dickey-Fuller = -3.647, Lag order = 3, p-value = 0.04603
# alternative hypothesis: stationary
# Because p is less than 0.05, we can reject the null hypothesis indicating
# that our data is stationary. Note that we take the log and difference to 
# data

isoTime <- forecast::tslm(timeSeriesCount ~ time(timeSeriesCount))

# Breusch-Godfrey test for serial correlation of order up to 6
# 
# data:  Residuals from Linear regression model
# LM test = 14.318, df = 6, p-value = 0.02628

fit <- forecast::tslm(timeSeriesCount ~ trend )
forecast::checkresiduals(fit)
fitFC <- forecast::forecast(fit, h=10)

fitFC$mean <- exp(fitFC$mean)
fitFC$upper <- exp(fitFC$upper)
fitFC$lower <- exp(fitFC$lower)
plot(fitFC)


