library(dplyr)
library(ggplot2)
library(magrittr)
library(here)
library(grid)
library(feasts)
library(fable )
library(astsa)
library(fpp3)
library(forecast) #v8.23.0

## Process of producing a forecast
## Following updated forecast book - 3rd version
## https://otexts.com/fpp3/
## Another nice site to help with interpretation
## https://online.stat.psu.edu/stat501/lesson/1

## 1. Data preparation (Tidy)
## 2. Plot the data (Visualize)
## 3. Define a model (Specify)
## 4. Train the model (Estimate) using a mable
## 5. Check model performance (Evaluation)
## 6. Producing forecasts (Forecast) using fable

################
#    1. Tidy
################
## read in the data
isoLoc <- read.csv(here::here('02_Manuscript/inputFiles',
                              'supplementalTable2Metadata-ST222.csv'),
                   header = T, stringsAsFactors = F)

## summarize case counts by years with some isolates excluded due to unknown
## or incomplete meta data
isoCount <- isoLoc %>% 
  dplyr::filter(Source != 'Environmental' & Source != 'Unknown') %>%
  dplyr::select(Year, Type, State) %>%      ## select these two columns 
  na.omit %>%                               ## remove NA
  dplyr::filter(State != "Canada") %>%
  dplyr::rename(year = Year) %>%
  dplyr::group_by(year) %>%          
  dplyr::summarize(caseCount = n())         ## sum count & rename column header

## remove years which are on either side of the years w/zero cases (1992-1996)
isoCount <- isoCount %>% 
  dplyr::slice(-c(1:2))       

## remove years which are on either side of the years w/zero cases (1998-2000)
isoCount <- isoCount %>%
  dplyr::slice(-c(2:3))       

## convert year to numeric to bind future prediction rows
isoCount$year <- as.numeric(isoCount$year)

## 1992 - 1996
isoCount <- rbind(isoCount, c(1992, 2/5))
isoCount <- rbind(isoCount, c(1993, 2/5))
isoCount <- rbind(isoCount, c(1994, 2/5))
isoCount <- rbind(isoCount, c(1995, 2/5))
isoCount <- rbind(isoCount, c(1996, 2/5))

## 1998-2000
isoCount <- rbind(isoCount, c(1998, 2/3))
isoCount <- rbind(isoCount, c(1999, 2/3))
isoCount <- rbind(isoCount, c(2000, 2/3))

isoCount <- isoCount[order(isoCount$year),]

## data frame which is rates of cases per year rather than counts; for plotting
#isoCountRate222 <- isoCount

## convert to tsibble to use with feasts and fable packages. The dependent 
## variable is considered a rate because we have case counts by year and include
## an estimate rate value for years with zero isolates. 
iso222 <- as_tsibble(isoCount, index = year)

#########################
#    2. Visualization
#########################
## Scatter plot of yearly rate changes with the fitted regression line (blue)
## First, we can see there is a slight trend (upward) over the time span. Added
## the mean line of rate of case counts (6.9), we can see that the series tends
## to stay on the same side of the mean (above or below) for a while and than
## wanders to the other side, with points flipping back and forth across the
## mean

## Calculate mean of caseCount for plot: 6.9
mean222 <- mean(iso222$caseCount)

iso222 |>
  ggplot(aes(x = year, y = caseCount)) +
  labs(x = "1992 - 2020, yearly",
       y = "rate of case counts") +
  geom_point() +
  geom_line()+
  geom_smooth(method = "lm", se = FALSE) +         ## doesn't include 95% CI
  geom_hline(yintercept = mean222, color = 'red')

## auto_plot() function from feasts package will automatically produce an 
## appropriate plot of whatever you pass to the first argument. This plot
## is the same as above but without the mean line and the smoothed lm line
iso222 |>
  autoplot(caseCount)

######################################
#    3 and 4. Specify and Estimate
######################################
## Linear model using time trend: time series linear model (TSLM). 
## y variable = caseCount rate (response or dependent variable) and
## x variable = year (predictor or independent variable) 
## Assumptions of linear models 
## 1) Linearity: relationship between x and the mean of Y is linear
## 2) Homoscedasticity: variance of the residuals is the same for any value of x
## 3) Independence: observations are independent, no diagnostic test available
## 4) Normality: for any fixed value of X, Y is normally distributed

## Plot of lag 1 values to determine if an autoregressive model of order 1 (AR1)  
## also known as linear regression is a good choice...
## First make a time series out of data and than plot, in which we see a 
## moderately strong positive linear association so AR(1) model will be useful
x = ts(isoCount, frequency = 1)
astsa::lag1.plot(x, 1)

## Plot visualizations after just running lm 
## Plot 1 - residuals vs fitted values; where residuals are measured as 
## residual = observed y - modelPredicted y. This plot can be used to check the
## assumptions of linearity. 
## Linearity (1) is confirmed by finding equally spread residuals around a 
## horizontal line without distinct patterns and is a good indication you have
## a linear relationship. Homoescedascity (2) can also be checked w/this plot.
## Plot 2 - QQ-plot used to evaluate the normality assumption (4) by comparing
## residuals to the 'ideal' normal observations. Expect observations to lie
## along the 45-degree line. 
## Plot 3 - Scale-Location plot is used to show if the residuals are spread
## equally along the ranges of predictors and allows you to check the assumption
## of equal variance (homoescedascity, assumption 2). Expect to see a horizontal
## line w/equally and randomly spread points. 

## regmod plot interpretation and evaluation
## Plot 1 interpretation - mostly horizontal band of points with no distinct
## pattern though possibly 3 outliers (labeled 22, 28, and 29)
## Plot 2 interpretation - Follows a mostly straight line though again three 
## observations may be a little off (22, 28, and 29)
## Plot 3 interpretation - The red line is mostly horizontal, though it does
## show a slight positive angle 
## Note that I have NOT taken the log of the data in the above interpretation 

regmod <- lm(isoCount$caseCount ~ isoCount$year)
layout(matrix(1:4, byrow = F, ncol = 2))
plot(regmod)

## transmod plot interpretation and evaluation
## Plot 1 interpretation - better horizontal band of points with no distinct
## pattern though possibly 1 outlier (labeled as 29, corresponding to year 2020)
## Plot 2 interpretation - Follows a mostly straight line and are thus more or 
## less normally distributed, though again one observation (labeled as 29, 
## corresponding to year 2020) may be an outlier
## Plot 3 interpretation - The red line is mostly horizontal
## Note that I have taken the log of the data for this interpretation 

transmod <- lm(log(isoCount$caseCount) ~ isoCount$year)
layout(matrix(1:4, byrow = F, ncol = 2))
plot(transmod)
dev.off()

mean(transmod$residuals) 
## -7.8e-18 mean of residuals is close to zero

summary(transmod)
# Call:
#   lm(formula = log(isoCount$caseCount) ~ isoCount$year)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.37180 -0.35523 -0.05116  0.40558  1.31763 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   -289.91627   30.62340  -9.467 4.54e-10 ***
#   isoCount$year    0.14504    0.01527   9.501 4.21e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6878 on 27 degrees of freedom
# Multiple R-squared:  0.7698,	Adjusted R-squared:  0.7612 
# F-statistic: 90.27 on 1 and 27 DF,  p-value: 4.208e-10

## using fable TSLM; REQUIRES fpp3 package! We use year here instead of trend
## but gives the same answer when using trend or year as the y variable. 
## Additionally, same values in report generated from fable::TSLM as in 
## regular time series above. 
fit <- iso222 %>%
  model(time_series = fable::TSLM(log(caseCount) ~ year )) |>
  report()

## mean of innovation residuals is close to zero and matches
## transmod mean above (~ L. 164)
fitTB <- augment(fit) %>% select(.innov) %>% as_tibble() 
mean(fitTB$.innov) #  5.278943e-18

## formal test whether the mean of the innovations is significantly different
## from 0, by using one-sample t-test based on the null hypothesis of E(et) = 0
## p-value for t-test of H0: E(innov) = 0; p-value >> 0.05 so we cannot reject
## the null hypothesis that E(et) = 0

t.test(fitTB$.innov, mu = 0)$p.value ## 1 

## The y variable time of fit2, is the years in the time series. Note, that
## I do not know how to do this with the updated TSLM function from fable
fit2 <- ts(log(iso222$caseCount), start = 1992, end = 2020, frequency = 1)
tslm(fit2 ~ time(fit2))

# Call:
#   tslm(formula = fit2 ~ time(fit2))
# 
# Coefficients:
#   (Intercept)   time(fit2)  
# -289.916        0.145  

tslm(fit2 ~ trend)
# 
# Call:
#   tslm(formula = fit2 ~ trend)
# 
# Coefficients:
#   (Intercept)        trend  
# -1.141        0.145  

## This plot is the case counts by year plot with a regression line like the 
## above ggplot (~ L.89) generated with the geom_smooth function. However, this 
## plot also includes the 95% confidence interval (gray area).

tslm(fit2 ~ trend(fit2))

# Call:
#   tslm(formula = fit2 ~ trend(fit2))
# 
# Coefficients:
#   (Intercept)  trend(fit2)fit  trend(fit2)lwr  trend(fit2)upr  
# 0.7562         -1.0906          2.0906              NA 

## Compare the residual RMSE with the RMSE obtained with the time series
## cross-validation
## First calculation implements a one-step time series cross-validation forecast
## errors using a rolling forecast origin. The second calculation estimates the
## error once for the whole data set and than computes the RMSE from the 
## one-step forecasts. The RMSE from the residuals is smaller (n = 0.663) as the 
## corresponding 'forecasts' are based on a model fitted to the entire data set,
## rather than being true forecasts. 

fun1 <- function(x,h){forecast(tslm(x ~ trend ), h = h)}
#fun2 <- function(x){tslm(x ~ time(x) )}

test<-fun1(fit2, h = 5)

e <- tsCV(fit2, fun1, h = 1)

sqrt(mean(e^2, na.rm=TRUE))
#0.7537545

sqrt(mean(residuals(fun1(fit2, h = 1))^2, na.rm=TRUE))
# 0.6636662

## A good way to choose the best forecasting model is to find the model with the
## smallest RMSE computed using time series cross-validation but already are 
## using the best model 


####################
#    5. Evaluate
####################
## After fitting a regression model, one should plot residuals to check the 
## assumptions of the model have been satisfied via ACF plots using the 
## gg_tsresiduals() function from the feast package. For an ACF to make sense,
## the series must be at least a weakly stationary series. This means that the 
## autocorrelation of any lag is the same regardless of where we are in time or
## that the values of the time series fluctuate around a constant mean or with
## a constant variance or one whose properties do not depend on the time at 
## which the series is observed or means that the statistical properties of a
## process generating a time series do not change over time.

## Common to find autocorrelation in the residuals when fitting a regression
## model to time series data. If the estimated model violates the assumption of
## no autocorrelation in errors, than the forecasts may be inefficient, meaning 
## there is information left over which should be accounted for in the model.
## Forecasts from a model with autocorrelated errors are still unbiased but
## will have larger prediction intervals than needed. 

## Forecast Assumptions
## (1), we can use a Q-Q plot to see whether the innovations are normally
## distributed with a mean of zero.
## 2) we can use the sample autocorrelation function (ACF) to examine whether
## the innovations covary with a time-lagged version of themselves.

## ACF plots - 
## * Time plot - If some changing variation over time (heteroscedasticity)
## will make prediction interval coverage inaccurate. 
## * Autocorrelation plot - correlation between an observation and previous
## observation at a defined 'lag' or at some value in the past
## * Histogram of residuals - check whether residuals are normally distributed.
## Not essential for forecasting but makes calculation of prediction intervals
## much easier.

fit |> feasts::gg_tsresiduals()

## ACF plot intepretation
## * Time plot - variation of the residuals stays much the same across the 
## historical data, apart from the 2020 outlier and therefore the residual
## variance can be treated as constant (constant variance)
## * Autocorrelation plot - lack of correlation b/c no lines outside blue lines
## so this is considered a white noise series (independent or uncorrelated)
## * Histogram of residuals plot - normal residuals but one outlier (mean of 
## residuals is 0 and approximately a normal distribution)

## More evaluation by using a testing and train set
## Select case count data up to 2009
train <- iso222 |>
  filter(year <= 2009)

## train the model using the subset data
count_fit <- train |>
  model(trend_model = fable::TSLM(log(caseCount) ~ trend())) |>
  report()

# Series: caseCount 
# Model: TSLM 
# Transformation: log(caseCount) 
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.5920 -0.2911 -0.1011  0.3194  0.7380 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     -1.3807     0.2035  -6.784 4.39e-06 ***
#   year           0.1741     0.0188   9.262 7.90e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4138 on 16 degrees of freedom
# Multiple R-squared: 0.8428,	Adjusted R-squared: 0.833
# F-statistic: 85.78 on 1 and 16 DF, p-value: 7.896e-08

## forecast into future for 10 years starting at year 2010
count_fc <- count_fit |>
  forecast(h = 10) 

## Plot point forecasts are both lower (2013) and higher (2014) for years but
## the estimated point distribution from 2015 - 2019 is much better. 
count_fc |>
  autoplot(iso222) 

## look at test set, would be informative if we were testing multiple models
## but we are not, providing as just example (see below section: fitting and
## plotting multiple models for an example)
fabletools::accuracy(count_fc, iso222) |>
  arrange(.model) |>
  select(.model, .type, RMSE, MAE, MAPE, MASE, RMSSE)

## Plot of fitted values compared to actual values for caseCounts, which shows
## that the fitted values follow the actual data 'fairly' closely. 
augment(fit) |>
  ggplot(aes(x = year)) +
  geom_line(aes(y = caseCount, colour = "Data")) +
  geom_line(aes(y = .fitted, colour = "Fitted")) +
  scale_colour_manual(values=c(Data="red",Fitted="blue")) +
  guides(colour = guide_legend(title = NULL))

## same plot as above but with points instead of lines
augment(fit) |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = caseCount, colour = "Data")) +
  geom_point(aes(y = .fitted, colour = "Fitted")) +
  scale_colour_manual(values=c(Data="red",Fitted="blue")) +
  guides(colour = guide_legend(title = NULL))

## fitted values plotted against caseCount, in which we see a positive
## relationship, that verifies fitted values follow the actual data as well 
## as the previous plot 
augment(fit) |>
  ggplot(aes(x = caseCount, y = .fitted)) + 
  geom_point() + 
  geom_abline(intercept =  0, slope = 1)

## same as plot 1 from plot() function call above
augment(fit) |>
  ggplot(aes(x = .fitted, y = .resid))+
  geom_point()

#################
#    Forecast
#################
## next step is evaluation. This makes a fable with is a forecast table w/point
## forecasts (.mean) and distributions. 
## Remember a forecast is the whole 
## distribution and a point forecast is some point selected from the 
## distribution, usually the mean

## without approx_normal as FALSE, than prediction intervals are different than
## from predict

fit <- iso222 %>%
  model(time_series = fable::TSLM(log(caseCount) ~ year )) 
# |>
#   fabletools::forecast(h = "10 year", approx_normal=FALSE,
#                        point_forecast = list(.mean = median)) %>%
#   hilo(level = c (95))

fitFc <- fit |> forecast(h = "10 year" , approx_normal=FALSE,
                         point_forecast = list(.mean = median))
  
## just plot the future forecast
autoplot(fitFc)

fitFc %>%
  hilo(level = c (80, 95))
 
# A tsibble: 10 x 6 [1Y]
# Key:       .model [1]
#.model       year           caseCount .mean                    `80%`                   `95%`
#<chr>       <dbl>              <dist> <dbl>                   <hilo>                  <hilo>
# 1 time_series  2021 t(t(27, 3.2, 0.74))  24.8 [ 9.421422,  65.16941]80 [ 5.472078, 112.2039]95
# 2 time_series  2022 t(t(27, 3.4, 0.74))  28.6 [10.822214,  75.82710]80 [ 6.263024, 131.0257]95
# 3 time_series  2023 t(t(27, 3.5, 0.75))  33.1 [12.426701,  88.26022]80 [ 7.164173, 153.0928]95
# 4 time_series  2024 t(t(27, 3.6, 0.75))  38.3 [14.263923, 102.76902]80 [ 8.190369, 178.9772]95
# 5 time_series  2025 t(t(27, 3.8, 0.76))  44.3 [16.366991, 119.70511]80 [ 9.358397, 209.3534]95
# 6 time_series  2026 t(t(27, 3.9, 0.76))  51.2 [18.773657, 139.48033]80 [10.687238, 245.0171]95
# 7 time_series  2027 t(t(27, 4.1, 0.77))  59.2 [21.526955, 162.57717]80 [12.198346, 286.9070]95
# 8 time_series  2028 t(t(27, 4.2, 0.78))  68.4 [24.675931, 189.56098]80 [13.915968, 336.1314]95
# 9 time_series  2029 t(t(27, 4.4, 0.78))  79.1 [28.276473, 221.09431]80 [15.867496, 393.9984]95

## plot of historical and future 
fit |> forecast(h = "10 year") |>
  autoplot(iso222) 

## The ljung_box test statistic tests whether any of a group of autocorrelations
## of a time series are different than zero. The test is Q = 6.13 and the
## p-value of the test is 0.804, which is much larger than 0.05. Thus, we 
## fail to reject the null hypothesis of the test and conclude that the data
## values are independent.

augment(fit) |>  features(.innov, ljung_box, lag =10) 
# .model      lb_stat lb_pvalue
# <chr>         <dbl>     <dbl>
#   1 trend_model    6.13     0.804

## fitting and plotting multiple models
fitTest <- iso222 %>% 
  model(
    ets = ETS(log(caseCount)),
    arima = ARIMA(log(caseCount)),
    lm = TSLM(log(caseCount) ~ trend())
  )
fitTest

fitOut <-fitTest %>% 
  forecast(h = "10 years")  

fitOut %>%   
  hilo(level = c (80, 95))

# # A tsibble: 20 x 6 [1Y]
# # Key:       .model [2]
# .model  year       caseCount .mean                    `80%`                   `95%`
# <chr>  <dbl>          <dist> <dbl>                   <hilo>                  <hilo>
# 1 ets     2021 t(N(3.2, 0.51))  31.1 [ 9.902264,  61.88804]80 [ 6.096490, 100.5221]95
# 2 ets     2022 t(N(3.4, 0.51))  35.9 [11.446542,  71.53961]80 [ 7.047250, 116.1987]95
# 3 ets     2023 t(N(3.5, 0.51))  41.5 [13.231653,  82.69636]80 [ 8.146282, 134.3201]95
# 4 ets     2024 t(N(3.6, 0.51))  48.0 [15.295155,  95.59303]80 [ 9.416710, 155.2676]95
# 5 ets     2025 t(N(3.8, 0.51))  55.5 [17.680463, 110.50097]80 [10.885263, 179.4820]95
# 6 ets     2026 t(N(3.9, 0.51))  64.2 [20.437764, 127.73384]80 [12.582839, 207.4726]95
# 7 ets     2027 t(N(4.1, 0.51))  74.2 [23.625069, 147.65422]80 [14.545153, 239.8284]95
# 8 ets     2028 t(N(4.2, 0.51))  85.7 [27.309439, 170.68125]80 [16.813491, 277.2303]95
# 9 ets     2029 t(N(4.4, 0.51))  99.1 [31.568389, 197.29941]80 [19.435578, 320.4651]95
# 10 ets     2030 t(N(4.5, 0.51)) 115.  [36.491526, 228.06875]80 [22.466580, 370.4425]95
# 11 lm      2021 t(N(3.2, 0.54))  31.5 [ 9.647051,  63.64520]80 [ 5.854909, 104.8673]95
# 12 lm      2022 t(N(3.4, 0.55))  36.5 [11.083133,  74.04198]80 [ 6.704202, 122.4034]95
# 13 lm      2023 t(N(3.5, 0.56))  42.3 [12.728420,  86.16808]80 [ 7.672473, 142.9504]95
# 14 lm      2024 t(N(3.6, 0.56))  49.1 [14.612807, 100.31539]80 [ 8.775867, 167.0364]95
# 15 lm      2025 t(N(3.8, 0.57))  57.0 [16.770395, 116.82566]80 [10.032658, 195.2835]95
# 16 lm      2026 t(N(3.9, 0.58))  66.1 [19.240077, 136.09904]80 [11.463533, 228.4249]95
# 17 lm      2027 t(N(4.1, 0.59))  76.7 [22.066201, 158.60417]80 [13.091900, 267.3249]95
# 18 lm      2028  t(N(4.2, 0.6))  89.0 [25.299331, 184.89001]80 [14.944240, 313.0031]95
# 19 lm      2029 t(N(4.4, 0.61)) 103.  [28.997107, 215.59969]80 [17.050504, 366.6617]95
# 20 lm      2030 t(N(4.5, 0.62)) 120.  [33.225226, 251.48680]80 [19.444563, 429.7194]95

fitOut |>
  autoplot(iso222, level = 80, alpha = 0.5)

fabletools::accuracy(fitTest) |>
  arrange(.model) |>
  select(.model, .type, ME, RMSE, MAE, MAPE, MASE, RMSSE)

# # A tibble: 3 × 7
# .model .type     RMSE   MAE  MAPE  MASE RMSSE
# <chr>  <chr>    <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 arima  Training 1.79  8.22  4.42  81.0 0.895 0.851
# 2 ets    Training 1.51  7.13  3.74  71.4 0.757 0.738
# 3 lm     Training   1.50  7.13  3.74  71.4 0.757 0.738
