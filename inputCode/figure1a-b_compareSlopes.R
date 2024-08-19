###compare regression coefficients across two groups
##this uses data from st222Emergence and ST1Emergence R scripts

require(emmeans)

st222 <- isoCountLog[1:29,]
#add years 1982 - 1993
time <- seq(1982, 1993, by = 1)
count <- rep(0, times = 12)
timeCount <- cbind(time, count)
timeCount <- as.data.frame(timeCount)

timeCount <- timeCount %>%
  dplyr::rename(year = time,
                caseCount = count)

st222Com <- rbind(timeCount, st222)
st222Com <- st222Com %>%
  dplyr::mutate("ST" = "ST222")

st1 <- isoBioFullCountLog[1:37,]
st1 <- st1 %>%
  dplyr::mutate("ST" = "ST1")

st222_1 <- rbind(st1, st222Com)

model <- lm(caseCount ~ year*ST, st222_1)

st222_2 <-anova(model)
# Analysis of Variance Table
# s
# Response: caseCount
# Df Sum Sq Mean Sq F value    Pr(>F)    
# year       1 38.128  38.128 71.5933 2.333e-12 ***
#   ST         1  9.084   9.084 17.0564 9.769e-05 ***
#   year:ST    1  5.069   5.069  9.5182  0.002899 ** 
#   Residuals 71 37.812   0.533                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

slopes <- emmeans::emtrends(model, 'ST', var = 'year') #gets each slope
pairs(slopes)
# contrast    estimate    SE df t.ratio p.value
# ST1 - ST222  -0.0462 0.015 71  -3.085  0.0029

model <- lm(caseCount~year*ST,data=st222_1)
summary(model)

# Call:
# lm(formula = caseCount ~ year * ST, data = st222_1)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -1.79734 -0.46067  0.04531  0.59502  1.48106 
# 
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -74.54025   21.65610  -3.442 0.000972 ***
# year           0.03805    0.01082   3.517 0.000765 ***
# STST222      -93.08311   29.94580  -3.108 0.002706 ** 
# year:STST222   0.04616    0.01496   3.085 0.002899 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7298 on 71 degrees of freedom
#   (3 observations deleted due to missingness)
# Multiple R-squared:  0.5803,	Adjusted R-squared:  0.5626 
# F-statistic: 32.72 on 3 and 71 DF,  p-value: 2.151e-13

anova(lm(caseCount~year,data=st222_1),lm(caseCount~year*ST,data=st222_1))
# Analysis of Variance Table
# 
# Model 1: caseCount ~ year
# Model 2: caseCount ~ year * ST
#   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
# 1     73 51.965                                  
# 2     71 37.812  2    14.153 13.287 1.254e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1