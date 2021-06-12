---
  title: "Analysis of Heart Failure based on Follow UP period"
author: "Ruchil Barya"
date: "06/06/2020"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


Including all libraries used in R 

```{r rlibrary}
library(WVPlots)
library(leaps)
library(dplyr)
library(ggplot2)
library(pwr)
library(readr)
library(reticulate)
library(aod)
library(broom)
library(infer)
library(pROC)
library(sigr)
library(WVPlots)
library(ranger)
library(boot)
```

# Objective
To analyze various contributing factors related to heart failure like age, red blood cells, high blood pressure, level of the creatinine phosphokinase in the blood, if the patient has diabetes or not, percentage of blood leaving the heart, number of platelets, sex, level of serum creatinine, serum sodium and if the patient is smoking.
The goal of this research is to use regression models to predict the odds of death event during follow up period, time of follow up period based on the various factors analyzed.


# Load dataset

```{r heartDataLoad}
heartFailData = read.csv("heartFailData.csv")
```

Data source:  University of California
Link: https://archive.ics.uci.edu/ml/datasets/Heart+failure+clinical+records      
File Path: https://archive.ics.uci.edu/ml/machine-learning-databases/00519/heart_failure_clinical_records_dataset.csv
Size of data set: 299 observations 
Years - 2020
Total Variables - 13
Total Observations - 299
Dependent Variable - Death event(binary), time(num)
Independent Variable considered - Age(num), anaemia(binary), high blood pressure(binary), creatinine phosphokinase(num), Diabetes(categorical), Ejection Fraction(num), Platelets(num), sex(binary), serum creatinine (num), serum sodium (num), smoking (binary)


# Data Cleaning - Outliers 

```{r Outliers}

#leverage
modeltime = lm(time~.,data=heartFailData)
summary(modeltime)
k=12
leverage = hatvalues(modeltime)
cutleverage = (2*k+2) / nrow(heartFailData)
badleverage = as.numeric(leverage > cutleverage)
table(badleverage)

#cooks distance 
cooks = cooks.distance(modeltime)
cutcooks = 4 / (nrow(heartFailData) - k - 1)
badcooks = as.numeric(cooks > cutcooks)
table(badcooks)

#mahalanobis 

datanew <- heartFailData[,c(1,3,5,7,8,9,12)]
mahal = mahalanobis(datanew, 
                    colMeans(datanew), 
                    cov(datanew))
cutmahal = qchisq(1-.001, ncol(datanew))
badmahal = as.numeric(mahal > cutmahal) ##note the direction of the > 
table(badmahal)

#total outlers 

totalout = badmahal + badleverage + badcooks
table(totalout)

heartFailureData = subset(heartFailData, totalout < 3)
#total outliers = 5 (which satisfies all the three conditions)


```

Outliers were removed from Leverage distance(0.08), Mahalanobis distance(24.32), Cooks' distance(0.01)
All the assumptionsLinearity, Additivity, Multicollinearity, Homogeneity, Homoscedasticity were met except normality, for which we used boot strapping.

# Data Screening 

```{r heartDataScreeningClean}
heartFailureData$anaemia = as.factor(heartFailureData$anaemia)
heartFailureData$high_blood_pressure = as.factor(heartFailureData$high_blood_pressure)
heartFailureData$diabetes = as.factor(heartFailureData$diabetes)
heartFailureData$sex = as.factor(heartFailureData$sex)
heartFailureData$smoking = as.factor(heartFailureData$smoking)
heartFailureData$DEATH_EVENT = as.factor(heartFailureData$DEATH_EVENT)
summary(heartFailureData)
```

# Exploratory Data analysis 

```{r heartDataEDA1, echo=FALSE}

heartFailData %>%
  ggplot(aes(x=creatinine_phosphokinase, y=DEATH_EVENT)) +
  geom_jitter(width=0, height=0.05) +
  labs(x = "Creatinine Phosphokinase", y= "death event") +
  ggtitle("Death vs CPK Levels") +
  geom_smooth(method = 'glm', method.args = list(family="binomial"))
```
death odds seem to increase with cpk levels.

```{r heartDataEDA2, echo=FALSE}
heartFailData %>%
  ggplot(aes(x=serum_creatinine, y=DEATH_EVENT)) +
  geom_jitter(width=0, height=0.05) +
  labs(x = "Serum Creatinine", y= "death event") +
  ggtitle("Death vs serum creatinine test results") +
  geom_smooth(method = 'glm', method.args = list(family="binomial"))
```
death event odds are less when serum creatinine levels are low

```{r heartDataEDA3, echo=FALSE}
heartFailData %>%
  ggplot(aes(x=serum_sodium, y=DEATH_EVENT)) +
  geom_jitter(width=0, height=0.05) +
  labs(x = "serum sodium", y= "death event") +
  ggtitle("Death vs serum sodium test results") +
  geom_smooth(method = 'glm', method.args = list(family="binomial"))
```

```{r heartDataEDA4, echo=FALSE}

heartFailData %>%
  ggplot(aes(x=time, y=DEATH_EVENT)) +
  facet_grid(~ sex) + 
  geom_jitter(width=0, height=0.05) +
  geom_smooth(method = 'glm', method.args = list(family="binomial"))
```

```{r heartDataEDA5, echo=FALSE}
heartFailData %>%
  ggplot(aes(x=age, y=DEATH_EVENT)) +
  geom_jitter(width=0, height=0.05) +
  labs(x = "age", y= "death event") +
  ggtitle("Death vs age") +
  geom_smooth(method = 'glm', method.args = list(family="binomial"))

heartFailData %>%
  ggplot(aes(x=ejection_fraction, y=DEATH_EVENT)) +
  geom_jitter(width=0, height=0.05) +
  geom_smooth(method = 'glm', method.args = list(family="binomial"))


heartFailDataEda = heartFailureData
levels(heartFailDataEda$diabetes) = c('Without Diabetes', 'With Diabetes')
levels(heartFailDataEda$smoking) = c('Doesnt smoke', 'Does Smoke')
levels(heartFailDataEda$DEATH_EVENT) = c('Alive','Passed Away')
levels(heartFailDataEda$high_blood_pressure) = c('Without high BP','With high BP')
heartFailDataEda %>%
  ggplot(aes(x=ejection_fraction, y=time)) +
  facet_wrap(~diabetes) +
  labs(x = "ejection fraction", y= "time") +
  ggtitle("Time vs Ejection Fraction") +
  geom_jitter(width=0, height=0.05) +
  geom_smooth(method = 'lm')

heartFailDataEda %>%
  ggplot(aes(x=serum_creatinine, y=time)) +
  labs(x = "serum creatinine", y= "time") +
  ggtitle("Time vs serum creatinine") +
  geom_jitter(width=0, height=0.05) +
  geom_smooth(method = 'lm')

heartFailDataEda %>%
  ggplot(aes(x=age, y=time)) +
  facet_wrap(~smoking) +
  labs(x = "age", y= "time") +
  ggtitle("Time vs Age by smoking") +
  geom_jitter(width=0, height=0.05) +
  geom_smooth(method = 'lm')

heartFailDataEda %>%
  ggplot(aes(x=age, y=time)) +
  facet_wrap(~DEATH_EVENT) +
  labs(x = "age", y= "time") +
  ggtitle("Time vs Age by smoking") +
  geom_jitter(width=0, height=0.05) +
  geom_smooth(method = 'lm')

heartFailDataEda %>%
  ggplot(aes(x=creatinine_phosphokinase, y=time)) +
  labs(x = "creatinine phosphokinase", y= "time") +
  ggtitle("Time vs cpk levels") +
  geom_jitter(width=0, height=0.05) +
  geom_smooth(method = 'lm')

heartFailDataEda %>%
  ggplot(aes(x=age, y=time)) +
  facet_grid(~DEATH_EVENT + smoking) +
  labs(x = "age", y= "time") +
  ggtitle("Time vs Age by smoking") +
  geom_jitter(width=0, height=0.05) +
  geom_smooth(method = 'lm')


heartFailDataEda %>%
  ggplot(aes(x=ejection_fraction, y=time)) +
  facet_grid(~DEATH_EVENT + smoking + high_blood_pressure) +
  labs(x = "ejection fraction", y= "time") +
  ggtitle("Time vs Ejection Fraction by death event, smoking, high blood pressure") +
  geom_jitter(width=0, height=0.05) +
  geom_smooth(method = 'lm')


heartFailDataEda %>%
  ggplot(aes(x=platelets, y=time)) +
  facet_grid(~DEATH_EVENT + smoking ) +
  labs(x = "platelets", y= "time") +
  ggtitle("Time vs platelets by death event, smoking") +
  geom_jitter(width=0, height=0.05) +
  geom_smooth(method = 'lm')

heartFailDataEda %>%
  ggplot(aes(x=DEATH_EVENT, y=ejection_fraction)) +
  geom_boxplot()

heartFailDataEda %>%
  ggplot(aes(x=DEATH_EVENT, y=age)) +
  geom_boxplot()

heartFailDataEda %>%
  ggplot(aes(x=DEATH_EVENT, y=time)) +
  geom_boxplot()

heartFailDataEda %>%
  ggplot(aes(x=DEATH_EVENT, y=serum_creatinine)) +
  geom_boxplot()

heartFailDataEda %>%
  ggplot(aes(x=DEATH_EVENT, y=serum_creatinine)) +
  facet_grid(~smoking) +
  geom_boxplot()

heartFailDataEda %>%
  ggplot(aes(x=high_blood_pressure, y=time)) +
  facet_grid(~DEATH_EVENT) +
  geom_boxplot()

```
There seems to be a difference in median time of the follow up period between people with high BP and people without high blood pressure. the difference of time is seen in when the patient was alive , patient passed away.

# Inferencial Analysis 

```{r edaAndInferences2, echo=FALSE, warning= FALSE}

heartFailDataInf = heartFailDataEda

heartFailDataInf %>%
  ggplot(aes(x=DEATH_EVENT, fill=diabetes)) +
  geom_bar()

heartFailDataInf %>%
  ggplot(aes(x=high_blood_pressure, fill=DEATH_EVENT)) +
  geom_bar()

tab = heartFailDataInf %>%
  select(high_blood_pressure, DEATH_EVENT) %>%
  table()

summary(heartFailDataInf$DEATH_EVENT)

null_spec = heartFailDataInf %>%
  specify(response = DEATH_EVENT , explanatory = high_blood_pressure, success = "Passed Away") %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1, type = "permute") %>%
  calculate(stat = "Chisq")


p_hats = heartFailDataInf %>%
  group_by(high_blood_pressure) %>%
  summarize(mean(DEATH_EVENT == "Passed Away")) %>%
  pull()

d_hat = diff(p_hats)
d_hat

null = heartFailDataInf %>%
  specify( DEATH_EVENT ~ high_blood_pressure, success = "Passed Away") %>%
  hypothesize(null = "independence") %>%
  generate(reps = 500, type = "permute") %>%
  calculate(stat="diff in props", order = c("Without high BP", "With high BP"))

ggplot(null, aes(x=stat)) +
  geom_density() +
  geom_vline(xintercept = d_hat, color="red")

null %>%
  summarize(pval = 2 * mean(stat <= d_hat))
```
The graph suggests there is no difference between two groups one with high blood pressure and one with not high blood pressure in the occurrence of death


```{r inferenceHighBPDeathEvent}
heartFailDataInf %>%
  ggplot(aes(x=time)) +
  geom_density(adjust=2)
  
aov(time ~ high_blood_pressure, heartFailDataInf) %>%
  tidy()

74442.97 / (74442.97 + 1675451.91)
```

high blood pressure seem to be statistically significant with p < 0.05 and explains 4% of variance

By plotting observed statistic against null distribution, we learned that observed difference of 0.07 is not as greater then sort of difference we would see if there was no association between these two variables. This suggests there is no difference between two groups , one with high blood pressure and one with not high blood pressure in the occurence of death during the follow up period.

# Train and Test

```{r trainAndTest}
smp_size <- floor(0.75 * nrow(heartFailureData))

set.seed(299)
train_ind <- sample(seq_len(nrow(heartFailureData)), size = smp_size)
train <- heartFailureData[train_ind, ]
test <- heartFailureData[-train_ind, ]

```

# Binary Logistic regression
```{r deathEventModel}
m=glm(DEATH_EVENT ~ time + platelets + high_blood_pressure , data = train, family="binomial")
prob = predict(m, test, type="response")
mean(as.numeric(as.character(train$DEATH_EVENT)))
test$pred = ifelse(prob > 0.31, 1, 0)
mean(test$DEATH_EVENT == test$pred)

ROC = roc(test$DEATH_EVENT, prob)
plot(ROC, col="blue")
auc(ROC)


```


#Stepwise Binary Logistic regression 

```{r}
null_model = glm(DEATH_EVENT ~ 1 , data = train, family="binomial")

full_model = glm(DEATH_EVENT ~., data = train, family="binomial")

step_model = step(null_model, scope=list(lower=null_model, upper = full_model), directoin="forward")

step_prob = predict(step_model, test, type= "response")
ROC3 = roc(test$DEATH_EVENT, step_prob)
plot(ROC3, col="red")
auc(ROC3)
```
Evaluation of AUC:

From the graph we can see that area under the curve is greater area  - 
Our curve is not close to the 45-degree diagonal of the ROC space, which suggest the model has higher accuracy. 
Which suggest we have good rate for positive test where people are correctly identified with the heart failure at the time when they visit the hospital.

# Assumption for logistic regression

## Additivity
```{r Additivity}
correl = cor(heartFailData[,c(1,3,5,7,8,9,12)])
symnum(correl)
#it does met the assumption of additivity as correlation values are less than 0.7 and there is no multicollinearity

#Additivity - we already meet the assumption of additivity (verified for logistic regression)

#new model after removing outliers(to check assumptions)
modeltime = lm(time~.,heartFailureData)
summary(modeltime)
```

## Linearity
```{r linearity}
#linearity 
standardized = rstudent(modeltime)
fitted = scale(modeltime$fitted.values)
qqnorm(standardized)
abline(0,1)
```
Yes , it met the assumption of linearity  as the dots are close to line

## Normality
```{r Normality}
hist(standardized)

library(moments)
skewness(datanew)
kurtosis(datanew)
shapiro.test(modeltime$residuals)

```
creatinine_phosphokinase and serum_creatinine has very high skew value does not met the assumption of normality as the p value is less 0.05

#Heterocadasticity
```{r Heterocadasticity}
plot(fitted, standardized)
abline(0,0)
abline(v = 0)
library(car)
vif(modeltime)
```
It meet the assumption of Homogeneity and Homoscedasticity as the variance is spread evenly across y axis from -2 to +2,also VIF is less than 5


# Linear regression 

We are having a null model with no predictors and full model using all of the predictors

```{r timeModel}
myTimeModel1 = lm(time~ejection_fraction+high_blood_pressure+DEATH_EVENT, data = train)
wrapFTest(myTimeModel1)
test$timePred1 = predict(myTimeModel1, newdata = test)
test$residuals = test$time - test$timePred1
summary(test$time)
summary(test$timePred1)

test %>%
  ggplot(aes(x=timePred1, y=time)) +
  geom_point() +
  geom_abline(color = "blue")

test %>%
  ggplot(aes(x=timePred1, y=residuals)) +
  geom_pointrange(aes(ymin=0, ymax=residuals)) +
  geom_hline(yintercept = 0, linetype = 3) +
  ggtitle("linear model prediction for time vs ejection_fraction, high bp, death event")

GainCurvePlot(test, "timePred1", "time", "time model")

rmse = sqrt(mean(test$residuals ^ 2))
print("rmse")
rmse

print("standard deviation")
sd(test$time)

print("r square")
glance(myTimeModel1)$r.squared
```

diagnol line shows the if the time is random
green is the curve where perfect model would trace out
blue is the curve our model traces out

rmse of 67.71 is better than standard deviation of 76.04

# Random Forest

```{r randomForestTimeModel}
myTimeModel2 = ranger(time~ejection_fraction+high_blood_pressure+DEATH_EVENT, train, num.trees = 500, respect.unordered.factors = "order", seed = 299)

test$timePred2 = predict(myTimeModel2, test)$predictions
test$residuals2 = test$time - test$timePred2
summary(test$time)
summary(test$timePred2)

test %>%
  ggplot(aes(x=timePred2, y=time)) +
  geom_point() +
  geom_abline(color = "blue")

test %>%
  ggplot(aes(x=timePred2, y=residuals2)) +
  geom_pointrange(aes(ymin=0, ymax=residuals2)) +
  geom_hline(yintercept = 0, linetype = 3) +
  ggtitle("linear model prediction for time vs ejection_fraction, high bp, death event")

GainCurvePlot(test, "timePred2", "time", "time model")

print("rmse")
rmse = sqrt(mean(test$residuals2 ^ 2))
rmse

print("standard deviation")
sd(test$time)

print("r square")
glance(myTimeModel1)$r.squared
```
rmse improved from 67.71 to 66.12 utilizing random forest algorithm for test data.

```{r anovamodel1}
anova(myTimeModel1, myTimeModel2)
```
To evaluate the linear regression model 
We are using train and test data - with Sample size of 0.75
Train data is for running the models and test data is for the predictions

We are having a null model with no predictors and full model using all of the predictors
We also perfomed anova test on the models
These are the results obtained

#Prediction for Linear regression 

```{r}
mymodel3 = lm( time ~ ejection_fraction + high_blood_pressure + DEATH_EVENT , data = train)
summary(mymodel3)

predicted_model = predict(mymodel3, newdata = test)
test$predictedTime = predicted_model
#head(test)

 
# time vs predicted time
ggplot(test, aes(predictedTime, time)) + 
  geom_point() + 
  geom_abline() +
  ggtitle("Prediction vs. Real values")

 
residuals_model = test$time - predicted_model
test$residuals =  residuals_model

 
# Predicted values vs residuals
ggplot(test, aes(predictedTime, residuals)) +
  geom_pointrange(aes(ymin = 0, ymax = residuals)) +
  geom_hline(yintercept = 0) +
  labs(x="Predicted Charges", y="Residuals") 

# Histogram of Residuals
ggplot(test, aes(residuals)) + 
  geom_histogram(bins = 10) +
  labs(x="Residuals", y="Frequency") 
 

# The gain curve plot measures how well the model score sorts the data compared to the true outcome value.
GainCurvePlot(test, "predictedTime", "time", "Model")

act_pred <- data.frame(cbind(actuals=test$time, predicteds=test$predictedTime)) # actuals_predicteds 
cor(act_pred) # correlation_accuracy

min_max <- mean(apply(act_pred, 1, min) / apply(act_pred, 1, max))  
print(min_max) # show the result

```


# Bootstrapping

Boot Co-ef
```{r bootstrapCoef}
bootcoef = function(formula, data, indices){
  d = data[indices, ] 
  model = lm(formula, data = d)
  return(coef(model))
}
```

```{r mymodel4}
mymodel4 = lm( time ~ . , data = heartFailureData)
summary(mymodel4)
AIC(mymodel4)
plot(mymodel4)
```


```{r modelbootstrap1}
mymodel4.boot = boot(formula = time ~ .,,
                  data = heartFailureData,
                  statistic = bootcoef,
                  R = 1000)
mymodel4.boot
plot(mymodel4.boot)
```

```{r modelbootstrap2}
mymodel4.boot2 = boot(formula = time ~ .,
                  data = heartFailureData[sample(1:nrow(heartFailureData), 50, replace=FALSE),],
                  statistic = bootcoef,
                  R = 1000)
mymodel4.boot2
plot(mymodel4.boot2)
```

From Bootstrapping:-
Standard error are much higher for the boot2 - in comparision with original model
Co-effient values are much higher for the subset of data - in comparision with original model

```{r ciEjectionfac}
length(heartFailureData$ejection_fraction)
mean(heartFailureData$ejection_fraction)
sd(heartFailureData$ejection_fraction)

error <- qt(0.95,df=length(heartFailureData$ejection_fraction)-1)*sd(heartFailureData$ejection_fraction)/sqrt(length(heartFailureData$ejection_fraction))

left <- mean(heartFailureData$ejection_fraction)-error
right <- mean(heartFailureData$ejection_fraction)+error
left
right
```
CONFIDENCE INTERVAL for ejection_fraction is between 36 and 39

```{r ciEjectionfacboot}
boot.ci(mymodel4.boot2, index = 5)
```
CONFIDENCE INTERVAL for variable 5 that is in the range which suggest that our data is consistent


# Regression Subset
```{r subseteval}
# Additional evaluations
subsets = leaps::regsubsets(time~.,data=train, nvmax=12)
summary(subsets)
names(summary(subsets))
subsets_measures = data.frame(model=1:length(summary(subsets)$cp),
                              rsq=summary(subsets)$rsq,
                              bic=summary(subsets)$bic,
                              adjr2=summary(subsets)$adjr2)
subsets_measures

which.min(summary(subsets)$bic)
coef(subsets,which.min(summary(subsets)$bic))

 
library(tidyr)
subsets_measures %>%
  gather(key = type, value=value, 2:4)%>%
  ggplot(aes(x=model,y=value))+
  geom_line()+
  geom_point()+
  facet_grid(type~.,scales='free_y')
```

Evaluation 
We are creating the - regression subset of the predictors from all the variables 
Where predictors are changing for each model to get the best one.
which suggested Best Model is Model3 with 
predictors:  were selected
also these are the values obtained for R2
.	R2 0.368 - Also some of the studies suggest that for human behavior and for physical process we generally have R2 value less than 50%. Which is true in our case.
.	BIC -0.79
.	SD - 76.05
.	RMSE - 66.12 - rmse improved from 67.71 to 66.12 utilizing random forest algorithm for test data.
We are using the same in random forest 




# CONCLUSION

Follow Up Time:

There seem to be a difference in median time of the follow up period between people with high BP and people without high BP. The difference of time is seen in when the patient was alive, patient passed away.

We fitted a random forest model to predict the follow up time with an RMSE of 67.71 is better than standard deviation of 76.04 with ejection fraction, high blood pressure, death event as statistically significant.



Death Event

We used forward step wise algorithm we found the best model with AIC of 156, with time, ejection_fraction, age, serum_creatininge, serum_sodium , diabetes being statistically significant.

The Area under curve was 0.7984 showing 79% accuracy when we predicted the values against test data. 



