# Echocardiogram-Survival-Analysis
**STA 592 Analysis of Time to Event Data Final Project**

In 2019 alone, heart disease was responsible for 840,768 deaths in America. For years researchers have been studying this disease, desperately trying to establish an understanding of it and finding methods to prevent or ease the effects of this disease. In order to establish a better understanding of the timeline of this disease, we have conducted an analysis of 132 heart attack patients in an attempt to better understand the factors that contribute to the length of the patient’s
recovery or eventual demise. It is our hope that we will be able to accurately model the timeline of heart attack patients from heart attack occurrence to death.

First we fit a Weibull regression with only the intercept and use the Kaplan-Meier estimator to check the goodness of fit. Next we investigated the relationship between the survival time and the other meaningful variables in our data set and estimated the months that a given heart attack patient will survive by using a log-location-scale regression model. Likelihood ratio tests, smallest AIC and BIC, and the distribution of the error term were used as the critieron of model selection.

The dataset was retrieved from UCI Machine Learning Repository: https://archive.ics.uci.edu/ml/datasets/echocardiogram
