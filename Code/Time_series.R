# Script for the Example
# You can copy and paste this script into RStudio to reproduce the example. The dataset was provided with the video.
library() # zoo # handle irregular time series , forecast # , lubridate # preprocessing packages for date and times, 
# tseries # modelling, formatting, tatistical test (stationary)
### Example of Time Series Analysis Workflow

# Data import, 2 alternatives

myts <- scan()

# Conversion to a time series
mycounts <- ts(myts2$X2, start = 1,
               frequency = 12)

# Alternative with myts vector
mycounts_check <- ts(myts,start = 1,
                     frequency = 12)

# Visualization
plot(mycounts, ylab = "Customer Counts",
     xlab = "Weeks")

library(forecast)
monthplot(mycounts, labels = 1:12,
          xlab = "Bidaily Units")
seasonplot(mycounts, season.labels = F,
           xlab = "")

# Model forecast
plot(forecast(auto.arima(mycounts)))

