# data.csv : csv file including training data
# data_prediction1.csv : csv file including test data with Y
# data_prediction2.csv : csv file including test data without Y

# X: X of training data
# y: Y of training data
# X_prediction1: X of test data with Y
# y_prediction1: Y of test data with Y
# X_prediction2: X of test data without Y

data <- read.csv("data.csv", row.names = 1)
y <- as.matrix(data[1])
X <- as.matrix(data[c(2:ncol(data))])
data_prediction1 <- read.csv("data_prediction1.csv", row.names = 1)
y_prediction1 <- as.matrix(data_prediction1[1])
X_prediction1 <- as.matrix(data_prediction1[c(2:ncol(data_prediction1))])
data_prediction2 <- read.csv("data_prediction2.csv", row.names = 1)
X_prediction2 <- as.matrix(data_prediction2[c(1:ncol(data_prediction2))])