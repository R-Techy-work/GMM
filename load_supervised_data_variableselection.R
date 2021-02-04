# data.csv : csv file including training data

# X: X of training data
# y: Y of training data

data = read.csv("data.csv", row.names = 1)
y = as.matrix(data[1])
X = as.matrix(data[c(2:ncol(data))])
