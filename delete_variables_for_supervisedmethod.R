# Delete variables whose numbers are 'Var0Variable' from X, X_prediction1 and X_prediction2
if (length(Var0Variable) != 0) {
  X = X[,-Var0Variable]
  X_prediction1 = X_prediction1[,-Var0Variable]
  X_prediction2 = X_prediction2[,-Var0Variable]
}
