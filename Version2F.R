View(He_Kelly_Manela_Factors_And_Test_Assets_monthly_65272913)
data<-He_Kelly_Manela_Factors_And_Test_Assets_monthly_65272913
head(data)
class=data[,3]
table(class)
X=data[,3]
head=(X)


dat = as.matrix(data[, -ncol(data)])

dat = center_scale(dat)

gmm = GMM(dat, 2, "maha_dist", "random_subset", 10, 10)