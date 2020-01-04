usps.train <- read.table("zip.train.txt", header=F, sep="")
usps.test <- read.table("zip.test.txt", header=F, sep="")



# find features in every digit when k=100
# k <- 150
success_rate <- list()
for (k in seq(5, 100, 5)){ # test k from 5 to 100 
digits_feature <- list()
for (i in c(0:9)){
  temp<-subset(usps.train,usps.train[,1]==i)
  A <- as.matrix(temp[,2:257])
  A <- t(A)
  A.k <- svds(A,k) #k=100
  digits_feature <- c(digits_feature,list(A.k$u)) # 先把 matrix convert to list , then list in list
}



# 解 list to matrix 
# output <- matrix(unlist(digits_feature[1]), ncol = 100, byrow = FALSE)
success = 0
for (j in c(1:2006)){
  compare_set <- list()
  for (d in c(1:10)){
    z<-t(as.matrix(usps.test[j,2:257])) 
    # residual for digit 0~9
    digit.matrix <- matrix(unlist(digits_feature[d]), ncol = k, byrow = FALSE)
    res.digit <- norm((diag(256)-(digit.matrix%*%t(digit.matrix)))%*%z,type="2")/norm(z,type="2") 
    compare_set <- c(compare_set,  res.digit )
  }
  min.index <- which.min(compare_set) # 找最小的 residual 的 index
  guess <- min.index-1
  if (guess == usps.test[j,1]) { success = success +1 }  
  
}
success_r = success/2006
success_rate <- c(success_rate,  success_r )

}

plot(result.csv$k,result.csv$success_r,type='l',col="green",xlab = "k value",
     +      ylab = "success rate", main = "Success rate under different k")
