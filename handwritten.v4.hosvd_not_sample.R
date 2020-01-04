library(rTensor)
library(matrixcalc) # for frobenius.prod
usps_train <- read.table("zip.train.txt", header=F, sep="")
usps_test <- read.table("zip.test.txt", header=F, sep="")

### function ###
rotate <- function(x){x[,nrow(x):1]}

# extract feature of every digit from training set into a list
digits_u1 <- list()
digits_u2 <- list()
digits_z <- list()
digits_dim <- list()

for (i in c(0:9)){
  temp <- subset(usps_train,usps_train[,1]==i) 
  temp_t <-t(temp[,2:257]) # 256x658
  tnsr <- as.tensor(array(temp_t, dim = c(16,16,nrow(temp)))) # 利用 array 將 dataset 摺疊成 16x16xnrow
  hosv_decomp <- hosvd(tnsr,ranks=c(16,16,50)) # truncated k=50
  digits_u1 <- c(digits_u1,list(hosv_decomp$U[[1]]))
  digits_u2 <- c(digits_u2,list(hosv_decomp$U[[2]]))
  core_unfold <- k_unfold(hosv_decomp$Z, 3)
  digits_z <- c(digits_z, list(core_unfold@data)) # 16x16x658 => 658x256
  digits_dim <- c(digits_dim,list(hosv_decomp$Z@modes))
}


start_time <- Sys.time()
success = 0
aj = 30 # aj is the number of Aj for linear combination
f = 5 # f is the number of digits to be tesetd in usps_test

sample_list =c(1:nrow(usps_test))
draw_list <- sample(sample_list,f,replace = F)

pb <- txtProgressBar(min = 0, max = f, style = 3)
for (t in c(1:f)){
  setTxtProgressBar(pb, t)
  # test_index <- draw_list[t]
  test_index <- draw_list[t]
  Z <- rotate(array(unlist(usps_test[test_index,2:257]),dim=c(16,16))) # Z is the digit to be tested by algorithm
  compare_set <- list()
  for (n in c(1:10)){
    digit.u1 <- array(unlist(digits_u1[n]),dim=c(16,16))  
    digit.u2 <- array(unlist(digits_u2[n]),dim=c(16,16))  
    digit.z.m <- matrix(unlist(digits_z[n]), ncol = 256, byrow = FALSE)
    digit.z <- k_fold(digit.z.m,m=3,modes=c(digits_dim[[n]]))
    
    # A.j.=S(:,:,j)x1U(1)x2U(2), S is core tensor
    # A.1 = S(:,:,1)x1U(1)x2U(2), A.2 = S(:,:,2)x1U(1)x2U(2)
    # digit.z is core tensor
    A.j <- function(j){ 
      rotate(
        ttl (
          tnsr=digit.z[,,j],
          list_mat = list(digit.u1,digit.u2),
          ms=c(1,2)
        )
      )
    }
    
    # z.j=<Z,Aj>/<Aj,Aj>
    z.j <- function(j){ 
      frobenius.prod(Z,A.j(j)@data)/frobenius.prod(A.j(j)@data,A.j(j)@data)
    } 
    
    P1=0
    # summation of zj<Z,Aj>, j=1,2,3...,aj
    for (j in c(1:aj)){
      P1 <- P1+ z.j(j)*frobenius.prod(Z,A.j(j)@data)
    }
    
    P2=0
    for (j in c(1:aj)){
      P2 <- P2+ (z.j(j))^2*frobenius.prod(A.j(j)@data,A.j(j)@data)
    }
    
    G.z <- 0.5*(frobenius.prod(Z,Z))-P1+0.5*P2
    
    compare_set <- c(compare_set,G.z) 
  }
  min.index <- which.min(compare_set) # 找最小的 G.z 的 index
  guess <- min.index-1
  if (guess == usps_test[test_index,1]) { success = success +1 }  
}
close(pb)
success
success_r = success/f
success_r
end_time <- Sys.time()
dur <- end_time - start_time
dur
