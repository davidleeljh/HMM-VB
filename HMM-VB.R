library(mvtnorm)

#Forward probabiliy calculation of all states in all blocks
#x is the imput vector, mu and sigma are lists, pi is prior, a is transition matrix
forw = function(x, mu, sigma, pi, a){
  alpha = vector("list", Ti)
  alpha[[1]] = rep(0, S[1])
    for(k in 1:S[1])
    alpha[[1]][k]=pi[k]*dmvnorm(x[1:B[1]],mu[[1]][[k]], sigma[[1]][[k]])
  
  for(t in 2:Ti){
    alpha[[t]] = rep(0, S[t])
    for(k in 1:S[t])
    alpha[[t]][k] = dmvnorm(x[(B[t-1]+1):B[t]],mu[[t]][[k]], sigma[[t]][[k]])*sum(alpha[[t-1]]*a[[t-1]][,k])
  }
  return(alpha)
}

#Backward probability
backw = function(x, mu, sigma,a){
  beta = vector("list", Ti)
  beta[[Ti]]=rep(1, S[Ti])
  for(t in Ti-1:1){
    beta[[t]] = rep(0, S[t])
    for(k in 1:S[t])
      for(l in 1:S[t+1])
      beta[[t]][k] = beta[[t]][k] + (dmvnorm(x[(B[t]+1):B[t+1]], mu[[t+1]][[l]], sigma[[t+1]][[l]])*beta[[t+1]][l]*a[[t]][k,l])
  }
  return(beta)
}

#Compute likelihood
lik = function(x, mu, sigma, pi,a){
    lik = vector("list", Ti)
    alpha = forw(x, mu, sigma, pi,a)
    beta = backw(x, mu, sigma,a)
    for(t in 1:Ti){
      lik[[t]] = rep(0, S[t])
      for (k in 1:S[t]){
        lik[[t]][k] = alpha[[t]][k]*beta[[t]][k]/sum(alpha[[t]]*beta[[t]])
      }
    }
  return (lik)
}

#Compute joint probability
h = function(x, mu, sigma, pi,a){
    H = vector("list", Ti-1)
    alpha = forw(x, mu, sigma, pi,a)
    beta = backw(x, mu, sigma,a)
    for(t in 1:(Ti-1)){
      H[[t]] = matrix(0, S[t], S[t+1])
      for (k in 1: S[t]){
        for(l in 1:S[t+1]){
          H[[t]][k,l] = alpha[[t]][k]*a[[t]][k,l]*dmvnorm(x[(B[t]+1):B[t+1]], mu[[t+1]][[l]], sigma[[t+1]][[l]])*beta[[t+1]][l]/sum(alpha[[t]]*beta[[t]])
        }
      }
    }
  return(H)
}

#BW algorithm
BW = function(x, mu0, sigma0, pi, a,tol = 1e-4, maxit = 10){
  mu1 = mu0
  sigma1 = sigma0
  a1 = a
  pi1 = pi
  n=dim(x)[1]
#calculate likelihood for each observation in the data set
  for(j in 1:maxit){
    l_temp = vector("list", Ti)
    h_temp = vector("list", Ti-1)
    for(t in 1:Ti){
      l_temp[[t]] = vector("list", S[t])
      if(t < Ti)
        h_temp[[t]] = matrix(0, S[t], S[t+1])
        for( k in 1:S[t]){
          l_temp[[t]][[k]] = rep(0, n)
          for( i in 1:n){
            l1 = lik(x[i,], mu1, sigma1, pi1,a1)
            h1 = h(x[i,], mu1, sigma1, pi1,a1)
            l_temp[[t]][[k]][i] = l1[[t]][k]
            if(t < Ti && k ==1 )
              h_temp[[t]] = h1[[t]] + h_temp[[t]]
          }
        }
    }
#Use the likelihood to iteratively estimate the parameters
    for(t in 1:Ti){
      for(k in 1:S[t]){
        if(t==1) x_temp = x[,1:B[1]]
        else x_temp = x[,(B[t-1]+1):B[t]]
        mu1[[t]][[k]] = apply(l_temp[[t]][[k]]*x_temp, 2, sum)/sum(l_temp[[t]][[k]])
        sigma1[[t]][[k]] = (t(x_temp - mu1[[t]][[k]]))%*%(l_temp[[t]][[k]]*(x_temp - mu1[[t]][[k]]))/sum(l_temp[[t]][[k]])
        if(t ==1) pi1[k] = sum(l_temp[[1]][[k]])
        if(t < Ti)
          a1[[t]][k,] = h_temp[[t]][k,]/sum(l_temp[[t]][[k]])
      }
    }
    pi1 = pi1/sum(pi1)
  }
  return(list(mu = mu1, sigma = sigma1, pi = pi1, a = a1))
}

#Modal BW Algorithm
MBW = function(x0, mu, sigma, pi, a, tol = 1e-4, maxit = 100){
  x1 = x0
  x2 = x0
  for (i in 1:maxit){
    l1 = lik(x1, mu, sigma, pi,a)
    for(t in 1:Ti){
      if (t == 1) {
        num = rep(0, B[1])
        den = matrix(0, B[1], B[1])
      }
      else {
        num = rep(0, B[t] - B[t-1])
        den = matrix(0, B[t] - B[t-1],B[t] - B[t-1])
      }
      for(k in 1:S[t]){
        den = den + l1[[t]][k]*solve(sigma[[t]][[k]])
        num = num + l1[[t]][k]*solve(sigma[[t]][[k]])%*%mu[[t]][[k]]
      }
      if(t==1)
        x2[0:B[1]] = solve(den)%*%num
      else x1[(B[t-1]+1):B[t]] = solve(den)%*%num
    }
#    l2 = lik(x2, mu, sigma, pi,a)
#    if (abs(l1-l2)<tol)
#      break;
    x1 = x2
  }
  return(l1)
}

#simulation set up

x1_1 = rmvnorm(100, mean = c(0,5.5,5.5,0,0), sigma = diag(nrow = 5, ncol = 5))
x1_2 = rmvnorm(300, mean = c(10,10,5,20,0), sigma = diag(nrow = 5, ncol = 5))
#x1_3 = rmvnorm(100, mean = c(10,5,40,0,15), sigma = diag(nrow = 5, ncol = 5))
#x1_4 = rmvnorm(100, mean = c(0,10,20,5,5), sigma = diag(nrow = 5, ncol = 5))

x2_1 = rmvnorm(100, mean = c(0,10,20), sigma = diag(nrow = 3, ncol = 3))
x2_2 = rmvnorm(100, mean = c(10,20,0), sigma = diag(nrow = 3, ncol = 3))
x2_3 = rmvnorm(200, mean = c(20,5,10), sigma = diag(nrow = 3, ncol = 3))
#x2_4 = rmvnorm(60, mean = c(80,20,30), sigma = diag(nrow = 3, ncol = 3))
#x2_5 = rmvnorm(70, mean = c(10,40,40), sigma = diag(nrow = 3, ncol = 3))
#x2_6 = rmvnorm(60, mean = c(5,30,80), sigma = diag(nrow = 3, ncol = 3))

#S is number of states in each block, Ti is block number, B is position of last 
S = c(2,3)
Ti = 2
pi = rep(1,S[1])/S[1]
B = c(5, 8)

data1 = rbind(x1_1,x1_2)
data2 = rbind(x2_1,x2_2,x2_3)
data = cbind(rbind(x1_1,x1_2), rbind(x2_1,x2_2,x2_3))

init1 = kmeans(data1, S[1])
init2 = kmeans(data2, S[2])

mu1 = list(apply(data1[init1$cluster==1,],2,mean),apply(data1[init1$cluster==2,],2,mean))#,apply(data1[init1$cluster==3,],2,mean),apply(data1[init1$cluster==4,],2,mean))
sigma1 = list(var(data1[init1$cluster==1,]),var(data1[init1$cluster==2,]))#,var(data1[init1$cluster==3,]),var(data1[init1$cluster==4,]))

mu2 = list(apply(data2[init2$cluster==1,],2,mean),apply(data2[init2$cluster==2,],2,mean),apply(data2[init2$cluster==3,],2,mean))#,apply(data2[init2$cluster==4,],2,mean),apply(data2[init2$cluster==5,],2,mean),apply(data2[init2$cluster==6,],2,mean))
sigma2 = list(var(data2[init2$cluster==1,]),var(data2[init2$cluster==2,]),var(data2[init2$cluster==3,]))#,var(data2[init2$cluster==4,]),var(data2[init2$cluster==5,]),var(data2[init2$cluster==6,]))

mu0 = list(mu1, mu2)
sigma = list(sigma1, sigma2)

a = list(matrix(rep(1,prod(S))/S[2], nrow = S[1]))

#Use BW to estimate parameters
para = BW(data, mu0, sigma, pi, a, maxit =4)

#Get the group index of each observation with MBW
group = matrix(0, dim(data)[1], 2)
for( i in 1:dim(data)[1]){
  p = MBW(data[i,], para$mu, para$sigma, para$pi, para$a)
  group[i, ] = c(which.max(p[[1]]), which.max(p[[2]]))
}
group





