# Designing to detect heteroscedasticity in a regression model
# by Lanteri, A., Leorato, S., Lopez-Fidalgo, J. and Tommasi C.

# Code by Lanteri, A.

# this function returns the log likelihood ratio statistics
# under the settings of the manuscript

ddhrm_llr_test=function(mean_function,               # mean function, as a function of beta and X
                        h_function,                  # function h(), as a function of gamma and X
                        beta,                        # parameters for the mean function, as a vector
                        sigma_sq=1,                  # constant par of the variance function, as a number
                        lambda,                      # parameter lambda of local alternative. Remember that gamma=lambda/sqrt(n)
                        xi=rbind(c(0,1),c(0.5,0.5)), # design of experiment, as a 2 by k matrix: first row design points, second row design weights
                        n=100                        # sample size n, as an integer
){
  require(optimx)
  
  gamma=lambda/sqrt(n)
  
  X = rep(xi[1,],floor(n*xi[2,]))
  y = rnorm(n = n,
            mean = mean_function(beta=beta,X=X),
            sd = sqrt(sigma_sq*h_function(gamma=gamma,X=X)))
 
  m=length(beta)
  s=length(gamma)
  
  #log likelihood under H1
  ll_1=function(par1,y,X,mean_function,h_function){
       beta1=par1[1:m]
       gamma1=par1[-(1:(m+1))]
       sigma_sq1=par1[m+1]
    -sum(dnorm(y,
               mean = mean_function(beta=beta1,X=X),
               sd=sqrt(sigma_sq1*h_function(gamma=gamma1,X=X)),
               log=T))}

  #log likelihood under H0
  ll_0=function(par0,y,X,mean_function){
    beta0=par0[1:m]
    sigma_sq0=par0[m+1]
    -sum(dnorm(y,
               mean = mean_function(beta=beta0,X=X),
               sd=sqrt(sigma_sq0),
               log=T))}

# a single vector for all starting parameters in H1, for optimization
# (if 1 is not a viable value for some paramenter, modify below)
par1=rep(1,m+s+1)

# a single vector for all starting parameters in H0, for optimization
# (if 1 is not a viable value for some paramenter, modify below)
par0=rep(1,m+1) 
  
opt1=optimx(par=par1, fn=ll_1, method = "Nelder-Mead",
            y=y, X=X, mean_function=mean_function, h_function=h_function)
opt0=optimx(par=par0, fn=ll_0, method = "Nelder-Mead",
            y=y, X=X, mean_function=mean_function)

LR_test=2*(opt0$value-opt1$value)  

return(LR_test)
}

# an example of application on Case 1 with Ds design, n=100 and lambda=5:
ddhrm_llr_test(mean_function=function(beta,X){beta[1]+beta[2]*X}, 
               h_function=function(gamma,X){exp(gamma*X)},        
               beta=c(1,1),                                       
               sigma_sq=1,                                        
               lambda=c(5),                                       
               xi=rbind(c(0,1),c(0.5,0.5)),                      
               n=100                                             
)

# # commented below an example on how to reiterate the function to obtain one
# # entry of one of the tables in the manuscript.
# # WARNING: it may take more than 30minutes for a single entry
#
# MCMC_ITER=10000
# LR_tests=rep(NA,MCMC_ITER)
# for (i in 1:MCMC_ITER) {
#   LR_tests[i]=ddhrm_llr_test(mean_function=function(beta,X){beta[1]+beta[2]*X}, 
#                              h_function=function(gamma,X){exp(gamma*X)},        
#                              beta=c(1,1),                                       
#                              sigma_sq=1,                                        
#                              lambda=c(5),                                       
#                              xi=rbind(c(0,1),c(0.5,0.5)),                      
#                              n=25  
#   )
#   }
# mean(pchisq(LR_tests, df = 1)>(1-0.05),na.rm = T)
