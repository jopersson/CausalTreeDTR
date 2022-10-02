data_generation <- function(n, model = "linear_interaction"){
  expit <- function(x) { 1 / (1 + exp(-x))}
  
  # baseline covariates 
  C1 <- rnorm(n, mean = 2)
  C2 <- rnorm(n, mean = 1)
  
  # Stage 1
  X11 <- rnorm(n, mean = 1,  sd = 1)
  X12 <- rnorm(n, mean =  3,  sd = 1) 
  X13 <- rnorm(n, mean = -1.5, sd = 1)
  X14 <- rnorm(n, mean = 9,   sd = 1)
  X15 <- rnorm(n, mean = -5,  sd = 1)
  X16 <- rnorm(n, mean =  3,  sd = 1) 
  X17 <- rnorm(n, mean = 1.5, sd = 1)
  X18 <- rnorm(n, mean = 10,   sd = 1)
  X19 <- rbinom(n, size = 3, prob = c(0.3, 0.2, 0.5))
  X19 <- as.factor(X19)
  
  A1 <- rbinom(n, 1, 1-expit(2*C1 - C2 + 0.5*X11 + X12 + X13 - as.numeric(X19)))
  
  # Stage 2
  X21 <- X11 + 0.5*A1 + rnorm(n, sd = 2) 
  X22 <- X12 + 2*A1 + rnorm(n, sd = 2)  
  X23 <- X13 - 0.5*A1 + rnorm(n, sd = 2)
  X24 <- X14 - 2*A1 + rnorm(n, sd = 2)
  X25 <- X15 - 0.5*A1 + rnorm(n, sd = 2) 
  X26 <- X16 + 2*A1 + rnorm(n, sd = 2)  
  X27 <- X17 + 0.5*A1 + rnorm(n, sd = 2)
  X28 <- as.numeric(X18) + A1*rbinom(n, size = 2, prob = c(0.3, 0.2, 0.5))
  X28 <- as.factor(X28)
  X29 <- rbinom(n, size = 3, prob = c(0.3, 0.2, 0.5))
  X29 <- as.factor(X29)
  
  A2 <- rbinom(n, 1, 1-expit(2*C1 - C2 + 3*A1*(0.5*X11 + X12 + X13 + as.numeric(X19)) -0.5*X21 - X22 - X23 + as.numeric(X29)))
  
  # Stage 3
  X31 <- X21 + 0.5*A2 + rnorm(n, sd = 2) 
  X32 <- X22 + 2*A2 + rnorm(n, sd = 2)  
  X33 <- X23 - 0.5*A2 + rnorm(n, sd = 2)
  X34 <- X24 - 2*A2 + rnorm(n, sd = 2)
  X35 <- X25 - 0.5*A2 + rnorm(n, sd = 2) 
  X36 <- X26 + 2*A2 + rnorm(n, sd = 2)  
  X37 <- X27 + 0.5*A2 + rnorm(n, sd = 2)
  X38 <- as.numeric(X28) + A2*rbinom(n, size = 2, prob = c(0.3, 0.4, 0.1, 0.1,0.1))
  X38 <- as.factor(X38)
  X39 <- rbinom(n, size = 3, prob = c(0.3, 0.4, 0.3))
  X39 <- as.factor(X39)
  
  A3 <- rbinom(n, 1, 1-expit(C1 + C2 + 2*A1*(0.5*X11 + X12 + X13 + as.numeric(X19)) + 4*A2*(X21 - X22 - 2*X23 + as.numeric(X29)) - 0.5*X31 - X32 - X33 + as.numeric(X39)))
  
  
  if(model == "linear_interaction"){
    # outcome model (Chakraborty, p. 64)
    tau1 <- (C1 -5*C2 + 0.5*X12 - 20*X13*as.numeric(X19 == "0")*as.numeric(X18 != "0") + 3*X14*as.numeric(X18 == "3") + X15*as.numeric(X19 == "3")  + as.numeric(X19 == "2")*X13 - as.numeric(X19 == "1")*(X15) + as.numeric(X19 == "0")*(X11)) + as.numeric(X18 == "3")*as.numeric(X19 == "0")*X15  - as.numeric(X19 == "3")*as.numeric(X18 == "2")*X11 + as.numeric(X18 == "1")*(X14) + as.numeric(X18 == "0")*(X15)
    tau2 <- -(C1 + C2 + .5*X21 + 5*X22*as.numeric(X19 == "3") - 5*X24 + 10*as.numeric(X29 == "0")*X24 + as.numeric(X29 == "1")*X25 - as.numeric(X19 == "0")*as.numeric(X29 == "2")*(X22) + as.numeric(X28 == "3")*(X24) + A1*(X12 - as.numeric(X19 == "0")*X13 + as.numeric(X19 == "3")*X11  - as.numeric(X19 == "2")*X13 - as.numeric(X19 == "1")*(X12) + as.numeric(X19 == "0")*(X11)))
    tau3 <- (C1 + C2 - 5*X34*as.numeric(X29 != "1")*as.numeric(X39 == "0") + 5*X34*as.numeric(X39 == "1") +10*X33*as.numeric(X28 == "2") + X34*as.numeric(A2)*as.numeric(X29 == "3") + 3*as.numeric(X39 == "0")*X31 + as.numeric(X38 == "1")*as.numeric(X39 == "1")*10 - as.numeric(X29 == "2")*(X32) + as.numeric(X19 == "0")*as.numeric(X39 == "3")*(X31) + A2*(as.numeric(X19 == "0")*X22 - X23 + as.numeric(X29 == "3")*X21  - as.numeric(X19 == "0")*as.numeric(X29 == "2")*X23 - as.numeric(X29 == "1")*(X22) + as.numeric(X29 == "0")*(X25))) -2*(- as.numeric(X19 == "3")*as.numeric(X18 == "2")*X11 + as.numeric(X18 == "1")*(X14) + as.numeric(X18 == "0")*(X15))
    
    tf1 <- X11 + X12 + X13 + as.numeric(X19)
    tf2 <- X21 + X22 - X23 + as.numeric(X29)
    tf3 <- X31 + X32 - X33 + as.numeric(X39)
  }
  
  if(model == "nonlinear_interaction"){
    tau1 <- C1^2*cos(X11)*as.numeric(X18 == "0")*as.numeric(X19 != "3")*20 - 5*X12*X14*as.numeric(X18 != "1")*as.numeric(X19 != "2") - 2*X11*X15*as.numeric(X18 == "0")*as.numeric(X19 == "0") - 3*cos(C2 - 5)*exp(X15/X14)*as.numeric(X19 != "1")*as.numeric(X18 == "0") + as.numeric(X19 == "1")*(X13/X12)*5 + as.numeric(X18)^2
    tau2 <- C1^2*cos(X21 + 5*A1)*as.numeric(X29 == "0")*5 - 0.5*X22*X23 + 0.2*X21*X23*as.numeric(X29 == "3") - cos(3*C2 - 5*A1)*X21*as.numeric(X29 == "1")*5 + as.numeric(X29 == "2")*(X23/X22 + 5*A1)
    tau3 <- C1^2*cos(X31 + 5*A2)*as.numeric(X39 == "0")*5 - 0.5*X32*X33 + 0.2*X31*X33*as.numeric(X39 == "3") - cos(3*C2 - 5*A2)*X21*as.numeric(X39 == "1")*5 + as.numeric(X39 == "2")*(X33/X32 + 5*A2) + X11*as.numeric(X19 == "1")*3 + as.numeric(X19 == "2")*(X13/X12)*5
    
    tf1 <- 5 + (as.numeric(X19))^2 - 3*sin(C2^2 + X11 - X12 + X23)
    tf2 <- (as.numeric(X29))^2 - 3*sin(X21 - X22 + X23)
    tf3 <- (as.numeric(X39))^2 - 3*sin(X31 - X32 + X33)
  }
  
  # treatment assignment
  opt1 <- as.numeric(tau1 > 0)
  opt2 <- as.numeric(tau2 > 0)
  opt3 <- as.numeric(tau3 > 0)
  
  Y <-  C1 + C2 + tf1 + A1*tau1 + tf2 + (A2)*tau2 + tf3 + A3*tau3 + rnorm(n)
  Y.opt <- C1 + C2 + tf1 + opt1*tau1 + tf2 + opt2 * tau2  + tf3 + opt3*tau3
  
  # data frames 
  data <- as.data.frame(cbind(C1, C2, 
                              X11, X12, X13, X14, X15, X16, X17, X18, X19, A1,
                              X21, X22, X23, X24, X25, X26, X27, X28, X29, A2,
                              X31, X32, X33, X34, X35, X36, X37, X38, X39, A3,
                              Y, Y.opt))
  
  # true optimal decisions and ITEs
  optimum <- as.data.frame(cbind(opt1, opt2, opt3))
  ite <- as.data.frame(cbind(tau1, tau2, tau3))
  
  # outputs
  out <- list(data = data,       # observational patient data
              optimum = optimum, # optimal decisions (ground truth)
              ite = ite)         # true individual treatment effects (ground truth)
  return(out)
}




