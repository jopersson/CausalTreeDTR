# Doubly-robust dynamic treatment regimen estimation using causal trees and causal forests
# Structure of function below adapted from https://github.com/cran/DTRreg/blob/master/R/DTRreg.R 
DTR <- function (outcome, 
                 outcome.mod, 
                 treatment.mod, 
                 data = NULL, 
                 optimum = NULL, 
                 ite = NULL,
                 approach = "causaltree", 
                 weight = "default", 
                 missing = "default",
                 type = "DTR", 
                 simulation = T,
                 minsize = 30) 
{
  # sanity checks
  try(match.arg(approach, c("qlearn", "dwols", "gestimation", "knn", "cart", "causaltree", "causalforest")))
  
  if (length(treatment.mod) != length(outcome.mod)) {
    stop("Treatment and outcome models must be of the same length.")
  }
  if (any(is.na(data$Y))) {
    stop("Missing values in the outcome not allowed.")
  }
  
  # initialize
  obj <- list()
  obj$approach <- approach
  K <- length(treatment.mod)
  obj$K <- K
  obj$minsize <- minsize
  
  # Q-Learning as special case of DWOLS
  qlearn <- 0
  if (approach == "qlearn") {
    approach <- "dwols"
    qlearn <- 1
  }
  
  if (class(treatment.mod) == "formula") {
    treatment.mod <- list(treatment.mod)
    outcome.mod <- list(outcome.mod)
  }
  
  
  # simulated data and ground truth optimal decisions
  data <- data.frame(data)
  rownames(data) <- seq(1:nrow(data))
  N <- length(data$Y)
  
  # assign inputs  
  obj$data <- data
  obj$optimum <- optimum
  obj$ite <- ite
  
  # assign models
  obj$outcome.mod <- outcome.mod
  obj$treatment.mod <- treatment.mod
  Y <- data$Y
  obj$opt.value <- as.numeric(0)
  
  # initialize outcomes
  obj$Y.obs <- Y           ## observed outcome
  obj$Y.hat <- list()      ## pseudo outcome
  obj$Y.hat[[K]] <- Y      ## pseudo outcome in final stage = observed outcome
  
  
  # initialize outputs
  obj$tau <- obj$beta <- obj$psi <- list()
  obj$trees <- obj$forests <- list()
  
  
  # initialize indices to cope with missing data and weighting
  dismiss <- c()
  keep <- c(1:N)
  remove <- rep(0, N)
  propscore <- rep(1,N)
  
  
  
  # BACKWARD RECURSION
  for (j in K:1) {
    Y.orig <- Y ## intermediately save pseudo outcome as factual outcome
    
    # response generation
    if(simulation == T){
      Y.hat <- obj$Y.hat[[j]]
    }
    
    # treatment model
    treat.var <- all.vars(treatment.mod[[j]])[1]
    A <- model.response(model.frame(treatment.mod[[j]], data, na.action = "na.pass"))
    Halpha <- model.matrix(treatment.mod[[j]], model.frame(treatment.mod[[j]], data, na.action = "na.pass"))
    
    # outcome model
    Hbeta <- model.matrix(outcome.mod[[j]], model.frame(outcome.mod[[j]], data, na.action = "na.pass"))
    Hpsi <- model.matrix(outcome.mod[[j]], model.frame(outcome.mod[[j]], data, na.action = "na.pass"))
    
    
    # missing data 
    if(missing == "keep"){
      keep <- c(1,N)
      propscore <- rep(1, N)
    }
    if (missing == "default") {
      propscore <- rep(1, N)
      new.dismiss <- which(apply(is.na(cbind(A, Hpsi, Halpha)), 1, any))
      if (j == K) 
        remove[new.dismiss] <- j
      if (j != K) {
        dismissk <- new.dismiss[!(new.dismiss %in% dismiss)]
        remove[dismissk] <- j
      }
      dismiss <- unique(new.dismiss, dismiss)
      
      if (length(new.dismiss) > 0 & missing == "ipcw" & j > 1) {
        Miss <- rep(1, N)
        Miss[new.dismiss] <- 0
        Miss <- Miss[keep]
        H.list <- unique(c(unlist(sapply(outcome.mod[(j - 1):1], all.vars)), unlist(sapply(treatment.mod[(j - 1):1], all.vars))))
        H <- data.matrix(data[, which(colnames(data) %in% H.list)])
        H <- H[keep, ]
        propscore[which(!apply(is.na(H), 1, any))] <- 1 - fitted(glm(Miss ~ ., data = as.data.frame(H), family = "binomial"))
        obj$ipcw[[j]] <- 1/propscore
      }
      
      if (missing == "ipcw") {
        keep <- c(1:N)
        keep <- keep[!(keep %in% new.dismiss)]
      }
      else {
        keep <- keep[!(keep %in% dismiss)]
      }
    }   
    propscore <- propscore[keep]
    A <- A[keep]
    Hpsi <- Hpsi[keep, ]
    Hbeta <- Hbeta[keep, ]
    Halpha <- Halpha[keep, ]
    Y <- Y[keep]
    
    
    Hpsi <- as.matrix(Hpsi)
    Halpha <- as.matrix(Halpha)
    obj$n[[j]] <- length(keep)
    n <- length(Y)
    
    # Treatment Model
    alpha <- glm(A ~ Halpha, binomial) 
    obj$treatment.mod.fitted[[j]] <- alpha
    Ahat <- fitted(alpha) # fitted propensity score
    
    #Q-Learning as special case of DWOLS
    qlearn <- 0
    if (obj$approach == "qlearn") {
      approach <- "dwols"
      qlearn <- 1
    }
    
    if (is.function(weight) == FALSE) {
      w <- abs(A - Ahat)
    }
    else {
      w <- weight(Ahat)
      w[A == 0] <- w[A == 0] * Ahat[A == 0]/(1 - Ahat[A == 0])
    }
    if (qlearn == 1) {
      w <- rep(1, n)
    }
    w <- w * propscore
    
    
    # approachS
    if (approach == "dwols") {
      beta <- solve(t(Hbeta) %*% (w * Hbeta)) %*% t(Hbeta) %*% (w * Y)
      tf <- Hbeta%*%beta
      
      psi <- solve(t(A * Hpsi) %*% (w * A * Hpsi)) %*% t(A * Hpsi) %*% (w * Y)
      tau <- Hpsi %*% psi
      
      obj$beta[[j]] <- beta
      obj$psi[[j]] <- psi
      
      Y.fit <- tf + as.numeric(tau > 0)*(tau)
      
    }
    
    if (approach == "gestimation") {
      Hd <- cbind(Hbeta,A*Hpsi)
      Hw <- cbind(Hbeta,Hpsi*(A - Ahat)*propscore)
      beta <- solve(t(Hbeta) %*% Hbeta) %*% t(Hbeta) %*% Y
      tf <- Hbeta%*%beta
      
      Hd.psi <- A*Hpsi
      Hw.psi <- Hpsi * (A - Ahat) * propscore
      psi <- solve(t(Hw.psi) %*% Hd.psi) %*% t(Hw.psi) %*% Y
      
      obj$beta[[j]] <- beta
      obj$psi[[j]] <- psi
      
      tau <- Hpsi %*% psi
      Y.fit <- tf + as.numeric(tau > 0)*(tau)
    }
    
    if (approach == "causaltree") {
      vars <- all.vars(outcome.mod[[j]])
      
      df <- data %>% select(c(vars, Y, treatment = all.vars(treatment.mod[[j]])[1]))
      df$Y <- Y
      df$weights <- w
      
      formula <- paste("Y", "~", paste(vars, collapse = "+"))
      formula <- as.formula(formula)
      
      # use causal tree as tau estimator
      tree <- causalTree(formula, data = df, treatment = df$treatment, weights = df$weights,
                         split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = T,
                         xval = 5, cp = 0, 
                         minsize = minsize, bucketMax = 20)
      
      
      obj$trees[[j]] <- tree
      tau <- predict(tree, type = "vector")
      Y.fit <- Y + as.numeric(tau > 0)*(tau)
    }
    if (approach == "causalforest"){
      vars <- all.vars(outcome.mod[[j]])
      
      cfformulaX <- paste(vars, collapse = "+")
      cfformula <- as.formula(paste("Y ~ -1 + ", cfformulaX))
      
      covariatematrix <- model.matrix(lm(cfformula, data = data))
      outcome <- as.matrix(Y)
      treatment <- as.matrix(data %>% select(treatment = all.vars(treatment.mod[[j]])[1]))
      
      # use causal forest as tau estimator
      forest <- causal_forest(covariatematrix, 
                              outcome, 
                              treatment, 
                              sample.weights = w,
                              num.trees = 500,
                              min.node.size = minsize)
      
      
      obj$forests[[j]] <- forest
      tau <- forest$predictions
      Y.fit <- Y + as.numeric(tau > 0)*(tau)
    }
    if(approach == "cart"){
      vars <- all.vars(outcome.mod[[j]])
      treatment <- all.vars(treatment.mod[[j]])[1]
      
      treeformulaX <- paste(vars, collapse = "+")
      treeformula <- as.formula(paste("Y ~ -1 +", treeformulaX))
      
      treematrix <- data %>% select(all_of(c(vars, treatment = all.vars(treatment.mod[[j]])[1]))) 
      treematrix$Y <- Y
      treematrix$weights <- w
      control <- treematrix %>% filter(.[["treatment"]] == 0)
      treated <- treematrix %>% filter(.[["treatment"]] == 1)
      
      cart0 <- rpart(treeformula,
                     data = control,
                     weights = control$weights,
                     minbucket = minsize/2) 
      cart1 <- rpart(treeformula,
                     data = treated,
                     weights = treated$weights,
                     minbucket = minsize/2) 
      
      cart0hat <- predict(cart0, newdata = treematrix)
      cart1hat <- predict(cart1, newdata = treematrix)
      
      tau <- cart1hat - cart0hat
      
      obj$trees[[j]] <- list(cart0,cart1)
      
      Y.fit <- Y + as.numeric(tau > 0)*(tau)
    }
    if(approach == "knn"){
      vars <- all.vars(outcome.mod[[j]])
      
      knnformulaX <- paste(vars, collapse = "+")
      knnformula <- as.formula(paste("Y ~ -1 + ", knnformulaX))
      
      covariatematrix <- model.matrix(lm(knnformula, data = data))
      
      treatment <- as.matrix(data %>% select(treatment = all.vars(treatment.mod[[j]])[1]))
      outcome <- as.matrix(data %>% select(Y))
      
      knn0 <- knn.reg(train = covariatematrix[treatment==0,], 
                      test = covariatematrix,
                      y = outcome[treatment==0], 
                      k = minsize/2)$pred 
      knn1 <- knn.reg(train = covariatematrix[treatment==1,],
                      test = covariatematrix,
                      y = outcome[treatment==1], 
                      k = minsize/2)$pred 
      
      tau <- knn1 - knn0
      Y.fit <- Y + as.numeric(tau > 0)*(tau)
    }
    
    
    ite.true <- ite[[j]]
    obj$tau[[j]] <- tau
    
    # SUBSEQUENT: treatment recommendations based on estimated DTR (optimization stage)
    optimum.true <- as.numeric(optimum[[j]])   ## true optimal decision
    optimum.true <- as.vector(optimum.true)
    obj$opt.true[[j]] <- optimum[[j]]          ## save true optimal
    
    opt <- as.numeric(tau > 0)                 ## estimated optimal decision
    opt <- as.vector(opt)
    obj$opt.hat[[j]] <- opt                    ## save estimated optimal 
    
    
    # Outcome Construction
    if(simulation == T){
      if(j > 1){
        obj$Y.hat[[j-1]] <- obj$Y.hat[[j]] + (opt - A) * (tau) 
      }
    }
    
    # final counterfactual outcome under optimal DTR
    Y.orig[keep] <- Y.orig[keep] + (as.numeric(type == "DTR") * (opt) - A) * (tau) # - A * tau 
    
    # Regret and reward (do we have to include true ITE somewhere?)
    if(simulation == T){
      obj$regret[[j]] <- (optimum.true - opt) * (ite.true) 
      ## optimum = opt: zero regret
      ## optimum = 1, opt = 0: regret is true ITE (>0*>0 > 0)
      ## optimum = 0, opt = 1: regret is true ITE (<0*<0 > 0)
    }
    
    obj$reward[[j]] <- (opt - A) * (tau) 
    ## opt = A: zero reward
    ## opt = 1, A = 0: reward is estimated ITE (>0*>0 > 0)
    ## opt = 0, A = 1: reward is estimated ITE (>0*>0 > 0)
    
    k <- j
    while (k < K) {
      # do nothing in j = K final stage
      # j = K - 1: one step back: construct pseudo outcome by substracting reward from outcome
      # j = K - 2: another step back: construct pseudo outcome from Y.hat(1) to Y.hat(2): iterate forwards over time points
      k <- k + 1
      
      # substract all gains to create a new dependant variable for the next iteration 
      Y.fit[which(rownames(obj$reward[[k]]) %in% keep)] <- Y.fit[which(rownames(obj$reward[[k]]) %in% keep)] - obj$reward[[k]][which(rownames(obj$reward[[k]]) %in% keep)]
    }
    
    # save fitted outcomes
    obj$fitted[[j]] <- Y.fit
    
    # update Y for next iteration: increased Y by pseudo outcome 
    Y <- Y.orig 
    
  }
  
  obj$remove <- remove
  obj$Y.opt <- Y 
  obj$type <- type
  
  class(obj) <- "DTR"
  obj
}


# optimal treatment and outcome prediction for new data
predict.DTR <- function(object, newdata) {
  x <- object
  approach <- x$approach 
  
  newdata <- as.data.frame(newdata)
  n <- nrow(newdata)
  Y <- newdata$Y
  
  optimum <- list()
  
  for(j in x$K:1){
    Y.orig <- Y
    
    A <- model.response(model.frame(x$treatment.mod[[j]], newdata, na.action = "na.pass"))
    
    
    if(approach == "qlearn" | approach == "dwols" | approach == "gestimation"){
      Hbeta <- model.matrix(x$outcome.mod[[j]], model.frame(x$outcome.mod[[j]], newdata, na.action = "na.pass"))
      tf <- Hbeta%*%x$beta[[j]]
      
      Hpsi <- model.matrix(x$outcome.mod[[j]], model.frame(x$outcome.mod[[j]], newdata, na.action = "na.pass"))
      tau <- Hpsi %*% as.numeric(x$psi[[j]])
    }
    if(approach == "causaltree"){
      vars <- all.vars(x$outcome.mod[[j]])
      df <- newdata %>% select(all_of(vars), Y, treatment = all.vars(x$treatment.mod[[j]])[1])
      
      tau <- predict(x$trees[[j]], newdata = df, type = "vector")
      
    }
    if(approach == "causalforest"){
      vars <- all.vars(x$outcome.mod[[j]])
      cfformulaX <- paste(vars, collapse = "+")
      cfformula <- as.formula(paste("Y ~ -1 + ", cfformulaX))
      newdatamatrix <- model.matrix(cfformula, data = newdata)
      
      tau <- predict(x$forests[[j]], newdatamatrix)
      
    }
    if(approach == "cart"){
      vars <- all.vars(x$outcome.mod[[j]])
      treeformulaX <- paste(vars, collapse = "+")
      treeformula <- as.formula(paste("Y ~ -1 + ", treeformulaX))
      
      newtreematrix <- newdata %>% select(all_of(vars), treatment = all.vars(x$treatment.mod[[j]])[1])
      
      cart0hat <- predict(x$trees[[j]][[1]], newdata = newtreematrix)
      cart1hat <- predict(x$trees[[j]][[2]], newdata = newtreematrix)
      
      tau <- cart1hat - cart0hat
    }
    if(approach == "knn"){
      treatment <- as.matrix(x$data %>% select(treatment = all.vars(x$treatment.mod[[j]])[1]))
      outcome <- as.matrix(x$data %>% select(Y))
      
      vars <- all.vars(x$outcome.mod[[j]])
      knnformulaX <- paste(vars, collapse = "+")
      knnformula <- as.formula(paste("Y ~ -1 + ", knnformulaX))
      
      covariatematrix <- model.matrix(lm(knnformula, data = x$data))
      newdatamatrix <- model.matrix(knnformula, data = newdata)
      
      
      knn0 <- knn.reg(train = covariatematrix[treatment==0,], 
                      test = newdatamatrix,
                      y = outcome[treatment==0], 
                      k = x$minsize*2)$pred 
      knn1 <- knn.reg(train = covariatematrix[treatment==1,],
                      test = newdatamatrix,
                      y = outcome[treatment==1], 
                      k = x$minsize*2)$pred 
      
      tau <- knn1 - knn0
    }
    
    
    opt <- as.numeric(tau > 0)
    opt <- as.vector(opt)
    optimum[[j]] <- opt  
    
    Y.orig <- Y.orig + (opt - A) * (tau) 
    Y <- Y.orig
  }
  
  
  value <- Y
  
  out <- list(value = value, opt = optimum)
  return(out)
}

#' Evaluation Metric
decision_accuracy.DTR <- function(fit, optimum){
  # decision accuracy: estimated versus optimal, mean over stages
  acc <- mean(sapply(fit$opt.hat == optimum, mean))
  return(mean(acc))
}

#' Evaluation Run
evaluate.DTR <- function(scenario="linear_interaction", approach = "causaltree", n = 5000, minsize = 20){
  
  data1 <- data_generation(n, model = scenario)
  ind <- sample(n, size = 0.75*n, replace = F)
  
  # train test split
  data1_train <- data1$data[ind,]
  data1_test <- data1$data[-ind,]
  
  optimum_train <- data1$optimum[ind,]
  optimum_test <- data1$optimum[-ind,]
  
  ite_train <- data1$ite[ind,]
  ite_test <- data1$ite[-ind,]
  
  # Models
  outcome.mod <- list(~C1 + C2 + X11 + X12 + X13 + X14 + X15 + X18 + X19,
                      ~C1 + C2 + X21 + X22 + X23 + X24 + X25 + X28 + X29 + A1 + (X11 + X12 + X13 + X14 + X15 + X18 + X19),
                      ~C1 + C2 + X31 + X32 + X33 + X34 + X35 + X38 + X39 + A2 + (X21 + X22 + X23 + X24 + X25 + X28 + X29) + A1 + (X11 + X12 + X13 + X14 + X15 + X18 + X19))
  treatment.mod <- list(A1 ~ C1 + C2 + X11 + X12 + X13, 
                        A2 ~ C1 + C2 + X21 + X22 + X23 + A1 + (X11 + X12 + X13 + X19),
                        A3 ~ C1 + C2 + X31 + X32 + X33 + A1 + (X11 + X12 + X13 + X19) + A2 + (X21 + X22 + X23 + X29))
  
  # Scenario of a linear outcome function
  ## fit models to standard case
  f1 <- DTR(outcome = Y, 
            outcome.mod = outcome.mod, 
            treatment.mod = treatment.mod, 
            approach = approach, 
            data = data1_train, 
            optimum = optimum_train, 
            ite = ite_train,
            minsize = minsize)
  
  # Evaluation
  ## Decision Accuracy
  accuracy.train <- decision_accuracy.DTR(f1, optimum_train)
  pred.opt <- predict.DTR(f1, newdata = data1_test)$opt
  accuracy.test <- mean(apply(pred.opt == optimum_test, 2, mean))
  
  ## Cumulative regret: sum over stages, mean over patients
  regret.train <- sum(sapply(f1$regret, mean))
  regret.test <- sum(sapply((optimum_test - pred.opt) * (ite_test), mean))
  
  # Evaluation metrics for training and test set
  train <- as.data.frame(cbind(accuracy = accuracy.train, 
                               regret = regret.train))
  test <- as.data.frame(cbind(accuracy = accuracy.test, 
                              regret = regret.test))
  
  return(list(train = train,
              test = test))
}

