library(PLNmodels)
library(nloptr)
library(parallel)
library(progress)
# Fonction pour sécuriser l'accès aux éléments
safe_get <- function(x, index) {
  tryCatch({
    return(x[index])  # Essayer d'accéder à l'élément
  }, error = function(e) {
    return(NA)  # Retourner NA si une erreur survient
  })
}
logfactorial <- function(n) { # Ramanujan's formula
  n[n == 0] <- 1 ## 0! = 1!
  n*log(n) - n + log(8*n^3 + 4*n^2 + n + 1/30)/6 + log(pi)/2
}
Compute_Sigma<-function(X,Y,O,M,S,B){
  n<-nrow(X)
  S2_bar<-colSums(S^2)
  MtM<-t(M)%*%M
  (MtM+diag(S2_bar))/n
}

PLN_vloglik <- function(X,Y,O,M,S,B,Omega) {
  n <- nrow(Y)
  p <- ncol(Y)
  d <-ncol(X)
  params<-c(M,S)
  Y <- as.matrix(Y) # réponses (n,p)
  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]
  O <- as.matrix(O)#matrix(0,n,p) # offsets (n,p)

  w <- rep(1,n) # poids (n)
  M <- matrix(params[1:(n*p)],n,p)  # (n,p)
  wbar<-sum(w)
  S <-matrix(params[(n*p+1):length(params)],n,p)

  # M <- matrix(params[((d*p)+(p*p)+1):(((d*p)+(p*p))+(n*p))],n,p)  # (n,p)
  # S <-  matrix(params[(((d*p)+(p*p))+(n*p)+1):(((d*p)+(p*p))+(n*p)+(n*p))],n,p)  # (n,p)
  # B<-matrix(params[1:(d*p)],d,p)
  B<-B
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)

  # Omega<-solve(as.matrix(Sigma))
  # M <- matrix(prams_vec_E_step[1:(n * p)], n, p)
  # S <- matrix(prams_vec_E_step[(n * p + 1):(n * p + n * p)], n, p)

  S2 <- S^2
  Z <- O + X %*% B + M
  A <- exp(Z + 0.5 * S2)
  # Sigma <- Compute_Sigma(X,Y,O=matrix(0,nrow(Y),ncol(Y)),M,S,B)#(1/n)* (t(M) %*% (M * w) + diag(w %*% S2))

  #nSigma <- t(M) %*% (M * w) + diag(w %*% S2)
  # Omega<-solve(Sigma)
  Omega<-Omega#MASS::ginv(Sigma)
  objective <- 0.5*log(det(Omega)) +rowSums(Y * Z- A + 0.5 * log(S2)) - 0.5 * rowSums((M%*%Omega)*M+S2*diag(Omega))

  ji<--rowSums(logfactorial(Y))+objective
  ji<-0.5*p+ji
  objective<- sum(ji)

  # sum(Y*Z-A+(1/2)*log(S2)) +(n/2)*(log(det(Omega)))-(n/2)*sum(diag(((1/n)*t(M)%*%M+(colSums(S2)))%*%Omega))+(n*p)/2-sum(.logfactorial(Y))
  #
  # sum(Y*Z-A+(1/2)*log(S2)) +(n/2)*(log(det(Omega)))-sum(diag((1/1)*nSigma%*%Omega))+(n*p)/2-sum(.logfactorial(Y))
  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  return(objective)
}

#-------------------------------- debut Fonctions à supprimer ----------------
PLN_elbo_Estep<- function(X,Y,O,M,S,B,Omega){
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  params<-c(M,S)
  Y <- as.matrix(Y) # réponses (n,p)
  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]
  O <- as.matrix(O)#matrix(0,n,p) # offsets (n,p)

  w <- rep(1,n) # poids (n)
  M <- matrix(params[1:(n*p)],n,p)  # (n,p)
  wbar<-sum(w)
  S <-matrix(params[(n*p+1):length(params)],n,p)

  # M <- matrix(params[((d*p)+(p*p)+1):(((d*p)+(p*p))+(n*p))],n,p)  # (n,p)
  # S <-  matrix(params[(((d*p)+(p*p))+(n*p)+1):(((d*p)+(p*p))+(n*p)+(n*p))],n,p)  # (n,p)
  # B<-matrix(params[1:(d*p)],d,p)
  B<-B
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)

  # Omega<-solve(as.matrix(Sigma))
  # M <- matrix(prams_vec_E_step[1:(n * p)], n, p)
  # S <- matrix(prams_vec_E_step[(n * p + 1):(n * p + n * p)], n, p)

  S2 <- S^2
  Z <- O + X %*% B + M
  A <- exp(Z + 0.5 * S2)
  # Sigma <- Compute_Sigma(X,Y,O=matrix(0,nrow(Y),ncol(Y)),M,S,B)#(1/n)* (t(M) %*% (M * w) + diag(w %*% S2))

  #nSigma <- t(M) %*% (M * w) + diag(w %*% S2)
  # Omega<-solve(Sigma)
  Omega<-Omega#MASS::ginv(Sigma)
  objective <- 0.5*log(det(Omega)) +rowSums(Y * Z- A + 0.5 * log(S2)) - 0.5 * rowSums((M%*%Omega)*M+S2*diag(Omega))

  ji<-objective
  ji<-ji
  objective<- sum(ji)

  # sum(Y*Z-A+(1/2)*log(S2)) +(n/2)*(log(det(Omega)))-(n/2)*sum(diag(((1/n)*t(M)%*%M+(colSums(S2)))%*%Omega))+(n*p)/2-sum(.logfactorial(Y))
  #
  # sum(Y*Z-A+(1/2)*log(S2)) +(n/2)*(log(det(Omega)))-sum(diag((1/1)*nSigma%*%Omega))+(n*p)/2-sum(.logfactorial(Y))
  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  return(objective)
}
PLN_elbo_Mstep<- function(X,Y,O,M,S,B,Omega) {
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  params<-c(M,S)
  Y <- as.matrix(Y) # réponses (n,p)
  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]
  O <- as.matrix(O)#matrix(0,n,p) # offsets (n,p)

  w <- rep(1,n) # poids (n)
  M <- matrix(params[1:(n*p)],n,p)  # (n,p)
  wbar<-sum(w)
  S <-matrix(params[(n*p+1):length(params)],n,p)


  B<-B


  S2 <- S^2
  Z <- O + X %*% B + M
  A <- exp(Z + 0.5 * S2)
  # Sigma <- Compute_Sigma(X,Y,O=matrix(0,nrow(Y),ncol(Y)),M,S,B)#(1/n)* (t(M) %*% (M * w) + diag(w %*% S2))

  #nSigma <- t(M) %*% (M * w) + diag(w %*% S2)
  # Omega<-solve(Sigma)
  Omega<-Omega#MASS::ginv(Sigma)
  objective <- 0.5*log(det(Omega)) +rowSums(Y * Z- A + 0.5 * log(S2)) - 0.5 * rowSums((M%*%Omega)*M+S2*diag(Omega))

  ji<-objective
  ji<-ji
  objective<- sum(ji)

  # sum(Y*Z-A+(1/2)*log(S2)) +(n/2)*(log(det(Omega)))-(n/2)*sum(diag(((1/n)*t(M)%*%M+(colSums(S2)))%*%Omega))+(n*p)/2-sum(.logfactorial(Y))
  #
  # sum(Y*Z-A+(1/2)*log(S2)) +(n/2)*(log(det(Omega)))-sum(diag((1/1)*nSigma%*%Omega))+(n*p)/2-sum(.logfactorial(Y))
  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  return(objective)
}
objective_E_step_SICPLN_old <- function(X,Y,O,params,B,Omega,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  nb_params<-(lambda/2)*((((p+1)*p)/2)+(p))
  M <- matrix(params[1:(n*p)],n,p)  # (n,p)
  S <-matrix(params[(n*p+1):length(params)],n,p)
  objective<-PLN_vloglik(X,Y,O,M,S,B,Omega)#PLN_elbo_Estep(X,Y,O,M,S,B,Omega)

  SIC_penalty<-(lambda/2) *(B^2 / (B^2 + epsilon^2))
  SIC_penalty[1,]<-0
  # entropy<-0.5 * (n * p * log(2*pi*exp(1)) + sum(log(S^2)))
  objective <- objective-sum(SIC_penalty)-nb_params#-(lambda/2)*entropy
  # sum(w * (A - Y * Z - 0.5 * log(S2))) - 0.5 *wbar* sum(det(Omega))
  # grad_M <- (diag(w) %*% (M %*% Omega + A - Y))
  # grad_S <- (diag(w) %*% (S * diag(Omega) + S * A - S^-1))
  objective<- objective

  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  return(-objective)
}
grad_E_step_SICPLN_old<- function(X,Y,O,params,B,Omega,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  Y <- as.matrix(Y) # réponses (n,p)
  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]
  O <- as.matrix(O) # offsets (n,p)
  w <- rep(1,n) # poids (n)
  # M <- matrix(params[((d*p)+(p*p)+1):(((d*p)+(p*p))+(n*p))],n,p)  # (n,p)
  # S <-  matrix(params[(((d*p)+(p*p))+(n*p)+1):(((d*p)+(p*p))+(n*p)+(n*p))],n,p)  # (n,p)
  M <-matrix(params[1:(n*p)],n,p)  # (n,p)
  S <-matrix(params[(n*p+1):length(params)],n,p)  # (n,p)
  S2<-S^2
  # B<-matrix(params[1:(d*p)],d,p)
  # matrix(params_init_E_step[(n*p+1):length(params_init_E_step)],n,p)
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  B<-B
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  # Sigma <- Compute_Sigma(X,Y,O=O,M,S,B)
  # Omega<-solve(Sigma)
  # Omega<-MASS::ginv(Sigma)
  # Omega<-Omega
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  # M <- matrix(prams_vec_E_step[1:(n * p)], n, p)
  # S <- matrix(prams_vec_E_step[(n * p + 1):(n * p + n * p)], n, p)
  S2 <- S^2
  Z <- O + X %*% B + M
  A <- exp(Z + 0.5 * S2)
  # Omega<-solve(Sigma)
  Omega<-Omega#MASS::ginv(Sigma)
  # nSigma <- t(M) %*% (M * w) + diag(w %*% S2)
  # SIC_penalty<-lambda * sum(beta^2 / (beta^2 + a^2))
  # objective <- sum(w * (A - Y * Z - 0.5 * log(S2))) + 0.5 * sum(diag(Omega %*% nSigma))+SIC_penalty

  grad_M <- (diag(w) %*% (M %*% Omega + A - Y))

  grad_S <- (diag(w) %*% (S * diag(Omega) + S * A - (S^-1)/2))

  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  gradient <- c(as.vector(grad_M), as.vector(grad_S))
  return(gradient)
}
objective_M_step_SICPLN_old <- function(X,Y,O,paramsM,M,S,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  nb_params<-(lambda/2)*((((p+1)*p)/2)+(p))#nb_params<-log(n)*((((p+1)*p)/2)+(0))/2
  B <- matrix(paramsM[1:(d*p)],d,p)  # (n,p)
  Omega <-  matrix(paramsM[(d*p+1):length(paramsM)],p,p)

  objective<-PLN_elbo_Mstep(X,Y,O,M,S,B,Omega)
  # entropy<-0.5 * (n * p * log(2*pi*exp(1)) + sum(log(S^2)))
  SIC_penalty<-(lambda/2) *(B^2 / (B^2 + epsilon^2))
  SIC_penalty[1,]<-0
  objective <- objective-sum(SIC_penalty)-nb_params#-(lambda/2)*entropy
  # sum(w * (A - Y * Z - 0.5 * log(S2))) - 0.5 *wbar* sum(det(Omega))
  # grad_M <- (diag(w) %*% (M %*% Omega + A - Y))
  # grad_S <- (diag(w) %*% (S * diag(Omega) + S * A - S^-1))
  objective<- objective
  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  return(-objective)
}
grad_M_step_SICPLN_old<- function(X,Y,O,paramsM,M,S,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)

  Y <- as.matrix(Y) # réponses (n,p)

  X <- as.matrix(X) # covariables (n,d)
  d<-ncol(X)
  # X<-X[,-2]
  O <- as.matrix(O) # offsets (n,p)
  w <- rep(1,n) # poids (n)
  # M <- matrix(params[((d*p)+(p*p)+1):(((d*p)+(p*p))+(n*p))],n,p)  # (n,p)
  # S <-  matrix(params[(((d*p)+(p*p))+(n*p)+1):(((d*p)+(p*p))+(n*p)+(n*p))],n,p)  # (n,p)
  B <- matrix(paramsM[1:(d*p)],d,p)  # (n,p)
  Omega <-  matrix(paramsM[(d*p+1):length(paramsM)],p,p)
  # B<-matrix(params[1:(d*p)],d,p)
  # matrix(params_init_E_step[(n*p+1):length(params_init_E_step)],n,p)
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  M<-M
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  S<-S
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  # M <- matrix(prams_vec_E_step[1:(n * p)], n, p)
  # S <- matrix(prams_vec_E_step[(n * p + 1):(n * p + n * p)], n, p)

  S2 <- S ^ 2
  Z <- O + X %*% B + M
  A <- exp(Z + 0.5 * S2)
  nSigma <- n* Compute_Sigma(X,Y,O=O,M,S,B)
  Sigma <- 1/n*(t(M) %*% (M * w) + diag(w %*% S2))

  # SIC_penalty<-lambda * sum(beta^2 / (beta^2 + a^2))
  # objective <- sum(w * (A - Y * Z - 0.5 * log(S2))) + 0.5 * sum(diag(Omega %*% nSigma))+SIC_penalty
  SIC_deriv<-(lambda/2) * (2*(B*epsilon^2) / (B^2 + epsilon^2)^2)
  SIC_deriv[1,]<-0
  grad_B <-t(X) %*% ((A - Y)) + SIC_deriv #(t(X) %*% (w * (A - Y)) - (lambda/2) * (B^3) / (B^2 + epsilon^2)^2)
  # grad_B[1,]<-0
  grad_Omega <-(-nSigma+MASS::ginv(nSigma))/2#solve((t(M) %*% (M * w) + diag(w %*% S2))%*%solve(nSigma))#
  # grad_Omega <-diag(0,p)#(-nSigma+solve(nSigma))/2
  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  gradient_M <- c(as.vector(grad_B), as.vector(grad_Omega))
  return(gradient_M)
}
SICPLN_old<-function(X,Y,offset=FALSE,lambda_fixed=(log(nrow(X))*ncol(X)),optim_method="BFGS",max_it=200,length_lambda=100,grid_search=FALSE){
  prep_pln<-prepare_data(Y,X)
  # Y<-prep_pln$Abundance
  # X<-prep_pln[,-c(1,length(prep_pln))]
  # d<-ncol(X)
  p<-ncol(prep_pln$Abundance)
  n<-nrow(prep_pln$Abundance)
  if(offset==TRUE){
    O<-matrix(rep(log(prep_pln$Offset),p),ncol=p,byrow = FALSE)
  }
  if(offset==FALSE){
    O<-matrix(0,nrow=n,ncol=p)
  }
  solution_pln<-PLNmodels::PLN(Abundance~.-Offset+offset((O[,1])),data=prep_pln,control=PLN_param(backend = "nlopt",config_optim=list(algorithm="LBFGS",maxeval=1000)))
  BIC_PLN<-solution_pln$BIC
  if(grid_search==FALSE){
    solution<-SICPLN_optim(X=X,Y=Y,offset=offset,lambda_fixed=lambda_fixed,optim_method=optim_method,max_it=200)
    solution<-solution
    return(solution)
  }
  if(grid_search==TRUE){
    # Creation de la grille
    length_lambda=length_lambda
    # lambda_min=0#(log(nrow(X))*(ncol(X)*ncol(Y))+ncol(Y))/10
    # lambda_max=3*(log(nrow(X))*(ncol(X)))
    # Charger le package parallel
    # library(parallel)

    # Calculer la valeur log(n) une seule fois pour l'utiliser plus tard
    log_n <- log(nrow(X))

    # Définir les valeurs spécifiques de lambda
    lambda_min <- 0
    lambda_max <- 5 * log_n * ncol(X)#*ncol(Y)
    lambda_sic1 <- log_n
    lambda_sic2 <- log_n * ncol(X)
    lambda_sic3 <- log_n * ncol(X) / 2
    lambda_sic4 <- log_n * ncol(X) * ncol(Y)

    # Créer la séquence de valeurs lambda
    # seq_lambda <- c(0, lambda_sic1, lambda_sic2, lambda_sic3, lambda_sic4,
    #                 10^(seq(0, log(lambda_max), length.out = length_lambda)))
    seq_lambda <- c(0, lambda_sic1,
                    10^(seq(0, log(lambda_max), length.out = length_lambda)))

    # Créer un cluster avec un nombre de cœurs spécifié
    # Créez le cluster
    n_cores <- detectCores() - 2  # Utiliser tous les cœurs sauf deux pour ne pas saturer le système
    cl <- makeCluster(n_cores)

    # Exporter les variables nécessaires et charger les bibliothèques sur les nœuds
    clusterExport(cl, varlist = ls(envir = .GlobalEnv))
    clusterEvalQ(cl, library(PLNmodels))

    # Envelopper l'exécution du code dans un tryCatch pour garantir l'arrêt propre du cluster
    tryCatch({
      res_fs <- list()
      BIC_lambda <- numeric(length(seq_lambda))

      # Paralléliser l'exécution de la boucle avec parLapply
      res_fs <- parLapply(cl, seq_lambda, function(lambda) {
        # Exécuter l'optimisation pour chaque valeur de lambda avec gestion des erreurs
        result <- tryCatch({
          SICPLN_optim(X = X, Y = Y, offset = offset, lambda_fixed = lambda,
                       optim_method = optim_method, max_it = 200)
        }, error = function(e) {
          message("Erreur lors de l'optimisation avec lambda = ", lambda, ": ", e$message)
          return(NULL)  # Retourner NULL en cas d'erreur
        })

        # Vérification si le résultat est NULL avant d'extraire les valeurs du BIC
        if (is.null(result)) {
          message("Aucun résultat pour lambda = ", lambda, ". BIC non disponible.")
          return(list(result = NULL, BIC = NA, BIC_approximer = NA))  # Retourner des valeurs manquantes ou NA
        } else {
          # Assurez-vous que le résultat a bien les valeurs attendues pour le BIC
          BIC_value <- ifelse(length(result) >= 11, result[10], NA)
          BIC_approx_value <- ifelse(length(result) >= 11, result[11], NA)

          return(list(result = result, BIC = BIC_value, BIC_approximer = BIC_approx_value))
        }
      })

    }, error = function(e) {
      # Gestion des erreurs (si nécessaire, afficher le message d'erreur)
      message("Une erreur est survenue: ", e$message)
    }#,
    # finally = {
    #   # Ce bloc est toujours exécuté, même en cas d'erreur ou d'arrêt manuel
    #   stopCluster(cl)
    #   message("Le cluster a été arrêté.")
    # }
    )

    stopCluster(cl)
    message("Le cluster a été arrêté.")
    # Extraire les BIC et les résultats
    BIC_lambda <- sapply(res_fs, function(x) x$BIC)
    BIC_lambda_approximer <- sapply(res_fs, function(x) x$BIC_approximer)
    res_fs <- lapply(res_fs, function(x) x$result)

    # Nommer les éléments dans la liste des résultats
    names(res_fs) <- paste0("SICPLN_lambda_", round(1:length(seq_lambda), 0))

    # Trouver le meilleur lambda basé sur le BIC (plus bas BIC)
    best_lambda_indice <- which(BIC_PLN == min(unlist(BIC_lambda)))
    best_lambda_indice_approximer <- which(BIC_PLN == min(unlist(BIC_lambda_approximer)))

    # Créer la solution finale
    solution <- list(solution = res_fs, seq_lambda = seq_lambda,
                     BIC_all_lambda = unlist(BIC_lambda), best_lambda_indice = best_lambda_indice,best_lambda_indice_approximer=best_lambda_indice_approximer,BIC_lambda_approximer=unlist(BIC_lambda_approximer))

    # Arrêter le cluster après l'exécution
    stopCluster(cl)

  }
  return(solution)
}

#-------------------------------- Fin Fonctions à supprimer ----------------

objective_E_step_SICPLN <- function(X,Y,O,params,B,Omega,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  nb_params<-(lambda/2)*((((p+1)*p)/2)+(p))
  M <- matrix(params[1:(n*p)],n,p)  # (n,p)
  S <-matrix(params[(n*p+1):length(params)],n,p)
  objective<-PLN_vloglik(X,Y,O,M,S,B,Omega)

  SIC_penalty<-(lambda/2) *(B^2 / (B^2 + epsilon^2))
  SIC_penalty[1,]<-0
  objective <- objective-sum(SIC_penalty)-nb_params
  # sum(w * (A - Y * Z - 0.5 * log(S2))) - 0.5 *wbar* sum(det(Omega))
  # grad_M <- (diag(w) %*% (M %*% Omega + A - Y))
  # grad_S <- (diag(w) %*% (S * diag(Omega) + S * A - S^-1))
  objective<- objective

  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  return(-objective)
}


grad_E_step_SICPLN<- function(X,Y,O,params,B,Omega,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d <-ncol(X)
  Y <- as.matrix(Y) # réponses (n,p)
  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]
  O <- as.matrix(O) # offsets (n,p)
  # w <- rep(1,n) # poids (n)
  # M <- matrix(params[((d*p)+(p*p)+1):(((d*p)+(p*p))+(n*p))],n,p)  # (n,p)
  # S <-  matrix(params[(((d*p)+(p*p))+(n*p)+1):(((d*p)+(p*p))+(n*p)+(n*p))],n,p)  # (n,p)
  M <-matrix(params[1:(n*p)],n,p)  # (n,p)
  S <-matrix(params[(n*p+1):length(params)],n,p)  # (n,p)
  S2<-S^2
  # B<-matrix(params[1:(d*p)],d,p)
  # matrix(params_init_E_step[(n*p+1):length(params_init_E_step)],n,p)
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  B<-B
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  # Sigma <- Compute_Sigma(X,Y,O=O,M,S,B)
  # Omega<-solve(Sigma)

  # Omega<-Omega
  # n <- nrow(Y)
  # p <- ncol(Y)
  # d<-ncol(X)
  # M <- matrix(prams_vec_E_step[1:(n * p)], n, p)
  # S <- matrix(prams_vec_E_step[(n * p + 1):(n * p + n * p)], n, p)
  # S2 <- S^2
  Z <- O + X %*% B + M
  A <- exp(Z + 0.5 * S2)
  # pas nessaire car omega est fixer ?
  # S_hat<-diag(colSums(S^2))
  # Omega_cal<-chol2inv(chol((1/n)*(S_hat+t(M)%*%(M))))
  # Omega<-Omega_cal
  Omega<-Omega
  # nSigma <- t(M) %*% (M * w) + diag(w %*% S2)
  # SIC_penalty<-lambda * sum(beta^2 / (beta^2 + a^2))
  # objective <- sum(w * (A - Y * Z - 0.5 * log(S2))) + 0.5 * sum(diag(Omega %*% nSigma))+SIC_penalty

  grad_M <-((M %*% Omega + A - Y))

  grad_S <- ( (S * diag(Omega) + S * A - (S^-1)/2))#(S^-1)/2)
  # grad_S <--(S*A+(1/S)-S%*%diag((diag(Omega))))
  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  gradient <- c(as.vector(grad_M), as.vector(grad_S))
  return(gradient)
}

objective_M_step_SICPLN <- function(X,Y,O,paramsM,Omega,M,S,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  nb_params<-(lambda/2)*((((p+1)*p)/2)+(p))#nb_params<-log(n)*((((p+1)*p)/2)+(0))/2
  B <- matrix(paramsM[1:(d*p)],d,p)  # (n,p)
  Omega<-Omega
  #Omega <-  matrix(paramsM[(d*p+1):length(paramsM)],p,p)

  objective<-PLN_vloglik(X,Y,O,M,S,B,Omega)

  SIC_penalty<-(lambda/2) *(B^2 / (B^2 + epsilon^2))
  SIC_penalty[1,]<-0
  objective <- objective-sum(SIC_penalty)-nb_params
  # sum(w * (A - Y * Z - 0.5 * log(S2))) - 0.5 *wbar* sum(det(Omega))
  # grad_M <- (diag(w) %*% (M %*% Omega + A - Y))
  # grad_S <- (diag(w) %*% (S * diag(Omega) + S * A - S^-1))
  objective<- objective
  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  return(-objective)
}

grad_M_step_SICPLN <- function(X,Y,O,paramsM,Omega,M,S,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)

  Y <- as.matrix(Y) # réponses (n,p)

  X <- as.matrix(X) # covariables (n,d)
  d<-ncol(X)
  # X<-X[,-2]
  O <- as.matrix(O) # offsets (n,p)
  w <- rep(1,n) # poids (n)
  # M <- matrix(params[((d*p)+(p*p)+1):(((d*p)+(p*p))+(n*p))],n,p)  # (n,p)
  # S <-  matrix(params[(((d*p)+(p*p))+(n*p)+1):(((d*p)+(p*p))+(n*p)+(n*p))],n,p)  # (n,p)
  B <- matrix(paramsM[1:(d*p)],d,p)  # (n,p)
  #Omega <-  matrix(paramsM[(d*p+1):length(paramsM)],p,p)
  # B<-matrix(params[1:(d*p)],d,p)
  # matrix(params_init_E_step[(n*p+1):length(params_init_E_step)],n,p)
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  M<-M
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  S<-S
  # n <- nrow(Y)
  # p <- ncol(Y)
  # d<-ncol(X)
  # M <- matrix(prams_vec_E_step[1:(n * p)], n, p)
  # S <- matrix(prams_vec_E_step[(n * p + 1):(n * p + n * p)], n, p)

  S2 <- S ^ 2
  Z <- O + X %*% B + M
  A <- exp(Z + 0.5 * S2)
  # nSigma <- n* Compute_Sigma(X,Y,O=O,M,S,B)
  # Sigma <- 1/n*(t(M) %*% (M * w) + diag(w %*% S2))

  # SIC_penalty<-lambda * sum(beta^2 / (beta^2 + a^2))
  # objective <- sum(w * (A - Y * Z - 0.5 * log(S2))) + 0.5 * sum(diag(Omega %*% nSigma))+SIC_penalty
  SIC_deriv<-(lambda/2) * (2*(B*epsilon^2) / (B^2 + epsilon^2)^2)
  SIC_deriv[1,]<-0
  grad_B <-t(X) %*% ((A - Y)) + SIC_deriv #(t(X) %*% (w * (A - Y)) - (lambda/2) * (B^3) / (B^2 + epsilon^2)^2)
  # grad_B[1,]<-0
  # grad_Omega <-(-nSigma+MASS::ginv(nSigma))/2#solve((t(M) %*% (M * w) + diag(w %*% S2))%*%solve(nSigma))#
  # grad_Omega <-diag(0,p)#(-nSigma+solve(nSigma))/2
  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  gradient_M <- (as.vector(grad_B))
  return(gradient_M)
}
SICPLN_optim<-function(X,Y,offset=FALSE,lambda_fixed=(log(nrow(X))*ncol(X)),optim_method="BFGS",max_it=200){
  prep_pln<-prepare_data(Y,X)
  Y<-prep_pln$Abundance
  X<-prep_pln[,-c(1,length(prep_pln))]
  d<-ncol(X)
  p<-ncol(Y)
  n<-nrow(X)
  lambda<-lambda_fixed
  # if(lambda_fixed==TRUE){lambda<-log(nrow(X))*(ncol(X))}
  # if(lambda_fixed!=TRUE){lambda<-log(nrow(X))*(ncol(X))}
  if(offset==TRUE){
    O<-matrix(rep(log(prep_pln$Offset),p),ncol=p,byrow = FALSE)
  }
  if(offset==FALSE){
    O<-matrix(0,nrow=n,ncol=p)
  }

  X<-model.matrix(~.,as.data.frame(X))

  # Offset_mat<-matrix(0,n,p)
  res_PLN<-PLNmodels::PLN(Abundance~.-Offset+offset((O[,1])),data=prep_pln,control=PLN_param(backend = "nlopt",config_optim=list(algorithm="LBFGS",maxeval=1000)))
  # params_init_E_step<-c(as.vector(res_PLN$var_par$M),as.vector(res_PLN$var_par$S))
  # params_init_M_step<-c(as.vector(res_PLN$model_par$B),as.vector(res_PLN$model_par$Omega))
  B<-res_PLN$model_par$B
  Omega<-res_PLN$model_par$Omega
  S<-res_PLN$var_par$S
  M<-res_PLN$var_par$M
  X<-as.matrix(X)

  Y<-as.matrix(Y)
  d<-ncol(X)
  # param_initpln=param_init(variables =X[,-1],abondance =Y,offset= matrix(0,ncol=ncol(Y),nrow=nrow(Y)))
  # formatp=format_plnnew(abondance=Y,covariables=X,B_reg=param_initpln$B,Sigma =diag(1,p) ,M_reg=param_initpln$M,S_reg=param_initpln$S,O = matrix(0,ncol=ncol(Y),nrow=nrow(Y)),w =rep(1,nrow(Y)))
  # params_init_E_step<-c(as.vector(param_initpln$M),as.vector(param_initpln$S))
  params_init_E_step<-c(as.vector(M),as.vector(S))
  # v <- res_PLN$model_par$B
  # #c(0.44844978, -0.08463290, -0.31860861, 0.50380695, 0.96334493, -0.18933897,
  #        0.50788630, -0.02754789, 3.14563421)
  # v<-rep(3,9)
  params_init_M_step<-c(as.vector(B),as.vector((Omega)))
  # param_initpln$B
  # B<-v#matrix(v,3,3)
  # Omega<-diag(1,p)

  length_epsilon<-100
  E<-c()
  e1<-10
  E[1]<-e1
  for(t in 2:length_epsilon){
    E[t]=e1*(0.87)^(t-1)
  }
  iter<-0
  for(epsilon_val in E) {
    iter<-iter+1
    cat("Epsilon value for telescoping :\t",iter,epsilon_val,"\t","loglik :\t",PLN_vloglik(X,Y,O=O,M,S,B,Omega),"\n")
    for(i in 1:max_it){
      result_E_step <- optim(
        fn = objective_E_step_SICPLN,
        gr = grad_E_step_SICPLN,
        Y=Y,
        X=(X),
        O=O,
        par = params_init_E_step,
        B=B,
        Omega=Omega,
        lambda=lambda,
        epsilon=epsilon_val,
        method = optim_method,
        control = list(
          maxit = 1000
          # gradtol = 1e-8
        )
      )
      M <- (matrix(result_E_step$par[1:(n*p)],n,p))  # (n,p)
      S <-  matrix(result_E_step$par[(n*p+1):length(result_E_step$par)],n,p)
      #cat("E_step : ",result_E_step$value,"\n")
      # result_E_step$value
      params_init_E_step<-c(as.vector(M),as.vector(S))

      result_M_step <- optim(
        fn = objective_M_step_SICPLN,
        gr = grad_M_step_SICPLN,
        Y=Y,
        X=(X),
        O=O,
        par = params_init_M_step,
        M=M,
        S=S,
        lambda=lambda,
        epsilon=epsilon_val,
        method = optim_method,
        control = list(
          maxit = 1000
          # trace=TRUE
          # gradtol = 1e-8
        ))

      B <- matrix(result_M_step$par[1:(d*p)],d,p)  # (n,p)
      #Compute_Sigma(X,Y,O=0,M,S,B)
      Sigma<-Compute_Sigma(X,Y,O=0,M,S,B)
      Omega<-MASS::ginv(Sigma)
      # Omega<- solve( Sigma)#res_PLN$model_par$Omega#matrix(result_M_step$par[(d*p+1):length(result_M_step$par)],p,p)
      params_init_M_step<-c(as.vector(B),as.vector(Omega))
      #cat("M_step : ",result_M_step$value,"\n")
      convergence<-abs(result_E_step$value-result_M_step$value)
      if(convergence<=1e-8){
        cat("epsilon telescoping : ",epsilon_val,"Convergence avant maxit\n")
        break}
    }

  }
  # Post traitement
  # Y
  if (!is.null(colnames(Y))) {
    colnames(Omega)<-colnames(Y)
    rownames(Omega)<-colnames(Y)
    colnames(Sigma)<-colnames(Y)
    rownames(Sigma)<-colnames(Y)
    colnames(M)<-colnames(Y)
    colnames(S)<-colnames(Y)
    colnames(B)<-colnames(Y)
  }
  if(!is.null(rownames(Y))) {
    rownames(M)<-rownames(Y)
    rownames(S)<-rownames(Y)
  }
  # X
  if (!is.null(colnames(X))) {
    rownames(B)<-colnames(X)#c("Intercept",colnames(X))
    # rownames(Omega)<-colnames(Y)
    # colnames(Sigma)<-colnames(Y)
    # rownames(Sigma)<-colnames(Y)
    # colnames(M)<-colnames(Y)
    # colnames(S)<-colnames(Y)
  }
  # if(!is.null(rownames(X))) {
  #   rownames(B)<-rownames(X)
  #   # rownames(S)<-rownames(X)
  # }

  loglik<-PLN_vloglik(X,Y,O=O,M,S,B,Omega)
  B[abs(B)<=1e-5]<-0

  BIC_SICPLN<-loglik-0.5*log(n)*(res_PLN$nb_param)
  approxi_BIC_SICPLN<-loglik-0.5*log(n)*(res_PLN$nb_param-length(which(c(B[-1,])==0)))
  model_par<-list(B=B,Omega=Omega,Sigma=Sigma)
  var_par<-list(M=M,S=S,S2=S^2)
  data<-list(X=X,Y=Y,O=O)
  SICPLN_prediction<-(exp(O+X%*%B+M+S^2/2))
  q<-p
  entropy<-0.5 * (n * q * log(2*pi*exp(1)) + sum(log(S^2)))
  ICL<-BIC_SICPLN-entropy
  ICL_approx<-approxi_BIC_SICPLN-entropy
  return(list(model_par=model_par,var_par=var_par,data=data,elbo=result_M_step$value,loglik=loglik,BIC_SICPLN=BIC_SICPLN,approxi_BIC_SICPLN=approxi_BIC_SICPLN,SICPLN_prediction=SICPLN_prediction,ICL=ICL,ICL_approx=ICL_approx,res_PLN=res_PLN))
}
SICPLN_NLOPTR<-function(X,Y,offset=FALSE,lambda_fixed=(log(nrow(X))*ncol(X)),optim_method="BFGS",max_it=200,initialisation_method="LM"){
  # if (is.null(w)) w <- rep(1.0, nrow(Y))
  prep_pln<-PLNmodels::prepare_data(Y,X)
  Y<-prep_pln$Abundance
  X<-prep_pln[,-c(1,length(prep_pln))]
  d<-ncol(X)
  p<-ncol(Y)
  n<-nrow(X)
  w <- rep(1.0, nrow(Y))
  lambda<-lambda_fixed
  # if(lambda_fixed==TRUE){lambda<-log(nrow(X))*(ncol(X))}
  # if(lambda_fixed!=TRUE){lambda<-log(nrow(X))*(ncol(X))}
  if(offset==TRUE){
    O<-matrix(rep(log(prep_pln$Offset),p),ncol=p,byrow = FALSE)
  }
  if(offset==FALSE){
    O<-matrix(0,nrow=n,ncol=p)
  }

  X<-model.matrix(~.,as.data.frame(X))
  PLN_algo_name<- macth_nloptr_algo_name_with_algo_name_PLN(optim_method)
  res_PLN<-PLNmodels::PLN(Abundance~.-Offset+offset((O[,1])),data=prep_pln,control=PLN_param(backend = "nlopt",config_optim=list(algorithm=PLN_algo_name,maxeval=1000)))

  # Offset_mat<-matrix(0,n,p)
  # Initialisation avec LM
  ## PLN component: fast multivariate LM (identical to compute_PLN_starting_point "LM")
  if(initialisation_method=="LM"){
    cat("***** Initialisation with LM and ",optim_method," algorithm for the optimizer *****", "\n")
  Starting_LM <- PLNmodels::compute_PLN_starting_point(Y=Y, X=X, O=O, w)
  B<-Starting_LM$B
  S<-Starting_LM$S
  M<-Starting_LM$M
  # Sigma<-Compute_Sigma(X,Y,O,M,S,B)
  # Omega<-chol2inv(chol(Sigma))
  # # Calcul init Omega
  S_hat<-diag(colSums(S^2))
  Omega_init<-chol2inv(chol((1/n)*(S_hat+t(M)%*%(M))))
  Omega<-Omega_init#res_PLN$model_par$Omega

  cat("*** Initialization is complete ***", "\n")

  } else  if(initialisation_method=="PLN"){
    cat("***** Initialisation with PLN and ",optim_method," algorithm for the optimizer *****", "\n")
  # Initialisation avec PLN
  # params_init_E_step<-c(as.vector(res_PLN$var_par$M),as.vector(res_PLN$var_par$S))
  # params_init_M_step<-c(as.vector(res_PLN$model_par$B),as.vector(res_PLN$model_par$Omega))
  B<-res_PLN$model_par$B
  Omega<-res_PLN$model_par$Omega
  S<-res_PLN$var_par$S
  M<-res_PLN$var_par$M

  cat("*** Initialization is complete ***", "\n")
  }
  else{
    cat("***** Choose initialisation method between LM or PLN *****", "\n")
  }
  X<-as.matrix(X)

  Y<-as.matrix(Y)
  d<-ncol(X)
  # param_initpln=param_init(variables =X[,-1],abondance =Y,offset= matrix(0,ncol=ncol(Y),nrow=nrow(Y)))
  # formatp=format_plnnew(abondance=Y,covariables=X,B_reg=param_initpln$B,Sigma =diag(1,p) ,M_reg=param_initpln$M,S_reg=param_initpln$S,O = matrix(0,ncol=ncol(Y),nrow=nrow(Y)),w =rep(1,nrow(Y)))
  # params_init_E_step<-c(as.vector(param_initpln$M),as.vector(param_initpln$S))
  params_init_E_step<-c(as.vector(M),as.vector(S))
  # v <- res_PLN$model_par$B
  # #c(0.44844978, -0.08463290, -0.31860861, 0.50380695, 0.96334493, -0.18933897,
  #        0.50788630, -0.02754789, 3.14563421)
  # v<-rep(3,9)
  params_init_M_step<-(as.vector(B))#params_init_M_step<-c(as.vector(B),as.vector((Omega)))
  # param_initpln$B
  # B<-v#matrix(v,3,3)
  # Omega<-diag(1,p)

  length_epsilon<-100
  E<-c()
  e1<-10
  E[1]<-e1
  for(t in 2:length_epsilon){
    E[t]=e1*(0.87)^(t-1)
  }
  iter<-0
  ELBO_epsilon <- rep(NA, length_epsilon)
  # Barre de progression
  pb <- progress_bar$new(
    format = "  Epsilon telescoping procedure [:bar] :percent",
    total = length_epsilon,
    clear = FALSE,
    width = 80
  )
  for(epsilon_val in E){
    # epsilon_val<-0
    # lambda<-0
    iter<-iter+1
    v_loglik_tmp<-PLN_vloglik(X,Y,O=O,M,S,B=B,Omega=Omega)
    # cat("Epsilon value for telescoping :\t",iter,epsilon_val,"\t","loglik :\t",v_loglik_tmp,"\n")
    ELBO_epsilon[iter]<-v_loglik_tmp
    ELBO <- c(v_loglik_tmp, rep(NA, max_it))
    for(i in 1:max_it){
      # function(X,Y,O,params,B,Omega,lambda, epsilon)
      objective_E_step_SICPLN_Wrap <- function(par) {
        objective_E_step_SICPLN(par, Y = Y, X = X, O = O, B = B, Omega = Omega, lambda = lambda, epsilon = epsilon_val)
      }
      # function(X,Y,O,params,B,Omega,lambda, epsilon)
      grad_E_step_SICPLN_Wrap <- function(par) {
        grad_E_step_SICPLN(par, Y = Y, X = X, O = O, B = B, Omega = Omega,lambda = lambda, epsilon = epsilon_val)
      }
      result_E_step <- nloptr::nloptr(
        x0 = params_init_E_step,
        eval_f = objective_E_step_SICPLN_Wrap,
        eval_grad_f = grad_E_step_SICPLN_Wrap,
        opts = list(
          algorithm = optim_method,
          maxeval   = 1000,
          maxtime  = 200,
          xtol_rel  = 1e-5,
          xtol_abs  =0,
          ftol_rel  = 1e-6,
          ftol_abs  = 1e-6

        )
      )
      # result_E_step <- optim(
      #   fn = objective_E_step_SICPLN,
      #   gr = grad_E_step_SICPLN,
      #   Y=Y,
      #   X=(X),
      #   O=O,
      #   par = params_init_E_step,
      #   B=B,
      #   Omega=Omega,
      #   lambda=lambda,
      #   epsilon=epsilon_val,
      #   method = optim_method,
      #   control = list(
      #     maxit = 1000
      #     # gradtol = 1e-8
      #   )
      # )
      M <- (matrix(result_E_step$solution[1:(n*p)],n,p))  # (n,p)
      S <-  matrix(result_E_step$solution[(n*p+1):length(result_E_step$solution)],n,p)
      #cat("E_step : ",result_E_step$value,"\n")
      # result_E_step$value
      params_init_E_step<-c(as.vector(M),as.vector(S))

      # M step
      # function(X,Y,O,paramsM,M,S,lambda, epsilon)
      objective_M_step_SICPLN_Wrap <- function(par) {
        objective_M_step_SICPLN(par, Y = Y, X = X, O = O, M = M, S = S, Omega=Omega,lambda = lambda, epsilon = epsilon_val)
      }

      # Avec mise à jour B non analytique
      # function(X,Y,O,paramsM,M,S,lambda, epsilon)
      grad_M_step_SICPLN_Wrap <- function(par) {
        grad_M_step_SICPLN(par, Y = Y, X = X, O = O, M = M, S = S,Omega=Omega,
                           lambda = lambda, epsilon = epsilon_val)
      }
      result_M_step <- nloptr::nloptr(
        x0 = params_init_M_step,
        eval_f = objective_M_step_SICPLN_Wrap,
        eval_grad_f = grad_M_step_SICPLN_Wrap,
        opts = list(
          algorithm = optim_method,
          maxeval   = 1000,
          maxtime  = 200,
          xtol_rel  = 1e-5,
          xtol_abs  =0,
          ftol_rel  =1e-6,
          ftol_abs  =1e-6
        )
      )
      # result_M_step <- optim(
      #   fn = objective_M_step_SICPLN,
      #   gr = grad_M_step_SICPLN,
      #   Y=Y,
      #   X=(X),
      #   O=O,
      #   par = params_init_M_step,
      #   M=M,
      #   S=S,
      #   lambda=lambda,
      #   epsilon=epsilon_val,
      #   method = optim_method,
      #   control = list(
      #     maxit = 1000
      #   ))

      B <- matrix(result_M_step$solution[1:(d*p)],d,p)  # (n,p)
      #Compute_Sigma(X,Y,O=0,M,S,B)
      # Sigma<-Compute_Sigma(X,Y,O=0,M,S,B)
      # Omega<-MASS::ginv(Sigma)
      S_hat<-diag(colSums(S^2))
      Omega_cal<-chol2inv(chol((1/n)*(S_hat+t(M)%*%(M))))
      Omega<-Omega_cal
      Sigma<-chol2inv(chol(Omega))
      # Omega<- solve( Sigma)#res_PLN$model_par$Omega#matrix(result_M_step$par[(d*p+1):length(result_M_step$par)],p,p)
      params_init_M_step<-c(as.vector(B))#params_init_M_step<-c(as.vector(B),as.vector(Omega))
      #cat("M_step : ",result_M_step$value,"\n")
      ELBO[i + 1]<-PLN_vloglik(X,Y,O=O,M,S,B=B,Omega=Omega)
      diff_ELBO   <- abs(ELBO[i] - ELBO[i + 1])/(abs(ELBO[i]) + 1e-8)
      convergence<-abs(result_E_step$objective-result_M_step$objective)
      # if(convergence<=1e-5){
        if(diff_ELBO<=1e-8){
        # cat("epsilon telescoping : ",epsilon_val,"Convergence avant maxit\n")
        break
          }
    }
    # Affichage de la barre de progression
    # Sys.sleep(0.05)
    pb$tick()

  }
  # Post traitement
  # Y
  if (!is.null(colnames(Y))) {
    colnames(Omega)<-colnames(Y)
    rownames(Omega)<-colnames(Y)
    colnames(Sigma)<-colnames(Y)
    rownames(Sigma)<-colnames(Y)
    colnames(M)<-colnames(Y)
    colnames(S)<-colnames(Y)
    colnames(B)<-colnames(Y)
  }
  if(!is.null(rownames(Y))) {
    rownames(M)<-rownames(Y)
    rownames(S)<-rownames(Y)
  }
  # X
  if (!is.null(colnames(X))) {
    rownames(B)<-colnames(X)#c("Intercept",colnames(X))
    # rownames(Omega)<-colnames(Y)
    # colnames(Sigma)<-colnames(Y)
    # rownames(Sigma)<-colnames(Y)
    # colnames(M)<-colnames(Y)
    # colnames(S)<-colnames(Y)
  }
  # if(!is.null(rownames(X))) {
  #   rownames(B)<-rownames(X)
  #   # rownames(S)<-rownames(X)
  # }

  loglik<-PLN_vloglik(X,Y,O=O,M,S,B,Omega)
  B[abs(B)<=1e-5]<-0

  BIC_SICPLN<-loglik-0.5*log(n)*(res_PLN$nb_param)
  approxi_BIC_SICPLN<-loglik-0.5*log(n)*(res_PLN$nb_param-length(which(c(B[-1,])==0)))
  model_par<-list(B=B,Omega=Omega,Sigma=Sigma)
  var_par<-list(M=M,S=S,S2=S^2)
  data<-list(X=X,Y=Y,O=O)
  SICPLN_prediction<-(exp(O+X%*%B+M+S^2/2))
  q<-p
  entropy<-0.5 * (n * q * log(2*pi*exp(1)) + sum(log(S^2)))
  ICL<-BIC_SICPLN-entropy
  ICL_approx<-approxi_BIC_SICPLN-entropy
  return(list(ELBO_epsilon=ELBO_epsilon,ELBO=na.omit(ELBO), model_par=model_par,var_par=var_par,data=data,elbo=loglik,loglik=loglik,BIC_SICPLN=BIC_SICPLN,approxi_BIC_SICPLN=approxi_BIC_SICPLN,SICPLN_prediction=SICPLN_prediction,ICL=ICL,ICL_approx=ICL_approx,res_PLN=res_PLN))
}
SICPLN<-function(X=X,Y=Y,offset=FALSE,lambda_fixed=(log(nrow(X))*ncol(X)),optim_method="BFGS",max_it=200,length_lambda=100,grid_search=FALSE,initialisation_method="PLN"){
  cat(lambda_fixed, "\n")
  lambda_fixed<-lambda_fixed
  prep_pln<-PLNmodels::prepare_data(Y,X)
  p<-ncol(prep_pln$Abundance)
  n<-nrow(prep_pln$Abundance)
  if(offset==TRUE){
    O<-matrix(rep(log(prep_pln$Offset),p),ncol=p,byrow = FALSE)
  }
  if(offset==FALSE){
    O<-matrix(0,nrow=n,ncol=p)
  }

  if(grid_search==FALSE){
    solution<-SICPLN_NLOPTR(X=X,Y=Y,offset=offset,lambda_fixed=lambda_fixed,optim_method=optim_method,max_it=200,initialisation_method=initialisation_method)
    solution<-solution
    return(solution)
  }
  if(grid_search==TRUE){
    PLN_algo_name<- macth_nloptr_algo_name_with_algo_name_PLN(optim_method)
    solution_pln<-PLNmodels::PLN(Abundance~.-Offset+offset((O[,1])),data=prep_pln,control=PLN_param(backend = "nlopt",config_optim=list(algorithm=PLN_algo_name,maxeval=1000)))
    BIC_PLN<-solution_pln$BIC
    # Creation de la grille
    length_lambda=length_lambda
    # lambda_min=0#(log(nrow(X))*(ncol(X)*ncol(Y))+ncol(Y))/10
    # lambda_max=3*(log(nrow(X))*(ncol(X)))
    # Charger le package parallel
    # library(parallel)

    # Calculer la valeur log(n) une seule fois pour l'utiliser plus tard
    log_n <- log(nrow(X))

    # Définir les valeurs spécifiques de lambda
    lambda_min <- 0
    lambda_max <- log_n * ncol(X)#*ncol(Y)
    lambda_sic1 <- log_n/2
    lambda_sic2 <- log_n
    lambda_sic3 <- log_n * ncol(X) / 2
    lambda_sic4 <- log_n * ncol(X) * ncol(Y)

    # Créer la séquence de valeurs lambda
    # creation de la grille
    # Paramètres
    # lambda_sic4       # Valeur initiale (à t = 0)
    k <- 0.05        # Taux de décroissance (plus k est grand, plus la décroissance est rapide)
    t <- seq(0, length_lambda, by = 1)  # Séquence de temps allant de 0 à 100

    # Calcul de la décroissance exponentielle
    seq_lambda <-c(lambda_sic4 * exp(-k * t),lambda_min)

    ############################## old
    # seq_lambda <- c(0, lambda_sic1, lambda_sic2, lambda_sic3, lambda_sic4,
    #                 10^(seq(0, log(lambda_max), length.out = length_lambda)))

    # Créer un cluster avec un nombre de cœurs spécifié
    n_cores <- detectCores() - 2  # Utiliser tous les cœurs sauf un pour ne pas saturer le système
    cl <- makeCluster(n_cores)
    # clusterExport(cl, list("SICPLN_optim", "X", "Y"))
    # Exporter toutes les fonctions et objets de l'environnement global vers le cluster
    clusterExport(cl, varlist = ls(envir = .GlobalEnv))  #
    clusterEvalQ(cl, library(PLNmodels))
    # Initialiser une liste pour stocker les résultats et un vecteur pour le BIC
    res_fs <- list()
    BIC_lambda <- numeric(length(seq_lambda))

    # Paralléliser l'exécution de la boucle avec parLapply
    res_fs <- parLapply(cl, seq_lambda, function(lambda) {
      # Exécuter l'optimisation pour chaque valeur de lambda
      # result <- SICPLN_optim(X = X, Y = Y, offset = offset, lambda_fixed = lambda,
      #                        optim_method = optim_method, max_it = 200)
      result <- tryCatch({
        SICPLN_NLOPTR(X = X, Y = Y, offset = offset, lambda_fixed = lambda,
                     optim_method = optim_method, max_it = 200,initialisation_method=initialisation_method)
      }, error = function(e) {
        message("Erreur lors de l'optimisation avec lambda = ", lambda)
        return(NULL)  # Retourner NULL en cas d'erreur
      })
      # Retourner le résultat et la valeur du BIC
      list(result = result, BIC = safe_get(result[10]),BIC_approximer = safe_get(result[11]))  # Assumons que le BIC est à la 10e position
    })

    # Extraire les BIC et les résultats
    BIC_lambda <- sapply(res_fs, function(x) x$BIC)
    BIC_lambda_approximer <- sapply(res_fs, function(x) x$BIC_approximer)
    res_fs <- lapply(res_fs, function(x) x$result)

    # Nommer les éléments dans la liste des résultats
    names(res_fs) <- paste0("SICPLN_lambda_", round(1:length(seq_lambda), 0))

    # Trouver le meilleur lambda basé sur le BIC (plus bas BIC)
    best_lambda_indice <- which(BIC_PLN <= min(unlist(BIC_lambda)))
    best_lambda_indice_approximer <- which(BIC_PLN <= min(unlist(BIC_lambda_approximer)))

    # Créer la solution finale
    solution <- list(solution = res_fs, seq_lambda = seq_lambda,
                     BIC_all_lambda = unlist(BIC_lambda), best_lambda_indice = best_lambda_indice,best_lambda_indice_approximer=best_lambda_indice_approximer,BIC_lambda_approximer=unlist(BIC_lambda_approximer))

    # Arrêter le cluster après l'exécution
    stopCluster(cl)

  }
  return(solution)
}

macth_nloptr_algo_name_with_algo_name_PLN<-function(algo){
  Supported_algo<-c("NLOPT_LD_LBFGS","NLOPT_LD_MMA","NLOPT_LD_CCSAQ","NLOPT_LD_VAR1","NLOPT_LD_VAR2","NLOPT_LD_TNEWTON","NLOPT_LD_TNEWTON_RESTART","NLOPT_LD_TNEWTON_PRECOND","NLOPT_LD_TNEWTON_PRECOND_RESTART")
  if(algo=="NLOPT_LD_LBFGS"){algo_PLN<-"LBFGS"}
  if(algo=="NLOPT_LD_VAR1"){algo_PLN<-"VAR1"}
  if(algo=="NLOPT_LD_VAR2"){algo_PLN<-"VAR2"}
  if(algo=="NLOPT_LD_TNEWTON"){algo_PLN<-"TNEWTON"}
  if(algo=="NLOPT_LD_TNEWTON_RESTART"){algo_PLN<-"TNEWTON_RESTART"}
  if(algo=="NLOPT_LD_TNEWTON_PRECOND"){algo_PLN<-"TNEWTON_PRECOND"}
  if(algo=="NLOPT_LD_TNEWTON_PRECOND_RESTART"){algo_PLN<-"TNEWTON_PRECOND_RESTART"}
  if(algo=="NLOPT_LD_MMA"){algo_PLN<-"MMA"}
  if(algo=="NLOPT_LD_CCSAQ"){algo_PLN<-"CCSAQ"}
  if(!(algo%in%Supported_algo)){cat("Choose an algorithm name in : \n ",Supported_algo)}
  return(algo_PLN)
}
