#' Safely access an element by index
#'
#' Attempts to extract the element at position \code{index} from \code{x},
#' returning \code{NA} instead of raising an error if the access fails
#' (e.g. an out-of-bounds index, or \code{x} of an incompatible type).
#' Useful for extracting elements from a list of possibly-failed model
#' fits or heterogeneous results without wrapping every call site in its
#' own \code{tryCatch()}.
#'
#' @details
#'
#'
#' @param x An object to index, typically a vector or list.
#' @param index Index used to access an element of \code{x} (e.g. a
#'   single integer position, but anything valid for \code{x[index]} is
#'   accepted).
#'
#' @return The element \code{x[index]}, or \code{NA} if accessing it
#'   raises an error.
#'
#' @examples
#' safe_get(c(10, 20, 30), 2)
#' safe_get(c(10, 20, 30), 10)   # out-of-bounds: returns NA (base R default)
#' safe_get(list(a = 1, b = 2), "c")
#'
#' @export
safe_get <- function(x, index) {
  tryCatch({
    return(x[index])
  }, error = function(e) {
    return(NA)
  })
}
####################################
#' Ramanujan's approximation of the log-factorial function
#'
#' Computes an accurate approximation of \eqn{\log(n!)} using Ramanujan's
#' formula, which is much faster than exact computation for large
#' \code{n} and more numerically stable than computing \eqn{\log(n!)} via
#' \code{log(factorial(n))} (which overflows for even moderately large
#' \code{n}). Typically used for evaluating count-data log-likelihoods
#' (e.g. Poisson, zero-inflated Poisson) where a
#' \eqn{\sum_{i,j}\log(Y_{ij}!)} term is required.
#'
#' @details
#' The approximation is:
#' \deqn{\log(n!) \approx n\log(n) - n + \frac{1}{6}\log\!\left(8n^3+4n^2+n+\frac{1}{30}\right) + \frac{1}{2}\log(\pi)}
#'
#'
#' @param n Numeric vector of non-negative integers (or values coercible
#'   to non-negative integers) for which to compute the log-factorial.
#'   Vectorized: any number of elements is supported, including zeros.
#'
#' @return A numeric vector of the same length as \code{n}, giving
#'   \eqn{\log(n!)} (approximated) element-wise.
#'
#' @references
#' Ramanujan, S. (1988). The lost notebook and other unpublished papers.
#' Springer.
#'
#' @examples
#' logfactorial(0)              # 0
#' logfactorial(5)              # approx log(120) = 4.787492
#' log(factorial(5))            # exact, for comparison
#' logfactorial(c(0, 1, 5, 20)) # vectorized, handles zero
#'
#' @export
logfactorial <- function(n) {
  n[n == 0] <- 1 ## 0! = 1!
  n * log(n) - n + log(8 * n^3 + 4 * n^2 + n + 1 / 30) / 6 + log(pi) / 2
}

##################################################
#' Compute the estimated residual covariance matrix Sigma in PLN models
#'
#' \eqn{\hat\Sigma = n^{-1}(M^\top M + S_+)}, combining the variational
#' means \code{M} and the accumulated variational variances \code{S},
#' as used in the M-step of a variational EM algorithm for PLN models.
#'
#' @details
#' The estimator is:
#' \deqn{\hat\Sigma = \frac{1}{n}\left(M^\top M + \mathrm{diag}\left(\sum_{i=1}^n S_i^2\right)\right)}
#'
#' where \eqn{n} is the number of observations.
#'
#' @param M Numeric matrix of dimension \eqn{n \times p} of variational
#'   means.
#' @param S Numeric matrix of dimension \eqn{n \times p} of variational
#'   standard deviations (see Details for the squaring convention used
#'   here).
#'
#' @return A numeric matrix of dimension \eqn{p \times p}, the estimated
#'   covariance matrix \eqn{\hat\Sigma}.
#'
#' @examples
#' n <- 50; p <- 5
#' M <- matrix(runif(n * p, 0.05, 0.2), n, p)
#' S <- matrix(runif(n * p, 0.05, 0.2), n, p)
#' Compute_Sigma(X, Y, O, M, S, B)
#'
#' @export
Compute_Sigma <- function(X, Y, O, M, S, B) {
  n <- nrow(X)
  S2_bar <- colSums(S^2)
  MtM <- t(M) %*% M
  (MtM + diag(S2_bar)) / n
}

########################################################
#' Variational log-likelihood of the PLN model
#'
#' Computes the value of the approximated (variational) log-likelihood
#' \eqn{J(Y;\psi,\theta)} of a Poisson Log-Normal (PLN) model at a given
#' set of variational parameters (\code{M}, \code{S}) and model
#' parameters (\code{B}, \code{Omega}).
#'
#' @details
#'
#'
#' @param X Numeric matrix of dimension \eqn{n \times d} of covariates
#'   (design matrix), including the intercept if applicable.
#' @param Y Numeric matrix of dimension \eqn{n \times p} of observed
#'   counts.
#' @param O Numeric matrix of dimension \eqn{n \times p} of offsets.
#' @param M Numeric matrix of dimension \eqn{n \times p} of variational
#'   means.
#' @param S Numeric matrix of dimension \eqn{n \times p} of variational
#'   standard deviations (squared internally to obtain variances).
#' @param B Numeric matrix of dimension \eqn{d \times p} of regression
#'   coefficients.
#' @param Omega Numeric matrix of dimension \eqn{p \times p}, the current
#'   estimate of the precision matrix (inverse covariance).
#'
#' @return A single numeric value: the variational log-likelihood
#'   \eqn{J}.
#'
#' @seealso
#'
#' @examples
#' n <- 50; p <- 5; d <- 2
#' X <- cbind(1, rnorm(n))
#' B <- matrix(rnorm(d * p), nrow = d)
#' Y <- matrix(rpois(n * p, 5), n, p)
#' O <- matrix(0, n, p)
#' M <- X %*% B + matrix(rnorm(n * p, sd = 0.1), n, p)
#' S <- matrix(runif(n * p, 0.05, 0.2), n, p)
#' Omega <- diag(p)
#' PLN_vloglik(X, Y, O, M, S, B, Omega)
#'
#' @export
PLN_vloglik <- function(X, Y, O, M, S, B, Omega) {
  n <- nrow(Y)
  p <- ncol(Y)
  d <- ncol(X)
  params <- c(M, S)
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  O <- as.matrix(O)

  w <- rep(1, n)
  M <- matrix(params[1:(n * p)], n, p)
  wbar <- sum(w)
  S <- matrix(params[(n * p + 1):length(params)], n, p)

  B <- B

  S2 <- S^2
  Z <- O + X %*% B + M
  A <- exp(Z + 0.5 * S2)
  Omega <- Omega

  objective <- 0.5 * log(det(Omega)) +
    rowSums(Y * Z - A + 0.5 * log(S2)) -
    0.5 * rowSums((M %*% Omega) * M + S2 * diag(Omega))

  ji <- -rowSums(logfactorial(Y)) + objective
  ji <- 0.5 * p + ji
  objective <- sum(ji)

  return(objective)
}
###################################################
#' Penalized negative variational objective for the E-step of SICPLN
#'
#' Computes the negated, SIC-penalized variational log-likelihood used as
#' the objective minimized during the variational E-step of
#' \code{SICPLN}/\code{SICZIPLN}-type sparse PLN models. Combines the
#' variational log-likelihood (\code{\link{PLN_vloglik}}) with a smooth
#' approximation to an \eqn{\ell_0} penalty on the regression
#' coefficients \code{B}, plus a model-complexity correction term, and
#' returns the negative of the result so that the function can be handed
#' directly to a minimizer (e.g. \code{nloptr}).
#'
#' @details
#'
#' The function returns \code{-objective} (not \code{objective}) so it
#' can be passed directly to a minimization routine that expects to
#' minimize its input, such as \pkg{nloptr}-based optimizers used
#' elsewhere in the E-step.
#'
#' @param X Numeric matrix of dimension \eqn{n \times d} of covariates
#'   (design matrix), including the intercept if applicable.
#' @param Y Numeric matrix of dimension \eqn{n \times p} of observed
#'   counts.
#' @param O Numeric matrix of dimension \eqn{n \times p} of offsets.
#' @param params Numeric vector of length \eqn{2np}, the concatenation of
#'   the variational means \code{M} and variational standard deviations
#'   \code{S}, each flattened column-major and stacked (\code{M} first,
#'   then \code{S}), as expected by an optimizer operating on a flat
#'   parameter vector.
#' @param B Numeric matrix of dimension \eqn{d \times p} of regression
#'   coefficients (held fixed during the E-step; only \code{M} and
#'   \code{S}, packed in \code{params}, are being optimized over).
#' @param Omega Numeric matrix of dimension \eqn{p \times p}, the current
#'   estimate of the precision matrix (inverse covariance), held fixed
#'   during the E-step.
#' @param lambda Numeric scalar, the regularization (penalty) parameter
#'   controlling the strength of the SIC penalty on \code{B} and the
#'   magnitude of the \code{nb_params} correction term.
#' @param epsilon Numeric scalar, the smoothing parameter of the SIC
#'   penalty; smaller values make the penalty a closer (but less smooth)
#'   approximation to an \eqn{\ell_0} penalty on \code{B}.
#'
#' @return A single numeric value: the negative of the SIC-penalized
#'   variational objective, suitable for direct use with a minimizer.
#'
#' @seealso \code{\link{PLN_vloglik}} for the unpenalized variational
#'   log-likelihood used internally.
#'
#' @examples
#' n <- 50; p <- 5; d <- 2
#' X <- cbind(1, rnorm(n))
#' B <- matrix(rnorm(d * p), nrow = d)
#' Y <- matrix(rpois(n * p, 5), n, p)
#' O <- matrix(0, n, p)
#' M <- X %*% B + matrix(rnorm(n * p, sd = 0.1), n, p)
#' S <- matrix(runif(n * p, 0.05, 0.2), n, p)
#' Omega <- diag(p)
#' params <- c(as.vector(M), as.vector(S))
#' objective_E_step_SICPLN(X, Y, O, params, B, Omega, lambda = log(n), epsilon = 1e-4)
#'
#' @export
objective_E_step_SICPLN <- function(X, Y, O, params, B, Omega, lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d <- ncol(X)
  nb_params <- (lambda / 2) * ((((p + 1) * p) / 2) + (p))
  M <- matrix(params[1:(n * p)], n, p)
  S <- matrix(params[(n * p + 1):length(params)], n, p)
  objective <- PLN_vloglik(X, Y, O, M, S, B, Omega)

  SIC_penalty <- (lambda / 2) * (B^2 / (B^2 + epsilon^2))
  SIC_penalty[1, ] <- 0
  objective <- objective - sum(SIC_penalty) - nb_params

  return(-objective)
}
#####################################"
#' Gradient of the penalized negative variational objective (E-step, SICPLN)
#'
#' Computes the gradient, with respect to the variational
#' parameters \code{M} (means) and \code{S} (standard deviations), of the
#' negative variational objective minimized during the E-step of
#' \code{SICPLN}. Meant
#' to be passed as the \code{gr} argument to a gradient-based optimizer
#' (e.g. an \pkg{nloptr} algorithm such as \code{"NLOPT_LD_LBFGS"}) so
#' that the E-step does not fall back to numerical differentiation.
#'
#' @details
#'
#' @param X Numeric matrix of dimension \eqn{n \times d} of covariates
#'   (design matrix), including the intercept if applicable.
#' @param Y Numeric matrix of dimension \eqn{n \times p} of observed
#'   counts.
#' @param O Numeric matrix of dimension \eqn{n \times p} of offsets.
#' @param params Numeric vector of length \eqn{2np}, the concatenation of
#'   the variational means \code{M} and variational standard deviations
#'   \code{S}, each flattened column-major and stacked (\code{M} first,
#'   then \code{S}), matching the layout expected by
#'   \code{\link{objective_E_step_SICPLN}}.
#' @param B Numeric matrix of dimension \eqn{d \times p} of regression
#'   coefficients, held fixed.
#' @param Omega Numeric matrix of dimension \eqn{p \times p}, the current
#'   estimate of the precision matrix, held fixed.
#' @param lambda Numeric scalar, the regularization parameter. Unused in
#'   this gradient (the SIC penalty does not depend on \code{M} or
#'   \code{S}), but kept for a signature consistent with
#'   \code{\link{objective_E_step_SICPLN}} as required by optimizers that
#'   expect matching \code{fn}/\code{gr} argument lists.
#' @param epsilon Numeric scalar, the SIC penalty smoothing parameter.
#'   Unused in this gradient, for the same reason as \code{lambda}.
#'
#' @return A numeric vector of length \eqn{2np}: the gradient with
#'   respect to \code{M} (first \eqn{np} entries, column-major) followed
#'   by the gradient with respect to \code{S} (last \eqn{np} entries,
#'   column-major), matching the layout of \code{params}.
#'
# @seealso \code{\link{objective_E_step_SICPLN}} for the corresponding
#   objective function; \code{\link{PLN_vloglik}} for the unpenalized
#   variational log-likelihood.
#'
#' @examples
#' n <- 50; p <- 5; d <- 2
#' X <- cbind(1, rnorm(n))
#' B <- matrix(rnorm(d * p), nrow = d)
#' Y <- matrix(rpois(n * p, 5), n, p)
#' O <- matrix(0, n, p)
#' M <- X %*% B + matrix(rnorm(n * p, sd = 0.1), n, p)
#' S <- matrix(runif(n * p, 0.05, 0.2), n, p)
#' Omega <- diag(p)
#' params <- c(as.vector(M), as.vector(S))
#' grad_E_step_SICPLN(X, Y, O, params, B, Omega, lambda = log(n), epsilon = 1e-4)
#'
#' @export
grad_E_step_SICPLN <- function(X, Y, O, params, B, Omega, lambda, epsilon) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  n <- nrow(Y)
  p <- ncol(Y)
  d <- ncol(X)
  O <- as.matrix(O)

  M <- matrix(params[1:(n * p)], n, p)
  S <- matrix(params[(n * p + 1):length(params)], n, p)
  S2 <- S^2

  B <- B
  Z <- O + X %*% B + M
  A <- exp(Z + 0.5 * S2)
  Omega <- Omega

  grad_M <- (M %*% Omega + A - Y)
  grad_S <- -1/2 * (S^-1 - A - rep(1, n) %*% t(diag(Omega)))

  gradient <- c(as.vector(grad_M), as.vector(grad_S))
  return(gradient)
}
#####################################################################
grad_E_step_SICPLN_old_non_correct_gradient<- function(X,Y,O,params,B,Omega,lambda, epsilon) {
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
####################################################################
#' Objective Function for the Penalized M-Step in SICPLN
#'
#' Computes the penalized objective function used during the M-step of the
#' SICPLN estimation procedure. The function evaluates the variational
#' log-likelihood and applies a SIC penalty to the
#' regression coefficient matrix \code{B}, together with a model-complexity
#' penalty based on the number of estimated parameters.
#'
#' The sparsity penalty is defined as
#' \deqn{
#' \frac{\lambda}{2} \frac{B_{jk}^2}{B_{jk}^2 + \epsilon^2},
#' }
#' where \code{lambda} controls the strength of the penalty and
#' \code{epsilon} controls its smoothness. The first row of \code{B} is not
#' penalized, corresponding to an unpenalized intercept.
#'
#' The function returns the negative penalized objective, as required by
#' optimization routines that minimize their objective function.
#'
#' @param X A numeric matrix of covariates
#' @param Y A numeric matrix of multivariate count responses
#' @param O A matrix of offset
#' @param paramsM A numeric vector containing the current values of the
#'   regression parameters \code{B}.
#' @param Omega A matrix corresponding to the precision matrix
#' @param M A matrix of variational means.
#' @param S A matrix of variational scale parameters.
#' @param lambda A non-negative numeric value controlling the strength of the
#'   sparsity.
#' @param epsilon A positive numeric value controlling the smoothness of the
#'   sparsity penalty.
#'
#' @return A numeric scalar corresponding to the negative penalized objective
#'   function.
#'
#' @seealso
#' \code{\link{PLN_vloglik}}
#'
#' @keywords internal optimization penalization SICPLN
#'
objective_M_step_SICPLN <- function(X, Y, O, paramsM, Omega, M, S, lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d <- ncol(X)

  nb_params <- (lambda / 2) * ((((p + 1) * p) / 2) + p)

  B <- matrix(paramsM[1:(d * p)], d, p)

  objective <- PLN_vloglik(X, Y, O, M, S, B, Omega)

  SIC_penalty <- (lambda / 2) * (B^2 / (B^2 + epsilon^2))
  SIC_penalty[1, ] <- 0

  objective <- objective - sum(SIC_penalty) - nb_params

  return(-objective)
}

##############################################################################
#' Gradient of the Penalized M-Step Objective in SICPLN
#'
#' Computes the gradient of the penalized M-step objective function with
#' respect to the regression coefficient matrix \code{B} in the SICPLN
#' estimation procedure.
#'
#' @param X A numeric matrix of covariates
#' @param Y A numeric matrix of multivariate count responses
#' @param O A numeric matrix of offsets, with the same number of rows and
#'   columns as \code{Y}.
#' @param paramsM A numeric vector containing the current regression
#'   parameters \code{B}.
#' @param Omega A precision matrix.
#' @param M A numeric matrix of variational means.
#' @param S A numeric matrix of variational scale parameters.
#' @param lambda A non-negative numeric value controlling the strength of the
#'   SIC penalty.
#' @param epsilon A positive numeric value controlling the smoothness of the
#'   SIC penalty.
#'
#' @return A numeric vector containing the gradient of the penalized objective
#'   with respect to the regression coefficients in \code{B}.
#'
#' @seealso
#' \code{\link{objective_M_step_SICPLN}}
#'
#' @keywords internal optimization gradient penalization SICPLN
#'
grad_M_step_SICPLN <- function(X, Y, O, paramsM, Omega, M, S, lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)

  Y <- as.matrix(Y)
  X <- as.matrix(X)
  d <- ncol(X)
  O <- as.matrix(O)

  w <- rep(1, n)

  B <- matrix(paramsM[1:(d * p)], d, p)

  S2 <- S^2
  Z <- O + X %*% B + M
  A <- exp(Z + 0.5 * S2)

  SIC_deriv <- (lambda / 2) *
    (2 * (B * epsilon^2) / (B^2 + epsilon^2)^2)

  SIC_deriv[1, ] <- 0

  grad_B <- t(X) %*% (A - Y) + SIC_deriv

  return(as.vector(grad_B))
}

##############################################################################
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

##########################################################################
#' SICPLN Estimation Using NLOPTR Optimization
#'
#' Fits a Sparse Poisson Lognormal (SPLN) model using a penalized
#' variational inference procedure with SIC regularization.
#' The model parameters are estimated through alternating E- and M-steps,
#' where both steps are solved numerically using the \code{nloptr} package.
#'
#' \itemize{
#' \item an E-step, which updates the variational parameters \code{M} and
#' \code{S};
#' \item an M-step, which updates the regression coefficient matrix
#' \code{B}.
#' }
#'
#' The precision matrix \code{Omega} is updated after each E-step from the
#' current variational parameters. The optimization is repeated until the
#' relative change in the variational lower bound (ELBO) falls below the
#' convergence threshold or the maximum number of iterations is reached.
#'
#' Two initialization strategies are available. With
#' \code{initialisation_method = "LM"}, the starting values are obtained
#' using a fast multivariate linear-model initialization. With
#' \code{initialisation_method = "PLN"}, the parameters are initialized
#' from a standard Poisson Lognormal model fitted using
#' \pkg{PLNmodels}.
#'
#' The final model is evaluated using the variational log-likelihood,
#' a BIC-type criterion, and ICL-type criteria. Predictions of the expected
#' abundance are also returned.
#'
#' @param X A data frame or matrix containing the covariates. An intercept is
#' added automatically to the design matrix.
#' @param Y A data frame or matrix containing the multivariate count response
#' data, with observations in rows and response variables in columns.
#' @param offset Logical. If \code{TRUE}, an offset is included in the model.
#' The offset is obtained from the data prepared by
#' \code{\link[PLNmodels]{prepare_data}}. Defaults to \code{FALSE}.
#' @param lambda_fixed Numeric value specifying the regularization parameter
#' \eqn{\lambda}. By default, it is set to
#' \code{log(nrow(X))}.
#' @param optim_method Character string specifying the optimization algorithm
#' passed to \code{\link[nloptr]{nloptr}}. Defaults to \code{"NLOPT_LD_LBFGS"}.
#' The algorithm name is converted to the corresponding NLOPTR algorithm
#' using the internal function
#' \code{macth_nloptr_algo_name_with_algo_name_PLN}.
#' @param max_it Integer specifying the maximum number of EM iterations
#' performed for each value of the epsilon-telescoping sequence.
#' Defaults to \code{200}.
#' @param initialisation_method Character string specifying the initialization
#' strategy. Supported values are \code{"LM"} and \code{"PLN"}.
#' \code{"LM"} uses a fast multivariate linear-model initialization,
#' whereas \code{"PLN"} initializes the parameters from a standard
#' Poisson Lognormal model fitted with \pkg{PLNmodels}.
#'
#' @return A list containing the following elements:
#' \describe{
#' \item{ELBO_epsilon}{
#' The variational lower bound evaluated at the beginning of each
#' epsilon-telescoping step.
#' }
#' \item{ELBO}{
#' The ELBO values obtained during the final epsilon step, excluding
#' missing values.
#' }
#' \item{model_par}{
#' A list containing the estimated model parameters:
#' \code{B}, the regression coefficient matrix;
#' \code{Omega}, the estimated precision matrix; and
#' \code{Sigma}, the corresponding covariance matrix.
#' }
#' \item{var_par}{
#' A list containing the variational parameters \code{M}, \code{S}, and
#' \code{S2 = S^2}.
#' }
#' \item{data}{
#' A list containing the matrices used for model fitting:
#' \code{X}, \code{Y}, and \code{O}.
#' }
#' \item{elbo}{
#' The final variational log-likelihood (ELBO).
#' }
#' \item{loglik}{
#' The final variational log-likelihood.
#' }
#' \item{BIC_SICPLN}{
#' The BIC-type model selection criterion computed from the final
#' log-likelihood and the number of parameters of the fitted PLN model.
#' }
#' \item{approxi_BIC_SICPLN}{
#' An approximate BIC-type criterion in which the number of zero
#' regression coefficients is subtracted from the parameter count.
#' }
#' \item{SICPLN_prediction}{
#' Matrix of predicted expected abundances,
#' \eqn{\exp(O + XB + M + S^2/2)}.
#' }
#' \item{ICL}{
#' The ICL criterion obtained by subtracting the variational entropy
#' from \code{BIC_SICPLN}.
#' }
#' \item{ICL_approx}{
#' The approximate ICL criterion based on
#' \code{approxi_BIC_SICPLN}.
#' }
#' \item{res_PLN}{
#' The fitted standard PLN model returned by
#' \pkg{PLNmodels}, used for initialization and parameter
#' information.
#' }
#' }
#'
#' @details
#' The SIC penalty is progressively sharpened through the sequence
#' \code{epsilon = 10 * 0.87^(t-1)}, for \code{t = 1, ..., 100}.
#' This continuation strategy is intended to facilitate numerical
#' optimization while approaching the non-smooth sparsity-inducing penalty.
#'
#' At each epsilon value, the E-step and M-step are alternated until the
#' relative change in the ELBO is smaller than \code{1e-8}, or until
#' \code{max_it} iterations have been performed.
#'
#' The first row of the regression coefficient matrix \code{B}, corresponding
#' to the intercept, is not penalized by the SIC regularization.
#'
#' Regression coefficients with absolute values smaller than \code{1e-5}
#' are set to zero in the final model. These thresholded coefficients are
#' used when computing the approximate BIC criterion.
#'
#' @seealso
#' \code{\link{objective_M_step_SICPLN}},
#' \code{\link{grad_M_step_SICPLN}},
#' \code{\link{objective_E_step_SICPLN}},
#' \code{\link{grad_E_step_SICPLN}},
#' \code{\link{PLN_vloglik}}
#'
#' @importFrom nloptr nloptr
#' @importFrom progress progress_bar
#'
#' @keywords models optimization penalization sparse PoissonLognormal
#' @export
#'
SICPLN_NLOPTR<-function(X,Y,offset=FALSE,lambda_fixed=log(nrow(X)),optim_method="BFGS",max_it=200,initialisation_method="LM"){
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
  res_PLN<-PLNmodels::PLN(Abundance~X[,-1]+offset((O[,1])),data=prep_pln,control=PLN_param(backend = "nlopt",config_optim=list(algorithm=PLN_algo_name,maxeval=1000)))

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
  pb <- progress::progress_bar$new(
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

      M <- (matrix(result_E_step$solution[1:(n*p)],n,p))  # (n,p)
      S <-  matrix(result_E_step$solution[(n*p+1):length(result_E_step$solution)],n,p)
      # Mise à jour de Omega après le Estep car depend unique de M et S
      S_hat<-diag(colSums(S^2))
      Omega_cal<-chol2inv(chol((1/n)*(S_hat+t(M)%*%(M))))
      Omega<-Omega_cal
      Sigma<-chol2inv(chol(Omega))
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
#######################################################################
#' Sparse Poisson Lognormal Model with SIC Regularization
#'
#' Fits a Sparse Poisson Lognormal (SICPLN) model using variational
#' inference and a SIC penalty. The model is estimated
#' through an alternating E-step/M-step procedure implemented with numerical
#' optimization.
#'
#' The function can either fit a model for a user-specified value of the
#' regularization parameter \eqn{\lambda}, or perform a grid search over a
#' sequence of candidate values of \eqn{\lambda} and select the best model
#' according to a BIC-type criterion.
#'
#' When \code{grid_search = FALSE}, the function directly calls
#' \code{\link{SICPLN_NLOPTR}} using the specified value of \code{lambda_fixed}.
#'
#' When \code{grid_search = TRUE}, a sequence of candidate regularization
#' parameters is generated and a SICPLN model is fitted independently for
#' each value. The optimization is parallelized across available CPU cores.
#' The standard PLN model is also fitted and its BIC is used as a reference
#' for selecting the regularized model.
#'
#'
#' @param X A data frame or matrix containing the covariates.
#' @param Y A data frame or matrix containing multivariate count data.
#' @param offset Logical. If \code{TRUE}, an offset is included in the model.
#' The offset is extracted from the data using
#' \code{\link[PLNmodels]{prepare_data}}. Defaults to \code{FALSE}.
#' @param lambda_fixed Numeric value specifying the regularization parameter
#' \eqn{\lambda} when \code{grid_search = FALSE}. Defaults to
#' \code{log(nrow(X))}.
#' @param optim_method Character string specifying the optimization algorithm
#' used by \code{\link[nloptr]{nloptr}}. Defaults to \code{"NLOPT_LD_LBFGS"}.
#' @param max_it Integer specifying the maximum number of E/M iterations
#' performed for each value of the regularization parameter. Defaults to
#' \code{200}.
#' @param length_lambda Integer specifying the number of steps used to
#' construct the regularization-parameter grid when
#' \code{grid_search = TRUE}. Defaults to \code{100}.
#' @param grid_search Logical. If \code{FALSE}, fit a single SICPLN model
#' using \code{lambda_fixed}. If \code{TRUE}, fit a model for each value
#' in the automatically generated lambda grid. Defaults to \code{FALSE}.
#' @param initialisation_method Character string specifying the initialization
#' strategy passed to \code{\link{SICPLN_NLOPTR}}. Supported values are
#' \code{"LM"} and \code{"PLN"}. Defaults to \code{"PLN"}.
#'
#' @return
#' If \code{grid_search = FALSE}, returns the result of
#' \code{\link{SICPLN_NLOPTR}}, including the estimated model parameters,
#' variational parameters, ELBO, BIC and ICL criteria, and predicted
#' expected abundances.
#'
#' If \code{grid_search = TRUE}, returns a list containing:
#' \describe{
#' \item{solution}{
#' A named list containing the fitted SICPLN models for all values of
#' the lambda grid.
#' }
#' \item{seq_lambda}{
#' The sequence of regularization parameters considered during the
#' grid search.
#' }
#' \item{BIC_all_lambda}{
#' The BIC-type criterion obtained for each value of lambda.
#' }
#' \item{best_lambda_indice}{
#' The indices of the lambda values for which the SICPLN BIC is at least
#' as good as the BIC of the unpenalized PLN model.
#' }
#' \item{best_lambda_indice_approximer}{
#' The indices of the lambda values for which the approximate SICPLN
#' BIC is at least as good as the BIC of the unpenalized PLN model.
#' }
#' \item{BIC_lambda_approximer}{
#' The approximate BIC-type criterion obtained for each value of lambda.
#' }
#' }
#'
#' @details
#' The response and covariate data are first processed using
#' \code{\link[PLNmodels]{prepare_data}}. An intercept is subsequently added
#' to the covariate matrix through \code{\link[stats]{model.matrix}}.
#'
#' Each lambda value is optimized independently and the computations are
#' parallelized using the \pkg{parallel} package.
#'
#' @section Model:
#' The SICPLN model extends the Poisson Lognormal framework by introducing
#' sparsity in the regression coefficient matrix.
#'
#' @section Lambda selection:
#' When \code{grid_search = TRUE}, models fitted for different values of
#' \eqn{\lambda} can be compared using \code{BIC_all_lambda} or
#' \code{BIC_lambda_approximer}. The latter accounts for the sparsity of the
#' estimated regression coefficient matrix by reducing the effective number
#' of parameters according to the number of coefficients set to zero.
#'
#' @seealso
#' \code{\link{SICPLN_NLOPTR}},
#' \code{\link{objective_M_step_SICPLN}},
#' \code{\link{grad_M_step_SICPLN}},
#' \code{\link[PLNmodels]{PLN}}
#'
#' @importFrom parallel detectCores makeCluster clusterExport
#' clusterEvalQ parLapply stopCluster
#' @importFrom stats model.matrix
#' @importFrom PLNmodels prepare_data PLN
#'
#' @keywords models sparse penalization PoissonLognormal optimization
#' @export
#'
SICPLN<-function(X=X,Y=Y,offset=FALSE,lambda_fixed=(log(nrow(X))),optim_method="BFGS",max_it=200,length_lambda=100,grid_search=FALSE,initialisation_method="PLN"){
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
    # solution<-solution
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
############################################################################
#' Match an `nloptr` algorithm name to its `PLN` package equivalent
#'
#' Converts an optimization algorithm name as used by the \pkg{nloptr}
#' package (e.g. \code{"NLOPT_LD_LBFGS"}) into the corresponding short name
#' expected by the \pkg{PLN} package's optimizer configuration (e.g.
#' \code{"LBFGS"}). Useful when the same algorithm choice needs to be passed
#' to both \pkg{nloptr}-based routines and \pkg{PLN} routines within the
#' package's fitting functions.
#'
#' @details
#' Only gradient-based, derivative-requiring \pkg{nloptr} algorithms
#' (prefixed \code{NLOPT_LD_}) are supported, since these are the ones used
#' by the underlying optimization routines. If \code{algo} is not one of the
#' supported names, the list of supported algorithms is printed to the
#' console via \code{cat()}, and the function still returns \code{NULL} (as
#' \code{algo_PLN} is never assigned in that case) rather than raising an
#' error — callers should check the return value if input validation is not
#' guaranteed upstream.
#'
#' @param algo Character string giving the \pkg{nloptr} algorithm name. Must
#'   be one of \code{"NLOPT_LD_LBFGS"}, \code{"NLOPT_LD_MMA"},
#'   \code{"NLOPT_LD_CCSAQ"}, \code{"NLOPT_LD_VAR1"}, \code{"NLOPT_LD_VAR2"},
#'   \code{"NLOPT_LD_TNEWTON"}, \code{"NLOPT_LD_TNEWTON_RESTART"},
#'   \code{"NLOPT_LD_TNEWTON_PRECOND"}, or
#'   \code{"NLOPT_LD_TNEWTON_PRECOND_RESTART"}.
#'
#' @return A character string giving the corresponding \pkg{PLN} algorithm
#'   name (e.g. \code{"LBFGS"}). Returns \code{NULL} invisibly if \code{algo}
#'   is not among the supported names (a message listing the supported
#'   algorithms is printed instead).
#'
#' @examples
#' macth_nloptr_algo_name_with_algo_name_PLN("NLOPT_LD_LBFGS")
#' macth_nloptr_algo_name_with_algo_name_PLN("NLOPT_LD_TNEWTON_RESTART")
#'
#' @export
macth_nloptr_algo_name_with_algo_name_PLN <- function(algo) {
  Supported_algo <- c("NLOPT_LD_LBFGS", "NLOPT_LD_MMA", "NLOPT_LD_CCSAQ",
                      "NLOPT_LD_VAR1", "NLOPT_LD_VAR2", "NLOPT_LD_TNEWTON",
                      "NLOPT_LD_TNEWTON_RESTART", "NLOPT_LD_TNEWTON_PRECOND",
                      "NLOPT_LD_TNEWTON_PRECOND_RESTART")
  if (algo == "NLOPT_LD_LBFGS") { algo_PLN <- "LBFGS" }
  if (algo == "NLOPT_LD_VAR1") { algo_PLN <- "VAR1" }
  if (algo == "NLOPT_LD_VAR2") { algo_PLN <- "VAR2" }
  if (algo == "NLOPT_LD_TNEWTON") { algo_PLN <- "TNEWTON" }
  if (algo == "NLOPT_LD_TNEWTON_RESTART") { algo_PLN <- "TNEWTON_RESTART" }
  if (algo == "NLOPT_LD_TNEWTON_PRECOND") { algo_PLN <- "TNEWTON_PRECOND" }
  if (algo == "NLOPT_LD_TNEWTON_PRECOND_RESTART") { algo_PLN <- "TNEWTON_PRECOND_RESTART" }
  if (algo == "NLOPT_LD_MMA") { algo_PLN <- "MMA" }
  if (algo == "NLOPT_LD_CCSAQ") { algo_PLN <- "CCSAQ" }
  if (!(algo %in% Supported_algo)) { cat("Choose an algorithm name in : \n ", Supported_algo) }
  return(algo_PLN)
}

#' Heatmap of estimated coefficients with sparsity pattern
#'
#' Draws a tile heatmap of a coefficient matrix (e.g. regression
#' coefficients estimated by \code{SICZIPLN}/\code{SICPLN}), colored by
#' coefficient value and overlaid with a hatching pattern that flags
#' coefficients shrunk exactly to zero. Useful for visualizing, at a
#' glance, both the magnitude/sign of the estimated effects and the
#' sparsity pattern produced by the variable selection procedure.
#'
#' @details
#' The coefficient matrix is reshaped to long format internally via
#' \code{reshape2::melt()}: rows of \code{matrice_coefficient} are mapped
#' to the y-axis, columns to the x-axis. Each coefficient is classified as
#' \code{"Zero"} or \code{"Nonzero"} and rendered with
#' \code{ggpattern::geom_tile_pattern()}: zero coefficients are hatched
#' with a circle pattern, nonzero coefficients are left plain. Fill color
#' follows a diverging red-white-green gradient
#' (\code{ggplot2::scale_fill_gradient2()}) centered at zero, with
#' negative coefficients in shades of dark red and positive coefficients
#' in shades of dark green.
#'
#' This function requires the \pkg{ggpattern} package (for
#' \code{geom_tile_pattern()} and \code{scale_pattern_manual()}) in
#' addition to \pkg{ggplot2} and \pkg{reshape2}.
#'
#' @param matrice_coefficient Numeric matrix of coefficients to plot (e.g.
#'   a regression coefficient matrix, variables in rows and response
#'   columns in columns, or vice versa depending on \code{nom_axes}).
#'   Row and column names, if present, are used as tick labels.
#' @param nom_axes Character vector of exactly length 4 giving, in order:
#'   the y-axis label, the x-axis label, the color legend title (for the
#'   coefficient value gradient), and the plot title. Defaults to
#'   \code{c("Columns of Y", "Variables", "Coefficients value", "SICZIPLN")}.
#' @param grad_echelle Numeric vector of exactly length 2 giving the lower
#'   and upper bounds \code{c(min, max)} of the color scale for
#'   coefficient values. Defaults to \code{c(-5, 5)}. Values outside this
#'   range are clipped by \code{ggplot2::scale_fill_gradient2()}'s
#'   \code{limit} argument; adjust to the actual range of
#'   \code{matrice_coefficient} if coefficients fall outside
#'   \eqn{[-5, 5]}.
#'
#' @return A \code{ggplot} object representing the coefficient heatmap,
#'   which can be further customized (e.g. with additional \pkg{ggplot2}
#'   layers or themes) or printed/saved directly.
#'
#' @examples
#' B_hat <- matrix(c(0, 1.2, -0.8, 0, 2.1, 0, -1.5, 0, 0.6, 0),
#'                  nrow = 5, ncol = 2,
#'                  dimnames = list(paste0("Var", 1:5), paste0("Sp", 1:2)))
#' coef_plot_barre(B_hat)
#' coef_plot_barre(B_hat, grad_echelle = c(-2, 2))
#'
#' @export
coef_plot_barre <- function(matrice_coefficient,
                            nom_axes = c("Columns of Y", "Variables",
                                         "Coefficients value", "SICZIPLN"),
                            grad_echelle = c(-5, 5)) {

  stopifnot(
    "`nom_axes` must be a character vector of length 4" =
      length(nom_axes) == 4,
    "`grad_echelle` must be a numeric vector of length 2" =
      length(grad_echelle) == 2
  )

  grad_max <- grad_echelle[2]
  grad_min <- grad_echelle[1]
  nom_axe_y <- nom_axes[1]
  nom_axe_x <- nom_axes[2]
  nom_gradient_couleur <- nom_axes[3]
  titre <- nom_axes[4]

  metlcoef_sic_genus <- reshape2::melt(matrice_coefficient)
  metlcoef_sic_genus$sparsity <- ifelse(metlcoef_sic_genus$value == 0, "Zero", "Nonzero")

  plot_genus2_sicpln <- ggplot2::ggplot(
    data = metlcoef_sic_genus,
    ggplot2::aes(x = as.factor(.data$Var2), y = as.factor(.data$Var1),
                 pattern = .data$sparsity, fill = .data$value)
  ) +
    ggpattern::geom_tile_pattern(
      width = 0.9, height = 0.9,
      pattern_fill = "black",
      pattern_angle = 45,
      pattern_density = 0.015,
      pattern_spacing = 0.06,
      pattern_key_scale_factor = 1
    ) +
    ggpattern::scale_pattern_manual(values = c(Zero = "circle", Nonzero = "none"), name = "") +
    ggplot2::scale_fill_gradient2(
      low = "darkred", high = "darkgreen", mid = "white", midpoint = 0,
      limit = c(grad_min, grad_max), name = nom_gradient_couleur,
      guide = ggplot2::guide_colourbar(barwidth = 0.3, barheight = 2.5)
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = nom_axe_y, y = nom_axe_x, title = titre) +
    ggplot2::guides(pattern = ggplot2::guide_legend(override.aes = list(fill = "white"))) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      axis.title = ggplot2::element_text(size = 12, colour = "black"),
      legend.text = ggplot2::element_text(size = 8, colour = "black"),
      legend.title = ggplot2::element_text(color = "black", size = 8),
      axis.text.x = ggplot2::element_text(size = 5, colour = "black", angle = 45, hjust = 1),
      strip.text.x = ggplot2::element_text(size = 5, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 8, colour = "black")
    )

  return(plot_genus2_sicpln)
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
