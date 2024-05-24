#' Function defined to enhance the usability for users on IDEs.
#' @export
f <- function(smoothing_var, model = "iid", sd.prior = NULL, boundary.prior = NULL, initial_location = c("middle", "left", "right"), ...) {
  # Capture the full call
  mc <- match.call(expand.dots = TRUE)
  
  # Replace the first argument with the function name
  mc[[1]] <- as.name("f")
  
  # Ensure all default arguments are included
  mc$model <- model
  mc$sd.prior <- sd.prior
  mc$boundary.prior <- boundary.prior
  mc$initial_location <- initial_location[1]
  
  # Replace smoothing_var with its unevaluated form
  mc$smoothing_var <- substitute(smoothing_var)
  
  # Return the modified call
  return(mc)
}

parse_formula <- function(formula) {
  components <- as.list(attributes(terms(formula))$ variables)
  fixed_effects <- list()
  rand_effects <- list()
  offset_effects <- list()
  # Index starts as 3 since index 1 represents "list" and
  # index 2 represents the response variable
  for (i in 3:length(components)) {
    if (startsWith(toString(components[[i]]), "f,")) {
      rand_effects[[length(rand_effects) + 1]] <- components[[i]]
    } 
    else if(startsWith(toString(components[[i]]), "offset,")){
      offset_effects[[length(offset_effects) + 1]] <- components[[i]]
    }
    else {
      fixed_effects[[length(fixed_effects) + 1]] <- components[[i]]
    }
  }
  return(list(response = components[[2]], fixed_effects = fixed_effects, rand_effects = rand_effects, offset_effects = offset_effects))
}

# Create a class for iwp using S4
setClass("iwp", slots = list(
  response_var = "name", smoothing_var = "name", order = "numeric",
  knots = "numeric", observed_x = "numeric", sd.prior = "list",
  psd.prior = "list",
  boundary.prior = "list", data = "data.frame", X = "matrix",
  B = "matrix", P = "matrix", initial_location = "ANY"
))

# Create a class for sgp using S4
setClass("sgp", slots = list(
  response_var = "name", smoothing_var = "name", 
  a = "numeric", freq = "numeric", period = "numeric",
  m = "numeric", k = "numeric",
  knots = "numeric", observed_x = "numeric", sd.prior = "list",
  psd.prior = "list",
  boundary.prior = "list", data = "data.frame", X = "matrix",
  B = "matrix", P = "matrix", initial_location = "ANY", region = "numeric", 
  accuracy = "numeric", boundary = "logical"
))

# Create a class for iid using S4
setClass("iid", slots = list(
  response_var = "name", smoothing_var = "name", sd.prior = "list",
  data = "data.frame", B = "matrix", P = "matrix"
))

# Create a class for customized using S4
setClass("customized", slots = list(
  response_var = "name", smoothing_var = "name", sd.prior = "list",
  data = "data.frame", B = "matrix", P = "matrix",
  compute_B = "function", compute_P = "function"
))

#' @import Matrix
#' @importClassesFrom Matrix dgCMatrix
Compute_Q_sB <- function(a,k,region, accuracy = 5000, boundary = TRUE){
  ss <- function(M) {Matrix::forceSymmetric(M + Matrix::t(M))}
  if((accuracy %% 1) == 0){
    x <- seq(min(region),max(region), length.out = accuracy)
  }
  else{
    x <- seq(min(region),max(region), by = accuracy)
  }
  if(boundary == TRUE){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  B1matrix <-  fda::eval.basis(x, B_basis, Lfdobj=1, returnMatrix=TRUE)
  B2matrix <-  fda::eval.basis(x, B_basis, Lfdobj=2, returnMatrix=TRUE)
  cos_matrix <- cos(a*x)
  sin_matrix <- sin(a*x)
  Bcos <- as(apply(Bmatrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  B1cos <- as(apply(B1matrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  B2cos <- as(apply(B2matrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  Bsin <- as(apply(Bmatrix, 2, function(x) x*sin_matrix), "dgCMatrix")
  B1sin <- as(apply(B1matrix, 2, function(x) x*sin_matrix), "dgCMatrix")
  B2sin <- as(apply(B2matrix, 2, function(x) x*sin_matrix), "dgCMatrix")
  
  ### Compute I, L, T:
  Numerical_I <- as(diag(c(diff(c(0,x)))), "dgCMatrix")
  
  ### T
  T00 <- Matrix::t(Bcos) %*% Numerical_I %*% Bcos
  T10 <- Matrix::t(B1cos) %*% Numerical_I %*% Bcos
  T11 <- Matrix::t(B1cos) %*% Numerical_I %*% B1cos
  T20 <- Matrix::t(B2cos) %*% Numerical_I %*% Bcos
  T21 <- Matrix::t(B2cos) %*% Numerical_I %*% B1cos
  T22 <- Matrix::t(B2cos) %*% Numerical_I %*% B2cos
  
  ### L
  L00 <- Matrix::t(Bsin) %*% Numerical_I %*% Bsin
  L10 <- Matrix::t(B1sin) %*% Numerical_I %*% Bsin
  L11 <- Matrix::t(B1sin) %*% Numerical_I %*% B1sin
  L20 <- Matrix::t(B2sin) %*% Numerical_I %*% Bsin
  L21 <- Matrix::t(B2sin) %*% Numerical_I %*% B1sin
  L22 <- Matrix::t(B2sin) %*% Numerical_I %*% B2sin
  
  ### I
  I00 <- Matrix::t(Bsin) %*% Numerical_I %*% Bcos
  I10 <- Matrix::t(B1sin) %*% Numerical_I %*% Bcos
  I11 <- Matrix::t(B1sin) %*% Numerical_I %*% B1cos
  I20 <- Matrix::t(B2sin) %*% Numerical_I %*% Bcos
  I21 <- Matrix::t(B2sin) %*% Numerical_I %*% B1cos
  I22 <- Matrix::t(B2sin) %*% Numerical_I %*% B2cos
  
  
  ### Inner product involving B spline:
  Bmatrix <- as(Bmatrix, "dgCMatrix")
  B1matrix <- as(B1matrix, "dgCMatrix")
  B2matrix <- as(B2matrix, "dgCMatrix")
  
  BB <- Matrix::t(Bmatrix) %*% Numerical_I %*% Bmatrix
  B2B2 <- Matrix::t(B2matrix) %*% Numerical_I %*% B2matrix
  BB2 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B2matrix
  BS <- Matrix::t(Bmatrix) %*% Numerical_I %*% Bsin
  BC <- Matrix::t(Bmatrix) %*% Numerical_I %*% Bcos
  BS1 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B1sin
  BC1 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B1cos
  BS2 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B2sin
  BC2 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B2cos
  B2S <- Matrix::t(B2matrix) %*% Numerical_I %*% Bsin
  B2C <- Matrix::t(B2matrix) %*% Numerical_I %*% Bcos
  B2S1 <- Matrix::t(B2matrix) %*% Numerical_I %*% B1sin
  B2C1 <- Matrix::t(B2matrix) %*% Numerical_I %*% B1cos
  B2S2 <- Matrix::t(B2matrix) %*% Numerical_I %*% B2sin
  B2C2 <- Matrix::t(B2matrix) %*% Numerical_I %*% B2cos
  
  ## G = <phi,phj>
  G <- rbind(cbind(T00, Matrix::t(I00), Matrix::t(BC)), cbind(I00,L00, Matrix::t(BS)), cbind(BC,BS,BB))
  
  ## C = <D^2phi,D^2phj>
  C11 <- T22 - 2*a*ss(I21) - (a^2)*ss(T20) + 2*(a^3)*ss(I10) + 4 * (a^2) * L11 + (a^4)*T00
  C22 <- L22 + 2*a*ss(I21) - (a^2)*ss(L20) - 2*(a^3)*ss(I10) + 4 * (a^2) * T11 + (a^4)*L00
  C12 <- I22 + 2*a*T21 - (a^2)* ss(I20) - 2*a*Matrix::t(L21) - 4*(a^2)*I11 + 2*(a^3)*L10 - 2*(a^3)*Matrix::t(T10) + (a^4)*I00
  
  C13 <- Matrix::t(B2C2) - 2*a*Matrix::t(B2S1) - (a^2)*Matrix::t(B2C)
  C23 <- Matrix::t(B2S2) + 2*a*Matrix::t(B2C1) - (a^2)*Matrix::t(B2S)
  C33 <- B2B2
  C <- rbind(cbind(C11,C12,C13), cbind(Matrix::t(C12), C22, C23), cbind(Matrix::t(C13), Matrix::t(C23), C33))
  
  
  ## M = <phi,D^2phj>
  M11 <- Matrix::t(T20) - (2*a)*Matrix::t(I10) - (a^2)*T00
  M12 <- Matrix::t(I20) + (2*a)*Matrix::t(T10) - (a^2)*I00
  M21 <- Matrix::t(I20) - (2*a)*Matrix::t(L10) - (a^2)*I00
  M22 <- Matrix::t(L20) + (2*a)*Matrix::t(I10) - (a^2)*L00
  
  M13 <- Matrix::t(B2C)
  M23 <- Matrix::t(B2S)
  M31 <- BC2 - (2*a)*BS1 - (a^2)*BC
  M32 <- BS2 + (2*a)*BC1 - (a^2)*BS
  M33 <- BB2
  
  M <- rbind(cbind(M11,M12,M13), cbind(M21,M22,M23), cbind(M31,M32,M33))
  
  
  ### Compute the final precision matrix: Q
  Q <- (a^4)*G + C + (a^2)*ss(M)
  Matrix::forceSymmetric(Q)
}


Compute_B_sB <- function(x, a, k, region, boundary = TRUE){
  if(boundary){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  cos_matrix <- cos(a*x)
  sin_matrix <- sin(a*x)
  Bcos <- apply(Bmatrix, 2, function(x) x*cos_matrix)
  Bsin <- apply(Bmatrix, 2, function(x) x*sin_matrix)
  cbind(Bcos, Bsin,Bmatrix)
}


Compute_B_sB_helper <- function(refined_x, a, k, m, region, boundary = TRUE, initial_location = NULL){
  if(is.null(initial_location)){
    initial_location <- min(refined_x)
  }
  refined_x <- refined_x - initial_location
  region <- region - initial_location
  B <- NULL
  for (i in 1:m) {
    B <- cbind(B, Compute_B_sB(x = refined_x, a = (i*a), region = region, k = k, boundary = boundary))
  }
  B
}


setGeneric("compute_B", function(object) {
  standardGeneric("compute_B")
})
setMethod("compute_B", signature = "iid", function(object) {
  smoothing_var <- object@smoothing_var
  x <- as.factor((object@data)[[smoothing_var]])
  B <- model.matrix(~ -1 + x)
  B
})
setMethod("compute_B", signature = "customized", function(object) {
  smoothing_var <- object@smoothing_var
  object@compute_B((object@data)[[smoothing_var]])
})
setMethod("compute_B", signature = "sgp", function(object) {
  smoothing_var <- object@smoothing_var
  a <- object@a
  k <- object@k
  m <- object@m
  initial_location <- object@initial_location
  region <- object@region - initial_location
  accuracy <- object@accuracy
  boundary <- object@boundary
  refined_x <- (object@data)[[smoothing_var]] - initial_location
  B <- NULL
  for (i in 1:m) {
    B <- cbind(B, Compute_B_sB(x = refined_x, a = (i*a), region = region, k = k))
  }
  B
})


setGeneric("compute_P", function(object) {
  standardGeneric("compute_P")
})
setMethod("compute_P", signature = "iid", function(object) {
  smoothing_var <- object@smoothing_var
  x <- (object@data)[[smoothing_var]]
  num_factor <- length(unique(x))
  diag(nrow = num_factor, ncol = num_factor)
})
setMethod("compute_P", signature = "customized", function(object) {
  x <- (object@data)[[object@smoothing_var]]
  object@compute_P(x)
})
setMethod("compute_P", signature = "sgp", function(object) {
  smoothing_var <- object@smoothing_var
  initial_location <- object@initial_location
  a <- object@a
  k <- object@k
  m <- object@m
  region <- object@region - initial_location
  accuracy <- object@accuracy
  boundary <- object@boundary
  refined_x <- (object@data)[[smoothing_var]] - initial_location
  Q <- Compute_Q_sB(a = a, k = k, region = region, accuracy = accuracy)
  if(m >= 2){
    for (i in 2:m) {
      Q <- bdiag(Q, Compute_Q_sB(a = (i*a), k = k, region = region, accuracy = accuracy))
    }
  }
  as(Q, "matrix")
})


setGeneric("local_poly", function(object) {
  standardGeneric("local_poly")
})
setMethod("local_poly", signature = "iwp", function(object) {
  knots <- object@knots
  initial_location <- object@initial_location
  smoothing_var <- object@smoothing_var
  refined_x <- (object@data)[[smoothing_var]] - initial_location
  p <- object@order
  D <- local_poly_helper(knots = knots, refined_x = refined_x, p = p, neg_sign_order = 0)
  D # Local poly design matrix
})

setGeneric("global_poly", function(object) {
  standardGeneric("global_poly")
})
setMethod("global_poly", signature = "iwp", function(object) {
  smoothing_var <- object@smoothing_var
  x <- (object@data)[[smoothing_var]] - object@initial_location
  p <- object@order
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
})
setMethod("global_poly", signature = "sgp", function(object) {
  smoothing_var <- object@smoothing_var
  a <- object@a
  m <- object@m
  initial_location <- object@initial_location
  refined_x <- (object@data)[[smoothing_var]] - initial_location
  X <- NULL
  for (i in 1:m) {
    X <- cbind(X, cos(i*a*refined_x),sin(i*a*refined_x))
  }
  X 
})


#' Constructing the precision matrix given the knot sequence
#'
#' @param x A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @return A precision matrix of the corresponding basis function, should be diagonal matrix with
#' size (k-1) by (k-1).
#' @export
setGeneric("compute_weights_precision", function(object) {
  standardGeneric("compute_weights_precision")
})
setMethod("compute_weights_precision", signature = "iwp", function(object) {
  knots <- object@knots
  if (min(knots) >= 0) {
    as(diag(diff(knots)), "matrix")
  } else if (max(knots) < 0) {
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    as(diag(diff(knots_neg)), "matrix")
  } else {
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    knots_pos <- knots
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    d1 <- diff(knots_neg)
    d2 <- diff(knots_pos)
    Precweights1 <- diag(d1)
    Precweights2 <- diag(d2)
    as(Matrix::bdiag(Precweights1, Precweights2), "matrix")
  }
})

#' Constructing the precision matrix given the knot sequence (helper)
#'
#' @param x A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @return A precision matrix of the corresponding basis function, should be diagonal matrix with
#' size (k-1) by (k-1).
#' @examples
#' compute_weights_precision(x = c(0,0.2,0.4,0.6,0.8))
#' @export
compute_weights_precision_helper <- function(x){
  d <- diff(x)
  Precweights <- diag(d)
  Precweights
}



get_local_poly <- function(knots, refined_x, p) {
  dif <- diff(knots)
  nn <- length(refined_x)
  n <- length(knots)
  D <- matrix(0, nrow = nn, ncol = n - 1)
  for (j in 1:nn) {
    for (i in 1:(n - 1)) {
      if (refined_x[j] <= knots[i]) {
        D[j, i] <- 0
      } else if (refined_x[j] <= knots[i + 1] & refined_x[j] >= knots[i]) {
        D[j, i] <- (1 / factorial(p)) * (refined_x[j] - knots[i])^p
      } else {
        k <- 1:p
        D[j, i] <- sum((dif[i]^k) * ((refined_x[j] - knots[i + 1])^(p - k)) / (factorial(k) * factorial(p - k)))
      }
    }
  }
  D # Local poly design matrix
}

#' Constructing and evaluating the local O-spline basis (design matrix)
#'
#' @param knots A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @param refined_x A vector of locations to evaluate the O-spline basis
#' @param p An integer value indicates the order of smoothness
#' @param neg_sign_order An integer value N such that D = ((-1)^N)*D for the splines at negative knots. Default is 0.
#' @return A matrix with i,j component being the value of jth basis function
#' value at ith element of refined_x, the ncol should equal to number of knots minus 1, and nrow
#' should equal to the number of elements in refined_x.
#' @examples
#' local_poly(knots = c(0, 0.2, 0.4, 0.6, 0.8), refined_x = seq(0, 0.8, by = 0.1), p = 2)
#' @export
local_poly_helper <- function(knots, refined_x, p = 2, neg_sign_order = 0) {
  if (min(knots) >= 0) {
    # The case of all-positive knots
    D <- get_local_poly(knots, refined_x, p)
  } else if (max(knots) <= 0) {
    # Handle the negative part only
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    D <- get_local_poly(knots_neg, refined_x_neg, p)
    D <- D * ((-1)^neg_sign_order)
  } else {
    # Handle the negative part
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    D1 <- get_local_poly(knots_neg, refined_x_neg, p)
    D1 <- D1 * ((-1)^neg_sign_order)
    
    # Handle the positive part
    refined_x_pos <- ifelse(refined_x > 0, refined_x, 0)
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    D2 <- get_local_poly(knots_pos, refined_x_pos, p)
    
    D <- cbind(D1, D2)
  }
  D # Local poly design matrix
}

#' Constructing and evaluating the global polynomials, to account for boundary conditions (design matrix)
#'
#' @param x A vector of locations to evaluate the global polynomials
#' @param p An integer value indicates the order of smoothness
#' @return A matrix with i,j componet being the value of jth basis function
#' value at ith element of x, the ncol should equal to p, and nrow
#' should equal to the number of elements in x
#' @examples
#' global_poly(x = c(0, 0.2, 0.4, 0.6, 0.8), p = 2)
#' @export
global_poly_helper <- function(x, p = 2) {
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
}

#' Constructing and evaluating the global polynomials, to account for boundary conditions (design matrix) of sgp
#'
#' @param refined_x A vector of locations to evaluate the sB basis
#' @param a The frequency of sgp.
#' @param m The number of harmonics to consider
#' @return A matrix with i,j componet being the value of jth basis function
#' value at ith element of x, the ncol should equal to (2*m), and nrow
#' should equal to the number of elements in x
#' @export
global_poly_helper_sgp <- function(refined_x, a, m, initial_location = NULL) {
  if(is.null(initial_location)){
    initial_location <- min(refined_x)
  }
  refined_x <- refined_x - initial_location
  X <- NULL
  for (i in 1:m) {
    X <- cbind(X, cos(i*a*refined_x),sin(i*a*refined_x))
  }
  X 
}

#' Construct prior based on d-step prediction SD (for iwp)
#'
#' @param prior A list that contains alpha and u. This specifies the target prior on the d-step SD \eqn{\sigma(d)}, such that \eqn{P(\sigma(d) > u) = alpha}.
#' @param d A numeric value for the prediction step.
#' @param p An integer for the order of iwp.
#' @return A list that contains alpha and u. The prior for the smoothness parameter \eqn{\sigma} such that \eqn{P(\sigma > u) = alpha}, that yields the ideal prior on the d-step SD.
#' @export
prior_conversion_iwp <- function(d, prior, p) {
  Cp <- (d^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2))
  prior_q <- list(alpha = prior$alpha, u = (prior$u * (1 / sqrt(Cp))))
  prior_q
}


#' Compute the SD correction factor for sgp
#' @param d A numeric value for the prediction step.
#' @param a The frequency parameter of the sgp.
#' @return The correction factor c that should be used to compute the d-step PSD as c*SD.
compute_d_step_sgpsd <- function(d,a){
  sqrt((1/(a^2))*((d/2) - (sin(2*a*d)/(4*a))))
}


#' Construct prior based on d-step prediction SD (for sgp)
#'
#' @param prior A list that contains alpha and u. This specifies the target prior on the d-step SD \eqn{\sigma(d)}, such that \eqn{P(\sigma(d) > u) = alpha}.
#' @param d A numeric value for the prediction step.
#' @param a The frequency parameter of the sgp.
#' @param m The number of harmonics that should be considered, by default m = 1 represents only the sgp.
#' @return A list that contains alpha and u. The prior for the smoothness parameter \eqn{\sigma} such that \eqn{P(\sigma > u) = alpha}, that yields the ideal prior on the d-step SD.
#' @export
prior_conversion_sgp <- function(d, prior, a, m = 1) {
  correction_factor <- 0
  for (i in 1:m) {
    correction_factor <- correction_factor + compute_d_step_sgpsd(d = d, a = (i*a))
  }
  prior_SD <- list(u = prior$u/correction_factor, alpha = prior$alpha)
  prior_SD
}

#' @import Matrix
#' @importClassesFrom Matrix dgCMatrix
dgTMatrix_wrapper <- function(matrix) {
  # result <- as(as(as(matrix, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  result <- as(matrix, "dgCMatrix")
  result
}

#' @export
get_default_option_list_MCMC <- function(option_list = list()){
  default_options <- list(chains = 1, cores = 1, init = "random", seed = 123, warmup = 10000, silent = TRUE, laplace = FALSE)
  required_names <- names(default_options)
  for (name in required_names) {
    if(! name %in% names(option_list)){
      option_list[[name]] <- default_options[[name]]
    }
  }
  option_list
}


#' Custom Template Function
#'
#' This function allows for the dynamic modification of a C++ template
#' within the BayesGP package. Users can specify custom content for the 
#' log-likelihood as well as the log-prior of the variance parameter in the template.
#' 
#' 
#' @param SETUP A character string or vector containing the 
#'   lines of C++ code to be inserted before the computation of log-likelihood.
#' @param LOG_LIKELIHOOD A character string or vector containing the 
#'   lines of C++ code to be inserted in the log-likelihood section of 
#'   the template. Should be NULL if no changes are to be made to this section.
#' @param LOG_PRIOR A character string or vector containing the 
#'   lines of C++ code to be inserted in the log-prior (of the variance parameter) section of 
#'   the template. Should be NULL if no changes are to be made to this section.
#' @param compile_template A indicator of whether the new template should be compiled. default is FALSE.
#' @return A string representing the path to the temporary .so (or .dll) file
#'   containing the compiled custom C++ code.
#'
#' @export
custom_template <- function(SETUP = NULL, LOG_LIKELIHOOD = NULL, LOG_PRIOR = NULL, compile_template = FALSE){
  template_path <- system.file("extsrc", "BayesGP.cpp", package = "BayesGP")
  file_content <- readLines(template_path)
  if(!is.null(SETUP)){
    start_line = which(file_content == "  // START OF YOUR SETUP")
    end_line = which(file_content == "  // END OF YOUR SETUP")
    if (length(start_line) == 1 && length(end_line) == 1 && start_line < end_line) {
      # Replace the content
      file_content <- c(file_content[1:(start_line)], SETUP, file_content[(end_line+1):length(file_content)])
    } else {
      warning("Markers not found or in incorrect order")
    }
  }
  if(!is.null(LOG_LIKELIHOOD)){
    start_line = which(file_content == "    // START OF YOUR LOG_LIKELIHOOD: ll")
    end_line = which(file_content == "    // END OF YOUR LOG_LIKELIHOOD")
    if (length(start_line) == 1 && length(end_line) == 1 && start_line < end_line) {
      # Replace the content
      file_content <- c(file_content[1:(start_line)], LOG_LIKELIHOOD, file_content[(end_line+1):length(file_content)])
    } else {
      warning("Markers not found or in incorrect order")
    }
  }
  if(!is.null(LOG_PRIOR)){
    start_line = which(file_content == "  // START OF YOUR LOG_PRIOR: lpT")
    end_line = which(file_content == "  // END OF YOUR LOG_PRIOR")
    if (length(start_line) == 1 && length(end_line) == 1 && start_line < end_line) {
      # Replace the content
      file_content <- c(file_content[1:(start_line)], LOG_PRIOR, file_content[(end_line+1):length(file_content)])
    } else {
      warning("Markers not found or in incorrect order")
    }
  }
  
  # Create a temporary file
  temp_file <- tempfile(fileext = ".cpp") # .cpp extension is added for clarity
  
  # Write the modified content to the temporary file
  writeLines(file_content, temp_file)
  temp_so_file <- sub("\\.cpp$", "", temp_file)
  
  if(compile_template == TRUE){
    suppressWarnings(suppressMessages(invisible(TMB::compile(file = temp_file)) ))
    dyn.load(TMB::dynlib(temp_so_file))
  }
  
  return(temp_so_file)
}


#' Roxygen commands
#'
#' @useDynLib BayesGP
#'
dummy <- function() {
  return(NULL)
}
