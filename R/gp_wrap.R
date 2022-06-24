#' @useDynLib gprcpp
#' @importFrom Rcpp sourceCpp

#' @export
gp_main <- function(x,
                    x_star,
                    y,
                    gp_param = list(phi = 0.1,
                                    nugget = 0.1,
                                    nu = 1),
                    phi_sample = FALSE){
  # Checking the arguments
  if(!is.matrix(x)){
    stop("Insert a valid 'x' matrix")
  }

  if(!is.matrix(x_star)){
    stop("Insert a valid 'x_star' matrix")
  }

  if(!is.matrix(y)){
    stop("Insert a valid 'y' matrix")
  }

  # Getting the parameters values
  nugget <- gp_param[["nugget"]]
  nu <- gp_param[["nu"]]

  # Getting the phi sample if is required
  if(phi_sample){
    phi_post <- phi_post_sample(X = x,y = y,n_mcmc = 2000,n_burn = 500,
                                nu = nu,nugget = nugget,phi_init = 0.1)
    phi <- mean(phi_post)
  } else {
    phi <- gp_param[["phi"]]
  }

  if(any(is.null(phi),is.null(nugget),is.null(nu))){
    stop("Insert a valid parameters list")
  }

  # Getting the values for the function
  K_y <- k_y_nugget(A = x,phi = phi,nu = nu,nugget = nugget)

  # Getting the functions
  K_a_b <- k_A_B(A = x,B = x_star,phi = phi,nu = nu,nugget = nugget)

  # Getting the K_star_star
  K_star_star <- k_y_nugget(A = x_star,phi = phi,nu = nu,nugget = nugget)

  # Getting the GP-mean
  gp_mean <- gp_mean(K_y_nug = K_y,K_A_B = K_a_b,y = y)
  gp_cov_matrix <- gp_cov(K_y_nug = K_y,K_new = K_star_star,K_A_B = K_a_b)

  # Store the results
  mean_sd <- list(y_hat = c(gp_mean),
                  y_sd = sqrt(diag(gp_cov_matrix)))

  class(mean_sd) <- "gp.object"
  cat(paste0("The phi value is ",phi,"\n"))

  return(mean_sd)
}



