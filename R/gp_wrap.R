#' @useDynLib gprcpp
#' @importFrom Rcpp sourceCpp

#' @export
gp_main <- function(x,
                    x_star,
                    y,
                    gp_param = list(phi = 0.1,
                                    nugget = 0.1,
                                    nu = 1)){
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
  phi <- gp_param[["phi"]]
  nugget <- gp_param[["nugget"]]
  nu <- gp_param[["nu"]]

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

  return(mean_sd)
}


ggplot()+
  geom_point(data = tibble(x = x, y = y),mapping = aes(x = x, y = y))+
  geom_line(data = tibble (x = x_star, y = mean_sd$y_hat),
            mapping = aes(x = x, y = y),
            col = "blue")+
  geom_ribbon(data = tibble(x = x_star,
                            ymin = mean_sd$y_hat-1.96*mean_sd$y_sd,
                            ymax = mean_sd$y_hat+1.96*mean_sd$y_sd),
              mapping = aes(x = x, ymax = ymax, ymin = ymin),alpha = 0.2)
