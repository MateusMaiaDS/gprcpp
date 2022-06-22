#' Main GP-function
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
  # K_y <-
  # Getting the functions


}
