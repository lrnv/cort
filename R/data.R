#' Dataset funcdep_data
#'
#' This dependence structure is constructed by applying the function :
#' \deqn{h(u_1,u_2,u_3) = (u_{1},\sin(2\pi u_{1})-\frac{u_{2}}{\pi},(1+\frac{u_{3}}{\pi^{2}})(\frac{u_{3}}{2} I_{\frac{1}{4}\ge u_1}-\sin(\pi^{x_{1}}) I_{\frac{1}{4} < u_{1}}))}
#' to uniformly drawn 3-dimensional random vectors. The dataset is the ranks of \eqn{h(u)}.
#'
#' This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
#'
#' @format A matrix with 500 rows and 3 columns
#'
#' @references
#' \insertRef{laverny2020}{cort}
"funcdep_data"

#' Dataset impossible_data
#'
#' We simulate from a density inside the piecewise linear copula class, by applying the function:
#' \deqn{h(u) = (u_1,          \frac{u_2}{2} + \frac{1}{2}I_{u_1 \notin (\frac{1}{3}, \frac{2}{3})})}
#' to a 200x2 uniform sample, and taking ranks.
#'
#' This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
#'
#' @format A matrix with 200 rows and 2 columns
#'
#' @references
#' \insertRef{laverny2020}{cort}
"impossible_data"

#' Dataset recoveryourself_data
#'
#' This dataset is a simple test: we simulate random samples from a density inside the piecewise copula class,
#' and test whether or not the estimator can recover it. For that, we will use a 2-dimensional sample with 500
#' observations, uniform on the unit hypercube, and apply the following function:
#' \deqn{h(u) = (u_1, \frac{u_2 + I_{u_1 \le \frac{1}{4}} + 2I_{u_1 \le \frac{1}{2}} + I_{\frac{3}{4} \le u_1}}{4})}
#'
#' This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
#'
#'
#' @format A matrix with 500 rows and 2 columns
#'
#' @references
#' \insertRef{laverny2020}{cort}
"recoveryourself_data"

#' Dataset clayton_data
#'
#' This dataset is a simulation of 200 points from a 3-dimensional clayton copula with \eqn{\theta = 7},
#' hence highly dependent, for the first, third and fourth marginals. The second marginal is added
#' as independent uniform draws. Lastly, the third marginal is flipped,
#' inducing a negative dependence structure.
#'
#' This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
#'
#' @format A matrix with 200 rows and 4 columns
#'
#' @references
#' \insertRef{laverny2020}{cort}
"clayton_data"
