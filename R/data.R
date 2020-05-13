#' Dataset funcdep_data
#'
#' This dependence structure is constructed by applying the function :
#' \deqn{h(u_1,u_2,u_3) = (u_1,\sin(2\pi u_1)-\frac{u_2}{\pi},(1+\frac{u_3}{\pi^{2}})(\frac{u_3}{2}\mathbb{1}_{\frac{1}{4} \ge u_1}  - \sin(\pi^{x_1})    \mathbb{1}_{\frac{1}{4} < u_1}))}
#' to uniformly drawn 3-dimensional random vectors. The dataset is the ranks of \eqn{h(u)}.
#'
#' @format A matrix with 500 rows and 3 columns
"funcdep_data"

#' Dataset impossible_data
#'
#' We simulate from a density inside the piecewise linear copula class, by applying the function:
#' \deqn{h(u) = (u_1,          \frac{u_2}{2} + \frac{1}{2}\mathbb{1}_{u_1 \notin (\frac{1}{3}, \frac{2}{3})})}
#' to a 200x2 uniform sample, and taking ranks.
#'
#' @format A matrix with 200 rows and 2 columns
"impossible_data"

#' Dataset recoveryourself_data
#'
#' This dataset is a simple test: we simulate random samples from a density inside the piecewise copula class,
#' and test whether or not the estimator can recover it. For that, we will use a 2-dimensional sample with 500
#' observations, uniform on the unit hypercube, and apply the following function:
#' \deqn{h(u) = (u_1, \frac{u_2 + \mathbb{1}_{u_1 \le \frac{1}{4}} + 2\mathbb{1}_{u_1 \le \frac{1}{2}} + \mathbb{1}_{\frac{3}{4} \le u_1}}{4})}
#'
#'
#' @format A matrix with 500 rows and 2 columns
"recoveryourself_data"

#' Dataset clayton_data
#'
#' This dataset is a simulation of 200 points from a 3-dimensional clayton copula with \eqn{\theta = 7},
#' hence highly dependent, for the first, third and fourth marginals. The second marginal is added
#' as independent uniform draws. Lastly, the third marginal is flipped,
#' inducing a negative dependence structure.
#'
#' @format A matrix with 200 rows and 4 columns
"clayton_data"
