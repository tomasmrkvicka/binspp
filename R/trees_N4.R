#' Spanish oak trees
#'
#' The oak trees dataset sampled in 2009 in region consisting of 5 rectangles.
#'
#' The data contains point pattern of trees, 5 covariates
#'     (refor, reserve, slope, tdensity, tmi) and 4 vectors
#'     (x_left, x_right, y_bottom, y_top) of corners of rectangles forming
#'     the observation window.
#'
#' @source Jesús Fernández-Habas, Pilar Fernández-Rebollo,
#'         Mónica Rivas Casado, Alma María García Moreno,
#'         Begoña Abellanas. Spatio-temporal analysis of
#'         oak decline process in open woodlands: A case study
#'         in SW Spain, Journal of Environmental Management,
#'         248, 2019, 109308,
#'         \doi{10.1016/j.jenvman.2019.109308}.
#'
#' @format A list with columns:
#' \describe{
#'  \item{window}{A list of region window definition.}
#'  \item{n}{Number of oak trees in the region.}
#'  \item{x}{Array of x coordinates of oak trees.}
#'  \item{y}{Array of y coordinates of oak trees.}
#' }
#' @examples
#'  plot(trees_N4)
"trees_N4"
