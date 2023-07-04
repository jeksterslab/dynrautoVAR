#' Subgrouping Algorithm
#'
#' @author Jonathan Park
#'
#' This function takes a list of outputs from [dynr::dynr.cook()] or output from [dynr.var()] and
#' performs a subgrouping analysis. The output is a list of clusters.
#'
#' @details I encapsulated `dynrautoVAR` into a temporary package so I can render the documentation and do some tests.
#'
#' I made the following changes to the code.
#' 1. I changed the name to `dynr.sub()` to match the other functions in the package.
#' 1. I changed the argument `inputs` to `input`.
#' 1. I removed the argument `type` and based the type on the class of the input.
#' 1. I changed the argument `p.val` to `alpha`.
#' 1. I changed the argument `params` to `params.cook` and `var.opt`` to `params.var`.
#' 1. I used [styler::style_pkg()] to style the code for added readability.
#'
#' @references Add references to the subgrouping algorithms here.
#'
#' @param input list of [dynr::dynr.cook()] objects or output from [dynr.var()].
#' @param alpha significance level for testing of transition matrix coefficients.
#' @param params.cook parameters for clustering for [dynr::dynr.cook()] input.
#' @param params.var parameters for clustering for [dynr.var()] input.
#' @param method Subgrouping method.
#' @param k number of clusters.
#' @param m value of fuzzifier.
#' @importFrom fclust FKM
#' @importFrom igraph cluster_walktrap graph_from_adjacency_matrix
#' @export
dynr.sub <- function(input,
                     alpha = 0.05,
                     params.cook = NULL,
                     params.var = c(
                       "Betas", "ResidCov",
                       "PDC", "PCC", "Both"
                     ),
                     method = c("hard", "fuzz"),
                     k = 2,
                     m = 2) {
  if (inherits(input, "dynrVar")) {
    type <- "dynr.var"
  } else {
    if (is.list(input)) {
      lapply(
        X = input,
        FUN = function(x) {
          if (!inherits(x, "dynrCook")) {
            stop("Input must be a list of dynrCook objects or a dynrVar object.")
          }
        }
      )
      type <- "dynr.cook"
    } else {
      stop("Input must be a list of dynrCook objects or a dynrVar object.")
    }
  }
  if (type == "dynr.var") {
    comp.list <- NULL
    N <- length(input)
    if (params.var == "Both") {
      params.var <- c("PDC", "PCC")
    }
    for (i in 1:N) {
      comp.list <- cbind(comp.list, unlist(input[[i]][params.var]))
    }
    adj.mat <- matrix(NA, N, N)
    for (i in 1:N) {
      for (j in 1:N) {
        adj.mat[i, j] <- sum(abs(comp.list[, i]) > 0 &
          abs(comp.list[, j]) > 0 &
          sign(comp.list[, i]) == sign(comp.list[, j]))
      }
    }
    adj.mat <- adj.mat - min(adj.mat)
    diag(adj.mat) <- 0
    if (method == "hard") {
      grapp <- graph_from_adjacency_matrix(
        adj.mat,
        weighted = TRUE,
        mode = "undirected"
      )
      sub.sol <- cluster_walktrap(grapp)
    } else {
      if (method == "fuzz") {
        sub.sol <- FKM(adj.mat, k = k, m = m)
      }
    }
  }
  if (type == "dynr.cook") {
    N <- length(input)
    comp.list <- NULL
    for (i in 1:N) {
      temp <- summary(input[[i]])$Coefficients[, c(1, 6)]
      noms <- rownames(temp)
      comp.list <- cbind(comp.list, temp[, 1] * as.numeric(temp[, 2] < alpha))
    }
    row.names(comp.list) <- noms
    if (!is.null(params.cook)) {
      comp.list[which(row.names(comp.list) %in% params.cook), ]
    } else {
      message(
        "\'params.cook\' is NULL so we are using all parameters provided to cluster."
      )
    }

    adj.mat <- matrix(NA, N, N)
    for (i in 1:N) {
      for (j in 1:N) {
        adj.mat[i, j] <- sum(abs(comp.list[, i]) > 0 &
          abs(comp.list[, j]) > 0 &
          sign(comp.list[, i]) == sign(comp.list[, j]))
      }
    }
    adj.mat <- adj.mat - min(adj.mat)
    diag(adj.mat) <- 0
    if (method == "hard") {
      grapp <- graph_from_adjacency_matrix(
        adj.mat,
        weighted = TRUE,
        mode = "undirected"
      )
      sub.sol <- cluster_walktrap(grapp)
    } else {
      if (method == "fuzz") {
        sub.sol <- FKM(adj.mat, k = k, m = m)
      }
    }
  }
  return(sub.sol)
}
