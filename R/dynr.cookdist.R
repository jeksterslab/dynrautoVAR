#' Cook's Distance in Dynr
#'
#' @author Jonathan Park
#'
#' Maybe we can put a description here on how Cook's distance is calculated.
#' Also Jonathan can you give me the data set you used in the example?
#'
#' I made the following changes to the code.
#' 1. I changed the name to `dynr.cookdist` to match the other functions in the package.
#' 1. I changed the argument `data` to `dataframe` to match [dynr::dynr.data()].
#' 1. I changed the argument `ID` to `id` to match [dynr::dynr.data()].
#' 1. I set a default value for the argument `dir`.
#' 1. I used [styler::style_pkg()] to style the code for added readability.
#'
#' @param dataframe a data frame object of data that contain a column of subject ID numbers (i.e., an ID variable), a column indicating subject-specific measurement occasions (i.e., a TIME variable), at least one column of observed values, and any number of covariates. The TIME variable should contain subject-specific sequences of (subsets of) consecutively equally spaced numbers (e.g, 1, 2, 3, ...). That is, the program assumes that the input data.frame is equally spaced with potential missingness. If the measurement occasions for a subject are a subset of an arithmetic sequence but are not consecutive, NAs will be inserted automatically to create an equally spaced data set before estimation. Missing values in the observed variables shoud be indicated by NA.
#' @param time a character string of the name of the TIME variable in the data.
#' @param id a character string of the name of the ID variable in the data.
#' @param ini.beta the matrix of starting/fixed values for the transition matrix in the specified linear dynamic model.
#' @param dir path for output files.
#' @export
dynr.cookdist <- function(dataframe,
                          id,
                          time,
                          ini.beta,
                          dir = getwd()) {
  bet <- ini.beta
  nvar <- ncol(dataframe) - 2
  names(dataframe) <- c(paste0("X", 1:nvar))
  pseudo.data <- dataframe[, -c(id, time)]
  pseudo.data$id <- rep(1, nrow(pseudo.data))
  pseudo.data$time <- 1:nrow(pseudo.data)
  rawdata <- dynr.data(dataframe,
    id = "id", time = "time",
    observed = c(paste0("X", 1:nvar))
  )
  load.mat <- diag(1, nvar)
  meas <- prep.measurement(
    values.load = load.mat,
    params.load = matrix(rep("fixed", (nvar^2)), ncol = nvar),
    state.names = paste0("X", 1:nvar),
    obs.names = paste0("X", 1:nvar)
  )
  letters <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k")
  dynm <- prep.matrixDynamics(values.dyn = bet, params.dyn = matrix(levels(interaction(letters[1:nvar], 1:nvar, sep = "_")), nvar, nvar, byrow = TRUE), isContinuousTime = FALSE)
  means <- paste0("mu", 1:nvar)
  initial <- prep.initial(
    values.inistate = rep(0, nvar),
    params.inistate = rep("fixed", nvar),
    values.inicov = load.mat,
    params.inicov = matrix("fixed", nvar, nvar)
  )
  covs <- c("C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10")
  m1 <- matrix(levels(interaction(covs[1:nvar], 1:nvar, sep = "_")), nvar, nvar, byrow = TRUE)
  m1[lower.tri(m1)] <- t(m1)[lower.tri(m1)]
  m1 <- matrix("fixed", nvar, nvar)
  diag(m1) <- paste0("VARIANCE_", 1:nvar)
  mdcov <- prep.noise(
    # values.latent = cov(pseudo.data[,1:nvar]),
    # params.latent = m1,
    values.latent = diag(1, nvar),
    params.latent = m1,
    values.observed = diag(rep(0, nvar)),
    params.observed = matrix(("fixed"), nvar, nvar)
  )

  model <- dynr.model(
    dynamics = dynm, measurement = meas,
    noise = mdcov, initial = initial, data = rawdata,
    outfile = paste0(dir, "/OGModel10.c")
  )
  model$ub[1:nvar^2] <- 1.00
  model$lb[1:nvar^2] <- -1.00
  results <- dynr.cook(model)
  obj <- summary(results)
  effects <- matrix(obj$Coefficients[1:nvar^2, 1], nvar, nvar, byrow = FALSE)
  sigs <- matrix(as.numeric(obj$Coefficients[1:nvar^2, 6] < 0.05), nvar, nvar)
  OG.betas <- effects * sigs
  OG.inv.Hess <- results$"inv.hessian"[1:nvar^2, 1:nvar^2]
  # LOO: Leave One Out Loop
  B.i <- list()
  B.i.sig <- list()
  index <- 1
  for (i in unique(dataframe[, id])) {
    B.i.signs <- NULL
    temp.dat <- dataframe[dataframe[, id] != i, ]
    temp.dat <- temp.dat[, -c(id, time)]
    temp.dat$id <- rep(1, nrow(temp.dat))
    temp.dat$time <- 1:nrow(temp.dat)
    rawdata <- dynr.data(temp.dat,
      id = "id", time = "time",
      observed = c(paste0("X", 1:nvar))
    )
    model <- dynr.model(
      dynamics = dynm, measurement = meas,
      noise = mdcov, initial = initial, data = rawdata,
      outfile = paste0(dir, "/LOO1_", index, ".c")
    )
    model$ub[1:nvar^2] <- 1.00
    model$lb[1:nvar^2] <- -1.00
    results <- dynr.cook(model)
    obj <- summary(results)
    effects <- matrix(obj$Coefficients[1:nvar^2, 1], nvar, nvar, byrow = FALSE)
    sigs <- matrix(as.numeric(obj$Coefficients[1:nvar^2, 6] < 0.05), nvar, nvar)
    B.i[[index]] <- effects * sigs
    index <- index + 1
    message(paste0("I am done with Subject ", index - 1))
  }
  p <- ncol(dataframe[, -c(id, time)])
  ps <- p^2 + p + ((p)^2 - p) / 2
  CookD <- list()
  for (i in 1:length(B.i)) {
    CookD[[i]] <- t(c(OG.betas - B.i[[i]])) %*%
      solve(OG.inv.Hess) %*% c(OG.betas - B.i[[i]]) / ps
  }
  return(CookD)
}
