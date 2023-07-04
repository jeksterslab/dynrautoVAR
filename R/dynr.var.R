#' Automatically Fit VAR(1) Models Based on Input Data
#'
#' @author Jonathan Park
#'
#' This function takes in a `data.frame`, fits individual models for each subject indicated by the ID variable
#' and stores the fitted models as a list.
#' Note that if `ini.*` arguments are not provided, initial values are set to null matrices or null vectors.
#'
#' @return Returns a list of fitted VAR(1) models for each subject.
#'
#' @details I encapsulated `dynrautoVAR` into a temporary package so I can render the documentation and do some tests.
#'
#' I made the following changes to the code.
#' 1. I changed the argument `data` to `dataframe` to match [dynr::dynr.data()].
#' 1. I changed the argument `ID` to `id` to match [dynr::dynr.data()].
#' 1. I changed the argument `p.val` to `alpha`.
#' 1. I removed the default `NULL` value on arguments that require explicit values (`dataframe`, `nv`, and `time`).
#' 1. I made all initial condition arguments begin with `ini.`.
#' 1. I set a default value for the argument `dir` to remove it from the argument handling section.
#' 1. I cleaned up the argument handling section.
#' 1. I changed `result$estimation.result` to `results$estimation.result`.
#' 1. I gave the output a class of `dynrVar` for subsequent use of the output.
#' 1. I used [styler::style_pkg()] to style the code for added readability.
#'
#' Note that I mainly based the way I documented the arguments from functions in the `dynr` package used within this function. If you are amenable to the way I documented this function, I will proceed to documenting the rest.
#'
#' @references Add references to [dynr.var()] here.
#'
#' @param dataframe a data frame object of data that contain a column of subject ID numbers (i.e., an ID variable), a column indicating subject-specific measurement occasions (i.e., a TIME variable), at least one column of observed values, and any number of covariates. The TIME variable should contain subject-specific sequences of (subsets of) consecutively equally spaced numbers (e.g, 1, 2, 3, ...). That is, the program assumes that the input data.frame is equally spaced with potential missingness. If the measurement occasions for a subject are a subset of an arithmetic sequence but are not consecutive, NAs will be inserted automatically to create an equally spaced data set before estimation. Missing values in the observed variables shoud be indicated by NA.
#' @param nv number of variables.
#' @param time a character string of the name of the TIME variable in the data.
#' @param id a character string of the name of the ID variable in the data.
#' @param dir path for output files.
#' @param alpha significance level for testing of transition matrix coefficients.
#' @param ini.mu a vector of the starting or fixed values of the initial state vector.
#' @param ini.cov a positive definite matrix of the starting or fixed values of the initial error covariance structure.
#' @param ini.beta the matrix of starting/fixed values for the transition matrix in the specified linear dynamic model.
#' @param ini.loadings matrix of starting or fixed values for factor loadings.
#' @param use.mi if `use.mi = TRUE`, use [dynr::dynr.mi()] to address missing data in covariates.
#' @param aux names of the auxiliary variables used in the imputation model
#'   if `use.mi = TRUE`.
#' @import dynr
#' @export
dynr.var <- function(dataframe,
                     nv,
                     id = NULL,
                     time,
                     dir = getwd(),
                     alpha = 0.05,
                     ini.mu = NULL,
                     ini.cov = NULL,
                     ini.beta = NULL,
                     ini.loadings = NULL,
                     use.mi = FALSE,
                     aux = NULL) {
  # User Argument Handling
  if (is.null(id)) {
    message(
      "We are assuming you are only fitting to N = 1 subject.\n",
      "If this is a mistake FIX IT."
    )
  }
  if (is.null(ini.mu)) {
    ini.mu <- rep(0, nv)
  }
  if (is.null(ini.cov)) {
    ini.cov <- diag(1, nv)
  }
  if (is.null(ini.beta)) {
    ini.beta <- diag(0.00, nv)
  }
  if (is.null(ini.loadings)) {
    ini.loadings <- diag(1, nv)
  }

  # Preparing Measurement
  meas <- prep.measurement(
    values.load = ini.loadings,
    params.load = matrix(rep("fixed", (nv^2)), ncol = nv), # Flag for customization
    state.names = paste0("X", 1:nv),
    obs.names = paste0("X", 1:nv)
  )
  # Matrix Dynamics
  # By using 'letters' we are limited to 26 variables
  dynm <- prep.matrixDynamics(
    values.dyn = ini.beta,
    params.dyn = matrix(levels(interaction(letters[1:nv], 1:nv, sep = "_")),
      nv, nv,
      byrow = TRUE
    ),
    isContinuousTime = FALSE
  )
  # Prepping Initial Conditions
  initial <- prep.initial(
    values.inistate = ini.mu,
    params.inistate = rep("fixed", nv),
    values.inicov = ini.cov,
    params.inicov = matrix("fixed", nv, nv)
  )

  # Prepping Noise
  # Creating a covariance matrix to be estimated
  covs <- paste0("C_", 1:nv)
  m1 <- matrix(levels(interaction(covs[1:nv], 1:nv, sep = "_")), nv, nv, byrow = TRUE)
  m1[lower.tri(m1)] <- t(m1)[lower.tri(m1)]
  diag(m1) <- paste0("VAR_", 1:nv)
  mdcov <- prep.noise(
    values.latent = diag(1, nv),
    params.latent = m1, # Flag; change later post-testing for GVAR
    values.observed = diag(rep(0, nv)),
    params.observed = matrix(("fixed"), nv, nv)
  )
  # Creating Subsetted Dataset
  individualMods <- list()
  for (k in unique(dataframe[, id])) {
    pseudo.dat <- dataframe[dataframe[, id] == k, ]
    rawdata <- dynr.data(pseudo.dat, id = id, time = time, observed = c(paste0("X", 1:nv)))
    # Model Fitting
    model <- dynr.model(
      dynamics = dynm, measurement = meas,
      noise = mdcov, initial = initial, data = rawdata,
      outfile = paste0(dir, "ModSubj", k, ".c")
    )
    ##
    # Plot Formula Here
    ##
    # Setting upper/lower bounds for parameter estimates.
    # I've found that this helps speed things up.
    model$ub[1:nv^2] <- 1.00
    model$lb[1:nv^2] <- -1.00
    # Cooking the dynr model
    if (use.mi == TRUE) {
      results <- dynr.mi(model,
        which.aux = aux,
        which.lag = paste0("X", 1:nv), lag = 1,
        which.lead = NULL, lead = 0,
        m = 5, iter = 10,
        imp.obs = FALSE, imp.exo = TRUE,
        diag = TRUE, Rhat = 1.1,
        conf.level = 0.95,
        verbose = FALSE, seed = 1472
      )
      obj <- results$estimation.result
    } else {
      results <- dynr.cook(model)
      obj <- summary(results)$Coefficients
    }
    # Significance Testing of Transition Matrix Coefficients
    effects <- matrix(obj[1:nv^2, 1], nv, nv, byrow = FALSE)
    sigs <- matrix(as.numeric(obj[1:nv^2, 6] < alpha), nv, nv)
    OG.betas <- effects * sigs
    # Significance Testing of Residual Covariance Matrix
    sig.cov <- matrix(NA, nv, nv)
    sig.cov[lower.tri(sig.cov, diag = TRUE)] <- as.numeric(obj[(nv^2 + 1):nrow(obj), 6] < alpha)
    sig.cov[upper.tri(sig.cov)] <- t(sig.cov)[upper.tri(sig.cov)]
    res.cov <- matrix(NA, nv, nv)
    res.cov[lower.tri(res.cov, diag = TRUE)] <- obj[(nv^2 + 1):nrow(obj), 1]
    res.cov[upper.tri(res.cov)] <- t(res.cov)[upper.tri(res.cov)]
    OG.cov <- sig.cov * res.cov
    PDC <- matrix(NA, nv, nv)
    PCC <- matrix(NA, nv, nv)
    OG.kappa <- solve(OG.cov)
    # Calculating the PDC/PCCs
    # See Wild, Eichler, et al (2010)
    for (p in 1:nv) {
      for (q in 1:nv) {
        PCC[p, q] <- -((OG.kappa[p, q]) / (sqrt(OG.kappa[p, p] * OG.kappa[q, q])))
        PDC[p, q] <- (OG.betas[p, q]) / (sqrt(OG.betas[p, p] * OG.kappa[q, q] + (OG.betas[p, q]^2)))
      }
    }

    individualMods[[k]] <- list(
      Betas = OG.betas,
      ResidCov = OG.cov,
      PCC = PCC,
      PDC = PDC,
      Res = results
    )
  }
  class(individualMods) <- c(
    "dynrVar",
    class(individualMods)
  )
  return(individualMods)
}
