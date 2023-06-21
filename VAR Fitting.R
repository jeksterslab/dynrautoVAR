# Automatically Fit VAR(1) Models Based on Input Data
  # 1. Take in data.frame object
  # 2. Separate by ID
  # 3. Fit Individual Models and Store as a list object
    # 3a. Should interface with subgrouping code
    # 3b. Should be able to call Cook's Distance code
  # 4. Should be able to interface with dynr.mi()
# data is the data.frame
# nv = number of variables
# ID = name of id variable
# time = name of time variable
# dir = where you want shit saved
# mu = mean vector
# ini.cov = initial covariance matrix
# beta = initial beta matrix
# loadings = factor loadings of variables
# If anything is NULL, we just use NULL matrices or 0's throughout.
dynr.var = function(data = NULL, nv = NULL, ID = NULL, 
                    time = NULL, dir = NULL, p.val = 0.05,
                    ini.mu = NULL, ini.cov = NULL,
                    beta = NULL, loadings = NULL,
                    use.mi = FALSE, aux = NULL){
  # User Argument Handling
    if(is.null(ID)){
      message('We are assuming you are only fitting to N = 1 subject. If this is a mistake FIX IT.')
    }else{ID = ID}
    if(is.null(time)){
      message('Please specify a time variable.')
      break
    }else{time = time}
    if(is.null(dir)){
      message('We are using the default directory.')
      dir = getwd()
    }else{dir = dir}
    if(is.null(ini.mu)){
      ini.mu = rep(0, nv)
    }else{ini.mu = ini.mu}
    if(is.null(ini.cov)){
      ini.cov = diag(1, nv)
    }else{ini.cov = ini.cov}
    if(is.null(beta)){
      beta = diag(0.00, nv)
    }else{beta = beta}
    if(is.null(loadings)){
      loadings = diag(1, nv)
    }else{loadings = loadings}
  
  # Preparing Measurement
    meas <- prep.measurement(
      values.load = loadings, 
      params.load = matrix(rep("fixed", (nv^2)), ncol = nv), # Flag for customization
      state.names = paste0('X', 1:nv),
      obs.names = paste0('X', 1:nv)
    )
  # Matrix Dynamics
    # By using 'letters' we are limited to 26 variables
  dynm = prep.matrixDynamics(values.dyn = beta, 
                             params.dyn = matrix(levels(interaction(letters[1:nv],1:nv,sep='_')), 
                                                 nv, nv, byrow = TRUE), 
                             isContinuousTime = FALSE
                            )
  # Prepping Initial Conditions
  initial <- prep.initial(
    values.inistate = ini.mu,
    params.inistate = rep('fixed', nv),
    values.inicov = ini.cov,
    params.inicov = matrix('fixed', nv, nv)
  )
  
  # Prepping Noise
    # Creating a covariance matrix to be estimated 
    covs = paste0('C_', 1:nv)
    m1 = matrix(levels(interaction(covs[1:nv],1:nv,sep='_')), nv, nv, byrow = TRUE)
    m1[lower.tri(m1)] <- t(m1)[lower.tri(m1)]
    diag(m1) = paste0('VAR_',1:nv)
    mdcov <- prep.noise(
      values.latent = diag(1,nv),
      params.latent = m1, # Flag; change later post-testing for GVAR
      values.observed = diag(rep(0, nv)), 
      params.observed = matrix(('fixed'), nv, nv)
    )
  # Creating Subsetted Dataset
  individualMods = list()
  for(k in unique(data[,ID])){
    pseudo.dat = data[data[,ID] == k,]
    rawdata <- dynr.data(pseudo.dat, id = ID, time = time, observed = c(paste0('X', 1:nv)))
    # Model Fitting
    model <- dynr.model(dynamics = dynm, measurement = meas,
                      noise = mdcov, initial = initial, data = rawdata,
                      outfile = paste0(dir, "ModSubj",k,".c"))
    ##
    # Plot Formula Here
    ##
    # Setting upper/lower bounds for parameter estimates.
      # I've found that this helps speed things up.
      model$ub[1:nv^2] = 1.00
      model$lb[1:nv^2] = -1.00
    # Cooking the dynr model
    if(use.mi == TRUE){
      results <- dynr.mi(model, which.aux = aux, 
                        which.lag = paste0('X', 1:nv), lag = 1,
                        which.lead = NULL, lead = 0,
                        m = 5, iter = 10, 
                        imp.obs=FALSE, imp.exo=TRUE,
                        diag = TRUE, Rhat=1.1,
                        conf.level = 0.95,
                        verbose = FALSE, seed=1472)
      obj = result$estimation.result
    }else{results <- dynr.cook(model); obj = summary(results)$Coefficients}
    # Significance Testing of Transition Matrix Coefficients
      effects = matrix(obj[1:nv^2,1], nv, nv, byrow = FALSE)
      sigs = matrix(as.numeric(obj[1:nv^2,6]<p.val), nv, nv)
      OG.betas = effects * sigs
    # Significance Testing of Residual Covariance Matrix
      sig.cov = matrix(NA, nv, nv)
      sig.cov[lower.tri(sig.cov,diag=TRUE)] = as.numeric(obj[(nv^2+1):nrow(obj),6] < p.val)
      sig.cov[upper.tri(sig.cov)] = t(sig.cov)[upper.tri(sig.cov)]
      res.cov = matrix(NA, nv, nv)
      res.cov[lower.tri(res.cov,diag=TRUE)]=obj[(nv^2+1):nrow(obj),1]
      res.cov[upper.tri(res.cov)] = t(res.cov)[upper.tri(res.cov)]
      OG.cov = sig.cov*res.cov
      PDC = matrix(NA, nv, nv)
      PCC = matrix(NA, nv, nv)
      OG.kappa = solve(OG.cov)
      # Calculating the PDC/PCCs
      # See Wild, Eichler, et al (2010)
      for(p in 1:nv){
        for(q in 1:nv){
          PCC[p,q] = -((OG.kappa[p,q])/(sqrt(OG.kappa[p,p]*OG.kappa[q,q])))
          PDC[p,q] = (OG.betas[p,q])/(sqrt(OG.betas[p,p] * OG.kappa[q,q] + (OG.betas[p,q]^2)))
        }
      }
      
      individualMods[[k]] = list(Betas = OG.betas, 
                                 ResidCov = OG.cov,
                                 PCC = PCC,
                                 PDC = PDC,
                                 Res = results)
  }
    return(individualMods)
}