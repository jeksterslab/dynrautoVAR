  # Cook's Distance in Dynr
Cooks.Dynr = function(data = NULL, ID = NULL, time = NULL, ini.beta = NULL, dir = NULL) {
  bet = ini.beta
  nvar = ncol(data)-2
  names(data) = c(paste0('X', 1:nvar))
  pseudo.data = data[,-c(ID, time)]
  pseudo.data$id = rep(1, nrow(pseudo.data))
  pseudo.data$time = 1:nrow(pseudo.data)
  rawdata <- dynr.data(data, id = "id", time = "time", 
                       observed = c(paste0('X', 1:nvar)))
  load.mat = diag(1, nvar)
  meas <- prep.measurement(
    values.load = load.mat, 
    params.load = matrix(rep("fixed", (nvar^2)), ncol = nvar),
    state.names = paste0('X', 1:nvar),
    obs.names = paste0('X', 1:nvar)
  )
  letters = c('a','b','c','d','e','f','g','h','i','j','k')
  dynm = prep.matrixDynamics(values.dyn = bet, params.dyn = matrix(levels(interaction(letters[1:nvar],1:nvar,sep='_')), nvar, nvar, byrow = TRUE), isContinuousTime = FALSE)
  means = paste0('mu', 1:nvar)
  initial <- prep.initial(
    values.inistate = rep(0, nvar),
    params.inistate = rep('fixed', nvar),
    values.inicov = load.mat,
    params.inicov = matrix('fixed', nvar, nvar))
  covs = c('C_1','C_2','C_3','C_4','C_5','C_6','C_7','C_8','C_9','C_10')
  m1 = matrix(levels(interaction(covs[1:nvar],1:nvar,sep='_')), nvar, nvar, byrow = TRUE)
  m1[lower.tri(m1)] <- t(m1)[lower.tri(m1)]
  m1=matrix('fixed', nvar, nvar)
  diag(m1) = paste0('VARIANCE_',1:nvar)
  mdcov <- prep.noise(
    # values.latent = cov(pseudo.data[,1:nvar]),
    # params.latent = m1,
    values.latent = diag(1,nvar),
    params.latent = m1,
    values.observed = diag(rep(0, nvar)), 
    params.observed = matrix(('fixed'), nvar, nvar))
  
  model <- dynr.model(dynamics = dynm, measurement = meas,
                      noise = mdcov, initial = initial, data = rawdata,
                      outfile = paste0(dir, "/OGModel10.c"))
  model$ub[1:nvar^2] = 1.00
  model$lb[1:nvar^2] = -1.00
  results <- dynr.cook(model)
  obj = summary(results)
  effects = matrix(obj$Coefficients[1:nvar^2,1], nvar, nvar, byrow = FALSE)
  sigs = matrix(as.numeric(obj$Coefficients[1:nvar^2,6]<0.05), nvar, nvar)
  OG.betas = effects * sigs
  OG.inv.Hess = results$'inv.hessian'[1:nvar^2,1:nvar^2]
    # LOO: Leave One Out Loop
    B.i = list(); B.i.sig = list()
    index = 1
      for(i in unique(data[,ID])){
        B.i.signs = NULL
        temp.dat = data[data[,ID]!=i,]
        temp.dat = temp.dat[,-c(ID, time)]
        temp.dat$id = rep(1, nrow(temp.dat))
        temp.dat$time = 1:nrow(temp.dat)
        rawdata <- dynr.data(temp.dat, id = "id", time = "time", 
                     observed = c(paste0('X', 1:nvar)))
        model <- dynr.model(dynamics = dynm, measurement = meas,
                      noise = mdcov, initial = initial, data = rawdata,
                      outfile = paste0(dir, "/LOO1_",index,".c"))
        model$ub[1:nvar^2] = 1.00
        model$lb[1:nvar^2] = -1.00
        results <- dynr.cook(model)
        obj = summary(results)
        effects = matrix(obj$Coefficients[1:nvar^2,1], nvar, nvar, byrow = FALSE)
        sigs = matrix(as.numeric(obj$Coefficients[1:nvar^2,6]<0.05), nvar, nvar)
        B.i[[index]] = effects * sigs
        index = index + 1
        message(paste0('I am done with Subject ', index-1))
      }
    p = ncol(data[,-c(ID, time)])
    ps = p^2 + p + ((p)^2 - p)/2
    CookD = list()
      for(i in 1:length(B.i)){
        CookD[[i]] = t(c(OG.betas - B.i[[i]])) %*% 
          solve(OG.inv.Hess) %*% c(OG.betas - B.i[[i]])/ps
      }
    return(CookD)
}
attempt3 = Cooks.Dynr(data = rdat, ID = 6, time = 7, ini.beta = phi0, dir = '~/Desktop/CookD/')
# Inclusion of 1 WRONG person [Group A + 1 Group B]
# RESULTS = c(0.01797025, 0.1984687, 0.05177024, 0.1230428, 0.1211222,
# 0.03256747, 0.0368129, 0.03601584, 0.01800689, 0.04814378,
# 1.166811)
# NOTE:
  # `attempt` object will be the one with a fuzzy subject while the numbers above are
  # for the inclusion of a single "wrong" subject
# For Thursday's "talk", I will show the effects of introducing one of those biased people and if our algorithms can filter them out compared to 

plot(1:11, RESULTS, type = 'p')
