model{
  
  ####### Priors
  for(a in 1:2){
#    meanOmega ~ dunif(0,1)
#    mu <- log(meanOmega/(1-meanOmega))
#    sigma ~ dunif(0,10)
#    tau <- pow(sigma, -2)
#    sigma2 <- pow(sigma, 2)  
    meanOmega[a] ~ dunif(0.01, 0.5)
    muOmega[a] <- log(meanOmega[a]/(1-meanOmega[a]))
    #mu[a] ~ dnorm(0,0.667)#log(meanOmega[a]/(1-meanOmega[a]))
     sigmaOmega[a] ~ dunif(0,3)
     tauOmega[a] <- pow(sigmaOmega[a], -2)
     sigma2Omega[a] <- pow(sigmaOmega[a],2)

    muGamma[a] ~ dnorm(0, 0.01)
    
    beta1[a] ~ dnorm(0, 0.001) 
    beta2[a] ~ dnorm(0, 0.001)
  }
  
  
  ####### Survival Model
  for(s in 1:nSites){
    for(a in 1:2){
      #      logit(omega[s,a]) <- mu[a] + eps[s]  #+ COV[s,t]   # survival: see Kery & Schaub p. 186
      logit(omega[s,a]) <- muOmega[a] + beta1[a]*elev[s]  + epsOmega[s,a]
      epsOmega[s,a] ~ dnorm(0, tauOmega[a])
    }      
    #     logit(omega[s]) <- mu + eps[s]
    #     eps[s] ~ dnorm(0, tau)
    
  }
  
  ###### Gains Model
  #log(gamma) <- muGamma #+ betaGamma*COV[i,t]
  
  ###### Priors
  #muGamma~dnorm(0,0.001)

  for(a in 1:2){
    #lambda[a] ~ dunif(0,500)
    rInit[a] ~ dunif(0,2000)            # initial abundance                            
    pInit[a] ~ dunif(0,1)
      
    for(s in 1:nSites){
      log(gamma[s,a]) <- muGamma[a] + beta2[a]*elev[s] + epsGamma[s,a]
      epsGamma[s,a] ~ dnorm(0, tauGamma[a])
    }
    tauGamma[a] <- pow(sigmaGamma[a], -2)
    sigmaGamma[a] ~ dunif(0,3)
  }
  
  #    sigma[a] ~ dunif(0,3)
  #    tau[a] <- pow(sigma[a], -2)
  #    sigma2[a] <- pow(sigma[a],2)  
  
  
  #Detection model
  for(s in 1:nSites) {
    for(t in 1:nYears){
      for(a in 1:2){
        y[s,t,a,1] ~ dbin(Q[a], N[s,t,a])
        y[s,t,a,2] ~ dbin(Q[a],(N[s,t,a]-y[s,t,a,1]))
        y[s,t,a,3] ~ dbin(Q[a],(N[s,t,a]-y[s,t,a,1]-y[s,t,a,2]))
      }
    }
  }
  
  ###### Detection Priors
  for(a in 1:2){
    Q[a] ~ dunif(0,1)
  }
    
  #for(a in 1:2){
  #  for(p in 1:3){
  #    Q[a,p] ~ dunif(0,1)
  #  }
  #}
    
  ##### N section
  for(s in 1:nSites){
    for(a in 1:2){
      # Year 1 - initial abundance
      N[s,1,a] ~ dnegbin(pInit[a], rInit[a])
    }
  }
  
  # Year 2+ 
  for(s in 1:nSites){
    for(t in 2:nYears) {
      for( a in 1:2){
        #Estimate survival
        S[s,t-1,a] ~ dbin(omega[s,a], N[s,t-1,a])
        #S[s,t-1,a] ~ dbin(omega[s], N[s,t-1,a])
        #Estimate "gains"
        G[s,t-1,a] ~ dpois(gamma[s,a])
      }
      #Sum the two above to get total N at each site i in each year t
      N[s,t,1] <- G[s,t-1,1]  
      N[s,t,2] <- S[s,t-1,1] + S[s,t-1,2] + G[s,t-1,2]
    }
  }
  
  #sum up the number of individuals in all sites to estimate annual
  #total N
  for (t in 1:nYears){
    for(a in 1:2){ 
      Ntotal[t,a] <- sum(N[,t,a])
    } 
  }
}
