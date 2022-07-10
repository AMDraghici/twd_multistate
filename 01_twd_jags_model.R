model{
  ##### Likelihood #####
  for(i in 1:Nmax){
    
    ### Data augmentation
    ## exists[i]=1 if individual i from the superpopulation is
    ## part of the actual population
    exists[i] ~ dbern(xi)
    
    ### Recruitment
    ## born[i,j]=1 if individual i entered the population prior to
    ## (and including) year j
    ## By definition born[i,Nyears]=1 since individuals must enter before
    ## the last year to be part of the population. This should be set in the
    ## data.
    born[i,1] ~ dbern(beta[1])
    for(j in 2:(Nyears-1)){
      born[i,j] ~ dbern(born[i,j-1] + (1-born[i,j-1]) * beta[j])
    }
    
    ### Mortality
    ## dead[i,j]=1 if individual i died prior to year j
    ## By definition dead[i,1]=0 since no individual can die
    ## before the first year. This should be set in the data.
    ## Note that dead[i,j]=0 does not indicate that individual i
    ## is alive in year j (see below). It simply means that
    ## the individual did not die before year j.
    for(j in 1:(Nyears-1)){
      dead[i,j+1] ~ dbern(dead[i,j] + (1-dead[i,j]) * born[i,j] * (1-phi[j]))
    }
    
    ### Alive or dead
    ## alive[i,j]=1 if individual i is part of the actual population,
    ## was recruited in or before year j and has not yet died
    for(j in 1:Nyears){
      alive[i,j] <- exists[i] * born[i,j] * (1-dead[i,j])
    }
    
    ### Movement
    ## block[i,j,k] indicates which block individual i was in
    ## on the k-th day of year j
    for(j in 1:Nyears){
      ## Movement over remaining occasions in year j
      for(k in 1:(Ndays[j])){
        block[i,j,k] ~ dcat(psi[1:Nblock])
      }
    }

    ### Capture
    ## capture[i,j,k]=1 if individual i was detected on the k-th day of
    ## year j
    for(j in 1:Nyears){
      for(k in 1:Ndays[j]){
        capture[i,j,k] ~ dbern(p[j] * alive[i,j] * overlap[j,k,block[i,j,k]])
      }
    }
  }
  
  ##### Prior Distributions #####
  
  ## Data augmentation
  xi ~ dunif(0,1)
  
  ## Recruitment
  ## beta[j] is the probability that an individual is recruited in year j
  ## given that it was not recruited prior to year j
  ## beta[1] is the probability that an individual is alive in the
  ## first year
  for(j in 1:(Nyears-1)){
    beta[j] ~ dunif(0,1)
  }

  ## Mortality
  ## phi[j] is the probabiity that an individual survives between occasions
  ## j and j+1 given that it was alive on occasion j
  for(j in 1:(Nyears-1)){
    # Logit-Normal Prior
    logit(phi[j]) <- phi.norm[j]
    
    # Normal draw from informative prior 
    phi.norm[j] ~ dnorm(phi.xi,phi.tau)

  }
  
  #Hyperprior (assume prob ~ 0.9 with most density around 0.75 to 0.95)
  phi.xi ~ dnorm(2.2,1) #logit(0.9) ~ 2.2 and sig  pm 1
  phi.tau ~ dt(0,0.55,4)T(0,)
  
  ## Capture
  for(j in 1:Nyears){
      # Logit - Normal Prior
      logit(p[j]) <- p.norm[j]
      
      # Normal draw from informative prior 
      p.norm[j] ~ dnorm(p.xi, p.tau)
    
  }
  
  #Hyperprior (assume prob ~ 0.25 with most density around 0.1 to 0.5)
  p.xi ~ dnorm(-1.1, 1) #logit(0.9) ~ 2.2 and sig  pm 1
  p.tau ~ dt(0,0.55,4)T(0,)
  
  ## Movement
  psi ~ ddirich(alpha_psi[1:Nblock]) 
  
  ##### Derived parameters #####
  
  ## Population size
  for(j in 1:Nyears){
    abundance[j] <- sum(alive[1:Nmax,j])
  }
  
  ## Recruitment at time t
  for(j in 1:Nyears){
    recruited[j] <- sum(born[1:Nmax,j] * exists[1:Nmax])
  }
  
  
}
