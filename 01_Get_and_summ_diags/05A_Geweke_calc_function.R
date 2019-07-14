#### Function to go in and pull out Geweke's diagnostic from MCMC chains
####   using an initial window of 0.3 instead of 0.1 as in the main diagnostics
get.gew<-function(chains){
    if(length(chains)>2){
      chains<-chains[1:2]
    }
    chains<-check.chains(chains)
    # Get the burnin for each chain individually
    est_burn_all<-lapply(chains, FUN=function(x) getburnin(x))
    ## Beacuse these chains are from the same analysis, we should probably use the same burnin for them, here using the max of the estimated burnins
    est_burn<-max(unlist(est_burn_all)) 
    
    if(est_burn>=90){
        res<-list(NA, "Chains have not adequately passed burn-in")
    }else{

      # need to convert to coda's format
      chains_for_coda<-rwty_to_coda(chains, burn=est_burn, cutoff=100)
      
      ## Use coda to calculate the Geweke diagnostic - only gives z-score and not p-value
      coda_geweke<-lapply(chains_for_coda, geweke.diag, frac1=0.3, frac2=0.5)
      
      # Use Gewekes.z.to.p function to convert z scores to p values
      gewekes_p<-lapply(coda_geweke, function(x) Gewekes.z.to.p(x$z))
      
      # Identify parameters with a Geweke p < 0.05
      params_fail_Geweke<-lapply(gewekes_p, function(x) which(x<0.05))
      # Identify which chains have any parameters that are above 0.05
      chains_pass_Geweke<-which(lapply(params_fail_Geweke, length)==0)
      
      ## Make a list describing which tests were passed
      converge_summ<-list()
      for(i in 1:length(chains)){
        pass_Geweke<-i %in% chains_pass_Geweke
        names(pass_Geweke)<-"Pass Geweke"
        converge_summ[[i]]<-pass_Geweke
      }
      
      converge_summ_table<-do.call(rbind, lapply(converge_summ, function(x) x))
      
      res<-list(converge_summ_table, est_burn_all, coda_geweke, gewekes_p, params_fail_Geweke)
      names(res)<-c("Convergence_Summary", "Burnin",  "Geweke_z_scores", "Gewekes_p", "Params_fail_Geweke")
  } 
       return(res)
}




## Write a tryCatch version of this function so that any errors don't trip up everything
get.gew.catch<-function(chains){
  return(tryCatch(get.gew(chains), error=function(e) {
    fail<-list("FAILED", names(chains), e)
    return(fail)
  }))
}