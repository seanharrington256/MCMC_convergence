#### Function to extract a number of convergence diagnostics from chains from a single analysis
# The full pipeline is currently only set up to handle analyses that have exactly 2 runs, so we include a statement to prune down any analyses with more runs to only 2 near the beginning of this.
# requires the connection to the database to be open already
get.conv.summ<-function(chains, topo_ESS_reps=20){
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
      ## Get the final ESS values for and autocorrelation times for all parameters in the chain
      final_ess_ACT_vals<-lapply(chains, final_ess_ACT_all, burn=est_burn)
      
      # Split out just the ESS values and autocorrelation times separately for downstream fucntions
      final_ess<-lapply(final_ess_ACT_vals, function(x) x[[1]][1,])
      params_ACT<-lapply(final_ess_ACT_vals, function(x) x[[1]][2,])
      
      # Get out the values of each parameter through the chains
      vals<-lapply(final_ess_ACT_vals, function(x) x[[2]])
      # Run Shapiro-Wilk normality test 
      shap_test<-lapply(vals, function(x) lapply(x, shapiro.test))
      
      ## Get out the tree length through the chain
      TL<-get_TL.catch(chains, burnin=est_burn)
      
      
      ## See when all parameters get above ess of 200 and what the final ess value is for each chain
      ## make a loop to get the generations to ESS of 200 only for those chains that end up higher than 200
      all_params_200<-list()
      for(i in 1:length(chains)){
        failed_pre<-which(final_ess[[i]]<200) # find which chains have a final ESS that is less than 200
        failed<-rep(NA, length(failed_pre))  # repeat NA for each parameter that has ESS < 200
        names(failed)<-names(failed_pre)  # Give these the names of the failed parameters
        passed_params<-names(which(final_ess[[i]]>=200))  # Get the names of the parameters that have ESS values above 200
        gens_passed<-sapply(passed_params, FUN=function(x) gens_to_200(chains=chains[[1]], burn=est_burn, parameter=x))  # For the parameters up above 200, get the gens to ESS of 200
        all_params_200[[i]]<-c(gens_passed, failed)  # Concatenate the passed and failed parameters together
        all_params_200[[i]]<-all_params_200[[i]][names(final_ess[[i]])]  # Reorder the parameters to match the ESS value ordering
      }
      names(all_params_200)<-names(final_ess) # assign the names of the chains to the elements of the list
  
      ## Topo final pseudo ESS
      topo_final_pseudo_ess_ACT<-lapply(chains, topo.pseudo.ess.ACT.cutoff, burn = est_burn, cutoff=100, n = topo_ESS_reps)
      
      ## Split out into object of just topo pseudo ESS
      topo_final_pseudo_ess<-lapply(topo_final_pseudo_ess_ACT, function(x) x[[1]]) 

      
      # Get out the final approximate topo ESS
      approx_topo_ess_ACT<-lapply(chains, topological.approx.ess.cutoff.ACT, burn=est_burn)
      # Split into separate objects for topo ESS and topo autocorrelation
      approx_topo_ess<-lapply(approx_topo_ess_ACT, function(x) x[[1]])
      topo_ACT<-lapply(approx_topo_ess_ACT, function(x) x[[2]])
      
      
      ## Get the generations to approx topo ess of 200 only if the final approx topo ess is >= 200; otherwise, this is wasted time
      topo_200<-list()
      for(i in 1:length(chains)){
        if(approx_topo_ess[[i]]$approx.ess<200){
          topo_200[[i]]<-NA
        } else{
          topo_200[[i]]<-topo_gens_to_200(chains[[i]], burn=est_burn)
        }
      }
      #put names of "topo" onto each element of the topo_200 list
      for(i in 1:length(topo_200)){
          names(topo_200[[i]][[1]])<-"topo"
      }
      
  
      # For each chain, concatenate the parameter and topology gens to 200 ESS
      params_topo_200<-list()
      for(i in 1:length(all_params_200)){
        params_topo_200[[i]]<-c(unlist(all_params_200[[i]]), topo_200[[i]][[1]])
      }
      
      # For each chain, concatenate the parameter ess with the topology ess
      approx_topo_ess_val<-lapply(approx_topo_ess, function(x) x$approx.ess)
        for(i in 1:length(approx_topo_ess_val)){   # rename the elements of the median topo ess vals
            names(approx_topo_ess_val[[i]])<-"topo"
        }
      params_topo_ESS_final<-list()
      for(i in 1:length(final_ess)){
        params_topo_ESS_final[[i]]<-c(unlist(final_ess[[i]]), approx_topo_ess_val[[i]])
      }
      
      
      # Identify parameters that have ESS less than 200 in each chain - exclude topo ESS, because this is not a standard check
      params_less_200<-lapply(final_ess, function(x) which(x<200))
      # Identify which chains have no parameters that are less thatn 200 (i.e., pass the ESS 200 test)
      chains_pass_200<-which(lapply(params_less_200, length)==0)
      
      # Identify which chains have pseudo topo ess less than 200 in each chain
      med_topo_pseudo_ess<-lapply(topo_final_pseudo_ess, function(x) x$median.ess) #Pull out just the median topo ess vals
      chains_pass_pseudo_topo_ESS<-lapply(med_topo_pseudo_ess, function(x) x>=200)
      
      
      # Identify which chains have approx topo ess less than 200 in each chain
      chains_pass_approx_topo_ESS<-lapply(approx_topo_ess, function(x) x$approx.ess>=200)

      # Get asdsf and correlations among chains
      asdsf_sf_corrs<-get.sf.cors.asdsf(chains, burn=est_burn)
      
      # For PSRF, need to convert to coda's format
      chains_for_coda<-rwty_to_coda(chains, burn=est_burn, cutoff=100)
      
      ## Use coda function to calculate the PSRF/Gelman and Rubin diagnostic - defaults for all, except autoburnin
      PSRF<-gelman.diag(chains_for_coda, confidence = 0.95, transform=FALSE, autoburnin=FALSE, multivariate=TRUE)
      
      # Find out which PSRF are above 1.02
      PSRF_above_1.02<-which(PSRF$psrf[,1]>1.02)
      
      ## Use coda to calculate the Geweke diagnostic - only gives z-score and not p-value
      coda_geweke<-lapply(chains_for_coda, geweke.diag, frac1=0.1, frac2=0.5)
      
      # Use Gewekes.z.to.p function to convert z scores to p values
      gewekes_p<-lapply(coda_geweke, function(x) Gewekes.z.to.p(x$z))
      
      
      # Identify parameters with a Geweke p < 0.05
      params_fail_Geweke<-lapply(gewekes_p, function(x) which(x<0.05))
      # Identify which chains have only parameters that are above 0.05 (i.e., chains that pass Geweke's test for all parameters)
      chains_pass_Geweke<-which(lapply(params_fail_Geweke, length)==0)
      
      ## Get out the correlations among the parameters in the analysis
      param_cors<-lapply(chains, function(x) cor(x$ptable[,-1]))
      
      ## Make a list describing which tests were passed
      converge_summ<-list()
      for(i in 1:length(chains)){
        pass_ESS<-i %in% chains_pass_200
        names(pass_ESS)<-"Pass ESS"
        pass_ASDSF<-asdsf_sf_corrs$ASDSF<0.01
        names(pass_ASDSF)<-"Pass ASDSF"
        pass_sf_corr<-asdsf_sf_corrs$sf_corrs$cor>=0.9
        names(pass_sf_corr)<-"Pass SF corr"
        pass_PSRF<-length(PSRF_above_1.02)==0
        names(pass_PSRF)<-"Pass PSRF"
        pass_topo_pseudo_ESS<-chains_pass_pseudo_topo_ESS[[i]]
        names(pass_topo_pseudo_ESS)<-"Pass topo pseudo ESS"
        pass_approx_topo_ESS<-chains_pass_approx_topo_ESS[[i]]
        names(pass_approx_topo_ESS)<-"Pass Approx topo ESS"
        pass_Geweke<-i %in% chains_pass_Geweke
        names(pass_Geweke)<-"Pass Geweke"
        converge_summ[[i]]<-c(pass_ESS, pass_ASDSF, pass_sf_corr, pass_PSRF, pass_topo_pseudo_ESS, pass_approx_topo_ESS, pass_Geweke)
      }
      
      converge_summ_table<-do.call(rbind, lapply(converge_summ, function(x) x))
      
      res<-list(converge_summ_table, est_burn_all, params_topo_200, params_topo_ESS_final, params_ACT, topo_ACT, approx_topo_ess, params_less_200, asdsf_sf_corrs, PSRF, PSRF_above_1.02, coda_geweke, gewekes_p, params_fail_Geweke, topo_final_pseudo_ess, topo_200, param_cors, shap_test, TL)
      names(res)<-c("Convergence_Summary", "Burnin", "Generations_to_ESS_200", "ESS_values", "ACT params", "topo act", "Approx topo ESS","Params_ESS_less_200", "ASDSF_and_SF_corr", "PSRF", "Params_PSRF_less_1.02", "Geweke_z_scores", "Gewekes_p", "Params_fail_Geweke", "Topo pseudo ess confidence", "Approx_ESS_200_w_any_errors", "Param_corrs", "Shap_test", "TL")
  } 
       return(res)
}




## Write a tryCatch version of get.conv.summ so that any errors in this part don't trip up everything
conv.summ.catch<-function(chains, topo_ESS_reps=20){
  return(tryCatch(get.conv.summ(chains, topo_ESS_reps), error=function(e) {
    fail<-list("FAILED", names(chains), e)
    return(fail)
  }))
}

