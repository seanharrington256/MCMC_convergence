The scripts and files in this directory extract convergence diagnostics from MrBayes chains, summarize these diagnostics, set up reanalyses, and then get diagnostics and summaries for these reanalyses as well


Folder 01_Get_and_summ_diags - contains scripts to get diagnostics from chains and summarize these diagnostics afterwards
	01_R_functions_calculate_diagnostics_V2.R - contains functions to read in chains and calculate diagnostics
		read.compressed.chain -- reads in chains from a .tar.gz archive
		read.comp.catch -- a trycatch version of read.compressed.chain that will continue in the event of an error so that loops over this function can continue
		nums_from_txt -- extract numerical values from a txt file, e.g., numbers from a Nexus file
		getburnin -- identifies the burnin of a chain
		calc.ess.act -- extract ESS and autocorrelation time for a parameter - modified from Tracerer package
		param_vals_ess -- extract parameter values, ESS, and autocorrelation times - modified from rwty function makeplot.param
		param_vals -- get only parameter values
		final_ess_ACT_all -- get the final ESS and autocorrelation time - also spits out the parameter values at each generation through the chain
		gens_to_200 -- get the number of generations to an ESS value of 200, as multiples of 10% of the chain length
		tree.ess.tracer -- modification of rwty function tree.ess to use the tracerer function for calculating ESS and autocorrelation time
		tree.ess.tracer.multi -- modification of rwty tree.ess.multi to use tree.ess.tracer function
		topo.pseudo.ess.ACT.cutoff -- modification of rwty topological.pseudo.ess to take burnin as a percentage, include a cutoff for the maximum generation to go up to, and also to output autocorrelation times
		topo_pseudo_gens_to_200 -- get number of generations to a topological pseudo ESS of 200
		topological.autocorr.cutoff -- edits to rwty topological.autocorr to add in a cutoff
		topological.approx.ess.cutoff.ACT -- modification of rwty topological.approx.ess to take burnin as a percentage and work from burnin to cutoff value and return autocorrelation time as well
		approx.ess.cutoff.catch -- trycatch version of topological.approx.ess.cutoff.ACT that cases when there are not enough samples to evaluate approx ESS, the function will not stop
		topo_gens_to_200 -- get the number of generations to an approximate topological ESS value of 200
		get.sf.cors.asdsf -- get the correlations in split frequencies among chains and ASDSF - takes a few pieces from rwty makeplot.splitfreq.matrix and from https://stackoverflow.com/questions/45825685/correlations-for-pairs-of-combinations
		rwty_to_coda -- convert the parameters of chains from rwty to an mcmc class
		Gewekes.z.to.p -- convert Geweke's z scores to Geweke's p values - uses code from https://github.com/mikeryanmay/bonsai/blob/master/package/R/NumericalParameterRefClass.R
		get.acc.rates -- takes the lines of a log file and the indices of the lines in this object that contain the acceptance rates for a chain and returns the acceptance rates for each parameter in the chain
		get_TL -- get tree length
		get_TL.catch -- trycatch version of get_TL
		get_info_comp -- takes a compressed analysis that includes chains, nexus file, etc, and reads in the chains and various information (number of taxa, characters, generations, chains, runs, sampling frequency; acceptance rates; credible set size)
		read.comp.info.catch -- trycatch version of get_info_comp
	02_Wrapper_func_calc_diags.R - function and trycatch version to extract convergence diagnostics using functions from 01_R_functions_calculate_diagnostics_V2.R
	03_Calc_diags.R - script to uncompress folders containing analyses, read in chains, delete uncompressed folder, and get convergence diagnostics from a set of .tar.gz archives containing the output from MrBayes analyses
			utilizes the functions from 01_R_functions_calculate_diagnostics_V2.R and 02_Wrapper_func_calc_diags.R
			the object "storage" needs to be set to the location of the full set of MrBayes output stored as compressed .tar.gz archives
	04_Analyze_diagnostic_output.R - R script that summarizes the diagnostics output from 03_Calc_diags.R
			lines 4-7 specify directories that need to exist for output or to contain diagnostics from 03_Calc_diags.R
	05A-05C are scripts to calculate Geweke's diagnostic for parameters using an initial window of 0.3 rather than 0.1 done to this point
		05A_Geweke_calc_function.R - R script with function to extract Geweke's diagnostic from chains with 0.3 initial window
		05B_Get_geweke_wide.R - R script to use function from 01A_01_Geweke_calc_function.R and extract out Geweke's over compressed chains
			the object "storage" needs to be set to the location of the full set of MrBayes output stored as compressed .tar.gz archives
		05C_Process_Geweke_wide_results.R - R script to process the output from 01A_02_Get_geweke_wide.R
			working directory should be set to wherever output from 05B_Get_geweke_wide.R was stored
		
		
Folder 02_Reanalyses - contains the scripts to prep MrBayes reanalyses and summarize the output from these
	Folder Large_set -- the large set of analyses that were reanalyzed using different heating, nst=mixed substitution model, or +G rather than +I+G
		files ending _LARGE.csv contain the names of the files to reanalyze
		02_01_Get_chains_for_reanalysis.R - R script to pull out the chains that need to be reanalyzed into separate folders from the folder that contains all chains
			the object "storage" needs to be set to the location of the full set of original MrBayes output stored as compressed .tar.gz archives
		02_02_Prep_reestimate.R - R script to pull the necessary files out from the archives copied over in the previous script and get new Nexus files ready to run in MrBayes, including generating SLURM files to run this on a cluster
				lines 6-14 specify directories that need to exist and what they need to contain or will contain
				New Nexus files have a new heating, the substitution model changed to Nst=mixed, or +I+G in the substitution model changed to +G only
					heating is changed based on rules from 04_Analyze_diagnostic_output.R and read in as an R data file - if the acceptance rate for swaps between the cold and first heated chain was >0.5, heating is doubled, if not, heating is halved.
						For analyses that did not have log files to extract chain swap info from, decrease heating--of the analyses that had this info, decreasing was more common
		Diagnostics are then calculated for each of these sets of reanalyses using the same 01_ through 03_ R scripts described above from the 01_Get_and_summ_diags folder
		02_03_Compare_diags_LARGE.R - R script that compares the failed diagnostic tests of the original analyses to the reanalyses of the large set of datasets
			Additionally requires the .RData file of diagnostics calculated for the original chains from script 04_Analyze_diagnostic_output.R above
			lines 10-17 specify directories and what they need to contain
		02_04_RF_distances_LARGE.R - R script to get the RF distances between original trees and reestimated trees
			requires functions from 01_R_functions_calculate_diagnostics_V2.R as well as the original burnins for the sets of datasets to be reanalyzed (output from 04_Analyze_diagnostic_output.R) and the burnins for the reanalyses spit out by 02_03_Compare_diags_LARGE.R
			lines 19-26 specify directories and what they need to contain
		02_05_Summarize_RF_LARGE.R - R script to summarize the RF distances calculated in the previous script
			working directories are specified throughout script
		
	Folder Manual -- scripts for analyzing manually examined and reanalyzed datasets - much of this mirrors scripts from folder 02_Reanalyses with slight modifications 
		Manually_reanalyze.csv contains the names of the datasets to be reanalyzed
		These were reanalyzed following manual examination of chains and editing of Nexus files to change substitution models to HKY, change base frequencies from estimated to empirical, or to run analyses for 3x as long (30 M generations instead of 10)
		Following rerunning analyses, diagnostics were calculated for each set of manual reanalyses using the he same 01_ through 03_ R scripts described above from the 01_Get_and_summ_diags folder
		02M_01_Compare_diags_MANUAL.R - R script that compares the failed diagnostic tests of the original analyses to the reanalyses
			Additionally requires the .RData file of diagnostics calculated for the original chains from script 04_Analyze_diagnostic_output.R above
			lines 8-13 set up directories
		02M_02_Compare_topo_MANUAL.R - R script to compare RF distances between original and manually reanalyzed chains
			requires functions in 01_R_functions_calculate_diagnostics_V2.R, as well as the burnin for the orignal chains from script 04_Analyze_diagnostic_output.R and the burnin for each reanalysis from 02M_01_Compare_diags_MANUAL.R
			lines 7-12 set up directories
	
Additional notes:	
The csv file Study_names.csv indicates the study from which each analyzed dataset was obtained and is read into 04_Analyze_diagnostic_output.R to assign chains to the study they are from

The calculation of approximate ESS from the RWTY package incorporated into our coding pipeline allows for calculations using either all possible pairs of trees in a chain or a sample of distances. In initial tests, we found that when a subsample was used to calculate approximate topological ESS, the variation in estimates of the approximate topological ESS was high, and so we only calculated approximate topological ESS using all trees in the chain
Additionally, topological pseudo ESS calculations from RWTY use coda to calculate ESS values from tree distances. In initial tests, we identified differences between this implementation for calculating ESS values and the values returned by the package tracerer. Tracerer values are concordant with those output by the commonly used program Tracer, and so we modified topological pseudo ESS calculations to use tracerer's implementation of ESS rather than coda's




