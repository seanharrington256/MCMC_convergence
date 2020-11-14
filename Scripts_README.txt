More detailed description of scripts in this repository in directories 01_Get_and_summ_diags and 02_Reanalyses
Note that PhyLoTA pipeline has its own detailed readme and is not further detailed here

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
		get.compTable.asdsf -- this function gets the average standard deviation and a table of the clade frequencies from MCMC chains: modified version of rwty:::get.comparison.table. if there are 4 or more chains, get asdsf for 4 chains. if there are fewer than 4 chains, get asdsf for 2 chains
		get.sf.cors.asdsf -- get the correlations in split frequencies among chains and ASDSF - takes a few pieces from rwty makeplot.splitfreq.matrix and from https://stackoverflow.com/questions/45825685/correlations-for-pairs-of-combinations
		rwty_to_coda -- convert the parameters of chains from rwty to an mcmc class
		get.acc.rates -- takes the lines of a log file and the indices of the lines in this object that contain the acceptance rates for a chain and returns the acceptance rates for each parameter in the chain
		get_TL -- get tree length
		get_TL.catch -- trycatch version of get_TL
		get_info_comp -- takes a compressed analysis that includes chains, nexus file, etc, and reads in the chains and various information (number of taxa, characters, generations, chains, runs, sampling frequency; acceptance rates; credible set size)
		read.comp.info.catch -- trycatch version of get_info_comp
	02_Wrapper_func_calc_diags.R - function and trycatch version to extract convergence diagnostics using functions from 01_R_functions_calculate_diagnostics_V2.R
	03_Calc_diags.R - script to uncompress folders containing analyses, read in chains, delete uncompressed folder, and get convergence diagnostics from a set of .tar.gz archives containing the output from MrBayes analyses
			utilizes the functions from 01_R_functions_calculate_diagnostics_V2.R and 02_Wrapper_func_calc_diags.R
	04_Analyze_diagnostic_output.R - R script that summarizes the diagnostics output from 03_Calc_diags.R
	Study_names.csv indicates the study from which each analyzed dataset was obtained and is read into 04_Analyze_diagnostic_output.R to assign chains to the study they are from
		remaining csv files were output from 04_Analyze_diagnostic_output.R and specify the chains that we reanalyzed - used in 04_Analyze_diagnostic_output.R and in scripts described below
		
Folder 02_Reanalyses - contains the scripts to prep MrBayes reanalyses and summarize the output from these
	Folder Large_set -- the large set of analyses that were reanalyzed using different heating, nst=mixed substitution model, or +G rather than +I+G
		02_01_Get_chains_for_reanalysis.R - R script to pull out the chains that need to be reanalyzed into separate folders from the folder that contains all chains
		02_02_Prep_reestimate.R - R script to pull the necessary files out from the archives copied over in the previous script and get new Nexus files ready to run in MrBayes, including generating SLURM files to run this on a cluster
				New Nexus files have a new heating, the substitution model changed to Nst=mixed, or +I+G in the substitution model changed to +G only
					heating is changed based on rules from 04_Analyze_diagnostic_output.R and read in as an R data file - if the acceptance rate for swaps between the cold and first heated chain was >0.5, heating is doubled, if not, heating is halved.
						For analyses that did not have log files to extract chain swap info from, decrease heating--of the analyses that had this info, decreasing was more common
		Diagnostics are then calculated for each of these sets of reanalyses using the same 01_ through 03_ R scripts described above from the 01_Get_and_summ_diags folder
		02_03_Compare_diags_LARGE.R - R script that compares the failed diagnostic tests of the original analyses to the reanalyses of the large set of datasets
		02_04_RF_distances_LARGE.R - R script to get the RF distances between original trees and reestimated trees
		02_05_Summarize_RF_LARGE.R - R script to summarize the RF distances calculated in the previous script
		
	Folder Manual -- scripts for analyzing manually examined and reanalyzed datasets - much of this mirrors scripts from folder 02_Reanalyses with slight modifications 
		These were reanalyzed following manual examination of chains and editing of Nexus files to change substitution models to HKY, change base frequencies from estimated to empirical, or to run analyses for 3x as long (30 M generations instead of 10)
		02M_00_Calc_diags_cluster_MANUAL.R - This is a modified version of 03_Calc_diags.R to get diagnostics for chains that were reanalyzed after manual examination
		02M_01_Compare_diags_MANUAL.R - R script that compares the failed diagnostic tests of the original analyses to the reanalyses
		02M_02_Compare_topo_MANUAL.R - R script to compare RF distances between original and manually reanalyzed chains

	

The calculation of approximate ESS from the RWTY package incorporated into our coding pipeline allows for calculations using either all possible pairs of trees in a chain or a sample of distances. In initial tests, we found that when a subsample was used to calculate approximate topological ESS, the variation in estimates of the approximate topological ESS was high, and so we only calculated approximate topological ESS using all trees in the chain
Additionally, topological pseudo ESS calculations from RWTY use coda to calculate ESS values from tree distances. In initial tests, we identified differences between this implementation for calculating ESS values and the values returned by the package tracerer. Tracerer values are concordant with those output by the commonly used program Tracer, and so we modified topological pseudo ESS calculations to use tracerer's implementation of ESS rather than coda's




