###imported modules
import re
import sys
import os

###Functions 

#function for listing files in any path
def ls(path):
	files = os.listdir(path)
	#remove .DS Store files in the list of files in the folder
	if (files[0] == ".DS_Store"):
		files = files[1:] 
	return files

#function for changing directory	
def cd( cmd ):
	os.chdir(cmd)
	
def ls2(path):
	files = os.listdir(path)
	return files

#function for getting AICc information out of jmodeltest output files and writing a bayes block to the corresponding nexus file for input into mr bayes. 	
def getmrbayesblock(file):
	base=os.path.splitext(os.path.basename(file))[0]
	modeltest_output = open(base + "_modeltest.txt")
	
	#inside jmodel output file, search for the AICc best fit model results
	for line in modeltest_output:
		if line.startswith("AICc "):
			#return copy of list under AICc model selection block with the trailing space codes "/n" removed
			model_line = line.rstrip("\n")
			#split up the information given under AICc by whitespace
			elements = re.split('\s*', model_line)
			#grab only the model selected information from that list of elements (the first model mentioned)
			best_model = elements[1]
			#match that model to two sets of parameters
			parameters_first = {"JC":"lset nst=1 rates=equal", "JC+I":"lset nst=1 rates=propinv", "JC+G":"lset nst=1 rates=gamma", "JC+I+G":"lset nst=1 rates=invgamma", "F81":"lset nst=1 rates=equal", "F81+I":"lset nst=1 rates=propinv", "F81+G":"lset nst=1 rates=gamma", "F81+I+G":"lset nst=1 rates=invgamma", "K80":"lset nst=2 rates=equal", "K80+I":"lset nst=2 rates=propinv", "K80+G":"lset nst=2 rates=gamma", "K80+I+G":"lset nst=2 rates=invgamma", "HKY":"lset nst=2 rates=equal", "HKY+I":"lset nst=2 rates=propinv", "HKY+G":"lset nst=2 rates=gamma", "HKY+I+G":"lset nst=2 rates=invgamma","SYM":"lset nst=6 rates=equal", "SYM+I":"lset nst=6 rates=propinv", "SYM+G":"lset nst=6 rates=gamma", "SYM+I+G":"lset nst=6 rates=invgamma", "GTR":"lset nst=6 rates=equal", "GTR+I":"lset nst=6 rates=propinv", "GTR+G":"lset nst=6 rates=gamma", "GTR+I+G":"lset nst=6 rates=invgamma"}
			parameters_second = {"JC":"prset pinvarpr=fixed(0) statefreqpr=fixed(equal)", "JC+I":"prset statefreqpr=fixed(equal)", "JC+G":"prset pinvarpr=fixed(0) statefreqpr=fixed(equal)", "JC+I+G":"prset statefreqpr=fixed(equal)", "F81":"prset pinvarpr=fixed(0)", "F81+I":"NA", "F81+G":"prset pinvarpr=fixed(0)", "F81+I+G":"NA", "K80":"prset pinvarpr=fixed(0) statefreqpr=fixed(equal)", "K80+I":"prset statefreqpr=fixed(equal)", "K80+G":"prset pinvarpr=fixed(0) statefreqpr=fixed(equal)", "K80+I+G":"prset statefreqpr=fixed(equal)", "HKY":"prset pinvarpr=fixed(0)", "HKY+I":"NA", "HKY+G":"prset pinvarpr=fixed(0)", "HKY+I+G":"NA","SYM":"prset pinvarpr=fixed(0) statefreqpr=fixed(equal)", "SYM+I":"prset statefreqpr=fixed(equal)", "SYM+G":"prset pinvarpr=fixed(0) statefreqpr=fixed(equal)", "SYM+I+G":"prset statefreqpr=fixed(equal)", "GTR":"prset pinvarpr=fixed(0)", "GTR+I":"NA", "GTR+G":"prset pinvarpr=fixed(0)", "GTR+I+G":"NA"}
			lset_params = parameters_first[best_model]
			prset_params = parameters_second[best_model]
			#cd("..") #?
			mbb_file=open(file, "a")
			if prset_params == "NA":
				mbb_file.write("\n" + "begin mrbayes;" + "\n" + lset_params + ";" + "\n" + "log start filename= "+ base + ".log" + ";" + "\n" + "mcmc nruns=4 nchains=4 ngen=10000000 diagnfreq=1000 printfreq=1000 samplefreq=10000 file=" + file + ";" + "\n" + "sump burnin=250;" + "\n" + "sumt burnin=250;" + "\n" + "log stop;" + "\n" + "end;")
				mbb_file.close()
			else:
				mbb_file.write("\n" + "begin mrbayes;" + "\n" + lset_params + ";" + "\n" + prset_params + ";" + "\n" + "log start filename= "+ base + ".log" + ";" + "\n" + "mcmc nruns=4 nchains=4 ngen=10000000 diagnfreq=1000 printfreq=1000 samplefreq=10000 file=" + file + ";" + "\n" + "sump burnin=250;" + "\n" + "sumt burnin=250;" + "\n" + "log stop;" + "\n" + "end;")
				mbb_file.close()
	modeltest_output.close()
	
			
##loop
cd("../04_NEXUS_CONVERSION")
files = ls(".")
print files
for file in files:
	if file.endswith(".nex"):
		getmrbayesblock(file)	
cd("..")

