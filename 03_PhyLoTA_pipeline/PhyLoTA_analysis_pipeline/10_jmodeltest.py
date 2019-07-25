###imported modules
import os
import subprocess

###functions

#function for listing files in any path
def ls(path):
	files = os.listdir(path)
	#remove .DS Store files in the list of files in the folder
	if (files[0] == ".DS_Store"):
		files = files[1:] 
	return files

# function for changing directories
def cd( cmd ):
	os.chdir(cmd)
	
# function for making directories
def mkdir( cmd ):
	os.mkdir(cmd)

# function for calling shell commands within python script	
def shCall( cmd ):
	subprocess.call(cmd,shell=True)
	
#function for running jmodeltest on a nexus file with AICc to find best fit model and put the output in a separate folder
def jmodeltesting(file):
	
	outfile= file + "_modeltest.txt"
	jmodel_output = outfile.replace(".nex_modeltest.txt", "_modeltest.txt")
	os.system("java -jar /usr/local/bin/jmodeltest2/jModelTest.jar -d " + file + " -g 4 -i -f -t BIONJ -s 3 -AICc -o " + jmodel_output )
	
###For performing jmodeltest in a multiple folder loop using a list of all nexus files in the Clade directory
cd("../04_NEXUS_CONVERSION")
files = ls(".")
print files
for file in files:
	if file.endswith(".nex"):
		jmodeltesting(file)
cd("..")