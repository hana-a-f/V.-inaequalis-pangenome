##this script will filter the column predicted helix for values greater than zero, this is indicative of a transmembrane domain##
import sys

#input file is the tmhmm output file
input_file = sys.argv[1]
#output file will be a filtered version
output_file = sys.argv[2]

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
	for line in infile:
		columns = line.strip().split() #splits each line in the input file into columns
		pred_hel_value = int(columns[4].split("=")[1]) #extarct the value after predHel
		if pred_hel_value > 0:
			outfile.write(line)

