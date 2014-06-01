import sys

def get_dictionary(in_reads):
	#erstellt dictionary
	#key = readname, value = bitflag
	d={}
	identify = ""
	for line in in_reads:
		line = line.strip()
		if line[0] == '>':			
			identify = str(line[1:])
		else:
			d[identify] = line
	return d



if __name__ == "__main__":
	clusterFile = sys.argv[1]
	originFastaFile = sys.argv[2]	
	fastaFile_folder = sys.argv[3]
	outputName = sys.argv[4]
	
	cluster = open(clusterFile, "r")
	origin_fasta = open(originFastaFile, "r")
	#fq_file = open(fastaFile, "w")	
	
	dictionary = get_dictionary(origin_fasta)
	print(dictionary)
	
	number = 0
	sequence = ""
	for line in cluster:
		line = line.strip()
		clusterA = str.split(line, '\t') 
		StringName = str(fastaFile_folder) + str(outputName) + str(number) + ".fasta"
		fq_file = open(StringName, "w")
		for i in clusterA:								
			if str(i) in dictionary:
				header = ">" + str(i) + "\n"				
				sequence = dictionary[str(i)] + "\n"			
				fq_file.write(header)
				fq_file.write(sequence)
		number += 1
		 
