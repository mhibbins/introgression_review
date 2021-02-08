import sys 
trees = []

with open(sys.argv[1], 'r') as treesfile:

	for line in treesfile:
		if line.startswith("("):
			trees.append(line)

with open(sys.argv[2], 'w+') as nex_file:

	genetree_string = []

	nex_file.write("#NEXUS\n")
	nex_file.write("\n")
	nex_file.write("BEGIN TREES;\n")
	nex_file.write("\n")
		
	for i in range(len(trees)):
		nex_file.write("Tree geneTree" + str(i) + " = " + trees[i] + "\n")
		genetree_string.append("geneTree" + str(i))

	genetree_string = ','.join(genetree_string)

	nex_file.write("\n")
	nex_file.write("END;\n")
	nex_file.write("\n")
	nex_file.write("\n")
	nex_file.write("BEGIN PHYLONET;\n")
	nex_file.write("\n")
	nex_file.write("InferNetwork_ML (" + genetree_string + ") 1 -bl -di;\n")
	nex_file.write("\n")
	nex_file.write("END;")
	
		

	
