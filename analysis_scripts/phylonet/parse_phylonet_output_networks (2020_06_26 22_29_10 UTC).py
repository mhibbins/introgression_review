import sys

with open(sys.argv[1], 'r+') as output_file:
	for line in output_file:
		if "Inferred Network" in line:
			print(next(output_file))	
