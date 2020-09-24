import argparse
from utils import *

##  it prepares the molecules for translation: it canonicalizes SMILES and then tokenizes them. The data are stored in a txt file
## -input_file: the input molecules in a csv file (1st colums indicates the molecule ID / name and 2nd column indicates the SMILES representation)
## -output_file: the filename where the tokenised data will be saved 

def main(opt):
	input_file = opt.input_file
	output_file = opt.output_file
	count_invalid = 0
	outfile = open(output_file,'w')
	lines = open(input_file).read().split('\n')
	for i in range(1,len(lines)-1):
		_,smiles = lines[i].split(',')
		if not check_smile(smiles):
			print('invalid SMILES: ', smiles)
			count_invalid = count_invalid + 1
		smiles = canonicalise_smile(smiles)
		smiles_tok = smi_tokenizer(smiles)
		if i<len(lines)-2:
			outfile.write(smiles_tok + '\n')
		else:
			outfile.write(smiles_tok)
	outfile.close() 
	if count_invalid>0:
		print(count_invalid, 'invalid SMILES removed.')



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-input_file', type=str, default='data/test/test_molecules_input_file.csv',
                       help='Input File')
	parser.add_argument('-output_file', type=str, default='data/test/test_molecules_source.txt',
                       help='Output File')
	opt = parser.parse_args()
	main(opt)