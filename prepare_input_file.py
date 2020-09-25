import argparse
from utils import *

##  it prepares the molecules for translation: it canonicalizes SMILES and then tokenizes them. The data are stored in a txt file
## -input_file: the input molecules in a csv ot txt file
## -output_file: the filename where the tokenised data will be saved 
## -col: the column in the csv file that contains the SMILES representations of the input molecules

def main(opt):
	input_file = opt.input_file
	output_file = opt.output_file
	col = opt.col
	count_invalid = 0
	outfile = open(output_file,'w')
	lines = open(input_file).read().split('\n')
	for i in range(0,len(lines)):
		data = lines[i].split(',')
		if len(data)==1:
			smiles = data[0]
		else:
			smiles = data[col]
		if not check_smile(smiles):
			print('invalid SMILES: ', smiles)
			count_invalid = count_invalid + 1
			continue
		smiles = canonicalise_smile(smiles)
		smiles_tok = smi_tokenizer(smiles)
		if i<len(lines)-1:
			outfile.write(smiles_tok + '\n')
		else:
			outfile.write(smiles_tok)
	outfile.close() 
	if count_invalid>0:
		print(count_invalid, 'invalid SMILES removed.')



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-input_file', type=str,
                       help='Input File')
	parser.add_argument('-output_file', type=str, default='datasets/test/test_molecules_source.txt',
                       help='Output File')
	parser.add_argument('-col', type=int, default=1,
                       help='CSV column with SMILES')
	opt = parser.parse_args()
	main(opt)