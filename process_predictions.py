import argparse
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from utils import *
import os

## it reads the output of the models, un-tokenises the predicted sequences and filters out unlikely metabolites
## -input_file: the csv file that has the input molecules (molecule ID and SMILES representations)
## -output_file: the filename where the processed predictions will be saved. It's a csv file. 
## predictions_directory: the directory where the output of the models from the tranaslate_molecules script is saved
## -beam_size: the beam_size. It can be in [5,10,15,20]
## -visualise_molecules (boolean): it visualises all predicted metabolites if True. They are stored within the predictions directory.



def main(opt):
	input_file = opt.input_file
	output_file = opt.output_file
	predictions_directory = opt.predictions_directory
	figures_directory = predictions_directory + 'Figures/'
	models = ['591e8', '8743b', '84545', 'ebe91','c6631']
	beam = opt.beam_size
	if not beam in [5,10,15,20]:
		print('Beam size can be 5, 10, 15 or 20. Not ', beam)
		exit()

	pred_lines_5 = {}
	pred_lines_10 = {}
	pred_lines_15 = {}
	pred_lines_20 = {}

	for num in range(0,len(models)):
		predictions_file = predictions_directory+models[num]+'_'+'beam'+str(5)+'.txt'
		with open(predictions_file) as f_pred:  
			pred_lines_5[num] = [''.join(line.strip().split(' ')) for line in f_pred.readlines()]

	for num in range(0,len(models)):
		predictions_file = predictions_directory+models[num]+'_'+'beam'+str(10)+'.txt'
		with open(predictions_file) as f_pred: 
			pred_lines_10[num] = [''.join(line.strip().split(' ')) for line in f_pred.readlines()]

	for num in range(0,len(models)):
		predictions_file = predictions_directory+models[num]+'_'+'beam'+str(15)+'.txt'
		with open(predictions_file) as f_pred:  
			pred_lines_15[num] = [''.join(line.strip().split(' ')) for line in f_pred.readlines()]

	for num in range(0,len(models)):
		predictions_file = predictions_directory+models[num]+'_'+'beam'+str(20)+'.txt'
		with open(predictions_file) as f_pred:   
			pred_lines_20[num] = [''.join(line.strip().split(' ')) for line in f_pred.readlines()]

	models_count = len(pred_lines_5.keys())

	if opt.visualise_molecules:
		if not os.path.exists(figures_directory):
			os.makedirs(figures_directory)

	molID2smiles = {}
	molID2metabolites = {}
	index_5 = 0
	index_10 = 0
	index_15 = 0
	index_20 = 0
	drug_lines = open(input_file).read().split('\n')
	pred_counts = []
	for i in range(1,len(drug_lines)-1):
		mol_id,smiles = drug_lines[i].split(',')
		smiles = canonicalise_smile(smiles)
		molID2smiles[mol_id] = smiles
		predictions = set()
		for j in range(index_5,index_5+5):
			for num in range(0,models_count):
				predictions.add(pred_lines_5[num][j])
		index_5 = index_5 + 5
		if beam>5:
			for j in range(index_10,index_10+10):
				for num in range(0,models_count):
					predictions.add(pred_lines_10[num][j])
			index_10 = index_10 + 10
		if beam>10:
			for j in range(index_15,index_15+15):
				for num in range(0,models_count):
					predictions.add(pred_lines_15[num][j])
			index_15 = index_15 + 15
		if beam>15:
			for j in range(index_20,index_20+20):
				for num in range(0,models_count):
					predictions.add(pred_lines_20[num][j])
			index_20 = index_20 + 20
		processed, invalid, invalid_count = process_predictions(predictions,smiles,0.2,0.2,False,True)
		pred_counts.append(len(processed))
		molID2metabolites[mol_id] = processed
		drug = Chem.MolFromSmiles(smiles)
		preds = [Chem.MolFromSmiles(pred_smiles) for pred_smiles in processed]
		fig_dir = figures_directory + '/' + mol_id + '/'
		if not os.path.exists(fig_dir):
			os.makedirs(fig_dir)
		filename = fig_dir + mol_id + '.png'
		img = Draw.MolToFile(drug,filename,size=(500,500),wedgeBonds=False)
		prd_count = 1
		for prd in preds:
			filename = fig_dir + 'Metabolite' + str(prd_count) + '.png'
			img = Draw.MolToFile(prd,filename,size=(500,500),wedgeBonds=False)
			prd_count = prd_count  + 1

	table = ['Molecule ID', 'SMILES', 'Metabolites']
	for mol_id in molID2metabolites.keys():
		metabolites_str = ''
		smiles = molID2smiles[mol_id]
		metabolites = molID2metabolites[mol_id]
		for metabolite in metabolites:
			metabolites_str = metabolites_str + metabolite + ' '
		metabolites_str = metabolites_str[:-1]
		entry = [mol_id, smiles, metabolites_str]
		table = np.vstack((table,entry))

	with open(output_file,'wb') as f:
		np.savetxt(f,table, fmt='%s', delimiter=',')


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-input_file', type=str, default='data/test/test_molecules_input_file.csv',help='Input File')
	parser.add_argument('-output_file', type=str, default='data/test/predicted_metabolites_on_test_set.csv',help='Processed Predictions File')
	parser.add_argument('-predictions_directory', type=str, default='predictions/test/extended/predictions_extended_min5_max120/',help='Predictions Directory')
	parser.add_argument('-beam_size', type=int, default=15,help='Beam Size')
	parser.add_argument('-visualise_molecules', type=bool, default=False,help='Visualise predicted metabolites')
	opt = parser.parse_args()
	main(opt)