import pandas as pd
import numpy as np
from numpy import nan
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit import DataStructs
import math


def canonicalise_smile(smi):
    mol = Chem.MolFromSmiles(smi)
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    return canonical

def randomise_smile(smi):
    # unrestricted
    mol = Chem.MolFromSmiles(smi)
    random = Chem.MolToSmiles(mol, canonical=False, doRandom=True, isomericSmiles=False)
    return random


def smi_tokenizer(smi):
    """
    Tokenize a SMILES molecule or reaction
    """
    import re
    pattern =  "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    return ' '.join(tokens)


def check_smile(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
    	return False
    else:
    	return True

def count_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumAtoms()

def get_added_atoms(reactant,product):
    #returns the set of symbols of the atoms that have been added in the product and 
    # a vocabulary which for each symbol gives the number of added atoms
    reactant_atoms = {}
    product_atoms = {}
    for atom in reactant.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in reactant_atoms.keys():
            counter = reactant_atoms[symbol]
        else:
            counter = 0
        counter = counter + 1
        reactant_atoms[symbol] = counter
    for atom in product.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in product_atoms.keys():
            counter = product_atoms[symbol]
        else:
            counter = 0
        counter = counter + 1
        product_atoms[symbol] = counter
    added_atoms_counts = {}
    added_atoms_set = set()
    for symbol in product_atoms.keys():
        counter_prod = product_atoms[symbol]
        if not symbol in reactant_atoms.keys():
            counter_reac = 0
            added_atoms_set.add(symbol)
        else:
            counter_reac = reactant_atoms[symbol]
        if counter_prod>counter_reac:
            diff = counter_prod-counter_reac
            added_atoms_counts[symbol] = diff
    return added_atoms_set, added_atoms_counts


def check_added_atoms(reactant,product):
    pass_test = True
    accepted_atoms = set(['C','c','O','o','H'])
    added_atoms, counts = get_added_atoms(reactant,product)
    for atom in added_atoms:
        if not atom in accepted_atoms:
            if (atom == 'S') and ('O' in counts.keys()) and (counts['O']==3):
                pass_test = True
            else:
                pass_test = False
    return pass_test


def aldehyde_To_Carboxyl(smiles):
    try:
        if smiles.startswith('O=C'):
            smiles = smiles[3:len(smiles)]
            smiles = 'OC(=O)'+smiles
        if smiles.endswith('C=O'):
            smiles = smiles[0:len(smiles)-3]
            smiles = smiles+'C(=O)O'
    except:
        print('Molecule Destroyed!')
    return smiles


def process_targets(targets):
    processed = set()
    for trg in targets:
        # trg = trg.replace('[C@H]', 'C')
        # trg = trg.replace('[C@@H]', 'C')
        # trg = trg.replace('[O-]','O')
        # trg = trg.replace('[C@]','C')
        # trg = trg.replace('[C@@]','C')
        # trg = trg.replace('[P@]','P')
        # trg = trg.replace('[P@@]','P')
        # trg = trg.replace('/','')
        # trg = trg.replace('\\','')
        trg = canonicalise_smile(trg)
        processed.add(trg)
    return processed

def process_target(target):
    # trg = target.replace('[C@H]', 'C')
    # trg = trg.replace('[C@@H]', 'C')
    # trg = trg.replace('[O-]','O')
    # trg = trg.replace('[C@]','C')
    # trg = trg.replace('[C@@]','C')
    # trg = trg.replace('[P@]','P')
    # trg = trg.replace('[P@@]','P')
    # trg = trg.replace('/','')
    # trg = trg.replace('\\','')
    trg = canonicalise_smile(trg)
    return trg


def process_predictions(predictions,drug,size_diff_thresh,similarity_threshold,aldehyde,filtering):
    invalid = False
    invalid_count = 0
    canonicalised = set()
    for pred in predictions:
        if '.' in pred:
            continue
        try:
            canonical = canonicalise_smile(pred)
            canonicalised.add(canonical)
        except:
            invalid_count = invalid_count + 1
    processed = set()
    if invalid_count == len(predictions):
        invalid = True
    else:
        drug = Chem.MolFromSmiles(drug)
        drug_fgp = Chem.RDKFingerprint(drug)
        initial_atomcounts = drug.GetNumAtoms()
        for pred in canonicalised:
            if aldehyde:
                if pred.startswith('O=C') and (not (pred[3] in set(['(','1','2','3']))) or pred.endswith('C=O'):
                    pred = aldehyde_To_Carboxyl(pred)
            # pred = pred.replace('[C@H]', 'C')
            # pred = pred.replace('[C@@H]', 'C')
            # pred = pred.replace('[O-]','O')
            # pred = pred.replace('[C@]','C')
            # pred = pred.replace('[C@@]','C')
            # pred = pred.replace('[P@]','P')
            # pred = pred.replace('[P@@]','P')
            # pred = pred.replace('/','')
            # pred = pred.replace('\\','')
            mol = Chem.MolFromSmiles(pred)
            if mol is None:
                continue
            try:
                pred_fgp = Chem.RDKFingerprint(mol)
                sim = DataStructs.FingerprintSimilarity(drug_fgp,pred_fgp)
            except:
                sim = 0
                print('Exception in computing fingerprint similarity!')
                print(drug_smiles)
            pred = canonicalise_smile(pred)
            if filtering:
                if check_added_atoms(drug,mol):
                    if mol.GetNumAtoms()>math.ceil(size_diff_thresh*initial_atomcounts):
                        if sim>similarity_threshold and sim<1:
                            processed.add(pred)
            else:
                if mol.GetNumAtoms()>math.ceil(size_diff_thresh*initial_atomcounts):
                    if sim>similarity_threshold and sim<1:
                        processed.add(pred)
    return processed, invalid, invalid_count


def get_similarity(prediction_smiles, target_smiles,drug):
    #the targets_set includes only the metabolites that have not been identified
    similarities = []
    closest_preds = set()
    for target in target_smiles:
        max_similarity = -5
        target_mol = Chem.MolFromSmiles(target)
        target_fgp = Chem.RDKFingerprint(target_mol)
        closest = ''
        for prediction in prediction_smiles:
            prediction_mol = Chem.MolFromSmiles(prediction)
            prediction_fgp = Chem.RDKFingerprint(prediction_mol)
            sim = DataStructs.FingerprintSimilarity(target_fgp,prediction_fgp)
            if sim > max_similarity:
                max_similarity = sim
                closest = prediction
        similarities.append(max_similarity)
        closest_preds.add((drug,target,closest))
    return similarities, closest_preds



def countAtoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol is None:
        return mol.GetNumAtoms()
    else:
    	return 0