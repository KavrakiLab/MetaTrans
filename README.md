# Metabolite Translator (MetaTrans)

Metabolite Translator [MetaTrans] is a deep learning, Transformer-based, architecture for predicting metabolites of small molecules in humans. 

MetaTrans is developed using transfer learning: First, we trained a Transformer model for predcting the outcome of general chemical reactions. Subsequently, we fine-tuned it on human metabolic reactions. Finally we constructed an ensemble model, 

The specifications for pre-training the Transformer model were based on the [Molecular Transformer for reaction outcome prediction] (https://github.com/pschwllr/MolecularTransformer)


The implementation of the Transformer model is based on the [OpenNMT toolkit] (http://opennmt.net/OpenNMT-py/). 

## Installation

Create a conda environment:

```bash
conda create -n metatrans python=3.5
source activate metatrans
conda install rdkit -c rdkit
conda install future six tqdm pandas
conda install pytorch=0.4.1 torchvision -c pytorch
pip install torchtext==0.3.1
pip install -e .
```

## Prediction of human metabolites for small molecules

### Get Trained models
Step 1: Download trained models from this [link](https://rice.box.com/s/5jeb5pp0a3jjr3jvkakfmck4gi71opo0) and place thm inside the folder models.

### Prepare data
Step 2: Prepare a file (csv or txt) with the molecules in SMILES notation (Sample input files are given in datasets/test/input.csv and input.txt). 
Recommended use: store the files in a csv file where the 1st colum indicating the molecule ID/name and the second colum containing the SMILES representation.
Then prepare the data (canonicalise and tokenise SMILES) for translation:

```bash
python prepare_input_file.py -input_file ${infile} -output_file ${outfile} -col ${col}
```
`infile` the name of the input, csv or txt, file.
`outfile` (optional) the name of the output txt file which will contain the processed data. Default: processed_data.txt
`col` (optional) if the input file is in csv format, the user can specify the colum that contains the molecules SMILES. Default: 1 (2nd column)

### Translate

```bash
./translate_molecules
```
The default beam size is 5. The user can change the beam size (change the variable BEAM in the translate_molecules file) in order to get fewer or more predictions per molecule. A beam size of 5 approximatily gives top-10 predictions. A beam size of 2 gives top-5. A beam size of 10 gives top-20. 

### Get predictions

```bash
python process_predictions.py ${outfile}
```
`outfile` (optional) the name of the output csv file. Default: predicted_metabolites.csv 
This will generate a csv file with the predicted metabolites.


## Datasets

The datasets of human metabolic transformations we constructed for training, validating and testing MetaTrans are in the folder datasets. 
