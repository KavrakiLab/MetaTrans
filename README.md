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
Download trained models from this [link](https://rice.box.com/s/5jeb5pp0a3jjr3jvkakfmck4gi71opo0) and place thm inside the folder models.

### Prepare data
Prepare a txt file with the molecules in SMILES notation (A sample is given: input.txt). Then prepare (tokenise) the data for translation:

```bash
python process_input.py ${infile}
```
`infile` (optional) the name of the input txt file. Default: input.txt
### Translate

```bash
./metatrans 
```
### Get predictions

```bash
python process_predictions.py ${outfile}
```
`outfile` (optional) the name of the output csv file. Default: predicted_metabolites.csv 
This will generate a csv file with the predicted metabolites.


## Datasets

The datasets of human metabolic transformations we constructed for training, validating and testing MetaTrans are in the folder datasets. 
