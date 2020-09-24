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
Prepare a txt file with the molecules in SMILES notation (A sample is given in datasets/test/input.txt). Then prepare (tokenise) the data for translation:

```bash
python prepare_input_file.py ${infile}
```
`infile` (optional) the name of the input txt file. Default: input.txt
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
