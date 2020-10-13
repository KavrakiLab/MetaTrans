# Metabolite Translator (MetaTrans)

Metabolite Translator (MetaTrans) is a deep learning, Transformer-based, architecture for predicting metabolites of small molecules in humans. 

MetaTrans is developed using transfer learning: First, we trained a Transformer model for predcting the outcome of general chemical reactions. Subsequently, we fine-tuned it on human metabolic reactions. Finally, we constructed an ensemble model, combining multiple fine-tuned models. 
The methodology is described in detain in the paper [Prediction of drug metabolites using neural machine translation](https://pubs.rsc.org/en/content/articlelanding/2020/SC/D0SC02639E#!divAbstract)

The specifications for pre-training the Transformer model were based on the [Molecular Transformer](https://github.com/pschwllr/MolecularTransformer) for reaction outcome prediction.


The implementation of the Transformer model is based on the [OpenNMT toolkit](http://opennmt.net/OpenNMT-py/). 

Processing of the chemical data is based on the [RDKit](https://www.rdkit.org/) software.

## Installation

Create a conda environment:

```bash
conda create -n metatrans python=3.5
source activate metatrans
conda install rdkit -c rdkit
conda install future six tqdm pandas
conda install pytorch=1.1.0 torchvision -c pytorch
pip install torchtext==0.3.1
pip install -e .
```

## Prediction of human metabolites for small molecules

### Get Trained models
Step 1: Download trained models from this [link](https://rice.box.com/s/5jeb5pp0a3jjr3jvkakfmck4gi71opo0) and place thm inside the folder models.

### Prepare data
Step 2: Prepare a csv with the molecules in SMILES notation (Sample input files are given in datasets/test/input.csv). The csv should contain the molecule id in the 1st column and the molecule SMILES in the 2nd column.
Then prepare the data (canonicalise and tokenise SMILES) for translation:

```bash
python prepare_input_file.py -input_file ${infile} -output_file ${outfile}
```
`input_file` the name of the input csv file.
`output_file` (optional) the name of the output txt file which will contain the processed data. Default: processed_data.txt

### Translate
Step 3: Edit the bash file translate_molecules: define the source file (output of the previous step) and the beam size. The default beam size is 5 and the user can change it to get fewer or more predictions per molecule. A beam size of 5 approximatily gives top-10 predictions. A beam size of 2 gives top-5. A beam size of 10 gives top-20. The user can also define the min and max length of the predicted SMILES sequences.
Then, translate the input molecules into metabolites:

```bash
./translate_molecules
```
 

### Get predictions
Step 4: Obtain the predictions:

```bash
python process_predictions.py -input_file ${infile} -output_file ${outfile} -beam_size ${beam} -visualise_predictions ${bool}
```
`input_file` the name of the input csv file (same as the input file in step 2)
`output_file` (optional) the name of the output csv file where the predictions will be saved. Default: predicted_metabolites.csv 
`beam_size` (optional) the beam size which has to be the same as the beam size used in step 3. Default: 5
`visualise_predictions` (optional) if True then the predicted metabolites will be visualised using RDKit. Default: False.

This will generate a csv file with the predicted metabolites. 


## Datasets

The datasets of human metabolic transformations we constructed for training, validating and testing MetaTrans are in the folder datasets. 


## Citation

```
@Article{metatrans,
author ="Litsa, Eleni E. and Das, Payel and Kavraki, Lydia E.",
title  ="Prediction of drug metabolites using neural machine translation",
journal  ="Chem. Sci.",
year  ="2020",
pages  ="-",
publisher  ="The Royal Society of Chemistry",
doi  ="10.1039/D0SC02639E",
url  ="http://dx.doi.org/10.1039/D0SC02639E"
}
```

