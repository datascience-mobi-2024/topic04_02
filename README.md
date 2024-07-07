# Topic 4, Group 2: 
## Authors
Preet Shah, Maximilian Vassen, Marik Müller, Tobias Kreusel
## Supervisors
Prof. Dr. Dominik Niopek\
Dr. Jan Mathony\
Benedict Wolf


Tutor: Maximilian Fidlin

## Abstract
Being able to rationally improve protein thermostability has direct application to research and industry, for example, optimising enzyme catalysis or or improving  organisms growth at higher temperatures. This project, based on an extensive thermal proteome profiling data set, extensively analysed how factors and structural motifs of proteins influence, and potentially improve thermal stability. Essential proteins, limiting an organism's maximum growth temperature, were identified and key biological processes identified. These findings were used in a regression analysis, yielding a model that accurately predicts protein thermostability (r2 = 0.72). Using this regression, a novel in-silico model was developed, capable of predicting beneficial mutations for increased thermostability with minimal structural changes.

## Repository structure
To use any of our code, the notebook **[initialisation.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/initialisation.ipynb)** should be executed to create the required directories and download the dataframes required for our results. It also downloads two example PDB files. A zip file containing all the PDBs we used can be downloaded **[here](https://drive.google.com/file/d/1XFvu7OAfv0gtHU_4MM0vuoPFaZVmM7T2/view?usp=sharing)** (256 MB zipped, ~ 1 GB unzipped) ([Source](https://alphafold.ebi.ac.uk/)).

All of our data cleanup and processing has been collected in **[dataframe_creation.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/dataframe_creation.ipynb)**, this also includes the feature calculation for primary, secondary and tertiary structure. The resulting dataframe 'prokaryotes_348columns.csv' is also provided in the data folder created by running initialisation.\
Our results and plots for general correlations and PCA analysis as well as the training and evaluation of our regression model can be found in **[PCA_and_regression.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/PCA_and_regression.ipynb)**.\
Our analyses of essential proteins (especially in *E. coli* and *B. subtilis*) are collected in **[essential_proteins.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/essential_proteins.ipynb)**.\
Our main functions **[SPARC](https://github.com/datascience-mobi-2024/topic04_02/blob/main/SPARC.py)**, **[ThERMOS and ThERMless](https://github.com/datascience-mobi-2024/topic04_02/blob/main/ThERMOS.py)** are defined in python scripts to be usable from anywhere. They are also all combined into a user friendly class called Protein defined in **[proteinclass.py](https://github.com/datascience-mobi-2024/topic04_02/blob/main/proteinclass.py)**.\
Our functions are applied on the proteins identified as essential in **[Mutation_results_essential_proteins.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/Mutation_results_essential_proteins.ipynb)**.\
The file **[function.py](https://github.com/datascience-mobi-2024/topic04_02/blob/main/function.py)** contains miscellaneous functions but is mostly deprecated.\
Basic usage of the Protein class is demonstrated with two example proteins in **[example_usage.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/example_usage.ipynb)**.

## Requirements
#### Packages
The following packages are required for our code, versions may not have to be the exact same, but numpy should be the given version:

| Package         | Version   |
|-----------------|-----------|
| biopython       | 1.83      |
| joblib          | 1.4.2     |
| matplotlib      | 3.9.0     |
| numpy           | 1.26.4    |
| pandas          | 2.2.2     |
| pdb2pqr         | 3.6.2     |
| scikit-learn    | 1.5.0     |
| scipy           | 1.13.1    |
| seaborn         | 0.13.2    |
| torch           | 2.3.1     |
| statsmodels     | 0.14.2    |

#### s4pred
We used s4pred ([Moffat and Jones, 2021](https://doi.org/10.1093/bioinformatics/btab491)) for secondary structure prediction. It is not listed above as it's not available as a package. Instead it needs to be installed according to the instructions in its [GitHub](https://github.com/psipred/s4pred). The resulting s4pred directory should be placed inside the data directory. Alternatively the s4pred directory can be downloaded **[here](https://drive.google.com/drive/folders/1IRUzcyfX_V62fG6OP2qfKQpIi4CGnjsd?usp=sharing)** (822 MB) to avoid a tedious installation.
