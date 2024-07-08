# Topic 4, Group 2: Decoding Protein Stability - a Data-Driven Approach for Predictive Modeling of Enhanced Thermostability

## Authors
Preet Shah (preet.shah@stud.uni-heidelberg.de) \
Maximilian Vassen (maximilian.vassen@stud.uni-heidelberg.de) \
Marik Müller (marik.mueller@stud.uni-heidelberg.de) \
Tobias Kreusel (tobias.kreusel@stud.uni-heidelberg.de)
## Supervisors
Prof. Dr. Dominik Niopek (dominik.niopek@uni-heidelberg.de) \
Dr. Jan Mathony (jan.mathony@uni-heidelberg.de)\
Benedict Wolf (b.wolf@stud.uni-heidelberg.de) \
Tutor: Maximilian Fidlin (maximilian.fidlin@stud.uni-heidelberg.de)

## Abstract
Extending protein thermostability has significant applications in research and industry, such as extending shelf life, optimising enzyme catalysis, or improving growth of organisms at higher temperatures. This project utilised an extensive thermal proteome profiling dataset of mesophile and thermophile organisms to comprehensively analyse how various factors and structural motifs of proteins influence thermostability and to identify parameters for improving thermal stability. These findings were then used in a regression analysis, yielding a model that accurately predicts protein thermostability from protein sequence (r2 =0.72). Using this model, essential proteins that presumably act as a bottleneck, restricting the maximum growth temperature of an organism, were identified and classified according to their function. In a second step, a mutational screen was applied to these proteins, predicting beneficial mutations for increased thermostability with minimal structural changes. To facilitate automated use, a user friendly python class was developed combining **SPARC** (**S**equence-based **P**rotein **A**ttribute-de**R**ived (melting point) **C**alculator) and **ThERMOS** (**Th**ermal **E**nhancement by **R**apid **M**utational **O**ptimization **S**creen) as the main features to predict mutations enhancing thermostability while maintaining structural integrity.  \
\
The full report can be found **[here](https://github.com/datascience-mobi-2024/topic04_02/blob/main/Decoding_Protein_Stability_-%20_a_Data-Driven_Approach_for_Predictive_Modeling_of_Enhanced_Thermostability.pdf)**

## Repository structure
To use any of our code, a specific directory structure is required. The directories are automatically created upon importing the module [Protein]() for the first time. This also downloads our important dataframes and two example PDB files. A zip file containing all the PDBs we used can be downloaded **[here].
To use any of our code, the notebook **[initialisation.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/initialisation.ipynb)** should be executed to create the required directories and download the dataframes required for our results. It also downloads two example PDB files. A zip file containing all the PDBs we used can be downloaded **[here](https://drive.google.com/file/d/1XFvu7OAfv0gtHU_4MM0vuoPFaZVmM7T2/view?usp=sharing)** (256 MB zipped, ~ 1 GB unzipped) ([Source](https://alphafold.ebi.ac.uk/)).
All of our data cleanup and processing has been collected in **[dataframe_creation.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/dataframe_creation.ipynb)**, this also includes the feature calculation for primary, secondary and tertiary structure. The resulting dataframe 'prokaryotes_348columns.csv' is also provided in the data folder created by running initialisation.\
Our results and plots for general correlations and PCA analysis as well as the training and evaluation of our regression model can be found in **[PCA_and_regression.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/PCA_and_regression.ipynb)**.\
Our analyses of essential proteins (especially in *E. coli* and *B. subtilis*) are collected in **[essential_proteins.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/essential_proteins.ipynb)**.\
Our main functions **[SPARC](https://github.com/datascience-mobi-2024/topic04_02/blob/main/SPARC.py)**, **[ThERMOS and ThERMless](https://github.com/datascience-mobi-2024/topic04_02/blob/main/ThERMOS.py)** are defined in python scripts to be usable from anywhere. They are also all combined into a user friendly class called Protein defined in **[proteinclass.py](https://github.com/datascience-mobi-2024/topic04_02/blob/main/proteinclass.py)**.\
Our functions are applied on the proteins identified as essential in **[Mutation_results_essential_proteins.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/Mutation_results_essential_proteins.ipynb)**.\
The file **[function.py](https://github.com/datascience-mobi-2024/topic04_02/blob/main/function.py)** contains miscellaneous functions but is mostly deprecated.\
Basic usage of the Protein class is demonstrated with two example proteins in **[demonstration.ipynb](https://github.com/datascience-mobi-2024/topic04_02/blob/main/demonstration.ipynb)**.
## Requirements
#### Packages
All required packages can be activated by creating a virtual environment using **[requirements.txt](https://github.com/datascience-mobi-2024/topic04_02/blob/main/requirements.txt)**.\
The packages required for our code can be also installed manually. We did not encouter any issues with newer package versions, except for the numpy package, where not all versions work with S4pred:
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
