<html>
<section>
    <h1>Topic 4, Group 2:</h1>
</section>
<section>
    <h2>Authors </h2>
        <p>Preet Shah, Maximilian Vassen, Marik Müller, Tobias Kreusel</p>
    <h2>Supervisors</h2>
        <p>
        Prof. Dr. Dominik Niopek <br>
        Dr. Jan Mathony <br>
        Benedict Wolf
        <br><br>
        Tutor: Maximilian Fidlin
        </p>
</section>




<h2>Abstract</h2>
<p align="justify"> Extending protein thermostability has significant applications in research and industry, such as extending shelf life, optimising enzyme catalysis, or improving growth of organisms at higher temperatures. This project utilised an extensive thermal proteome profiling dataset of mesophile and thermophile organisms to comprehensively analyse how various factors and structural motifs of proteins influence thermostability and to identify parameters for improving thermal stability. These findings were then used in a regression analysis, yielding a model that accurately predicts protein thermostability from protein sequence (r2 =0.72). Using this model, essential proteins that presumably act as a bottleneck, restricting the maximum growth temperature of an organism, were identified and classified according to their function. In a second step, a mutational screen was applied to these proteins, predicting beneficial mutations for increased thermostability with minimal structural changes. To facilitate automated use, a user friendly python class was developed combining <b>SPARC</b> (<b>S</b>equence-based <b>P</b>rotein <b>A</b>ttribute-de<b>R</b>ived (melting point) <b>C</b>alculator) and <b>ThERMOS</b> (<b>Th</b>ermal <b>E</b>nhancement by <b>R</b>apid <b>M</b>utational <b>O</b>ptimization <b>S</b>creen) as the main features to predict mutations enhancing thermostability while maintaining structural integrity. </p>
 

<h2>Repository structure</h2>
<p align="justify"> To use any of our code, the notebook <a href = https://github.com/datascience-mobi-2024/topic04_02/blob/main/initialisation.ipynb> initialisation.ipynb</a> should be executed to create the required directories and download the dataframes required for our results. It also downloads two example PDB files. A zip file containing all the PDBs we used can be downloaded <a href=https://drive.google.com/file/d/1XFvu7OAfv0gtHU_4MM0vuoPFaZVmM7T2/view?usp=sharing> here </a> (256 MB zipped, ~ 1 GB unzipped) (<a href=https://alphafold.ebi.ac.uk/> Source </a>). </p>

<p align="justify"> All of our data cleanup and processing has been collected in <a href=https://github.com/datascience-mobi-2024/topic04_02/blob/main/dataframe_creation.ipynb> dataframe_creation.ipynb </a>, this also includes the feature calculation for primary, secondary and tertiary structure. The resulting dataframe 'prokaryotes_348columns.csv' is also provided in the data folder created by running initialisation.<br>
Our results and plots for general correlations and PCA analysis as well as the training and evaluation of our regression model can be found in <a href=https://github.com/datascience-mobi-2024/topic04_02/blob/main/PCA_and_regression.ipynb> PCA_and_regression.ipynb</a> <br>
Our analyses of essential proteins (especially in <i>E. coli</i> and <i>B. subtilis</i>) are collected in **[essential_proteins.ipynb](<a href= https://github.com/datascience-mobi-2024/topic04_02/blob/main/essential_proteins.ipynb> essential_proteins.ipynb </a>.<br>
Our main functions <a href= https://github.com/datascience-mobi-2024/topic04_02/blob/main/SPARC.py>SPARC</a>, <a href=https://github.com/datascience-mobi-2024/topic04_02/blob/main/ThERMOS.py>ThERMOS and ThERMless</a> are defined in python scripts to be usable from anywhere. They are also all combined into a user friendly class called Protein defined in <a href= https://github.com/datascience-mobi-2024/topic04_02/blob/main/proteinclass.py>proteinclass.py</a>.<br>
Our functions are applied on the proteins identified as essential in <a href= https://github.com/datascience-mobi-2024/topic04_02/blob/main/Mutation_results_essential_proteins.ipynb>Mutation_results_essential_proteins.ipynb</a>.<br>
The file <a href=https://github.com/datascience-mobi-2024/topic04_02/blob/main/function.py>function.py</a> contains miscellaneous functions but is mostly deprecated.<br>
Basic usage of the Protein class is demonstrated with two example proteins in <a href=https://github.com/datascience-mobi-2024/topic04_02/blob/main/example_usage.ipynb>example_usage.ipynb</a>. </p>

<h2>Requirements</h2>
<h3>Packages</h3>
<p>The following packages are required for our code, versions may not have to be the exact same, but numpy should be the given version:</p>
<table>
 <tr>  <th>Package</th>     <th>Version</th> </tr>
 <tr>  <td>biopython</td>   <td>1.83</td>    </tr>
 <tr>  <td>joblib</td>      <td>1.4.2</td>   </tr>
 <tr>  <td>numpy</td>       <td>1.26.4</td>  </tr>
 <tr>  <td>pandas</td>      <td>2.2.2</td>   </tr>
 <tr>  <td>pdb2pqr</td>     <td>3.6.2</td>   </tr>
 <tr>  <td>skipy</td>       <td>1.13.1</td>  </tr>
 <tr>  <td>torch</td>       <td>2.3.1</td>   </tr>
 <tr>  <td>statsmodels</td> <td>0.14.2</td>  </tr>
 
</table>

<h3>S4pred</h3> 
<p align="justify"> We used S4pred (<a href=https://doi.org/10.1093/bioinformatics/btab491>Moffat and jones, 2021</a>) for secondary structure prediction. It is not listed above as it's not available as a package. Instead it needs to be installed according to the instructions in its <a href=https://github.com/psipred/s4pred>GitHub</a>. The resulting s4pred directory should be placed inside the data directory. Alternatively the s4pred directory can be downloaded (<a href=https://drive.google.com/drive/folders/1IRUzcyfX_V62fG6OP2qfKQpIi4CGnjsd?usp=sharing>here</a> (822 MB) to avoid a tedious installation.</p>
</html>
