# A Categorical Naïve Bayes model based on microhaplotypes for ancestry inference
This is a model for predicting superpopulations using microhaplotypes.The traing set is derived from genomic data from 2504 samples of 1000 Genomes Project phase 3.There are 46 available microhaps in this model.
## Minimum version requirements
### Compulsory
Python: 3.6.0
Numpy: 1.8.0
Pandas: 1.1.0
Scikit-learn: 0.24
SciPy: 1.5.1
### Optional
R: 3.6.1
haplo.stats: 1.7.9
## Usage
### 1. Input file
The input file needs to be organized like the example file input.txt we provided.The first column of the file is the sample number and the rest of the columns store the genotype information of microhaps. The first line, in addition to "IDs", records the name of the microhaps corresponding to each column. The genotype of each MH is represented by the genotypes of the SNPs it includes connected by short horizontal lines.Each sample uses two adjacent rows to store data. The two genotypes in the same column of the same sample correspond to the same locus on the homologous chromosome, so their upper and lower orders can be exchanged.
### 2. Run
Run the scripts written in the file model.py in order and the results will be output to the output.txt file.
### 3. Output file
The columns of the output file are, in order, the sample number, the predicted superpopulation, the probabilities in each class, the likelihood ratio (LR), the area under the curve (AUC) in each class, kept numbers of markers for the 46 MHs and 12 MHs. The likelihood ratio is calculated by dividing the probability that this sample is in the most likely superpopulation by the sum of the probabilities that it is in the other superpopulations. AUCs were calculated using the samples in the training set by the leave-one-out cross validation to reflect the effect of the number of available microhaps on the accuracy of the results.
## Phasing from genotypes of SNPs
In addition to sequencing of the microhap locus, it's also feasible to obtain genotype of microhap by phasing from the genotypes of included SNPs. Here we provide a R script Phasing.R to do this. The plink text files (.ped and .map) storing the genotypes of the desired SNPs are required. Fill in the following code with the prefix of two plink files and run it.

`Rscript Phasing.R <prefix>`

The resulting output file <prefix>_phased.txt can be used as an input file for the model above.
