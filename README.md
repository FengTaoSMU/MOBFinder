# MOBFinder: a tool for mobilization typing for plasmid metagenomic fragments based on a language model

## Change logs
2025.12.09 The README.md was updated with a few modification.

## Introduction
Mobilization typing (MOB) is a classification scheme to classify plasmids based on their mobility or transferability, and it can help us understand the mechanism by which plasmids are transferred between bacterial cells. MOBFinder is developed for MOB typing for plasmid fragments and bins in metagenomic data.

Based on the natural language processing technique, MOBFinder uses the word vector language model to characterize the plasmid fragments from different MOB types. Based on the plasmid fragments represented by the word vector, several random forest classification models were trained and integrated for predicting plasmid fragments with different lengths，which can be manually downloaded from [here](https://zenodo.org/records/13368995).

MOBFinder takes one or multiple plasmid DNA sequences as input and outputs predicted MOB types for each fragment, including MOBB, MOBC, MOBF, MOBH, MOBL, MOBM, MOBP, MOBQ, MOBT, MOBV and non-mobilizable (non-mob). Plasmid bins in metagenomics can also be annotated using MOBFinder.

## Citations
Feng T, Wu S, Zhou H, Fang Z. MOBFinder: a tool for mobilization typing of plasmid metagenomic fragments based on a language model. Gigascience. 2024 Jan 2;13:giae047.

## Requires
+ Python >= 3.9.11
+ ete3 >= 3.1.2
+ biopython >= 1.78
+ pycurl >= 7.45.1
+ numpy >= 1.24.3
+ scipy >= 1.10.1
+ six >= 1.16.0

## Dependencies
+ blast >= 2.12.0
+ r-base >= 4.2.0
+ r-randomforest >= 4.7_1.1

## Preparation
MOBFinder is developed using some dependencies, and we recommend using conda to install them.
```
% conda config --add channels defaults
% conda config --add channels conda-forge
% conda config --add channels bioconda

% conda create -n mobfinder
% conda activate mobfinder
% conda install python
% conda install biopython
% conda conda install -c etetoolkit ete3
% conda install pycurl
% conda install r-base
% conda install r-randomforest
% conda install -c bioconda blast
```

## Installation
Clone this repository to your local linux PC.
```
% git clone https://github.com/FengTaoSMU/MOBFinder.git
% cd MOBFinder
```

Since the model trained by MOBFinder is large and exceeds GitHub's file upload limits, we provide an additional download link for the model.

Users need to download the model and place it in the /model folder. The download links are: https://zenodo.org/records/13368995.

## MOB typing for plasmid metagenomic fragments
For plasmid metagenomic fragments, MOBFinder takes a FASTA file as input. The output file consists of 13 columns. The first column represents the fragment ID, the second column displays the predicted MOB type, and columns three to thirteen represent the scores for each MOB class, namely MOBB, MOBC, MOBF, MOBH, MOBL, MOBM, MOBP, MOBQ, MOBT, MOBV, and non-mob. Here, the input fasta file is [test.fa](testdata/test.fa) .
```
% python3 MOBFinder.py -i ./testdata/test.fa -o ./
```
### Output files
| File | Description |
| :------------: | :------------: |
| test_mobfinder_result.tsv | This file contians the predicted MOB class of each plasmid fragment and scores for each MOB classes |
### Output file format
| field | Description |
| :---------: | :---------: | 
| id | Plasmid fragments' id |
| class | Predicted classes of each plasmid fragment |
| mobb_score | Score of MOBB |
| mobc_score | Score of MOBC |
| mobf_score | Score of MOBF |
| mobh_score | Score of MOBH |
| mobl_score | Score of MOBL |
| mobm_score | Score of MOBM |
| mobp_score | Score of MOBP |
| mobq_score | Score of MOBQ |
| mobt_score | Score of MOBT |
| mobv_score | Score of MOBV |
| non-mob_score | Score of Non-mob |

## MOB typing for plasmid metagenomic bins
For plasmid metagenomic bins, MOBFinder requires two input files: a FASTA file containing the plasmid fragments and a meta table that records the mapping between plasmid fragment IDs and bin IDs. The output results are similar to the output of plasmid fragments. The first column is plasmid bin’s id. The second column is the predicted MOB class of plasmid bins. The other columns were MOB scores of different MOB categories produced by MOBFinder. Here, the input fasta file is [test_bin.fa](testdata/test_bin.fa), and the input meta table is [test.tsv](testdata/test.tsv). 
```
% python3 MOBFinder_bin.py -i ./testdata/test_bin.fa -b ./testdata/test.tsv -o .
```
### Output files
| File | Description |
| :------------: | :------------: |
| test_bin_mobfinder_binning_result.tsv | This file contians the predicted MOB class of each plasmid bin and scores of each MOB classes |
### Output file format
| field | Description |
| :---------: | :---------: | 
| id | Plasmid bins' id |
| class | Predicted classes of each plasmid bin |
| mobb_score | Score of MOBB |
| mobc_score | Score of MOBC |
| mobf_score | Score of MOBF |
| mobh_score | Score of MOBH |
| mobl_score | Score of MOBL |
| mobm_score | Score of MOBM |
| mobp_score | Score of MOBP |
| mobq_score | Score of MOBQ |
| mobt_score | Score of MOBT |
| mobv_score | Score of MOBV |
| non-mob_score | Score of Non-mob |


## Contact
Feng Tao - fengtaosmu@foxmail.com

## License
MOBFinder is distributed under a GPL-3.0 license.
