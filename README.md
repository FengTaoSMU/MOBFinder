# MOBFinder: a tool for MOB typing for plasmid metagenomic fragments based on language model

## Introduction
Mobilization (MOB) typing is a classification scheme to classify plasmids based on their mobility or transferability, and it can help us understand the mechanism by which plasmids are transferred between bacterial cells. MOBFinder is developed for MOB typing for plasmid fragments and bins in metagenomic data

Based on the natural language processing technique, MOBFinder uses the word vector language model to characterize the plasmid fragments. Then MOBFinder predicts the MOB type of the plasmid metagenomic fragments using several random forest classification models, which can be manually downloaded from [here](https://www.jianguoyun.com/p/DWMTw1oQ5cbbCxjfk44FIAA).

## Citations
Tao F, Shufang W, Zhencheng Fang and Hongwei Z. "MOBFinder: a tool for MOB typing for plasmid metagenomic fragments based on language model."


## Preparation ##

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


## Installation

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
% git clone https://github.com/phac-nml/mob-suite.git
%
```
## Usage

## Using MOBFinder to perform MOB typing for plasmid metagenomic fragments

```
# input: Plasmid fragments fasta file with unique id
# output: test_mobfinder_result.tsv
% python3 MOBFinder.py -i ./testdata/test.fa -o ./
```

# Output files
| File | Description |
| ------------ | ------------ |
| test_mobfinder_result.tsv | This file contians the predicted MOB class of each plasmid fragments and scores of each MOB classes |

# Output file format
| field  | Description |
| --------- |  --------- | 
| id | Plasmid fragment id |
| class | Predicted classes of each plasmid fragments |
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


## Using MOBFinder to perform MOB typing for plasmid metagenomic bins

```
# input: a plasmid fragments fasta file with unique id, a meta file records the mapping relationship of plasmid fragment id and bin id
# output: test_bin_mobfinder_binning_result.tsv
% python3 MOBFinder_bin.py -i ./testdata/test_bin.fa -b ./testdata/test.tsv -o .
```

# Output file format
| field  | Description |
| --------- |  --------- | 
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
Tao Feng - fengtaosmu@foxmail.com

## License

MOBFinder is distributed under a GPL-3.0 license.
