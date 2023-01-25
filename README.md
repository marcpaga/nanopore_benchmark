# Nanopore Benchmark for Basecallers

This repository contains download links to several Nanopore sequencing (R9.4.1) datasets that can be used for benchmarking basecallers. It also contains code to evaluate the basecalls with the reference sequences to benchmark the accuracy and performance of the basecallers. 

**Who is this for?**
- You work on the development of nanopore sequencing basecallers
- You want to compare basecallers on read-level accuracy
- You want to use several evaluation metrics (per base accuracy, error profiles, homopolymer error rates, etc.)
- You want a defined training and testing datasets

These datasets and evaluation metrics have been used in our benchmark of the most recent basecalling neural network architectures: https://www.biorxiv.org/content/10.1101/2022.05.17.492272v2. If you want to compare your model against these models; you can use the same datasets (train/test splits) and the same metrics and forget about having to re-implement/train the already evaluated models.

For information regarind the architectures evaluated please see: https://github.com/marcpaga/basecalling_architectures

## Installation

```
$ git clone https://github.com/marcpaga/nanopore_benchmark.git 
$ cd nanopore_benchmark
$ python3 -m venv venv3
$ source venv3/bin/activate
(venv3) $ pip install --upgrade pip
(venv3) $ pip install -r requirements.txt
```

## Usage

### Evaluation

```
source venv3/bin/activate

python3 evaluate.py \
--basecalls-path demo/model1 \
--references-path demo/reference.fasta \
--model-name model1

# output -> data/model1/model1_evaluation.csv
```

DO NOT USE '_' (underscores) in model names please.

`--basecalls-path` can be a directory with `.fastq` or `.fasta` files or a single file. 

If `.fasta` files are provided, then there will not be an analysis on the PhredQ scores.

If `.fastq` files are provided, every 3rd line can be used to add a comment to that particular read in case it passes the evaluation, if it does not, then the comment is appended to the type of failure that the evaluation gives for that particular read. If there is no comment to provide just put a `+` or `-` sign in that line.

If `.fasta` files are provided, but one wants to use the comment functionality of `.fastq` files just create `.fastq` with fake phredq scores (e.g. fill with `!`).


### Report

To create a report of the evaluation of your model based on the reference and basecalls do the following:

```
source venv3/bin/activate

python3 report.py \
--evaluation-file demo/model1/evaluation.csv \
--output-dir demo/model1/reports

```

This will generate a bunch of report csv files in the output dir:
- `absoultecounts`: this contains the counts across all reads for different metrics (matched bases, mismatches bases, homopolymer correct bases, etc.)
- `auc`: this contains the values necessary to plot the AUC for the model. 
    - `fraction`: top fraction of best reads according to phredq.
    - `match_rate`: match rate of reads in that fraction.
    - `phredq_mean`: average PhredQ score of reads in that fraction.
- `event rates`: this contains the boxplot statistics for the main alignment events: match, mismatch, insertion and deletion.
- `homopolymerrates`: this containts the boxplot statistics for the homopolymer error rates per base or all together.
- `phredq`: this contains the boxplot statistics for the PhredQ scores of correctly and incorrectly basecalled bases.
- `readoutcomes`: this contains the number of reads that are successfully evaluated or that had some sort of error.
- `signatures`: this contains the rates and counts of different types of errors for each base in a 3-mer context. The 3-mer contexts are based on the basecalls, not the reference.
- `singlevalues`: this contains single summary values across all metrics based on the absolute counts, the read outcomes and the PhredQ scores distributions.

### Plot

```

source venv3/bin/activate

python3 plot.py \
--reports demo \
--output-dir demo/plots \
--depth 3

```

## Data download

To download the data check the bash scripts in `./download`.

To download all the data do `bash download_all.sh OUTPUTDIR`, otherwise use one of the other three scripts to download specific datasets: `download_wick.sh` (bacterial data), `download_verm.sh` (lambda phage), `download_jain.sh` (human data). 

WARNING: This is a large data download, about 1TB in disk space.

## Dataset tasks

Benchmarking tasks are divided between `global`, `human` and `cross-species`. Each task has its own set of data that can be used for training and testing. Lists with the reads that can be used for train and testing can be found in `static/tasks`.

The idea of dataset tasks are to evaluate the models in difference scenarios:

- The `global` task evaluates the performance of a general-purpose model by training and testing the model using data from all available species. In this task, models have access to most data. 
- The `human` task  evaluates the performance of a human specialized model by training and testing exclusively on human data. 
- The `cross-species` task evaluates the performance of a trained model on never before seen species. This allows us to evaluate the robustness of the model to overfit on genomic features from the training set, such as k-mer distribution or base modifications. This is achieved by training using a subset of bacterial species and testing on data from all species.
