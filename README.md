# Nanopore Benchmark for Basecallers

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

```
source venv3/bin/activate

python3 report.py \
--evaluation-file demo/model1/evaluation.csv \
--output-dir demo/model1/reports

```

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

To download all the data do `bash download_all.sh OUTPUTDIR`, otherwise use one of the other three scripts to download specific datasets: `download_wick.sh` (bacterial data), `download_verm.sh` (lambda phage), `download_jain.sh` (human data). NOTE: the lambda phage data from `download_verm.sh` has not been uploaded yet.

WARNING: This is a large data download, about 1TB in disk space.

## Todo

- Upload lambda phage data 