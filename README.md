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