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
--basecalls-path data/model1 \
--references-path data/reference.fasta \
--model-name model1

# output -> data/model1/model1_evaluation.csv
```

### Report

```
source venv3/bin/activate

python3 report.py \
--basecalls-path data/model1 \
--references-path data/reference.fasta \
--model-name model1

# output -> data/model1/model1_evaluation.csv
```


## Development

