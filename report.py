import argparse
import os
import shutil

import numpy as np
import pandas as pd

from src.report import EvaluationReport
from src.io import find_files


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluate-path", type=str, required=True, help='Path to evaluation csv file or dir where to search them')
    parser.add_argument("--model-name", type=str, required = True, help='Name of the model being evaluated')
    parser.add_argument("--depth", type=int, help='How deep to look for fastq or fasta files', default = 1)
    parser.add_argument("--verbose", action="store_true", help='Print read ids as they are being evaluated')
    parser.add_argument("--overwrite", action='store_true', help='Overwrite existing files')
    args = parser.parse_args()

    eval_files = list()
    if os.path.isfile(args.evaluate_path):
        eval_files.append(os.path.abspath(args.evaluate_path))
    else:
        for eval_file in find_files(args.evaluate_path, ['_evaluation.csv'], args.depth):
            eval_files.append(os.path.abspath(eval_file))
    eval_files = np.unique(eval_files)

    for eval_file in eval_files:

        report_path = os.path.join("/".join(eval_file.split('/')[:-1]), 'reports')
        if os.path.exists(report_path):
            if args.overwrite:
                shutil.rmtree(report_path)
                os.makedirs(report_path)
            else:
                raise OSError('Output directory ({0}) exists'.format(report_path))
        else:
            os.makedirs(report_path)

        df = pd.read_csv(eval_file)

        ereport = EvaluationReport(df, args.model_name, output_path = report_path, overwrite = args.overwrite)

        ereport.count_abs_counts()
        ereport.read_outcome_counts()
        ereport.event_rates()
        ereport.homopolymer_rates()
        ereport.phredq_distributions()
        ereport.calculate_signatures()
        ereport.calculate_auc()