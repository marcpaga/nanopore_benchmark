import argparse
import os
import shutil

import pandas as pd

from src.report import EvaluationReport


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluation-file", type=str, required=True, help='Evaluation file to be analyzed')
    parser.add_argument("--output-dir", type=str, help='Dir where to save all the report tables', default = None)
    parser.add_argument("--model-name", type=str, help='Model name evaluated', default = None)
    parser.add_argument("--depth", type=int, help='How deep to look for fastq or fasta files', default = 1)
    parser.add_argument("--verbose", action="store_true", help='Print read ids as they are being evaluated')
    parser.add_argument("--overwrite", action='store_true', help='Overwrite existing files')
    args = parser.parse_args()

    if args.output_dir is None:
        output_dir = os.path.abspath(os.path.join("/".join(args.evaluation_file.split('/')[:-1]), 'reports'))
    else:
        output_dir = os.path.abspath(args.output_dir)

    if os.path.exists(output_dir):
        if args.overwrite:
            shutil.rmtree(output_dir)
            os.makedirs(output_dir)
        else:
            raise OSError('Output directory ({0}) exists'.format(output_dir))
    else:
        os.makedirs(output_dir)
    print('Output dir: ' + str(output_dir))

    df = pd.read_csv(args.evaluation_file, comment = '#')
    if args.model_name is None:
        with open(args.evaluation_file, 'r') as f:
            for line in f:
                modelname = line[1:].strip('\n')
                break
    else:
        modelname = args.model_name

    ereport = EvaluationReport(df, modelname = modelname, output_path = output_dir, overwrite = args.overwrite)

    ereport.count_abs_counts()
    ereport.read_outcome_counts()
    ereport.event_rates()
    ereport.homopolymer_rates()
    ereport.phredq_distributions()
    ereport.calculate_signatures()
    ereport.calculate_auc()