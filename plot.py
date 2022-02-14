import os
import argparse
from inspect import getmembers, isfunction
from pathlib import Path
import shutil

import pandas as pd

from src import plot
from src.io import find_files

if __name__ == "__main__":

    available_plot_functions = [n[0] for n in getmembers(plot, isfunction)]

    parser = argparse.ArgumentParser()
    parser.add_argument("--reports", type=str, nargs='+', required=True, help='Path to a single report file, a list if report files or a dir to be seached for report files')
    parser.add_argument("--output-dir", type=str, help='Path to where the plots are saved', default = None)
    parser.add_argument("--plots", type=str, nargs='+', 
        choices= ["all"] + available_plot_functions, help = 'which plots to produce', default = "all"
    )
    parser.add_argument("--depth", type=int, help='How deep to look for report files', default = 1)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    plots_to_make = args.plots
    if "all" in plots_to_make:
        plots_to_make = available_plot_functions

    all_input_files = list()
    for inp in args.reports:
        if os.path.isfile(inp):
            all_input_files.append(os.path.abspath(inp))
        elif os.path.isdir(inp):
            for l in find_files(inp, starts = ['report'], maxdepth = args.depth):
                all_input_files.append(l)
    all_input_files = set(all_input_files)

    plots_and_files = dict()
    modelnames = list()
    for f in all_input_files:
        fname = Path(f).with_suffix('').stem
        pname = fname.split('_')[2]
        modelnames.append(fname.split('_')[1])
        if not pname in plots_and_files.keys():
            plots_and_files[pname] = list()
        plots_and_files[pname].append(f)
    modelnames = set(modelnames)

    output_dir = args.output_dir
    if os.path.exists(output_dir):
        if args.overwrite:
            shutil.rmtree(output_dir)
            os.makedirs(output_dir)
            os.makedirs(os.path.join(output_dir, 'combined'))
            for modelname in modelnames:
                os.makedirs(os.path.join(output_dir, modelname))
        else:
            raise OSError('Output directory ({0}) exists'.format(output_dir))
    else:
        os.makedirs(output_dir)
        os.makedirs(os.path.join(output_dir, 'combined'))
        for modelname in modelnames:
            os.makedirs(os.path.join(output_dir, modelname))

    for k in plots_to_make:
        df = list()
        for f in plots_and_files[k]:
            df.append(pd.read_csv(f))
        df = pd.concat(df)
        getattr(plot, k)(df = df, output_path = args.output_dir)
