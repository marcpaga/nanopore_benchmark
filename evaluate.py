import os
import shutil
import multiprocessing as mp
import argparse

import mappy
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError

from src.io import read_fast, find_files
from src.evaluation import eval_pair

def results_queue_writer(output_file, q):
    
    while True:
        df = q.get()
        if df == 'kill':
            break
        else:
            df = pd.DataFrame(df, index=[0])
        header = True
        with open(output_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    header = False
                    break

        df.to_csv(output_file, mode='a', header=header, index=False)
            

def eval_pair_wrapper(reads_queue, writer_queue, tmp_dir, homopolymer_min_length, verbose):
    """Wrapper evaluate a prediction in the queue
    Args:
        references (dict): dictionary with reference sequences
        read_queue (multiprocessing.Queue): queue from where to get the predictions
        writer_queue (multiprocessing.Queue): queue where to send the results
        verbose (bool): whether to print the read_id being processed
    """
    
    while not reads_queue.empty():

        data = reads_queue.get()
        read_id, reference, prediction = data

        if verbose:
            print(read_id)

        if isinstance(prediction, tuple):
            pred, comment, phredq = prediction
            if comment == '+' or comment == '-':
                comment = None
        else:
            pred = prediction
            phredq = None
            comment = None
            
        tmp_fasta = os.path.join(tmp_dir, read_id + '.fasta')
        with open(tmp_fasta, 'w') as f:
            f.write('>'+read_id+'\n')
            f.write(reference+'\n')

        ref = mappy.Aligner(tmp_fasta) 

        result = eval_pair(
            ref = ref, 
            que = pred, 
            read_id = read_id, 
            homopolymer_min_length = homopolymer_min_length,
            phredq = phredq, 
            comment = comment,
        )
        
        writer_queue.put(result)

        os.remove(tmp_fasta)

    return None


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--basecalls-path", type=str, required=True, help='Path to a fasta or fastq file or dir to be searched')
    parser.add_argument("--references-path", type=str, required=True, help='Path to a fasta reference file')
    parser.add_argument("--model-name", type=str, required=True, help='Name of the model')
    parser.add_argument("--homopolymer-length", type=int, default = 5, help='Minimum length of same consecutive bases to be considered a homopolymer')
    parser.add_argument("--output-file", type=str, help='csv output file', default = None)
    parser.add_argument("--depth", type=int, help='How deep to look for fastq or fasta files', default = 1)
    parser.add_argument("--processes", type=int, help='Number of parallel processes', default = 2)
    parser.add_argument("--verbose", action="store_true", help='Print read ids as they are being evaluated')
    parser.add_argument("--overwrite", action='store_true', help='Overwrite existing output file')
    args = parser.parse_args()

    # get all the basecall files
    fast_files = list()
    if os.path.isfile(args.basecalls_path):
        fast_files.append(os.path.abspath(args.basecalls_path))
    else:
        for fast_file in find_files(args.basecalls_path, endings = ['.fasta', '.fastq'], maxdepth = args.depth):
            fast_files.append(os.path.abspath(fast_file))
    fast_files = np.unique(fast_files)

    # create a tmp dir to write tmp files
    tmp_path = os.path.join("/".join(fast_files[0].split('/')[:-1]), 'tmp')
    if os.path.exists(tmp_path):
        shutil.rmtree(tmp_path)
    os.makedirs(tmp_path)

    # output file name
    if args.output_file is None:
        output_file = os.path.join("/".join(fast_files[0].split('/')[:-1]), 'evaluation.csv')
    else:
        output_file = args.output_file
    assert output_file.endswith('.csv'), "output file must end with .csv"

    # check if output file exists to skip evaluated reads
    processed_ids = set()
    if os.path.isfile(output_file):
        if args.overwrite:
            os.remove(output_file)
            with open(output_file, 'w') as f:
                f.write('#'+args.model_name+'\n')
        else:
            try:
                df = pd.read_csv(output_file, header = 0, index_col = False, comment = '#')
                processed_ids = set(df['read_id'])
                print('Output file already exists, {0} reads already evaluated'.format(len(processed_ids)))
            except EmptyDataError:
                pass
    else:
        with open(output_file, 'w') as f:
            f.write('#'+args.model_name+'\n')


    # read all the reference sequences into memory
    print('Reading references: ' + args.references_path)
    references = read_fast(os.path.abspath(args.references_path))

    # start multiprocessing manager and queues for reading and writing
    manager = mp.Manager() 
    writer_queue = manager.Queue()
    reads_queue = manager.Queue()
    watcher_writer = mp.Process(target = results_queue_writer, args = (output_file, writer_queue, ))
    watcher_writer.start()

    for basecalls_file in fast_files:
        for read_id, basecalls in read_fast(basecalls_file).items():
            if read_id in processed_ids:
                continue
            try:
                reads_queue.put((read_id, references[read_id], basecalls))
            except KeyError:
                print('Not found in references: ' + read_id)
                continue

  
    with mp.Pool(processes=args.processes-1) as pool:
       
       multiple_results = [pool.apply_async(eval_pair_wrapper, (reads_queue, writer_queue, tmp_path, args.homopolymer_length, args.verbose)) for _ in range(args.processes-1)]
       results = [res.get() for res in multiple_results]
                    
    writer_queue.put('kill')
    watcher_writer.join()

