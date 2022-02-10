import os
import shutil
import re
import multiprocessing as mp
import argparse

import mappy
import numpy as np
import pandas as pd

############################ FILE READING/MANAGING #############################

def find_fast_files(top, maxdepth = 1):
    """Find fastq or fasta files recursively
    """
    dirs, nondirs = [], []
    for name in os.listdir(top):
        (dirs if os.path.isdir(os.path.join(top, name)) else nondirs).append(name)
        for nondir in nondirs:
            if nondir.endswith('.fasta') or nondir.endswith('.fastq'):
                yield os.path.join(top, nondir)
    if maxdepth > 1:
        for name in dirs:
            for x in find_fast_files(os.path.join(top, name), maxdepth-1):
                yield x

def iter_fasta(fasta_file):
    """Read a fasta file iteratively
    """
    c = 0
    with open(fasta_file, 'r') as handle:
        for line in handle:
            if c == 0:
                read_id = line[1:].strip('\n')
                c += 1
            elif c == 1:
                seq = line.strip('\n')
                c = 0
                yield read_id, seq

def read_fasta(fasta_file):
    """Read a fasta file
    """
    fasta_dict = dict()
    for k, v in iter_fasta(fasta_file):
        fasta_dict[k] = v
    return fasta_dict

def iter_fastq(fastq_file):

    c = 0
    with open(fastq_file, 'r') as f:
        for line in f:
            if c == 0:
                c += 1
                read_id = line[1:].strip('\n')
            elif c == 1:
                c += 1
                seq = line.strip('\n')
            elif c == 2:
                c += 1
                direction = line.strip('\n')
            elif c == 3:
                c = 0
                phredq = line.strip('\n')
                if len(seq) != len(phredq):
                    raise ValueError('{}: seq ({}) and phredq ({}) lenghts are different'.format(read_id, len(seq), len(phredq)))
                yield read_id, (seq, direction, phredq)

def read_fastq(fastq_file):
    """Read a fastq file
    """
    
    fastq_dict = dict()
    for k, v in iter_fastq(fastq_file):
        fastq_dict[k] = v
                
    return fastq_dict

def read_fna(file):
    """Read a fna file, like fasta but sequences are split by \n 
    """
    d = dict()
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                k = line.strip('\n')
                d[k[1:]] = list()
            else:
                d[k[1:]].append(line.strip('\n'))
                
    for k, v in d.items():
        d[k] = "".join(v)
    return d

def read_fast(fast_file):

    if fast_file.endswith('.fastq'):
        return read_fastq(fast_file)
    if fast_file.endswith('.fasta'):
        return read_fasta(fast_file)
    if fast_file.endswith('.fna'):
        return read_fna(fast_file)

############################ FILE READING/MANAGING #############################

############################ EVALUATE QUEUES ###################################

def results_queue_writer(output_file, q):
    
    while True:
        df = q.get()
        if df == 'kill':
            break
        else:
            df = pd.DataFrame(df, index=[0])
        if not os.path.isfile(output_file):
            df.to_csv(output_file, header=True, index=False)
        else: # else it exists so append without writing the header
            df.to_csv(output_file, mode='a', header=False, index=False)
            

def eval_pair_wrapper(reads_queue, writer_queue, tmp_dir, verbose):
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
            phredq = phredq, 
            comment = comment,
        )
        
        writer_queue.put(result)

        os.remove(tmp_fasta)

    return None
                
############################ EVALUATE QUEUES ###################################

###################### EVALUATION FUNCTIONS AND PARAMS #########################

BASES = ['A', 'C', 'G', 'T']
REPORT_COLUMNS = [
    'read_id', # id of the read
    'len_reference', # length of the reference
    'len_basecalls', # length of the basecalls
    'que_start',
    'que_end',
    'ref_start', # number of insertions/deletions at the start of the alignment
    'ref_end', # number of insertions/deletions at the end of the alignment
    'comment'
] 

ERRORS = list()
for b1 in BASES:
    for b2 in BASES + ['-']:
        for b3 in BASES + ['-']:
            for b4 in BASES:
                if b2 == '-' and b3 == '-':
                    continue
                ERRORS.append(b1 + b2 + '>' + b3 + b4)
REPORT_COLUMNS += ERRORS

for b in BASES:
    REPORT_COLUMNS.append('homo_'+b+'_counts')
    REPORT_COLUMNS.append('homo_'+b+'_errors')

REPORT_COLUMNS.append('phred_mean')
REPORT_COLUMNS.append('phred_median')
REPORT_COLUMNS.append('phred_std')
REPORT_COLUMNS.append('phred_mean_correct')
REPORT_COLUMNS.append('phred_median_correct')
REPORT_COLUMNS.append('phred_std_correct')
REPORT_COLUMNS.append('phred_mean_error')
REPORT_COLUMNS.append('phred_median_error')
REPORT_COLUMNS.append('phred_std_error')

def find_runs(x):
    """Find runs of consecutive items in an array."""

    # ensure array
    x = np.asanyarray(x)
    if x.ndim != 1:
        raise ValueError('only 1D array supported')
    n = x.shape[0]

    # handle empty array
    if n == 0:
        return np.array([]), np.array([]), np.array([])

    else:
        # find run starts
        loc_run_start = np.empty(n, dtype=bool)
        loc_run_start[0] = True
        np.not_equal(x[:-1], x[1:], out=loc_run_start[1:])
        run_starts = np.nonzero(loc_run_start)[0]

        # find run values
        run_values = x[loc_run_start]

        # find run lengths
        run_lengths = np.diff(np.append(run_starts, n))

        return run_values, run_starts, run_lengths

def elongate_cigar(short_cigar):
    """Converts a short cigar with format 9M1I10X16M to a long string of 
    repeated MMMMMMMMMIXXXXXXXXXXMMMMMMMMMMMMMMMM

    Args:
        short_cigar (str): cigar string to be elongated
    """
    cigar_counts = re.split('H|X|=|I|D|N|S|P|M', short_cigar)
    cigar_strs = re.split('[0-9]', short_cigar)
    
    cigar_counts = [c for c in cigar_counts if c != '']
    cigar_strs = [c for c in cigar_strs if c != '']
    
    assert len(cigar_strs) == len(cigar_counts)
    
    longcigar = ''
    for c, s in zip(cigar_counts, cigar_strs):
        longcigar += s*int(c)
    return longcigar, cigar_counts, cigar_strs

def make_align_arr(long_cigar, truth_seq, pred_seq, phredq = None):
    """Makes an alignment array based on the long cigar
    
    Args:
        long_cigar (str): output from `elongate_cigar`
        truth_seq (str): sequence 1
        pred_seq (str): sequence 2
        phredq (str): quality scores for the predicted sequence

    Returns:
        A np:array of shape [3 or 4, alignment_length]. The first dimensions are the
        reference, alignment chars, predicted sequence and phredq if given.
    """

    if phredq is not None:
        if len(pred_seq) != len(phredq):
            raise ValueError('pred_seq ({}) and phredq ({}) lenghts are different'.format(len(pred_seq), len(phredq)))
    
    tc = 0
    pc = 0
    if phredq is None:
        align_arr = np.full((3, len(long_cigar)), '')
    else:
        align_arr = np.full((4, len(long_cigar)), '')
    for i, c in enumerate(long_cigar):
        if c == 'D':
            align_arr[0, i] = truth_seq[tc]
            align_arr[1, i] = ' '
            align_arr[2, i] = '-'
            if phredq is not None:
                align_arr[3, i] = ' '

            tc += 1
        elif c == 'I':
            align_arr[0, i] = '-'
            align_arr[1, i] = ' '
            align_arr[2, i] = pred_seq[pc]
            if phredq is not None:
                align_arr[3, i] = phredq[pc]

            pc += 1
        elif c == 'X':
            align_arr[0, i] = truth_seq[tc]
            align_arr[1, i] = '.'
            align_arr[2, i] = pred_seq[pc]
            if phredq is not None:
                align_arr[3, i] = phredq[pc]

            pc += 1
            tc += 1
        elif c == '=':
            align_arr[0, i] = truth_seq[tc]
            align_arr[1, i] = '|'
            align_arr[2, i] = pred_seq[pc]
            if phredq is not None:
                align_arr[3, i] = phredq[pc]

            pc += 1
            tc += 1
        elif c == 'M':
            align_arr[0, i] = truth_seq[tc]
            align_arr[2, i] = pred_seq[pc]
            if truth_seq[tc] == pred_seq[pc]:
                align_arr[1, i] = '|'
            else:
                align_arr[1, i] = '.'
            if phredq is not None:
                align_arr[3, i] = phredq[pc]

            pc += 1
            tc += 1
            
    return align_arr

def error_profile(align_arr):
    """Counts the different signatures in a local alingment array

    Args:
        arr (np.array): array with the alignment

    Returns:
        A dictionary with signatures as keys and counts as values
    """

    if align_arr[0, 0] == '-' or align_arr[0, -1] == '-':
        raise ValueError('The reference must start and end with bases, not insertions')

    # calculate the mutational signature style errors
    mut_dict = dict()
    for e in ERRORS:
        mut_dict[e] = 0

    # we iterate over the positions for which we can calculate a signature,
    # which are all but the first and last bases
    # for each of these positions we look for the closest base in the predictions
    # on the left side and right side
    # then we get which code should be based on the chunk of the array that 
    # we have
    r = np.array(align_arr[2, :])
    nogaps = r != '-'
    pos = np.arange(0, len(nogaps), 1)

    for i in np.arange(1, len(r) - 1, 1):
        st = pos[:i][np.where(nogaps[:i])[0]][-1]
        nd = pos[i+1:][np.where(nogaps[i+1:])[0]][0]

        code = align_arr[2, st] + align_arr[2, i] + '>' + align_arr[0, i] + align_arr[2, nd]

        mut_dict[code] += 1
    
    return mut_dict

def homopolymer_errors(align_arr):

    result = dict()
    homo_counts = dict()
    homo_errors = dict()
    for b in BASES:
        homo_counts[b] = 0
        homo_errors[b] = 0

    ref_arr = align_arr[0, :]
    for b in BASES:
        base_or_gap = (ref_arr == b) | (ref_arr == '-')
        sections = find_runs(base_or_gap)
        for t, st, l in zip(*sections):
            if not t:
                continue
            if l < 5:
                continue
            if np.sum(align_arr[0, st:st+l] == b) < 5:
                continue
            h_arr = align_arr[:, st:st+l]
            for j in range(h_arr.shape[1]):
                if h_arr[0, j] == '-' and h_arr[2, j] == b:
                    homo_errors[b] += 1
                elif h_arr[0, j] == b:
                    if h_arr[2, j] == b:
                        homo_counts[b] += 1
                    else:
                        homo_counts[b] += 1
                        homo_errors[b] += 1

    for b in BASES:
        result['homo_'+b+'_counts'] = homo_counts[b]
        result['homo_'+b+'_errors'] = homo_errors[b]

    return result

def eval_phredq_scores(align_arr):

    correct = list()
    error = list()
    result = dict()

    for i in range(align_arr.shape[1]):

        phred_symbol = align_arr[3, i]
        align_symbol = align_arr[1, i]
        
        if phred_symbol == ' ':
            continue
        
        score = ord(phred_symbol) - 33
        if align_symbol == '|':
            correct.append(score)
        elif align_symbol == '.':
            error.append(score)
        elif align_symbol == ' ':
            error.append(score)

    result['phred_mean'] = np.mean(correct + error)
    result['phred_median'] = np.median(correct + error)
    result['phred_std'] = np.std(correct + error)
    result['phred_mean_correct'] = np.mean(correct)
    result['phred_median_correct'] = np.median(correct)
    result['phred_std_correct'] = np.std(correct)
    result['phred_mean_error'] = np.mean(error)
    result['phred_median_error'] = np.median(error)
    result['phred_std_error'] = np.std(error)

    return result


def eval_pair(ref, que, read_id, phredq = None, comment = None):
    """Align two sequences and evaluate the alignment
    
    Args:
        ref (str): reference sequence or aligner if using minimap2
        que (str): predicted sequence
        read_id (str): uuid of the read
        phredq (str): string with predq symbols
        comment (str): comment in the direction line of a fastq file
        
    Returns:
        results (dict): dictionary with metrics and long confusion matrix
    """

    result = dict()
    for k in REPORT_COLUMNS:
        result[k] = None
    result['read_id'] = read_id

    # if there are no basecalls return empty results
    len_que = len(que)
    if len_que == 0:
        result['len_basecalls'] = len_que
        result['comment'] = 'no prediction'
        return result

    # if there is no correct alignment
    correct_match = False
    for alignment in ref.map(seq = que):
        if read_id == alignment.ctg:
            correct_match = True
            break
    if not correct_match:
        result['comment'] = 'failed mapping'
        return result
    

    # prepare an array with the aligned sequences to be evaluated
    que_st = alignment.q_st
    que_nd = alignment.q_en + 1
    ref_st = alignment.r_st
    ref_nd = alignment.r_en + 1
    len_ref = len(ref.seq(read_id))
    
    long_cigar, _, _ = elongate_cigar(alignment.cigar_str)

    if phredq is not None:
        phredq = phredq[que_st:que_nd]

    local_arr = make_align_arr(
        long_cigar = long_cigar, 
        truth_seq = ref.seq(read_id)[ref_st:ref_nd], 
        pred_seq = que[que_st:que_nd], 
        phredq = phredq
    )

    # fix the cigar as it does not have mismatches as X
    longcigar_arr = np.array(list(long_cigar))
    longcigar_arr[np.where(local_arr[1, :] == '.')[0]] = 'X'
    longcigar_fixed = "".join(longcigar_arr.tolist())

    # minimap2 does not produce local alignments so we have to cut the overhangs
    local = np.where((np.array(list(longcigar_fixed)) == 'X') | (np.array(list(longcigar_fixed)) == 'M'))[0]
    local_st = local[0]
    local_nd = len(longcigar_fixed) - local[-1]
    longcigar_fixed = longcigar_fixed[local_st:local_nd+1]
    que_st += local_st
    que_nd -= local_nd
    ref_st += local_st
    ref_nd -= local_nd
    local_arr = local_arr[:, local_st:local[-1]+1]
    
    # report results on length and alignment
    result['len_reference'] = len_ref
    result['len_basecalls'] = len_que
    result['ref_start'] = ref_st
    result['ref_end'] = ref_nd
    result['que_start'] = que_st
    result['que_end'] = que_nd
    
    # calculate error profiles
    signatures = error_profile(local_arr)
    result = {**result, **signatures}
    
    # count for each base the amount of bases in homopolymers
    # and how many errors in these regions
    homopolymer_results = homopolymer_errors(local_arr)
    result = {**result, **homopolymer_results}

    # calculate mean phredq scores for correct and incorrect bases
    if phredq is not None:
        phredq_results = eval_phredq_scores(local_arr)
        result = {**result, **phredq_results}

    if comment is None:
        result['comment'] = 'pass'
    else:
        result['comment'] = comment

    return result

###################### EVALUATION FUNCTIONS AND PARAMS #########################

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--basecalls-path", type=str, required=True, help='Path to a fasta or fastq file or dir to be searched')
    parser.add_argument("--references-path", type=str, required=True, help='Path to a fasta reference file')
    parser.add_argument("--model-name", type=str, required=True, help='Name of the model being evaluated')
    parser.add_argument("--output-file", type=str, help='csv output file', default = None)
    parser.add_argument("--depth", type=int, help='How deep to look for fastq or fasta files', default = 1)
    parser.add_argument("--processes", type=int, help='Number of parallel processes', default = 1)
    parser.add_argument("--verbose", action="store_true", help='Print read ids as they are being evaluated')
    args = parser.parse_args()

    # get all the basecall files
    fast_files = list()
    if os.path.isfile(args.basecalls_path):
        fast_files.append(os.path.abspath(args.basecalls_path))
    else:
        for fast_file in find_fast_files(args.basecalls_path, args.depth):
            fast_files.append(os.path.abspath(fast_file))
    fast_files = np.unique(fast_files)

    # create a tmp dir to write tmp files
    tmp_path = os.path.join("/".join(fast_files[0].split('/')[:-1]), 'tmp')
    if os.path.exists(tmp_path):
        shutil.rmtree(tmp_path)
    os.makedirs(tmp_path)

    # output file name
    if args.output_file is None:
        output_file = os.path.join("/".join(fast_files[0].split('/')[:-1]), args.model_name+'_evaluation.csv')
    else:
        output_file = args.output_file
    assert output_file.endswith('.csv'), "output file must end with .csv"


    # check if output file exists to skip evaluated reads
    processed_ids = set()
    if os.path.isfile(output_file):
        df = pd.read_csv(output_file, header = 0, index_col = False)
        processed_ids = set(df['read_id'])
        print('Output file already exists, {0} reads already evaluated'.format(len(processed_ids)))

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
            reads_queue.put((read_id, references[read_id], basecalls))

  
    with mp.Pool(processes=args.processes-1) as pool:
       
       multiple_results = [pool.apply_async(eval_pair_wrapper, (reads_queue, writer_queue, tmp_path, args.verbose)) for _ in range(args.processes-1)]
       results = [res.get() for res in multiple_results]
                    
    writer_queue.put('kill')
    watcher_writer.join()

