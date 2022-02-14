import os

def find_files(top, starts = None, endings = None, maxdepth = 1):
    """Find files recursively

    Args:
        top (str): dir to start searching
        formats (str): string or list of strings with file formats (e.g '.csv')
        maxdepth (int): how deep to search for files
    """
    if not isinstance(endings, list):
        endings = [endings]
    if not isinstance(starts, list):
        starts = [starts]

    dirs, nondirs = [], []
    for name in os.listdir(top):
        (dirs if os.path.isdir(os.path.join(top, name)) else nondirs).append(name)
        for nondir in nondirs:
            for ending in endings:
                if ending is None:
                    continue
                if nondir.endswith(ending):
                    yield os.path.join(top, nondir)
            for start in starts:
                if start is None:
                    continue
                if nondir.startswith(start):
                    yield os.path.join(top, nondir)
    if maxdepth > 1:
        for name in dirs:
            for x in find_files(os.path.join(top, name), starts, endings, maxdepth-1):
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