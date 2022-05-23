#!/bin/bash

output_dir=$1

# files with links to the data to be downloaded
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
data_links="${SCRIPT_DIR}/links_jain_data.txt"
references_links="${SCRIPT_DIR}/links_jain_references.txt"

mkdir -p $output_dir

# make a tmp dir for all the downloads so we can delete the dir afterwards
tmp_dir="${output_dir}/tmp"
mkdir -p $tmp_dir

# download the reference and unpack it
genomes_dir="${output_dir}/genomes"
mkdir -p $genomes_dir

while IFS= read -r line || [[ -n $line  ]]; do
    if [ ! -f "${tmp_dir}/reference.fasta.gz" ]
    then
        echo "Downloading reference data"
        wget $line -O "${tmp_dir}/reference.fasta.gz"
    fi

    if [ ! -f "${genomes_dir}/Homo_sapiens.fna" ]
    then
        echo "Decompressing reference data"
        gunzip -c "${tmp_dir}/reference.fasta.gz" > "${genomes_dir}/Homo_sapiens.fna"
    fi 
done < $references_links


# download the data
while IFS= read -r line || [[ -n $line  ]]; do
    files=($line)
    filename=${files[0]##*/}
    runid=${filename%-*}
    
    if [ ! -f "${tmp_dir}/${runid}_fastq.tar" ]
    then
        echo "Downloading basecalls data"
        wget ${files[1]} -O "${tmp_dir}/${runid}_fastq.tar"
    fi

    if [ ! -f "${tmp_dir}/${runid}_fast5.tar" ]
    then
        echo "Downloading fast5 data"
        wget ${files[0]} -O "${tmp_dir}/${runid}_fast5.tar"
    fi

    run_dir="${output_dir}/Homo_sapiens-${runid}"
    fastq_dir="${run_dir}/fastq"
    fast5_dir="${run_dir}/fast5"
    mkdir -p $run_dir $fastq_dir $fast5_dir

    if [ -z "$(ls -A ${fastq_dir})" ]
    then
        echo "Decompressing fastq data"
        tar -xf "${tmp_dir}/${runid}_fastq.tar" -C "${run_dir}/fastq"
    fi

    if [ -z "$(ls -A ${fast5_dir})" ]
    then
        echo "Decompressing fast5 data"
        tar -xf "${tmp_dir}/${runid}_fast5.tar" -C "${run_dir}/fast5"
    fi

done < $data_links

