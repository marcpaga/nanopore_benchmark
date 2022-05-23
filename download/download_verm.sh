#!/bin/bash

output_dir=$1

# files with links to the data to be downloaded
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
data_links="${SCRIPT_DIR}/links_verm_data.txt"
references_links="${SCRIPT_DIR}/links_verm_references.txt"

mkdir -p $output_dir

# make a tmp dir for all the downloads so we can delete the dir afterwards
tmp_dir="${output_dir}/tmp"
mkdir -p $tmp_dir

# download the reference and unpack it
genomes_dir="${output_dir}/genomes"
mkdir -p $genomes_dir

while IFS= read -r line || [[ -n $line  ]]; do
    if [ ! -f "${genomes_dir}/Lambda_phage.fna" ]
    then
        echo "Downloading reference data"
        curl -L $line > ${genomes_dir}/Lambda_phage.fna
    fi
done < $references_links


# download the data
while IFS= read -r line || [[ -n $line  ]]; do
    files=($line)
    link=${files[0]}
    spe=${files[1]}
    
    if [ ! -f "${tmp_dir}/${spe}.tar" ]
    then
        echo "Downloading data"
        wget $link -O "${tmp_dir}/${spe}.tar"
    fi

    run_dir="${output_dir}/Lambda_phage"
    fastq_dir="${run_dir}/fastq"
    fast5_dir="${run_dir}/fast5"
    tmptmp_dir="${run_dir}/tmp"
    mkdir -p $run_dir $fastq_dir $fast5_dir $tmp_dir

    if [ -z "$(ls -A ${tmptmp_dir})" ]
    then
        echo "Decompressing data"
        tar -xf "${tmp_dir}/${spe}.tar" -C $tmptmp_dir
    fi

done < $data_links