#!/bin/bash

output_dir=$1

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

train_fast5_links="${SCRIPT_DIR}/links_wick_data_train.txt"
test_fast5_links="${SCRIPT_DIR}/links_wick_data_test.txt"

basecalls_links="${SCRIPT_DIR}/links_wick_basecalls_test.txt"
references_links_general="${SCRIPT_DIR}/links_wick_references_general.txt"
references_links_specific="${SCRIPT_DIR}/links_wick_references_specific.txt"

mkdir -p $output_dir

# make a tmp dir for all the downloads so we can delete the dir afterwards
tmp_dir="${output_dir}/tmp"
mkdir -p $tmp_dir

# download the reference and unpack it
genomes_dir="${output_dir}/genomes"
mkdir -p $genomes_dir

echo "Downloading general genomes"
# download general genomes
while IFS= read -r line || [[ -n $line  ]]; do

    files=($line)
    link=${files[0]}
    spe=${files[1]}

    if [ ! -f "${tmp_dir}/${spe}.fna.gz" ]
    then
        echo "Downloading reference data for ${spe}"
        wget $link -O "${tmp_dir}/${spe}.fna.gz"
    fi

    if [ ! -f "${genomes_dir}/${spe}.fna" ]
    then
        echo "Decompressing reference data for ${spe}"
        gunzip -c "${tmp_dir}/${spe}.fna.gz" > "${genomes_dir}/${spe}.fna"
    fi 
done < $references_links_general

echo "Downloading wick specific genomes"
# download wick specific genomes
while IFS= read -r line || [[ -n $line  ]]; do

    files=($line)
    link=${files[0]}
    spe=${files[1]}

    if [ ! -f "${tmp_dir}/${spe}.fna.gz" ]
    then
        echo "Downloading reference data for ${spe}"
        wget $link -O "${tmp_dir}/${spe}.fna.gz"
    fi

    if [ ! -f "${genomes_dir}/${spe}.fna" ]
    then
        echo "Decompressing reference data for ${spe}"
        gunzip -c "${tmp_dir}/${spe}.fna.gz" > "${genomes_dir}/${spe}.fna"
    fi 
done < $references_links_specific

echo "Downloading train fast5 data"
# download train fast5 data
while IFS= read -r line || [[ -n $line  ]]; do

    files=($line)
    link=${files[0]}
    spe=${files[1]}

    echo "Checking status for ${spe}"

    if [ ! -f "${tmp_dir}/${spe}.tar.gz" ]
    then
        echo "Downloading fast5 data for ${spe}"
        wget $link -O "${tmp_dir}/${spe}.tar.gz"
    fi

    run_dir="${output_dir}/${spe}"
    mkdir -p ${run_dir} "${run_dir}/fast5" "${run_dir}/tmp"

    if [ ! -f "${run_dir}/read_references.fasta" ]
    then
        echo "Extracting fast5 and reference data for ${spe}"
        if tar -xzf "${tmp_dir}/${spe}.tar.gz" -C "${run_dir}/tmp";
        then
            echo "Moving fast5 data for ${spe}"
            find "${run_dir}/tmp" -type f -name "*.fast5" -exec mv {} "${run_dir}/fast5" \;
            find "${run_dir}/tmp" -type f -name "read_references.fasta" -exec mv {} "${run_dir}/read_references.fasta" \;
        else
            echo "Error extracting data for ${spe}"
            rm -r "${run_dir}/tmp"
        fi
    fi

done < $train_fast5_links

echo "Downloading test fast5 data"
# download test fast5 data
while IFS= read -r line || [[ -n $line  ]]; do

    files=($line)
    link=${files[0]}
    spe=${files[1]}

    echo "Checking status for ${spe}"

    if [ ! -f "${tmp_dir}/${spe}.tar.gz" ]
    then
        echo "Downloading fast5 data for ${spe}"
        wget $link -O "${tmp_dir}/${spe}.tar.gz"
    fi

    echo "Extracting data for ${spe}"
    run_dir="${output_dir}/${spe}"
    mkdir -p ${run_dir} "${run_dir}/fast5"
    tar -xzf "${tmp_dir}/${spe}.tar.gz" -C "${run_dir}/fast5"
    

done < $test_fast5_links

echo "Downloading test fastq data"
# download test fastq data
while IFS= read -r line || [[ -n $line  ]]; do

    files=($line)
    link=${files[0]}
    spe=${files[1]}

    if [ ! -f "${tmp_dir}/${spe}_basecalls.tar.gz" ]
    then
        echo "Downloading fast5 data"
        wget $link -O "${tmp_dir}/${spe}_basecalls.tar.gz"
    fi

    run_dir="${output_dir}/${spe}"
    mkdir "${run_dir}/fastq"
    tar -xvzf "${tmp_dir}/${spe}_basecalls.tar.gz" -C "${run_dir}/fastq"

    if [ -f "${run_dir}/fastq/guppy_v2.1.3-v2.2.3.fastq" ]
    then
        mv "${run_dir}/fastq/guppy_v2.1.3-v2.2.3.fastq" "${run_dir}/fastq/basecalls.fastq"
    fi
    
    if [ -f "${run_dir}/fastq/01_guppy_v2.1.3.fastq" ]
    then
        mv "${run_dir}/fastq/01_guppy_v2.1.3.fastq" "${run_dir}/fastq/basecalls.fastq"
    fi

done < $basecalls_links

echo "Done"

