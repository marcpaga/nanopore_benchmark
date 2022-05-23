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

# # download general genomes
# while IFS= read -r line || [[ -n $line  ]]; do

#     files=($line)
#     link=${files[0]}
#     spe=${files[1]}

#     if [ ! -f "${tmp_dir}/${spe}.fna.gz" ]
#     then
#         echo "Downloading reference data"
#         wget $link -O "${tmp_dir}/${spe}.fna.gz"
#     fi

#     if [ ! -f "${genomes_dir}/${spe}.fna" ]
#     then
#         echo "Decompressing reference data"
#         gunzip -c "${tmp_dir}/${spe}.fna.gz" > "${genomes_dir}/${spe}.fna"
#     fi 
# done < $references_links_general


# # download wick specific genomes
# while IFS= read -r line || [[ -n $line  ]]; do

#     files=($line)
#     link=${files[0]}
#     spe=${files[1]}

#     if [ ! -f "${tmp_dir}/${spe}.fna.gz" ]
#     then
#         echo "Downloading reference data"
#         wget $link -O "${tmp_dir}/${spe}.fna.gz"
#     fi

#     if [ ! -f "${genomes_dir}/${spe}.fna" ]
#     then
#         echo "Decompressing reference data"
#         gunzip -c "${tmp_dir}/${spe}.fna.gz" > "${genomes_dir}/${spe}.fna"
#     fi 
# done < $references_links_specific

# download train fast5 data
while IFS= read -r line || [[ -n $line  ]]; do

    files=($line)
    link=${files[0]}
    spe=${files[1]}

    echo $spe

    if [ ! -f "${tmp_dir}/${spe}.tar.gz" ]
    then
        echo "Downloading fast5 data"
        wget $link -O "${tmp_dir}/${spe}.tar.gz"
    fi

    run_dir="${output_dir}/${spe}"
    mkdir -p ${run_dir} "${run_dir}/fast5" "${run_dir}/fastq" "${run_dir}/tmp"
    tar -xvzf "${tmp_dir}/${spe}.tar.gz" -C "${run_dir}/tmp"

    if [ ! -f "${run_dir}/read_references.fasta" ]
    then
        echo "Moving data"
        find "${run_dir}/tmp" -type f -name "*.fast5" -exec mv {} "${run_dir}/fast5" \;
        find "${run_dir}/tmp" -type f -name "read_references.fasta" -exec mv {} "${run_dir}/read_references.fasta" \;
    fi

done < $train_fast5_links

# # download test fast5 data
# while IFS= read -r line || [[ -n $line  ]]; do

#     files=($line)
#     link=${files[0]}
#     spe=${files[1]}

#     if [ ! -f "${tmp_dir}/${spe}.tar.gz" ]
#     then
#         echo "Downloading fast5 data"
#         wget $link -O "${tmp_dir}/${spe}.tar.gz"
#     fi

#     run_dir="${output_dir}/${spe}"
#     mkdir -p ${run_dir} "${run_dir}/fast5"
#     tar -xvzf "${tmp_dir}/${spe}.tar.gz" -C "${run_dir}/fast5"
    

# done < $test_fast5_links

# # download test fastq data
# while IFS= read -r line || [[ -n $line  ]]; do

#     files=($line)
#     link=${files[0]}
#     spe=${files[1]}

#     if [ ! -f "${tmp_dir}/${spe}_basecalls.tar.gz" ]
#     then
#         echo "Downloading fast5 data"
#         wget $link -O "${tmp_dir}/${spe}_basecalls.tar.gz"
#     fi

#     run_dir="${output_dir}/${spe}"
#     mkdir "${run_dir}/fastq"
#     tar -xvzf "${tmp_dir}/${spe}_basecalls.tar.gz" -C "${run_dir}/fastq"

#     if [ -f "${run_dir}/fastq/guppy_v2.1.3-v2.2.3.fastq" ]
#     then
#         mv "${run_dir}/fastq/guppy_v2.1.3-v2.2.3.fastq" "${run_dir}/fastq/basecalls.fastq"
#     fi
    
#     if [ -f "${run_dir}/fastq/01_guppy_v2.1.3.fastq" ]
#     then
#         mv "${run_dir}/fastq/01_guppy_v2.1.3.fastq" "${run_dir}/fastq/basecalls.fastq"
#     fi

# done < $basecalls_links



