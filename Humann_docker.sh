#!/bin/bash


# transfer input data

cp /fastqfiles/$1_R1_microbiome.fq .
cp /fastqfiles/$1_R2_microbiome.fq .
cat $1_R1_microbiome.fq $1_R2_microbiome.fq > $1_FR_microbiome.fq

#transfer and untar database dirs
cp /my-bowtie2db.tar.gz .
tar -xzf my-bowtie2db.tar.gz
rm my-bowtie2db.tar.gz

cp /humann_databases.tar.gz .
tar -xzf humann_databases.tar.gz
rm humann_databases.tar.gz

# create output dir 
mkdir humann_out_$1

# run humann
humann --input $1_FR_microbiome.fq --output humann_out_$1 --nucleotide-database ./databases/chocophlan --protein-database ./databases/uniref --metaphlan-options "--bowtie2db ./my-bowtie2db" --threads 8 

echo "ls ./databases\n"
ls databases -alF

# Normalize RPKs to relative abundance
humann_renorm_table --input humann_out_$1/$1_FR_microbiome_genefamilies.tsv --output humann_out_$1/$1_genefamilies-cpm.tsv --units cpm --update-snames

# Regroup genes to other functional categories
humann_regroup_table --input humann_out_$1/$1_genefamilies-cpm.tsv --output humann_out_$1/$1_level4ec-cpm.tsv --groups uniref90_level4ec

# Attaching names to features
humann_rename_table --input humann_out_$1/$1_level4ec-cpm.tsv --output humann_out_$1/$1_level4ec-cpm-named.tsv --names ec

# tar output 
tar -czf humann_out_$1.tar.gz humann_out_$1

echo "ls ./"
ls -alF

# move output to staging/
mv humann_out_$1.tar.gz /staging/

# clear other data
mkdir ignore
mv *.fq* ignore

echo "remove items\t ls ./"
ls -alF
