#!/usr/bin/env bash
#
# (c) University of St Andrews 2023
# (c) University of Strathclyde 2023
# (c) James Hutton Institute 2023
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# annotate_genomes

# Annotate the genome sequences of genomes for which a proteome FASTA (.faa) file is not available 
# in the NCBI Assembly database

ACC_FILE="data/missing_genomes"
ACCESSIONS=$(cat $ACC_FILE)
GENOME_DIR=data/genomes/*
PROTEOME_DIR=data/proteomes

PROD_OUTPUT=data/prodigal
mkdir $PROD_OUTPUT

# make output dirs for prodigal
PROTEIN_DIR=$PROD_OUTPUT/proteins
CDS_DIR=$PROD_OUTPUT/cds
GBK_DIR=$PROD_OUTPUT/gbk

mkdir $PROTEIN_DIR
mkdir $CDS_DIR
mkdir $GBK_DIR

# get paths of download genome seq files
GENOME_PATHS=$(ls data/genomes)

for ACC in $ACCESSIONS
do
    echo "Parsing genome $ACC"
    for GENOME in $GENOME_DIR
    do
        if [[ "$(basename $GENOME)" = $ACC* ]]
        then
            # found genome of interest
            echo "Found the genome for $ACC at $GENOME"

            # run prodigal
            echo "Running Prodigal"
            prodigal \
                -a $PROTEIN_DIR/`basename ${GENOME%%fna}`faa \
                -d $CDS_DIR/`basename ${GENOME%%fna}`fasta \
                -i ${GENOME} \
                -o $GBK_DIR/`basename ${GENOME%%fna}`gbk

            # copy the proteome FASTA file (.faa) to the proteome dir
            echo "Copying protein FASTA file from $PROTEIN_DIR/`basename ${GENOME%%fna}`faa to $PROTEOME_DIR/`basename ${GENOME%%fna}`faa"
            cp $PROTEIN_DIR/`basename ${GENOME%%fna}`faa $PROTEOME_DIR/`basename ${GENOME%%fna}`faa
        fi
    done

done