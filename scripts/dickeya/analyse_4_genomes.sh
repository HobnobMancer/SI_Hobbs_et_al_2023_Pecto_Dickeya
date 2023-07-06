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

# analyse 4 genomes

# analyse the 4 Dickeya genomes that cluster together in the PCA analysis

GENOME_DIR=data/pecto_dic/4_cluster_genomes
OUTDIR=results/pecto_dic/pyani
mkdir -p $GENOME_DIR
mkdir -p $OUTDIR

ncbi-genome-download \
    bacteria \
    -s refseq \
    -F fasta \
    --flat-output \
    -A GCF_000400505.1,GCF_000023545.1,GCF_007858975.2,GCF_000406125.1,GCF_000365285.1,GCF_002846975.1,GCF_017656535.1,GCF_014642955.1,GCF_000365405.2,GCF_018361245.1,GCF_017897585.1,GCF_000406045.1,GCF_012278555.1,GCF_000406165.1,GCF_000406245.1,GCF_007210685.1,GCF_000406105.1 \
    -v \
    -o $GENOME_DIR

gunzip $GENOME_DIR/*.gz

PYANI_DIR=data/pecto_dic/4_cluster_pyani
mkdir -p $PYANI_DIR
DBPATH=data/pecto_dic/4_cluster_pyani/pyani.db

# make output dir
mkdir -p $PYANI_DIR/logs

# Create database
pyani createdb -l $PYANI_DIR/logs/pyani_01_createdb.log --dbpath $DBPATH

# Index genomes
pyani index \
  -i $GENOME_DIR \
  -v \
  -l $PYANI_DIR/logs/pyani_02_index.log

# Run ANIm analysis
pyani anim \
  -v \
  -i $GENOME_DIR \
  -o $PYANI_DIR/anim_output \
  -l $PYANI_DIR/logs/pyani_03_anim.log \
  --recovery \
  --name "pectobact_ANIm" \
  --classes $GENOME_DIR/classes.txt \
  --labels $GENOME_DIR/labels.txt \
  --dbpath $DBPATH

# Generate graphical anim output
pyani plot \
  -v \
  -l $PYANI_DIR/logs/pyani_04_plot.log \
  --formats pdf \
  --method seaborn \
  -o $OUTDIR \
  --run_id 1 \
  --dbpath $DBPATH

pyani plot \
  -v \
  -l $PYANI_DIR/logs/pyani_05_plot.log \
  --formats png \
  --method seaborn \
  -o $OUTDIR \
  --run_id 1 \
  --dbpath $DBPATH

# for GFILE in $GENOME_DIR/*.fna
# do
#     for gfile in $GENOME_DIR/*.fna
#     do
#         if [ "$GFILE" != "$gfile" ]; then
#                 base_name_1=$(basename ${GFILE})
#                 base_name_2=$(basename ${gfile})

#                 echo $OUTDIR/"${GFILE}--${gfile}"
#                 blastn \
#                     -query $GFILE \
#                     -subject $gfile \
#                     -out $OUTDIR/"${base_name_1}--${base_name_2}" \
#                     -outfmt "6 qseqid sseqid pident length slen qlen mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"

#                 echo $OUTDIR/"${gfile}--${GFILE}"
#                 blastn \
#                     -query $gfile \
#                     -subject $GFILE \
#                     -out $OUTDIR/"${base_name_2}--${base_name_1}" \
#                     -outfmt "6 qseqid sseqid pident length slen qlen mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"
#         fi
#     done
# done
