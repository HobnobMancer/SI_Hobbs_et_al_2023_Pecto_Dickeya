#!/usr/bin/env bash
#
# build_anim_tree.sh
#
# coordinate scripts to build a denodrgram from the ANIm analysis with leaves labelled with the genomic acccessions

# get matrix of comparisons
pyani report \
  -v \
  --run_matrices 4 \
  --genomes \
  -o data/pectobact/tree/ani_tree/ \
  -l data/pectobact/tree/ani_tree/logs/pyani_05_reports.log

# replace genome ids with accessions in the matrix
python3 scripts/tree/ani/parse_anim_tab.py

# run the R script to build a dendrogram from the analysis
scripts/tree/ani/build_anim_tree.R
