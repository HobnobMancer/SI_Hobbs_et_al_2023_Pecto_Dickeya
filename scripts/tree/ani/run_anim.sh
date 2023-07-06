#!/usr/bin/env bash
#
# run_ANIm.sh
#
# Run ANIm analysis (using pyani v0.3) on downloaded genomes

OUTPUT_DIR=data/ani_tree
GENOME_DIR=data/genomes
PLOT_DIR=results/anim_output

# make output dir
mkdir -p $OUTPUT_DIR/logs

# Create database
pyani createdb -l $OUTPUT_DIR/logs/pyani_01_createdb.log

# Index genomes
pyani index \
  -v \
  -i $GENOME_DIR \
  -l $OUTPUT_DIR/logs/pyani_02_index.log

# Run ANIm analysis
pyani anim \
  -v \
  -i $GENOME_DIR \
  -o $OUTPUT_DIR/anim_output \
  -l $OUTPUT_DIR/logs/pyani_03_anim.log \
  --recovery \
  --name "pectobact_ANIm" \
  --classes $GENOME_DIR/classes.txt \
  --labels $GENOME_DIR/labels.txt

# Generate graphical anim output
pyani plot \
  -v \
  -l $OUTPUT_DIR/logs/pyani_04_plot.log \
  --formats pdf \
  --method seaborn \
  -o $PLOT_DIR \
  --run_id 1

pyani plot \
  -v \
  -l $OUTPUT_DIR/logs/pyani_05_plot.log \
  --formats png \
  --method seaborn \
  -o $PLOT_DIR \
  --run_id 1

pyani report \
  -v \
  -l $OUTPUT_DIR/logs/pyani_06_plot.log \
  -o $PLOT_DIR \
  --run_matrices 1

pyani plot \
  -v \
  -l $OUTPUT_DIR/logs/pyani_07_plot.log \
  --formats svg \
  --method seaborn \
  -o $PLOT_DIR \
  --run_id 1

pyani report \
  -v \
  -l $OUTPUT_DIR/logs/pyani_08_plot.log \
  -o $PLOT_DIR \
  --genomes
