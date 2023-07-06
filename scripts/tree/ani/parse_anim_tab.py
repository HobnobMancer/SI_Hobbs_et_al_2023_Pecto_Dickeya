#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Replace genome IDs with genomic accessions in the pyani anim .tab file"""


import pandas as pd
import re

from tqdm import tqdm


TAB_PATH = "data/pectobact/tree/ani_tree/matrix_identity_4.tab"
GENOME_FILE = "data/pectobact/tree/ani_tree/genomes.tab"

print("Loading tab file")
df = pd.read_table(TAB_PATH, index_col="Unnamed: 0")

# columns and rows list the genomes in the same order
genome_ids = list(df.columns)

print("Getting genome ids and accessions")
# genome ID	description	path	MD5 hash	genome length
genome_df = pd.read_table(GENOME_FILE, index_col="Unnamed: 0")

genome_dict = {}  # {id: acc}
for ri in tqdm(range(len(genome_df)), desc="Getting genome IDs and accessions"):
    genome_id = str(genome_df.iloc[ri]['genome ID'])
    genome_acc = re.findall(r"GC[AF]_\d+\d.\d+_", genome_df.iloc[ri]['path'])[0][:-1]
    genome_dict[genome_id] = genome_acc

genome_list = [genome_dict[gid.replace("Genome_id:","")] for gid in genome_ids]
df.columns = genome_list
df['genomes'] = genome_list
df = df.set_index('genomes')

df.to_csv("data/pectobact/tree/ani_tree/anim_matrix.tab", sep="\t")

print("Identifying redundant (i.e. duplicate/identical) genomes")

identical_pairs = set()

for col_acc in tqdm(genome_list[:1], desc="Identifying groups of identical genomes"):
    column = df[col_acc]
    for row_acc in genome_list:
        if row_acc == col_acc:
            continue

        if column[row_acc] == 1:
            # check if identical in the other direction
            temp_col = df[row_acc]
            if temp_col[col_acc] == 1:
                if ((col_acc, row_acc) not in identical_pairs) and ((row_acc, col_acc) not in identical_pairs):
                    identical_pairs.add((col_acc,row_acc))

print(f"Find {len(identical_pairs)} identical pairs of genomes")
if len(identical_pairs) > 0:
    print("The pairs of identical genomes are (and written to data/pectobact/tree/ani_tree/identical_genome_pairs):")
    with open("data/pectobact/tree/ani_tree/identical_genome_pairs", "w") as fh:
        for pair in identical_pairs:
            print(pair)
            fh.write(identical_pairs)
