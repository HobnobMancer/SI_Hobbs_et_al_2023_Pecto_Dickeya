# /usr/bin/env python
# -*- coding: utf-8 -*-
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
"""Gather the protein sequences of all CAZymes into a single FASTA file"""

import pandas as pd

from saintBioutils.utilities.file_io.get_paths import get_file_paths


PROTEOME_DIR = "data/proteomes"
CAZYME_SEQ_FILE = "data/cazomes/all_cazyme_seqs.fasta"
ANNOTATION_FILE = "data/cazomes/pd_fam_genomes_proteins"


cazyme_accessions = pd.read_table(ANNOTATION_FILE,header=None)
cazyme_accessions.columns=['Fam', 'Genome', 'Protein']
cazyme_accessions = set(cazyme_accessions['Protein'])
print(f"Found {len(cazyme_accessions)} unique CAZyme protein accessions")

fasta_files_paths = get_file_paths(PROTEOME_DIR, suffixes=['.fasta', '.faa'])

if len(fasta_files_paths) == 0:
    print(
        f"Found 0 fasta files in {PROTEOME_DIR}\n"
        "Check the path is correct. Terminating program"
    )
    sys.exit(1)

print(f"Retrieved {len(fasta_files_paths)} FASTA files")

seqs_dict = {}  # {id: seq_record}
for fasta_path in tqdm(fasta_files_paths, desc="Loading all seqs into memory"):
    for record in SeqIO.parse(fasta_path, "fasta"):
        seqs_dict[record.id] = record

cazyme_seqs = []
parsed_ids = set()
for acc in tqdm(cazyme_accessions, desc="Gathering CAZyme seqs"):
    if acc not in parsed_ids:
        cazyme_seqs.append(seqs_dict[acc])
        parsed_ids.add(acc)

SeqIO.write(cazyme_seqs, CAZYME_SEQ_FILE, "fasta")
print(f"Written {len(cazyme_accessions)} protein sequences to {CAZYME_SEQ_FILE}")
