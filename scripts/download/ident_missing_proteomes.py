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
"""Identify the accessions of genomes for which a proteome FASTA file (.faa) was not downloaded"""


from pathlib import Path

from saintBioutils.utilities.file_io.get_paths import get_file_paths


GENOME_DIR = Path("data/proteomes")
ACC_FILE = "data/genomic_accessions/pecto_dic_accessions"
MISSING_G_FILE = "data/missing_genomes"

with open(ACC_FILE, "r") as fh:
    accessions = fh.read().splitlines()

print(f"File of accessions contained {len(accessions)} assembly accessions")

genome_paths = get_file_paths(GENOME_DIR, suffixes=['.faa.gz', '.faa'])  # tolerate being compressed or not

downloaded_genomes = [f'{genome.name.split("_")[0]}_{genome.name.split("_")[1]}' for genome in genome_paths]

print(f"Found {len(downloaded_genomes)} directories in {GENOME_DIR}")

missing_genomes = set(accessions).difference(set(downloaded_genomes))

with open(MISSING_G_FILE, "w") as fh:
    for genome in missing_genomes:
        fh.write(f"{genome}\n")

print(
    f"Missing {len(missing_genomes)} proteome FASTA (.faa) files\n"
    "Genomic accessions for these proteome FASTA files are written to:\n"
    f"{MISSING_G_FILE}"
)