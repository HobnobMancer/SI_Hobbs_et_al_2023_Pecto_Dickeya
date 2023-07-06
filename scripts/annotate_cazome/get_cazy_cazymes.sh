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

# get_cazy_cazymes

# Get CAZy family classifications for proteins in the proteomes that are also in the local CAZymes db
# Proteins in the proteomes not in the local CAZyme db are written to FASTA Files to be parsed by dbCAN
# Write the protein accessions, genomic accessions and CAZy families to tab delimited lists

# cazevolve_get_cazy_cazymes ars
# input dir - dir containing .faa proteome files
# cazy db
# output dir - write out fasta files
# fam-genome tab delimited file
# fam-genome-protein tab delimited file
# force (-f) writing to existing output dir but do not delete (-n) data already in the output dir


cazevolve_get_cazy_cazymes \
    data/proteomes \
    data/cazy/cazy_db \
    data/cazomes/dbcan_input \
    data/cazomes/pd_fam_genomes \
    data/cazomes/pd_fam_genomes_proteins \
    -f -n
