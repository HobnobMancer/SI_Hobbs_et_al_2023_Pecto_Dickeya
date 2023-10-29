# run coinfinder to pectobacteriaceae, and the pectobacterium and dickeya data sets

mkdir results/coinfinder

coinfinder \
    -i data/cazomes/fam_genomes_lists \
    -p data/tree/pecto_dic_bestTree_ \
    -a \
    -o results/coinfinder/pectodic_coinfinder_ \
    > results/coinfinder/pectodic_pectobact_run.log
