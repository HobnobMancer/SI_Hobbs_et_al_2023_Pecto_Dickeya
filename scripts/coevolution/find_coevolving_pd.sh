# run coinfinder to pectobacteriaceae, and the pectobacterium and dickeya data sets

mkdir results/pectobact/coinfinder

coinfinder \
    -i data/cazomes/pd_fam_genomes \
    -p data/tree/pecto_dic_bestTree \
    -a \
    -o results/coinfinder/pectodic_coinfinder_ \
    > results/coinfinder/pectodic_pectobact_run.log
