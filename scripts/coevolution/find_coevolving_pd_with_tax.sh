# run coinfinder to pectobacteriaceae, and the pectobacterium and dickeya data sets

mkdir results/pectobact/coinfinder

coinfinder \
    -i data/cazomes/pd_fam_genomes_taxs \
    -p data/tree/pecto_dic_bestTree_taxs \
    -a \
    -o results/coinfinder/pectodic_coinfinder_taxs_rectangular \
    > results/coinfinder/coinfinder_pectodic_taxs_run.log
