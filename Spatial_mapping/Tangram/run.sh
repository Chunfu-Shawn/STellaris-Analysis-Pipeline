bash src/Tangram/spatial_mapping.sh \
--sc_h5ad "/home/user/data2/rbase/spatial_annotate_scripts/data/STW-M-Brain-Stereo-seq-1/sc.h5ad" \
--dataset "STW-M-Brain-Stereo-seq-1" \
--section "coronal_1" \
--species "Mus musculus" \
--outDir "/home/user/data2/rbase/spatial_annotate_scripts/data/STW-M-Brain-Stereo-seq-1" \
>"/home/user/data2/rbase/spatial_annotate_scripts/data/STW-M-Brain-Stereo-seq-1/log/spatial_mapping.log" \
2>"/home/user/data2/rbase/spatial_annotate_scripts/data/STW-M-Brain-Stereo-seq-1/log/Error.log"

:'
## Run cellphonedb
echo -e "`date -u`\tDetect ligand-receptor interactions for each group of cell type pairs..."

source /home/user/BGM/uplee/anaconda3/bin/activate cellphonedb

outDir=/home/user/data2/rbase/spatial_annotate_scripts/data/STW-M-Brain-Stereo-seq-1
meta=$outDir/meta.txt
counts=$outDir/sc_registered.h5ad
output_path_cpdb=$outDir/out/table
n_threads=30
species="Mus musculus"

cellphonedb method statistical_analysis $meta $counts \
        --counts-data hgnc_symbol \
        --species "$species" \
        --threads $n_threads \
        --output-path $output_path_cpdb 1>$outDir/log/cellphonedb.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tDetecting cell-cell ligand-receptor interactions failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed detecting cell-cell ligand-receptor interactions!"
fi
'