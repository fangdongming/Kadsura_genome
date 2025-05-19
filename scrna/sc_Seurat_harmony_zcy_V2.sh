#
#DAP=25DAP
#DAP=17DAP
#DAP=13_15_17DAP_try
#DAP=13DAP
#ls -d $PWD/*| grep "${DAP}_1">file
#ls -d $PWD/*| grep "${DAP}-1">file
#ls -d $PWD/*| grep "_web_">file
#ls -d $PWD/*| grep "${DAP}.SC">file

#try=try2.soupx_r0.5

for minCG in 2000; do
for maxCG in 8000; do
for lambda in 1 ;do

mkdir harmony_${try}_gene${minCG}_${maxCG}_lambda${lambda}
cd harmony_${try}_gene${minCG}_${maxCG}_lambda${lambda}

#cp ../soupx.file ./
cat>"sc_Seurat_harmony.sh"<<EOF
export R_LIBS_USER=/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/lib/R/library
#/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/wangfang/02.software/miniconda3/envs/R_4.1/bin/Rscript \
#/hwfssz1/ST_EARTH/P20Z10200N0035/USER/nianting/script/sc_Seurat_harmony_zcy_V2.1.R \
#-i soupx.file -o ./ \
#-s maize -t stem \
#-d 30 -r 0.5 -m ${minCG} -n ${maxCG} -l ${lambda} \
#-p T -k gene_ID_V5 \
#-g /hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhaocaiyao/C3_C4/marker_genes/C4_maize.sobi.seit.sevi.rice.marker_genes.all.csv

/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript \
/hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC_flower_sigleCell/03.cc_soupX_cluster/sc_Seurat_cc_zcy_V1.1.change.R \
-i /hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC_flower_sigleCell/02.soupX_cluster/02.freeze/01.harmony/r_1/harmony.all_tissue.rds \
-e /hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC/single_cell/04.ann/05.ann/cc.genes.arth.2.maize1_add_copy_copy.change.csv \
-b KCOC_ID \
-d 20 -r 1 \
-p T -k Gene_ID \
-g /hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC/single_cell/04.ann/ID.2/Gene.csv


EOF

#qsub -clear -cwd -q st.q -P P20Z10200N0035 -l vf=240g,num_proc=2 -binding linear:2 sc_Seurat_harmony.sh
qsub -cwd -l vf=200g,num_proc=2 -P P20Z10200N0035_super -binding linear:2 -q st_supermem.q sc_Seurat_harmony.sh
cd ../
done
done
done
