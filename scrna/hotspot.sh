export R_LIBS_USER=/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/lib/R/library
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript cluster.r
/hwfssz1/ST_EARTH/P20Z10200N0035/USER/xiangsunhuan/adata/Miniconda3/envs/python3.10.6/bin/python hot.py
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript umap.R

qstat -j $JOB_ID |grep cpu
