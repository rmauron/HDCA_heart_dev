#!/bin/bash

conda activate stereoscope

export CUDA_VISIBLE_DEVICES=0

stfiles=data/hdca_heart_ST_data.tsv

SCCOUNT=data/HDCA_heart_sc_cnt_data.tsv
VARGENES=data/HDCA_heart_sc_var.features_5000.txt
OUTPUT=output/
METADATA=data/HDCA_heart_sc_annotation.tsv

stereoscope run \
        --sc_cnt $SCCOUNT \
        --sc_labels $METADATA \
        --label_colname bio_celltype \
        -gl $VARGENES \
        -sce 50000 \
        -o $OUTPUT \
        -ste 50000 \
        --gpu \
        -stb 2048 \
        -scb 2048 \
        --st_cnt $stfiles \
        --sc_upper_bound 1000