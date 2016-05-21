#!/bin/bash 

name=$1

sex -c shuang.sex $name".fits" -CHECKIMAGE_TYPE SEGMENTATION \
    -CHECKIMAGE_NAME $name"_seg.fits" \
    -CATALOG_NAME $name"_cat.fits"
